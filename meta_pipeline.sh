
cd ~/meta

time fastqc seq/*.gz -t 6 # 9s

multiqc -d seq/ -o result/qc

# remove contaminate sequence of host
time parallel -j 3 --xapply \
  'kneaddata -i {1} -i {2} \
  -o temp/qc -v -t 12 --remove-intermediate-output \
  --trimmomatic ~/SOFT/trimmomatic/ --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" \
  --bowtie2-options "--very-sensitive --dovetail " --no-discordant -db db/bowtie2/Homo_sapiens' \
 ::: seq/*_1.fq ::: seq/*_2.fq

# quality control
kneaddata_read_count_table --input temp/qc --output result/01kneaddata_sum.txt
cat result/01kneaddata_sum.txt

fastqc temp/qc/*_1_kneaddata_paired_* -t 6
multiqc -d temp/qc/ -o result/qc/

# Reference-based Metagenomic Pipeline (HUMAnN2)

mkdir -p temp/concat
for i in `tail -n+2 result/design.txt | cut -f 1`;do \
  cat temp/qc/${i}*data_paired* > temp/concat/${i}.fq; done
ll temp/concat/*.fq

mkdir -p temp/humann2
time parallel -j 3 \
  'humann2 --input {}  \
  --output temp/humann2/ ' \
  ::: temp/concat/*.fq > temp/log &
cat temp/log

mkdir -p result/humann2

humann2_join_tables --input temp/humann2 --file_name pathabundance --output result/humann2/pathabundance.tsv
sed -i 's/_Abundance//g' result/humann2/pathabundance.tsv

humann2_renorm_table --input result/humann2/pathabundance.tsv --units relab \
  --output result/humann2/pathabundance_relab.tsv

humann2_split_stratified_table --input result/humann2/pathabundance_relab.tsv \
  --output result/humann2/

## Taxonomic profiling

mkdir -p result/metaphlan2

merge_metaphlan_tables.py temp/humann2/*_humann2_temp/*_metaphlan_bugs_list.tsv | \
  sed 's/_metaphlan_bugs_list//g' > result/metaphlan2/taxonomy.tsv

metaphlan_to_stamp.pl result/metaphlan2/taxonomy.tsv > result/metaphlan2/taxonomy.spf

# metaphlan2 to graphlan
export2graphlan.py --skip_rows 1,2 -i result/metaphlan2/taxonomy.tsv\
  --tree temp/merged_abundance.tree.txt \
  --annotation temp/merged_abundance.annot.txt \
  --most_abundant 100 --abundance_threshold 1 --least_biomarkers 10 \
  --annotations 5,6 --external_annotations 7 --min_clade_size 1
# graphlan annotation
graphlan_annotate.py --annot temp/merged_abundance.annot.txt \
  temp/merged_abundance.tree.txt  temp/merged_abundance.xml
# output PDF figure, annoat and legend
graphlan.py temp/merged_abundance.xml result/metaphlan2/graphlan.pdf \
  --external_legends --dpi 300 

# LEfSe
sed '1 s/p[0-9]*//g' result/metaphlan2/taxonomy.tsv | grep -v '#' > result/metaphlan2/lefse.txt

lefse-format_input.py result/metaphlan2/lefse.txt temp/input.in -c 1 -o 1000000

run_lefse.py temp/input.in temp/input.res

lefse-plot_cladogram.py temp/input.res result/metaphlan2/lefse_cladogram.pdf --format pdf --dpi 300 

lefse-plot_res.py temp/input.res result/metaphlan2/lefse_res.pdf --format pdf --dpi 300

grep -v '-' temp/input.res | sort -k3,3n
lefse-plot_features.py -f one --feature_name "k__Bacteria.p__Firmicutes" --format pdf \
  temp/input.in temp/input.res result/metaphlan2/Firmicutes.pdf 

lefse-plot_features.py -f diff --archive none --format pdf \
  temp/input.in temp/input.res result/metaphlan2/features


## De novo Metagenomic Pipeline

mkdir -p temp/kraken2

time parallel -j 3 \
  'kraken2 --db /db/kraken2 --paired temp/qc/{1}_1_kneaddata_paired*.fastq \
  --threads 3 --use-names --use-mpa-style --report-zero-counts \
  --report temp/kraken2/{1}_report \
  --output temp/kraken2/{1}_output' \
  ::: `tail -n+2 result/design.txt | cut -f 1`

mkdir -p result/kraken2
parallel -j 6 \
  'cut -f 2 temp/kraken2/{1}_report | sed "1 s/^/{1}\n/" > temp/kraken2/{1}_count ' \
  ::: `tail -n+2 result/design.txt | cut -f 1`
header=`tail -n 1 result/design.txt | cut -f 1`
cut -f 1 temp/kraken2/${header}_report | sed "1 s/^/Taxonomy\n/" > temp/kraken2/0header_count
paste temp/kraken2/*count > result/kraken2/taxonomy_count.txt

# Megahit

rm -rf temp/megahit
time megahit -t 9 \
  -1 `ls temp/qc/*_1_kneaddata_paired_1.fastq|tr '\n' ','|sed 's/,$//'` \
  -2 `ls temp/qc/*_1_kneaddata_paired_2.fastq|tr '\n' ','|sed 's/,$//'` \
  -o temp/megahit 

## Gene annotation & quantitfy

# prokka

time prokka temp/megahit/final.contigs.fa --outdir temp/prokka \
  --prefix mg --metagenome --kingdom Archaea,Bacteria,Mitochondria,Viruses \
  --force --cpus 9

# cd-hit
mkdir -p temp/NR

time cd-hit-est -i temp/prokka/mg.ffn -o temp/NR/mg.ffn.nr \
  -aS 0.9 -c 0.95 -G 0 -M 0 -T 9 -n 8 -d 0 -g 1
  
mkdir -p result/NR
ln temp/NR/mg.ffn.nr result/NR/nucleotide.fa

transeq -sequence result/NR/nucleotide.fa -outseq result/NR/protein.fa

sed -i 's/_1 / /' result/NR/protein.fa

# salmon
mkdir -p temp/salmon

salmon index -t result/NR/nucleotide.fa -p 9 --type quasi -k 31 \
  -i temp/salmon/index 

time parallel -j 3 \
  'salmon quant -i temp/salmon/index -l A -p 3 --meta \
  -1 temp/qc/{1}_1_kneaddata_paired_1.fastq -2 temp/qc/{1}_1_kneaddata_paired_2.fastq \
  -o temp/salmon/{1}.quant' \
  ::: `tail -n+2 result/design.txt | cut -f 1`

mkdir -p result/salmon
salmon quantmerge --quants temp/salmon/*.quant -o result/salmon/gene.TPM
salmon quantmerge --quants temp/salmon/*.quant --column NumReads -o result/salmon/gene.count
sed -i '1 s/.quant//g' result/salmon/gene.*



## Functional profiling KEGG/eggNOG

# KEGG

mkdir -p  temp/kegg
time diamond blastp --db /db/kegg_v91/kegg --query result/NR/protein.fa \
  --outfmt 6 --threads 9 --max-target-seqs 1 --quiet \
  --out temp/kegg_v91/gene_diamond.f6
  
mkdir -p temp/kegg/keggresult

cut -f 1,2 temp/kegg_v91/gene_diamond.f6 | uniq | sed 's/|/\t/g' | cut -f 1,3 | \
	cut -f 1,2 -d '_' |sed '1 i Name\tKO' > temp/kegg/keggresult/gene_kegg.list

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' temp/kegg/keggresult/gene_kegg.list \
  result/salmon/gene.count | sed '/\t$/d' > temp/kegg/keggresult/gene_kegg.count

mat_gene2ko.R -i temp/kegg/keggresult/gene_kegg.count -o result/kegg/keggresult/kegg_tab


# eggNOG

mkdir -p temp/eggnog
time emapper.py -m diamond --no_annot --no_file_comments --data_dir /db/eggnog \
  --cpu 9 -i result/NR/protein.fa -o temp/eggnog/protein.fa --override

time emapper.py --annotate_hits_table temp/eggnog/protein.fa.emapper.seed_orthologs --no_file_comments \
		-o temp/eggnog/output --cpu 9 --data_dir /db/eggnog --override

mkdir -p result/eggnog
sed '1 i Name\teggNOG\tEvalue\tScore\tGeneName\tGO\tKO\tBiGG\tTax\tOG\tBestOG\tCOG\tAnnotation' \
  temp/eggnog/output.emapper.annotations > temp/eggnog/output

cut -f 1,12 temp/eggnog/output|cut -f 1 -d ','|grep -v -P '\t$'>temp/eggnog/1cog.list

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' temp/eggnog/1cog.list result/salmon/gene.count | \
	sed '/\t$/d' | sed '1 s/COG/KO/' > temp/eggnog/gene_cog.count

mat_gene2ko.R -i temp/eggnog/gene_cog.count -o result/eggnog/cogtab -n 1000000

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3} NR>FNR{print a[$1],$0}' \
  /db/eggnog/COG.anno result/eggnog/cogtab.count > result/eggnog/cogtab.count.spf

cut -f 1,7 temp/eggnog/output|cut -f 1 -d ','|grep -v -P '\t$' > temp/eggnog/2ko.list

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' temp/eggnog/2ko.list \
  result/salmon/gene.count | sed '/\t$/d' > temp/eggnog/gene_ko.count

mat_gene2ko.R -i temp/eggnog/gene_ko.count -o result/eggnog/kotab -n 1000000

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' /db/eggnog/KO.anno \
  result/eggnog/kotab.count | sed 's/^\t/Undescription\t/' > result/eggnog/kotab.count.spf
