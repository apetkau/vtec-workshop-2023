Data from <https://wwwnc.cdc.gov/eid/article/22/12/16-0017_article>. Specifically, the reads are from this BioProject <https://www.ncbi.nlm.nih.gov/bioproject/PRJEB11886/>. The reference file from the tree is <https://www.ncbi.nlm.nih.gov/nuccore/NC_002695.2?report=fasta>.

Files
* SRA and BioSample mapping: [name-mapping.txt](name-mapping.txt)
* Phylogenetic tree: [snvphyl-results/phylogeneticTree.newick](snvphyl-results/phylogeneticTree.newick)
* Other SNVPhyl results files: [snvphyl-results/](snvphyl-results)
* Table 1 from paper: [table1.csv](table1.csv)

The below information is the steps for downloading and pre-processing reads.

# 1. Download run info file

Went to <https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB11886>, selected **SRA Experiments** and exported to an `SraRunInfo.csv` file.

Went to <https://www.ncbi.nlm.nih.gov/bioproject/PRJEB11886> and downloaded assembly details table to `PRJEB11886_AssemblyDetails.txt` file.

# 2. Download reads

```
for i in `cut -f 1 -d ',' PRJEB11886-SraRunInfo.csv | tail -n+2`; do prefetch $i; done

for i in `cut -f 1 -d ',' PRJEB11886-SraRunInfo.csv | tail -n+2`; do fasterq-dump $i; done

mv *.fastq data/
ls data/*.fastq -1 | sed -e 's/^/gzip /' > commands.txt

parallel -j 12 -a commands.txt
```

# 3. Download reference

Downloaded reference genome <https://www.ncbi.nlm.nih.gov/nuccore/NC_002695.2?report=fasta> into `reference/Sakai.fasta`.

# 4. Rename genomes

```
for srr in `cut -f 1 -d ',' PRJEB11886-SraRunInfo.csv | tail -n+2`; do biosample=`grep $srr PRJEB11886-SraRunInfo.csv | cut -d ',' -f 26`; strain=`grep $biosample PRJEB11886_AssemblyDetails.txt | cut -f 5 -d $'\t'`; echo $srr,$biosample,$strain; done > name-mapping.txt

for i in data/*_1.fastq.gz; do n=`basename $i`; b=`basename $i _1.fastq.gz`; strain=`grep $b name-mapping.txt | cut -f 3 -d ','`; ln -s ../data/$n data-renamed/${strain}_1.fastq.gz; done
for i in data/*_2.fastq.gz; do n=`basename $i`; b=`basename $i _2.fastq.gz`; strain=`grep $b name-mapping.txt | cut -f 3 -d ','`; ln -s ../data/$n data-renamed/${strain}_2.fastq.gz; done
```

# 5. Clean data

```
for i in `cut -f 3 -d ',' name-mapping.txt`; do echo fastp --detect_adapter_for_pe --in1 data-renamed/${i}_1.fastq.gz --in2 data-renamed/${i}_2.fastq.gz --out1 data-cleaned/${i}_1.fastq --out2 data-cleaned/${i}_2.fastq --json data-cleaned/${i}.json --html data-cleaned/${i}.html --report_title 'fastp report: ${i}' \&\& gzip data-cleaned/${i}_1.fastq \&\& gzip data-cleaned/${i}_2.fastq; done > commands-fastp.txt

parallel -j 12 -a commands-fastp.txt
```

# 6. Downsample data

```
# Create file with coverages/downsample proportions
MAX_COV=40; REF_LENGTH=`conda run -n bioperl bp_seq_length reference/Sakai.fasta | cut -d ' ' -f 2| tr -d '\n'`; (for i in data-cleaned/*_1.fastq.gz; do name=`basename $i _1.fastq.gz`; forward=`zcat data-cleaned/${name}_1.fastq.gz | sed -n 2~4p|tr -d '\n'|wc -c`; reverse=`zcat data-cleaned/${name}_2.fastq.gz | sed -n 2~4p | tr -d '\n'|wc -c`; cov=`echo "($forward+$reverse)/${REF_LENGTH}"|bc -l`; ratio=`echo "(${MAX_COV}/$cov)"|bc -l`; echo -e "$name\t$forward\t$reverse\t${REF_LENGTH}\t$cov\t${MAX_COV}\t$ratio"; done) | tee /tmp/cov.txt
sort -k5,5nr /tmp/cov.txt > coverages.txt

# Create downsample commands
for i in `cat coverages.txt | cut -f 1 -d $'\t'`; do ratio=`grep $i coverages.txt | cut -f 7 -d $'\t'`; if (( $(echo "$ratio < 1" | bc -l) )); then echo seqtk sample -s 121 data-cleaned/${i}_1.fastq.gz $ratio \| gzip \> data-downsampled/${i}_1.fastq.gz; echo seqtk sample -s 121 data-cleaned/${i}_2.fastq.gz $ratio \| gzip \> data-downsampled/${i}_2.fastq.gz; else echo cp data-cleaned/${i}_1.fastq.gz data-downsampled\; cp data-cleaned/${i}_2.fastq.gz data-downsampled; fi done > downsample-commands.txt
# Downsample in parallel
parallel -j 16 -a downsample-commands.txt
```

# 7. Rename sequences for SNVPhyl

```
prename 's/_1.fastq/_R1.fastq/' data-downsampled/*.fastq.gz
prename 's/_2.fastq/_R2.fastq/' data-downsampled/*.fastq.gz
```

# 8. Run SNVPhyl

```
pushd SNVPhyl_Nextflow/

nextflow run snvphyl.nf --outdir ./results -c snvphyl.config --refgenome ../reference/Sakai.fasta --input_reads ../data-downsampled/ --window_size 10 --density_threshold 2
```

# 9. Download table of metadata

Downloaded `table1.csv` and `table2.csv` from <https://wwwnc.cdc.gov/eid/article/22/12/16-0017_article>. Split first column of Table 2 into "Source" and "Strain". Got read of unneeded columns in table2 (which caused issues in CSV format).

Concatenated tables into `table1_2.csv` using Python:

```python
import pandas as pd
df1 = pd.read_csv('table1.csv')
df2 = pd.read_csv('table2.csv')
df = pd.concat([df1,df2])
df.to_csv('table1_2.csv', index=False)
```

# 10. Upload data to IRIDA

I uploaded sequence reads to IRIDA and ran the AMR Detection and ECTyper pipelines. I exported the metadata table.
