# 1. Download run info file

Went to <https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB11886>, selected **SRA Experiments** and exported to an `SraRunInfo.csv` file.

# 2. Download reads

```
for i in `cut -f 1 -d ',' PRJEB11886-SraRunInfo.csv | tail -n+2`; do prefetch $i; done

for i in `cut -f 1 -d ',' PRJEB11886-SraRunInfo.csv | tail -n+2`; do fasterq-dump $i; done

mv *.fastq data/
ls data/*.fastq -1 | sed -e 's/^/gzip /' > commands.txt

ls *.fastq | parallel -j 12 ::: gzip {}
```

# 3. Download reference

Downloaded reference genome <https://www.ncbi.nlm.nih.gov/nuccore/NC_002695.2?report=fasta> into `reference/Sakai.fasta`.

# 4. Rename genomes


