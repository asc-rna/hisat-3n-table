# hisat-3n-table

single thread hisat-3n-table specificaly used in ASC25

## Usage

Use standard input, output:

```sh
./hisat-3n-table u /mnt/ramdisk/rna/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

Use SAM file as input, tsv file as output:

```sh
./hisat-3n-table m /mnt/ramdisk/rna/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa < /mnt/ramdisk/rna/output/SRR23538290.mRNA.genome.mapped.sorted.dedup.filtered.sam > /mnt/ramdisk/rna/output/SRR23538290.filtered_multi.tsv
```

Use BAM file as input, tsv file as output:

```sh
~/rna-workflow/samtools/samtools view -e "rlen<12000" -h /mnt/ramdisk/rna/output/SRR23538290.mRNA.genome.mapped.sorted.dedup.filtered.bam | ./hisat-3n-table m /mnt/ramdisk/rna/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa > /mnt/ramdisk/rna/output/SRR23538290.filtered_multi.tsv
```

