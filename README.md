# Phasing and Imputation Pipeline (GSA → SHAPEIT5 → Beagle)

This repository provides a chromosome-wise pipeline for **variant ID updating, QC, phasing (common and rare variants), and genotype imputation** starting from Illumina GSA PLINK files.

The workflow uses **PLINK**, **bcftools**, **SHAPEIT5**, and **Beagle**, and is designed for **GRCh38** reference coordinates.

---

## Environment Setup

```bash
conda activate phase_impute
```

---

## Step 1: Update GSA Variant IDs

```bash
plink --bfile input_bfile   --update-name /mnt/faruq2/lab_users/ankit/shapeit_analysis/GSA-24v2-0_A1_b150_rsids.txt   --make-bed   --out output_updated_ids
```

---

## Step 2: Convert PLINK to VCF

```bash
plink --bfile output_updated_ids   --recode vcf   --out output_updated_ids
```

```bash
bcftools view -Oz -o output_updated_ids.vcf.gz output_updated_ids.vcf
```

---

## Step 3: Remove Duplicate Variants

```bash
bcftools norm -d all -Oz   -o output_updated_ids_no_dups.vcf.gz   output_updated_ids.vcf.gz
```

---

## Step 4: Add AC, AN, and AF Tags

```bash
bcftools plugin fill-tags -Oz   -o output_updated_ids_no_dups_AC_AN_AF.vcf.gz   -- output_updated_ids_no_dups.vcf.gz   -t AC,AN,AF
```

---

## Step 5: Rename Chromosomes (1 → chr1)

```bash
bcftools annotate --rename-chrs update_chr_name.txt -Oz   -o output_updated_ids_no_dups_AC_AN_AF_update_chr_name.vcf.gz   output_updated_ids_no_dups_AC_AN_AF.vcf.gz
```

```bash
tabix -p vcf output_updated_ids_no_dups_AC_AN_AF_update_chr_name.vcf.gz
```

---

## Step 6: Split VCF by Chromosome

```bash
#!/bin/bash

INPUT_VCF=output_updated_ids_no_dups_AC_AN_AF_update_chr_name.vcf.gz
OUTPUT_DIR=desired_folder

for CHR in chr{1..22}; do
  bcftools view -r $CHR $INPUT_VCF -Oz -o $OUTPUT_DIR/${CHR}.vcf.gz
  tabix -p vcf $OUTPUT_DIR/${CHR}.vcf.gz
done
```

---

## Step 7: Phase Common Variants (SHAPEIT5)

```bash
#!/bin/bash

for CHR in chr{1..22}; do
  CHUNKS=shapeit5/resources/chunks/b38/20cM/chunks_${CHR}.txt
  VCF=chr_wise_vcf/${CHR}.vcf.gz
  MAP=shapeit5/resources/maps/b38/${CHR}.b38.gmap.gz
  REF=aggarwal/aggarwal_157/aggarwal_ref_phasing_unique_singlersid_firstrsid.vcf.gz

  while read LINE; do
    REG=$(echo $LINE | awk '{print $3}')
    CHUNK_NBR=$(echo $LINE | awk '{print $1}')
    OUT=chr_wise_phase_common/${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.bcf

    SHAPEIT5_phase_common       --input $VCF       --map $MAP       --reference $REF       --region $REG       --filter-maf 0.001       --thread 20       --output $OUT

    bcftools index -f $OUT
    bcftools view $OUT -Oz -o ${OUT%.bcf}.vcf.gz
    tabix -p vcf ${OUT%.bcf}.vcf.gz
  done < ${CHUNKS}
done
```

---

## Step 8: Ligate Common Variant Chunks

```bash
#!/bin/bash

for CHR in chr{1..22}; do
  ls -1v chr_wise_phase_common/${CHR}.chunk_*.shapeit5_common.bcf     > list_${CHR}.txt

  SHAPEIT5_ligate     --input list_${CHR}.txt     --output chr_wise_common_ligate/${CHR}.shapeit5_common_ligate.bcf     --thread 20

  bcftools index -f chr_wise_common_ligate/${CHR}.shapeit5_common_ligate.bcf
  bcftools view chr_wise_common_ligate/${CHR}.shapeit5_common_ligate.bcf -Oz     -o chr_wise_common_ligate/${CHR}.shapeit5_common_ligate.vcf.gz
  tabix -p vcf chr_wise_common_ligate/${CHR}.shapeit5_common_ligate.vcf.gz
done
```

---

## Step 9: Phase Rare Variants

```bash
#!/bin/bash

for CHR in chr{1..22}; do
  CHUNKS=shapeit5/resources/chunks/b38/4cM/chunks_${CHR}.txt
  VCF=chr_wise_vcf/${CHR}.vcf.gz
  MAP=shapeit5/resources/maps/b38/${CHR}.b38.gmap.gz
  SCAF=chr_wise_common_ligate/${CHR}.shapeit5_common_ligate.bcf

  while read LINE; do
    CHUNK_NBR=$(echo $LINE | awk '{print $1}')
    SCAFFOLD_REG=$(echo $LINE | awk '{print $3}')
    INPUT_REG=$(echo $LINE | awk '{print $4}')
    OUT=chr_wise_phase_rare/${CHR}.chunk_${CHUNK_NBR}.shapeit5_rare.bcf

    SHAPEIT5_phase_rare       --input $VCF       --map $MAP       --scaffold $SCAF       --scaffold-region $SCAFFOLD_REG       --input-region $INPUT_REG       --thread 20       --output $OUT

    bcftools index -f $OUT
    bcftools view $OUT -Oz -o ${OUT%.bcf}.vcf.gz
    tabix -p vcf ${OUT%.bcf}.vcf.gz
  done < ${CHUNKS}
done
```

---

## Step 10: Concatenate Rare Variant Chunks

```bash
#!/bin/bash

for CHR in chr{1..22}; do
  ls -1v chr_wise_phase_rare/${CHR}.chunk_*.shapeit5_rare.bcf     > concat_${CHR}.txt

  bcftools concat -n -f concat_${CHR}.txt --threads 20     -o chr_wise_rare_concatenate/${CHR}.full.shapeit5_rare.bcf

  bcftools index chr_wise_rare_concatenate/${CHR}.full.shapeit5_rare.bcf
  bcftools view chr_wise_rare_concatenate/${CHR}.full.shapeit5_rare.bcf -Oz     -o chr_wise_rare_concatenate/${CHR}.full.shapeit5_rare.vcf.gz
  tabix -p vcf chr_wise_rare_concatenate/${CHR}.full.shapeit5_rare.vcf.gz
done
```

---

## Step 11: Concatenate All Chromosomes

```bash
bcftools concat --threads 20 -W -Oz   -o chr_all_phased_(common/rare).vcf.gz   chr_wise_*/*.vcf.gz
```

---

## Step 12: Imputation with Beagle

```bash
#!/bin/bash

export _JAVA_OPTIONS="-Xmx32g"

for CHR in {1..22}; do
  beagle     gt=all_chr_phased_(common/rare).vcf.gz     ref=aggarwal_ref_phasing_unique_singlersid_firstrsid.vcf.gz     map=map/update_chr/plink_chr_update.chr${CHR}.GRCh38.map     chrom=chr${CHR}     out=beagle_imputation/imputation_(common/rare)/chr${CHR}.imputed     nthreads=20

  tabix -p vcf beagle_imputation/imputation_(common/rare)/chr${CHR}.imputed.vcf.gz
done
```

---

## Step 13: Concatenate Imputed Variants

```bash
bcftools concat --threads 20 -W -Oz   -o chr_all_(common/rare)_variant_imputation.vcf.gz   beagle_imputation/imputation_(common/rare)/chr*.imputed.vcf.gz
```

---

## Step 14: Filter by Imputation Quality

```bash
bcftools view -i 'INFO/DR2 > x' -Oz   -o high_quality_imputed.vcf.gz   chr_all_(common/rare)_variant_imputation.vcf.gz
```

---

## Reference Genome

- GRCh38

## Tools Used

- PLINK
- bcftools
- SHAPEIT5
- Beagle
