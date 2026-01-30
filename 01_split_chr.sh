#!/bin/bash

INPUT_VCF=/mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/gsa_validation_rel025_snpflipped_qcfiltered_update_ids_no_dups_AC_AN_AF_update_chr_name.vcf.gz
OUTPUT_DIR=/mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_vcf

for CHR in chr{1..22}; 
do
    OUTPUT_VCF=$OUTPUT_DIR/${CHR}.vcf.gz
    bcftools view -r $CHR $INPUT_VCF -Oz -o $OUTPUT_VCF
    tabix -p vcf $OUTPUT_VCF
done
