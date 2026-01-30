#!/bin/bash

export _JAVA_OPTIONS="-Xmx32g"

for CHR in {1..22};
do
	beagle gt=/mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/all_chr_phased_common.vcf.gz \
	ref=/mnt/faruq2/lab_users/ankit/shapeit_analysis/aggarwal/aggarwal_157/aggarwal_ref_phasing_unique_singlersid_firstrsid.vcf.gz \
	map=/mnt/faruq2/lab_users/ankit/shapeit_analysis/map/update_chr/plink_chr_update.chr${CHR}.GRCh38.map \
	chrom=chr${CHR} \
	out=/mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/beagle_imputation/imputation_common/chr${CHR}.imputed \
	nthreads=20

	tabix -p vcf /mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/beagle_imputation/imputation_common/chr${CHR}.imputed.vcf.gz
done
