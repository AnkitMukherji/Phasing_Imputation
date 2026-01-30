#!/bin/bash

for CHR in chr{1..22};
do
	CHUNKS=/mnt/faruq2/lab_users/ankit/shapeit_analysis/shapeit5/resources/chunks/b38/20cM/chunks_${CHR}.txt
	VCF=/mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_vcf/${CHR}.vcf.gz
	MAP=/mnt/faruq2/lab_users/ankit/shapeit_analysis/shapeit5/resources/maps/b38/${CHR}.b38.gmap.gz

	while read LINE;
	do
		REG=$(echo $LINE | awk '{print $3}')
		CHUNK_NBR=$(echo $LINE | awk '{print $1}')
		OUT=/mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_phase_common/${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.bcf
		LOG=/mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_phase_common/${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.log
		REF=/mnt/faruq2/lab_users/ankit/shapeit_analysis/aggarwal/aggarwal_157/aggarwal_ref_phasing_unique_singlersid_firstrsid.vcf.gz
		
		SHAPEIT5_phase_common --input $VCF --map $MAP --reference $REF --output $OUT --thread 20 --log $LOG --filter-maf 0.001 --region $REG
		bcftools index -f $OUT
		bcftools view $OUT -Oz -o ${OUT%.bcf}.vcf.gz
		tabix -p vcf ${OUT%.bcf}.vcf.gz
	done < ${CHUNKS}
done
