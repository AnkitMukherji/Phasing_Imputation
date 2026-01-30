#!/bin/bash

for CHR in chr{1..22};
do
	CHUNKS=/mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_phase_rare/${CHR}.chunk_*.shapeit5_rare.bcf
	OUT=/mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_rare_concatenate/${CHR}.full.shapeit5_rare.bcf
	
	ls -1v ${CHUNKS} > /mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_rare_concatenate/concat_list_${CHR}.txt
	bcftools concat -n -f /mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_rare_concatenate/concat_list_${CHR}.txt --threads 20 -o ${OUT}
	bcftools index ${OUT}
	bcftools view $OUT -Oz -o ${OUT%.bcf}.vcf.gz
	tabix -p vcf ${OUT%.bcf}.vcf.gz
done
