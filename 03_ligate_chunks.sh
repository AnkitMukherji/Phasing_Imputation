#!bin/bash

for CHR in chr{1..22};
do
	ls -1v /mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_phase_common/${CHR}.chunk_*.shapeit5_common.bcf > /mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_common_ligate/list_ligate.${CHR}.txt
	SHAPEIT5_ligate --input /mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_common_ligate/list_ligate.${CHR}.txt --output /mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_common_ligate/${CHR}.shapeit5_common_ligate.bcf --thread 20
	bcftools index -f /mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_common_ligate/${CHR}.shapeit5_common_ligate.bcf

	bcftools view /mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_common_ligate/${CHR}.shapeit5_common_ligate.bcf -Oz -o /mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_common_ligate/${CHR}.shapeit5_common_ligate.vcf.gz
	tabix -p vcf /mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_common_ligate/${CHR}.shapeit5_common_ligate.vcf.gz
done
