#!/bin/bash

for CHR in chr{1..22};
do
        CHUNKS=/mnt/faruq2/lab_users/ankit/shapeit_analysis/shapeit5/resources/chunks/b38/4cM/chunks_${CHR}.txt
        VCF=/mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_vcf/${CHR}.vcf.gz
        MAP=/mnt/faruq2/lab_users/ankit/shapeit_analysis/shapeit5/resources/maps/b38/${CHR}.b38.gmap.gz
	SCAF=/mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_common_ligate/${CHR}.shapeit5_common_ligate.bcf

        while read LINE;
        do
                CHUNK_NBR=$(echo $LINE | awk '{print $1;}')

		SCAFFOLD_REG=$(echo $LINE | awk '{print $3;}')
		SCAFFOLD_REG_START=$(echo ${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f1)
		SCAFFOLD_REG_END=$(echo ${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f2)
		SCAFFOLD_REG_NAME=${CHR}_${SCAFFOLD_REG_START}_${SCAFFOLD_REG_END}

		INPUT_REG=$(echo $LINE | awk '{print $4;}')
		INPUT_REG_START=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f1)
		INPUT_REG_END=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f2)
		INPUT_REG_NAME=${CHR}_${INPUT_REG_START}_${INPUT_REG_END}
			
                OUT=/mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_phase_rare/${CHR}.chunk_${CHUNK_NBR}.shapeit5_rare.bcf
                LOG=/mnt/faruq2/lab_users/ankit/shapeit_analysis/sca12_validation_30_1_26/chr_wise_phase_rare/${CHR}.chunk_${CHUNK_NBR}.shapeit5_rare.log

                SHAPEIT5_phase_rare --input $VCF --map $MAP --output $OUT --thread 20 --scaffold $SCAF --scaffold-region $SCAFFOLD_REG --input-region $INPUT_REG
                bcftools index -f $OUT
                bcftools view $OUT -Oz -o ${OUT%.bcf}.vcf.gz
                tabix -p vcf ${OUT%.bcf}.vcf.gz
        done < ${CHUNKS}
done
