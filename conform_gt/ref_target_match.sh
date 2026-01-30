#!/bin/bash

for i in chr{1..22};
do
	java -Xmx10g -jar conform-gt.24May16.cee.jar \
	ref=/staging/ankit/shapeit_analysis/aggarwal/aggarwal_ref_phasing_unique_singlersid_firstrsid.vcf.gz \
	gt=/staging/ankit/shapeit_analysis/aggarwal_cases/new_no_maf/aggarwal_case.vcf.gz \
	chrom=$i \
	out=mod.$i
done
