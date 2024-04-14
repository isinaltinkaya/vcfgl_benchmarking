

bcftools(){
	/maps/projects/lundbeck/scratch/pfs488/ibdgl/ibdgl_paper_analyses/tools/bcftools/bcftools "${@}"
}

# nhet=$(bcftools view ../sim/sim_vcfgl_2401/model_OutOfAfrica_3G09/contig_chr22/vcfgl_platform0_gl1_qs00/sim_vcfgl_2401-OutOfAfrica_3G09-chr22-rep0-d20-e0.002-vcfgl_platform0_gl1_qs00.truth.bcf | bcftools +fill-tags -- -t all | bcftools query -f'%AC_Het\n' | datamash sum 1 )
bcftools view ../sim/sim_vcfgl_2401/model_OutOfAfrica_3G09/contig_chr22/vcfgl_platform0_gl1_qs00/sim_vcfgl_2401-OutOfAfrica_3G09-chr22-rep0-d20-e0.002-vcfgl_platform0_gl1_qs00.truth.bcf | bcftools +fill-tags -- -t all | bcftools query -f'%AC_Het\n' | datamash sum 1 
# nhet=$(bcftools view sim/sim_vcfgl_2312/model_OutOfAfrica_3G09/contig_chr22/vcfgl_gl2/sim_vcfgl_2312-OutOfAfrica_3G09-chr22-rep0-d1-e0.002-qs0_0.truth.bcf | bcftools +fill-tags -- -t all | bcftools query -f'%AC_Het\n' | datamash sum 1 )

# nSites=$( bcftools view -H sim/sim_vcfgl_2312/model_OutOfAfrica_3G09/contig_chr22/vcfgl_gl2/sim_vcfgl_2312-OutOfAfrica_3G09-chr22-rep0-d1-e0.002-qs0_0.truth.bcf | wc -l )
nSites=$( bcftools view -H ../sim/sim_vcfgl_2401/model_OutOfAfrica_3G09/contig_chr22/vcfgl_platform0_gl1_qs00/sim_vcfgl_2401-OutOfAfrica_3G09-chr22-rep0-d20-e0.002-vcfgl_platform0_gl1_qs00.truth.bcf  | wc -l )
bcftools view -H ../sim/sim_vcfgl_2401/model_OutOfAfrica_3G09/contig_chr22/vcfgl_platform0_gl1_qs00/sim_vcfgl_2401-OutOfAfrica_3G09-chr22-rep0-d20-e0.002-vcfgl_platform0_gl1_qs00.truth.bcf  | wc -l 

printf "${nhet}\t${nSites}\n"
printf "${nhet}\t${nSites}\n" | awk '{print $1 / ($2 * 100)}'
printf "${nSites}\n" | awk '{print ($1 * 100)}'
# total heterozygosity rate:
# 0.10967
#
#
# # average per-site heterozygosity rate:
# avghet=$(bcftools view sim/sim_vcfgl_2312/model_OutOfAfrica_3G09/contig_chr22/vcfgl_gl2/sim_vcfgl_2312-OutOfAfrica_3G09-chr22-rep0-d1-e0.002-qs0_0.truth.bcf | bcftools +fill-tags -- -t all | bcftools query -f'%AC_Het\n' | awk '{print $1/100}'|datamash mean 1)
#
# echo $avghet
#
# 0.10967026464303



