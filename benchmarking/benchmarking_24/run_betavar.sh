


for dp in 0.1 0.5 1 2 10 20;do
	{ time /maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/tools/vcfgl/vcfgl -i sim_vcfgl_2312-OutOfAfrica_3G09-chr22.vcf --error-qs 2 --beta-variance 1e-5 -s 1 -d ${dp} -e 0.002  -o outd${dp}_qs2_beta5 2> outd${dp}_qs2_beta5.log ; } 2> time_outd${dp}_qs2_beta5
done


#
# for dp in 0.1 0.5 1 2 10 20;do
	# { time /maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/benchmarking/vcfgl/vcfgl -i sim_vcfgl_2312-OutOfAfrica_3G09-chr22.vcf --threads 4 --error-qs 2 --beta-variance 1e-5 -s 1 -d ${dp} -e 0.002  -o outd${dp}_qs2_beta5_t4 2> outd${dp}_qs2_beta5_t4.log ; } 2> time_outd${dp}_qs2_beta5_t4
# done
#
# for dp in 0.1 0.5 1 2 10 20;do
	# { time /maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/benchmarking/vcfgl/vcfgl -i sim_vcfgl_2312-OutOfAfrica_3G09-chr22.vcf --threads 2 --error-qs 2 --beta-variance 1e-5 -s 1 -d ${dp} -e 0.002  -o outd${dp}_qs2_beta5_t2 2> outd${dp}_qs2_beta5_t2.log ; } 2> time_outd${dp}_qs2_beta5_t2
# done
#
# for dp in 20;do
	# { time /maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/benchmarking/vcfgl/vcfgl -i sim_vcfgl_2312-OutOfAfrica_3G09-chr22.vcf --threads 4 --error-qs 2 --beta-variance 1e-5 -s 1 -d ${dp} -e 0.002  -o outd${dp}_qs2_beta5_t4 2> outd${dp}_qs2_beta5_t4.log ; } 2> time_outd${dp}_qs2_beta5_t4
# done
#
# threads=8
# for dp in 20;do
	# { time /maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/benchmarking/vcfgl/vcfgl -i sim_vcfgl_2312-OutOfAfrica_3G09-chr22.vcf --threads ${threads} --error-qs 2 --beta-variance 1e-5 -s 1 -d ${dp} -e 0.002  -o outd${dp}_qs2_beta5_t${threads} 2> outd${dp}_qs2_beta5_t${threads}.log ; } 2> time_outd${dp}_qs2_beta5_t${threads}
# done
#
#
# for threads in 1 8;do
	# for dp in 20;do
		# { time /maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/benchmarking/vcfgl/vcfgl -i sim_vcfgl_2312-OutOfAfrica_3G09-chr22.vcf --threads ${threads} -explode 1 -s 1 -d ${dp} -e 0.002  -o outd${dp}_explode1_t${threads} 2> outd${dp}_explode1_t${threads}.log ; } 2> time_outd${dp}_explode1_t${threads}
	# done
# done
#
#
#
# dp=1
# { time /maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/benchmarking/vcfgl/vcfgl -i sim_vcfgl_2312-OutOfAfrica_3G09-chr22.vcf --threads 1 -dogvcf 1 --gvcf-dps 10,20 -s 1 -d ${dp} -e 0.002  -explode 1 -o outd${dp}_explode1_gvcf1 2> outd${dp}_explode1_gvcf1.log ; } 2> time_outd${dp}_explode1_gvcf1
