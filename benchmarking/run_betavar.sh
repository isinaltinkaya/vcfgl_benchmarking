
benchmark(){
	/maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/tools/time-1.9/time -v "$@"
}

exec=/maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/tools/vcfgl/vcfgl

bmdir=benchmark_results
outdir=benchmark_outputs

mkdir -pv ${bmdir}
mkdir -pv ${outdir}


dp=10
for rep in 1 2 3 4 5;do
	outfile=out_d${dp}_explode1_gvcf1_rep${rep}_qs0
	{ benchmark ${exec} -i sim_vcfgl_2312-OutOfAfrica_3G09-chr22.vcf -dogvcf 1 --gvcf-dps 10,20 -s 1 -d ${dp} -e 0.002  -explode 1 -o ${outdir}/${outfile} 2> ${outdir}/${outfile}.log ; } 2> ${bmdir}/${outfile} &
done

for threads in 1 8;do
	for dp in 1;do
		for rep in 1 2 3 4 5;do
			outfile=out_d${dp}_explode1_t${threads}_rep${rep}_qs0
			{ benchmark ${exec} -i sim_vcfgl_2312-OutOfAfrica_3G09-chr22.vcf --threads ${threads} -explode 1 -s 1 -d ${dp} -e 0.002  -o ${outdir}/${outfile} 2> ${outdir}/${outfile}.log; } 2> ${bmdir}/${outfile} 
		done
	done
done


for threads in 1 8;do
	for dp in 1;do
		for rep in 1 2 3 4 5;do
			outfile=out_d${dp}_explode1_t${threads}_rep${rep}_qs25
			{ benchmark ${exec} -i sim_vcfgl_2312-OutOfAfrica_3G09-chr22.vcf --threads ${threads} --error-qs 2 --beta-variance 1e-5 -explode 1 -s 1 -d ${dp} -e 0.002  -o ${outdir}/${outfile} 2> ${outdir}/${outfile}.log; } 2> ${bmdir}/${outfile} &
		done
	done
done

wait



# for dp in 0.1 0.5 1 2 10 20;do
# 	{ benchmark ${exec} -i sim_vcfgl_2312-OutOfAfrica_3G09-chr22.vcf -s 1 --error-qs 0 -d ${dp} -e 0.002  -o outd${dp} 2> outd${dp}.log ; } 2> benchmark_results/time_outd${dp}
# done
