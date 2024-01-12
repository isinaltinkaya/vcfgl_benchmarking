#
# for dp in 0.1 0.5 1 2 10 20;do
	# { time /maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/tools/vcfgl/vcfgl -i sim_vcfgl_2312-OutOfAfrica_3G09-chr22.vcf -s 1 --error-qs 0 -d ${dp} -e 0.002  -o outd${dp} 2> outd${dp}.log ; } 2> time_outd${dp}
# done
for i in time*;do printf "$(echo $i|sed 's/time_outd//g'),$(cat $i|grep real|cut -f2)\n" ;done
