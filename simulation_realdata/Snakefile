# BCFTOOLS = "/maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/tools/bcftools-1.18/bcftools"
BCFTOOLS="/maps/projects/lundbeck/scratch/pfs488/ibdgl/ibdgl_paper_analyses/tools/bcftools/bcftools"
ANGSD = "/maps/projects/lundbeck/scratch/pfs488/Programs/angsd/angsd"
VCFGL = "/maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/tools/vcfgl/vcfgl"
GET_GT_DISCORDANCE = "/maps/projects/lundbeck/scratch/pfs488/vcfgl/vcfgl_paper_analyses/tools/vcfgl/misc/gtDiscordance"
SAMTOOLS="/maps/projects/lundbeck/scratch/pfs488/Programs/samtools-1.20/samtools"

# DEPTHS = [1, 2, 10, 20, 100]
DEPTHS = [1, 5, 10, 20, 100]
#DEPTHS = [0.1,0.5,1,2, 10, 20, 100]
REPS = range(20)
# GLS = [1, 2]
GLS=1

VCFGL_VERSIONS=["v1", "v2", "v3"]

rule all:
	input:
		"sim_v1/depth/max_avg_persite_depth.txt",
		expand(
			"sim_v1/simulations/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl{vcfglv}.bcf",
			depth=DEPTHS,
			rep=REPS,
			gl=GLS,
			vcfglv=VCFGL_VERSIONS
		),
		expand(
			"sim_v1/depth/HG00096_chr21_subsample_depth{depth}.bam_avg_persite_depth",
			depth=DEPTHS
		),
		expand(
			"sim_v1/depth/HG00096_chr21_subsample_depth{depth}_gl1_bcftools.bcf_avg_persite_depth",
			depth=DEPTHS
		),
		expand("sim_v1/depth/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl{vcfglv}.bcf_avg_persite_depth",
				 depth=DEPTHS,
				 rep=REPS,
				 gl=GLS,
				 vcfglv=VCFGL_VERSIONS
				 ),
		expand(
			 "sim_v1/results/HG00096_chr21_usecommon{usecommonv}_gt_tgt_discordance.tsv",
			 usecommonv=["v1","v2","v3"]
		),
		expand("sim_v1/results/HG00096_chr21_alldepths_rep{rep}_gl{gl}_methods-subsamplev1-vcfgl{vcfglv}_usecommon{usecommonv}_gt_tgt_discordance_summary.csv", 
				 rep=REPS,
				 gl=GLS,
				 vcfglv=VCFGL_VERSIONS,
				 usecommonv=["v1","v2","v3"]
		),
		expand(
			"sim_v1/results/HG00096_chr21_alldepths_usecommon{usecommonv}_gt_tgt_discordance.tsv",
			usecommonv=["v1","v2","v3"]
		),




# --------------------------------------------------------------------------- #
# PREPARE RESOURCES

rule download_reference:
	output:
		"resources/human_g1k_v37_decoy.fasta"
	shell:
		"""
		wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/b37/human_g1k_v37_decoy.fasta.gz
		gunzip -c {input}.gz > {output}
		"""

# for gatk error rate estimation
rule download_reference_dict:
	output:
		"resources/human_g1k_v37_decoy.dict.gz"
	shell:
		"""
		wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/b37/human_g1k_v37_decoy.dict.gz
		"""

rule extract_reference_dict:
	input:
		"resources/human_g1k_v37_decoy.dict.gz"
	output:
		"resources/human_g1k_v37_decoy.dict"
	shell:
		"""
		gunzip -c {input} > {output}
		"""

rule download_bam:
	output:
		"resources/HG00096.chr21.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam"
	shell:
		"""
		wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/high_coverage_alignment/HG00096.wgs.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam
		wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/high_coverage_alignment/HG00096.wgs.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam.bai
		{SAMTOOLS} view -b resources/HG00096.wgs.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam 21 {output}
		"""

rule filter_bam:
	input:
		"resources/HG00096.chr21.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam"
	output:
		"resources/HG00096.chr21.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.filtered.bam"
	shell:
		"""
		{SAMTOOLS} view -b -f 2 -F 1804 {input} > {output}
		"""

# for gatk
rule download_vcfgz:
	output:
		"resources/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
	shell:
		"""
		wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
		"""

rule index_reference:
	input:
		"resources/human_g1k_v37_decoy.fasta"
	output:
		"resources/human_g1k_v37_decoy.fasta.fai"
	shell:
		"""
		{SAMTOOLS} faidx {input}
		"""

rule get_max_avg_persite_depth:
	input:
		"resources/HG00096.chr21.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.filtered.bam",
	output:
		"sim_v1/depth/max_avg_persite_depth.txt"
	shell:
		"""
		{SAMTOOLS} depth {input} | awk '{{sum += $3; n++}} END {{print sum/n}}' > {output}
		"""



# --------------------------------------------------------------------------- #
# SUBSAMPLING FROM REAL DATA TO GET REAL DATA AT VARIOUS DEPTHS

rule subsample_bam:
	input:
		bam="resources/HG00096.chr21.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.filtered.bam",
		max_depth="sim_v1/depth/max_avg_persite_depth.txt"
	output:
		"sim_v1/subsamples/bam/HG00096_chr21_subsample_depth{depth}.bam",
	params:
		seed=42,
		fraction_decimal_part=lambda wildcards,input: "{:.10f}".format(float(wildcards.depth) / float(open(input.max_depth).read().strip())).split(".")[1],
	shell:
		"""
		{SAMTOOLS} view -s {params.seed}.{params.fraction_decimal_part} -b {input.bam} > {output}
		"""

#		 {BCFTOOLS} view -h {input.gt} | grep "^##contig" | grep -Fv '##contig=<ID=21,length=48129895>' > resources/rmcontigs.txt
rule calculate_gl1_bcftools_for_subsamples:
	input:
		ref="resources/human_g1k_v37_decoy.fasta",
		bam="sim_v1/subsamples/bam/HG00096_chr21_subsample_depth{depth}.bam",
		rmcontigs="resources/rmcontigs.txt"
	output:
		"sim_v1/subsamples/bcf/HG00096_chr21_subsample_depth{depth}_gl1_bcftools.bcf",
	log:
		"sim_v1/subsamples/bcf/HG00096_chr21_subsample_depth{depth}_gl1_bcftools.log",
	shell:
		"""
		( 
		{BCFTOOLS} mpileup --max-depth 500 --skip-indels -Q 0 --no-BAQ -Ob --fasta-ref {input.ref} {input.bam}  | {BCFTOOLS} +tag2tag -- --PL-to-GL | {BCFTOOLS} view -i 'REF!=\"N\"' | grep -v -f {input.rmcontigs} | {BCFTOOLS} view -Ob -o {output}
		) 2> {log}
		"""


rule call_genotypes_naive_caller_for_subsamples_bcftools:
	input:
		"sim_v1/subsamples/bcf/HG00096_chr21_subsample_depth{depth}_gl1_bcftools.bcf",
	output:
		"sim_v1/genotype_calling/HG00096_chr21_subsample_depth{depth}_gl1_bcftools_gt.bcf",
	shell:
		"""
		{BCFTOOLS} +tag2tag {input} -- --GL-to-GT --threshold 1 | {BCFTOOLS} view -Ob -o {output}
		"""

# --------------------------------------------------------------------------- #
# CALCULATE AVG PERSITE DEPTH

rule get_depth_from_subsamples_bam:
	input:
		"sim_v1/subsamples/bam/HG00096_chr21_subsample_depth{depth}.bam",
	output:
		"sim_v1/depth/HG00096_chr21_subsample_depth{depth}.bam_avg_persite_depth",
	shell:
		"""
		{SAMTOOLS} depth {input} | awk '{{sum += $3; n++}} END {{print sum/n}}' > {output}
		"""

# only expected difference: we exclude sites with REF=N in the bcf
rule get_depth_from_subsamples_bcf:
	input:
		"sim_v1/subsamples/bcf/HG00096_chr21_subsample_depth{depth}_gl1_bcftools.bcf",
	output:
		"sim_v1/depth/HG00096_chr21_subsample_depth{depth}_gl1_bcftools.bcf_avg_persite_depth",
	shell:
		"""
		{BCFTOOLS} query -f '%CHROM\\t%POS\\t%DP\\n' {input} | awk '{{sum += $3; n++}} END {{print sum/n}}' > {output}
		"""


rule get_depth_from_vcfgl_bcf:
	input:
		"sim_v1/simulations/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl{vcfglv}.bcf",
	output:
		"sim_v1/depth/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl{vcfglv}.bcf_avg_persite_depth"
	shell:
		"""
		{BCFTOOLS} query -f '%CHROM\\t%POS\\t%DP\\n' {input} | awk '{{sum += $3; n++}} END {{print sum/n}}' > {output}
		"""


# --------------------------------------------------------------------------- #

rule collect_error_metrics:
	input:
		bam="sim_v1/subsamples/bam/HG00096_chr21_subsample_depth{depth}.bam",
		vcf="resources/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
		ref="resources/human_g1k_v37_decoy.fasta",
		refdict="resources/human_g1k_v37_decoy.dict",
	params:
		outprefix="sim_v1/error_estimation/HG00096_chr21_subsample_depth{depth}_errorMetrics",
	output:
		"sim_v1/error_estimation/HG00096_chr21_subsample_depth{depth}_errorMetrics.error_by_all",
	log:
		"sim_v1/error_estimation/HG00096_chr21_subsample_depth{depth}_errorMetrics.log",
	shell:
		"""
		module load openjdk/17.0.8 python gatk/4.4.0.0
		gatk CollectSamErrorMetrics -I {input.bam} --OUTPUT {params.outprefix} --VCF {input.vcf} --REFERENCE_SEQUENCE {input.ref} 2> {log}
		"""

rule calculate_avg_error_rate:
	input:
		"sim_v1/error_estimation/HG00096_chr21_subsample_depth{depth}_errorMetrics.error_by_all",
	output:
		"sim_v1/error_estimation/HG00096_chr21_subsample_depth{depth}_errorMetrics.avg_error_rate.txt",
	shell:
		"""
		cat {input} | grep -A2 '^## METRICS' | tail -1 | awk '{{print $1/$4}}' > {output}
		"""



# --------------------------------------------------------------------------- #
# SIMULATION WITH VCFGL 

rule prepare_reference_bcf_depth100_for_vcfgl:
	input:
		"sim_v1/genotype_calling/HG00096_chr21_subsample_depth100_gl1_bcftools_gt.bcf",
	output:
		"sim_v1/simulations/HG00096_chr21_depth100_ref.bcf",
	shell:
		"""
		 {BCFTOOLS} annotate -x "INFO/DP,FORMAT/DP,INFO/I16,INFO/QS,INFO/MQ0F,FORMAT/GL,FORMAT/PL,INFO/SGB,INFO/RPBZ,INFO/MQBZ,INFO/BQBZ,INFO/SCBZ,INFO/INDEL,INFO/IDV,INFO/IMF,INFO/VDB,INFO/MQSBZ" {input} -O b -o {output}
		"""


# VERSION V1: 
# --error-rate: use depth100 error rate for all output depth{depth}
# --depth: use target depth from wildcards in simulation
# VERSION V2: 
# --error-rate: use depth{depth} error rate for each output depth{depth}
# --depth: use depth{depth} observed actual depth for simulation
# VERSION V3:
# --error-rate: use depth100 error rate for all output depth{depth}
# --depth: use depth{depth} observed actual depth for simulation

# VERSION V1: 
# --error-rate: use depth100 error rate for all output depth{depth}
# --depth: use target depth from wildcards in simulation
rule simulate_vcfgl_v1:
	input:
		bcf="sim_v1/simulations/HG00096_chr21_depth100_ref.bcf",
		error_rate_file="sim_v1/error_estimation/HG00096_chr21_subsample_depth100_errorMetrics.avg_error_rate.txt",
	output:
		"sim_v1/simulations/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfglv1.bcf",
	params:
		# read error rate from file and save into error_rate variable
		error_rate_par=lambda wildcards, input: open(input.error_rate_file).read().strip(),
		outprefix="sim_v1/simulations/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfglv1",
	log:
		"sim_v1/simulations/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfglv1.log",
	shell:
		"""
		({VCFGL} -i {input.bcf} -o {params.outprefix} -O b --depth {wildcards.depth} --error-rate {params.error_rate_par} -explode 0 -addInfoAD 1 -addFormatAD 1 -addInfoDP 1 -addFormatDP 1 --rm-invar-sites 0 --rm-empty-sites 1 -doUnobserved 1 --source 1 -addQS 1 -addGL 1 -GL {wildcards.gl} --seed {wildcards.rep}) 2> {log}
		"""


# VERSION V2: 
# --error-rate: use depth{depth} error rate for each output depth{depth}
# --depth: use depth{depth} observed actual depth for simulation
rule simulate_vcfgl_v2:
	input:
		bcf="sim_v1/simulations/HG00096_chr21_depth100_ref.bcf",
		error_rate_file="sim_v1/error_estimation/HG00096_chr21_subsample_depth{depth}_errorMetrics.avg_error_rate.txt",
		avg_persite_depth_file="sim_v1/depth/HG00096_chr21_subsample_depth{depth}_gl1_bcftools.bcf_avg_persite_depth",
	output:
		"sim_v1/simulations/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfglv2.bcf",
	params:
		error_rate_par=lambda wildcards, input: open(input.error_rate_file).read().strip(),
		avg_persite_depth_par=lambda wildcards, input: open(input.avg_persite_depth_file).read().strip(),
		outprefix="sim_v1/simulations/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfglv2",
	log:
		"sim_v1/simulations/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfglv2.log",
	shell:
		"""
		({VCFGL} -i {input.bcf} -o {params.outprefix} -O b --depth {params.avg_persite_depth_par} --error-rate {params.error_rate_par} -explode 0 -addInfoAD 1 -addFormatAD 1 -addInfoDP 1 -addFormatDP 1 --rm-invar-sites 0 --rm-empty-sites 1 -doUnobserved 1 --source 1 -addQS 1 -addGL 1 -GL {wildcards.gl} --seed {wildcards.rep}) 2> {log}
		"""

# VERSION V3:
# --error-rate: use depth100 error rate for all output depth{depth}
# --depth: use depth{depth} observed actual depth for simulation
rule simulate_vcfgl_v3:
	input:
		bcf="sim_v1/simulations/HG00096_chr21_depth100_ref.bcf",
		error_rate_file="sim_v1/error_estimation/HG00096_chr21_subsample_depth100_errorMetrics.avg_error_rate.txt",
		avg_persite_depth_file="sim_v1/depth/HG00096_chr21_subsample_depth{depth}_gl1_bcftools.bcf_avg_persite_depth",
	output:
		"sim_v1/simulations/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfglv3.bcf",
	params:
		# read error rate from file and save into error_rate variable
		error_rate_par=lambda wildcards, input: open(input.error_rate_file).read().strip(),
		avg_persite_depth_par=lambda wildcards, input: open(input.avg_persite_depth_file).read().strip(),
		outprefix="sim_v1/simulations/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfglv3",
	log:
		"sim_v1/simulations/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfglv3.log",
	shell:
		"""
		({VCFGL} -i {input.bcf} -o {params.outprefix} -O b --depth {params.avg_persite_depth_par} --error-rate {params.error_rate_par} -explode 0 -addInfoAD 1 -addFormatAD 1 -addInfoDP 1 -addFormatDP 1 --rm-invar-sites 0 --rm-empty-sites 1 -doUnobserved 1 --source 1 -addQS 1 -addGL 1 -GL {wildcards.gl} --seed {wildcards.rep}) 2> {log}
		"""
	

# --------------------------------------------------------------------------- #
# CALL GENOTYPES FOR SIMULATED DATA


rule call_genotypes_naive_caller_for_vcfgl:
	input:
		"sim_v1/simulations/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl{vcfglv}.bcf",
	output:
		"sim_v1/genotype_calling/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl{vcfglv}_gt.bcf",
	shell:
		"""
		{BCFTOOLS} +tag2tag -Ob -o {output} {input} -- --GL-to-GT --threshold 1
		"""


# --------------------------------------------------------------------------- #
# GET %TGT 


rule get_tgt_subsamples_bcftools:
	input:
		"sim_v1/genotype_calling/HG00096_chr21_subsample_depth{depth}_gl1_bcftools_gt.bcf",
	output:
		"sim_v1/results/HG00096_chr21_depth{depth}_gl1_subsamplev1_gt_tgt.tsv",
	shell:
		"""
		{BCFTOOLS} query -f '%CHROM\\t%POS[\\t%TGT]\\n' {input} |  tr -d '/'  > {output}
		"""

rule get_tgt_vcfgl:
	input:
		"sim_v1/genotype_calling/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl{vcfglv}_gt.bcf",
	output:
		"sim_v1/results/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl{vcfglv}_gt_tgt.tsv",
	shell:
		"""
		{BCFTOOLS} query -f '%CHROM\\t%POS[\\t%TGT]\\n' {input} |  tr -d '/'  > {output}
		"""

rule get_joint_tgt_discordance:
	input:
		truthtgt="sim_v1/results/HG00096_chr21_depth100_gl{gl}_subsamplev1_gt_tgt.tsv",
		calltgt1="sim_v1/results/HG00096_chr21_depth{depth}_gl{gl}_subsamplev1_gt_tgt.tsv",
		calltgt2="sim_v1/results/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl{vcfglv}_gt_tgt.tsv",
		depthfile="sim_v1/depth/HG00096_chr21_subsample_depth{depth}_gl1_bcftools.bcf_avg_persite_depth",
	output:
		persite="sim_v1/results/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_methods-subsamplev1-vcfgl{vcfglv}_usecommon{usecommonv}_gt_tgt_discordance_persite.csv",
		summary="sim_v1/results/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_methods-subsamplev1-vcfgl{vcfglv}_usecommon{usecommonv}_gt_tgt_discordance_summary.csv",
	shell:
		"""
		~/.conda/envs/snakemake/bin/Rscript get_joint_common_tgt_discordance.R {input.truthtgt} {input.calltgt1} {input.calltgt2} "subsamplev1" "vcfgl{wildcards.vcfglv}" {wildcards.rep} {wildcards.depth} {wildcards.gl} FALSE {output.persite} {output.summary} {wildcards.usecommonv} {input.depthfile}
		"""


rule collect_joint_tgt_discordance_summary:
	input:
		expand(
			"sim_v1/results/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_methods-subsamplev1-vcfgl{vcfglv}_usecommon{{usecommonv}}_gt_tgt_discordance_summary.csv",
			depth=DEPTHS,
			rep=REPS,
			gl=1,
			vcfglv=VCFGL_VERSIONS
		),
	output:
		"sim_v1/results/HG00096_chr21_usecommon{usecommonv}_gt_tgt_discordance.tsv",
	shell:
		"""
		cat {input} > {output}
		"""


rule get_alljoint_tgt_discordance:
	input:
		truthtgt="sim_v1/results/HG00096_chr21_depth100_gl{gl}_subsamplev1_gt_tgt.tsv",
		calltgt1=expand("sim_v1/results/HG00096_chr21_depth{depth}_gl{{gl}}_subsamplev1_gt_tgt.tsv",depth=DEPTHS),
		calltgt2=expand("sim_v1/results/HG00096_chr21_depth{depth}_rep{{rep}}_gl{{gl}}_vcfgl{{vcfglv}}_gt_tgt.tsv", depth=DEPTHS),
		depthfile=expand("sim_v1/depth/HG00096_chr21_subsample_depth{depth}_gl1_bcftools.bcf_avg_persite_depth", depth=DEPTHS),
	output:
		persite="sim_v1/results/HG00096_chr21_alldepths_rep{rep}_gl{gl}_methods-subsamplev1-vcfgl{vcfglv}_usecommon{usecommonv}_gt_tgt_discordance_persite.csv",
		summary="sim_v1/results/HG00096_chr21_alldepths_rep{rep}_gl{gl}_methods-subsamplev1-vcfgl{vcfglv}_usecommon{usecommonv}_gt_tgt_discordance_summary.csv",
	shell:
		"""
		~/.conda/envs/snakemake/bin/Rscript get_alldepths_joint_common_tgt_discordance.R {input.truthtgt} {input.calltgt1} {input.calltgt2} {wildcards.rep} {wildcards.gl} FALSE {output.persite} {output.summary} {wildcards.usecommonv} {input.depthfile}
		"""
	


rule collect_alljoint_tgt_discordance_summary:
	input:
		expand(
			"sim_v1/results/HG00096_chr21_alldepths_rep{rep}_gl{gl}_methods-subsamplev1-vcfgl{vcfglv}_usecommon{{usecommonv}}_gt_tgt_discordance_summary.csv",
			rep=REPS,
			gl=1,
			vcfglv=VCFGL_VERSIONS
		),
	output:
		"sim_v1/results/HG00096_chr21_alldepths_usecommon{usecommonv}_gt_tgt_discordance.tsv",
	shell:
		"""
		cat {input} > {output}
		"""
		



# --------------------------------------------------------------------------- #

# colnames<-c("Sample","nSitesTotal","nSitesRetained","nSitesCompared","nSitesCallMis","nDiscordantSites","nConcordantSites","nSitesinTrueNotCall","MissingnessRate","DiscordanceRate","ConcordanceRate","T_Hom","T_Het","F_HomToHom","F_HomToHet","F_HetToHom","F_HetToHet","T_Hom_rate","T_Het_rate","F_HomToHom_rate","F_HomToHet_rate","F_HetToHom_rate","F_HetToHet_rate","Rep","Depth","GL","Method")
# rule collect_all_tidy_genotype_discordance:
#	 input:
#		 expand(
#			 "sim_v1/results/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl_gt_discordance.tsv.tidy",
#			 depth=DEPTHS,
#			 rep=REPS,
#			 gl=1,
#		 ),
#		 expand(
#			 "sim_v1/results/HG00096_chr21_subsample_depth{depth}_gl1_bcftools_filtered_gt_discordance.tsv.tidy",
#			 depth=DEPTHS_NO100X,
#			 rep=REPS,
#		 ),
#	 output:
#		 "sim_v1/results/HG00096_chr21_gt_discordance.tsv",
#	 shell:
#		 """
#		 cat {input} > {output}
#		 """


# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# rule get_genotype_discordance_subsamples_angsd:
#	 input:
#		 # truthgt="sim_v1/data/HG00096_chr21_gt.bcf",
#		 callgt="sim_v1/genotype_calling/HG00096_chr21_subsample_depth{depth}_gl{gl}_angsd_filtered_gt.bcf",
#	 output:
#		 "sim_v1/results/HG00096_chr21_subsample_depth{depth}_gl{gl}_angsd_filtered_gt_discordance.tsv",
#	 log:
#		 "sim_v1/results/HG00096_chr21_subsample_depth{depth}_gl{gl}_angsd_filtered_gt_discordance.log",
#	 shell:
#		 """
#		 ({GET_GT_DISCORDANCE} -t {input.truthgt} -i {input.callgt} -o {output}) 2> {log}
#		 """


# rule tidy_genotype_discordance_subsamples_angsd:
#	 input:
#		 "sim_v1/results/HG00096_chr21_subsample_depth{depth}_gl{gl}_angsd_filtered_gt_discordance.tsv",
#	 output:
#		 "sim_v1/results/HG00096_chr21_subsample_depth{depth}_gl{gl}_angsd_filtered_gt_discordance.tsv.tidy",
#	 params:
#		 method="angsd",
#	 shell:
#		 """
#		 awk -v REP={wildcards.rep} -v DEPTH={wildcards.depth} -v GL={wildcards.gl} -v METHOD={params.method} 'BEGIN{{FS="\t";OFS="\t"}}{{print $0,REP,DEPTH,GL,METHOD}}' {input} > {output}
#		 """


# rule collect_tidy_genotype_discordance_vcfgl_angsd:
#	 input:
#		 expand(
#			 "sim_v1/results/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl_gt_discordance.tsv.tidy",
#			 depth=DEPTHS,
#			 rep=REPS,
#			 gl=GLS,
#		 ),
#		 expand(
#			 "sim_v1/results/HG00096_chr21_subsample_depth{depth}_gl{gl}_angsd_filtered_gt_discordance.tsv.tidy",
#			 depth=DEPTHS_NO100X,
#			 rep=REPS,
#			 gl=GLS,
#		 ),
#	 output:
#		 "sim_v1/results/HG00096_chr21_gt_discordance.tsv",
#	 shell:
#		 """
#		 cat {input} > {output}
#		 """


# rule get_genotype_discordance_doGQ_vcfgl:
#	 input:
#		 truthgt="sim_v1/data/HG00096_chr21_gt.bcf",
#		 callgt="sim_v1/genotype_calling/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl_gt.bcf",
#	 output:
#		 "sim_v1/results/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl_gt_discordance_doGQ.tsv"
#	 log:
#		 "sim_v1/results/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl_gt_discordance_doGQ.log"
#	 shell:
#		 """
#		 ({GET_GT_DISCORDANCE} -t {input.truthgt} -i {input.callgt} -o {output} -doGQ 8 ) 2> {log}
#		 """

# rule get_genotype_discordance_doGQ_subsamples:
#	 input:
#		 truthgt="sim_v1/data/HG00096_chr21_gt.bcf",
#		 callgt="sim_v1/genotype_calling/HG00096_chr21_subsample_depth{depth}_gl{gl}_angsd_filtered_gt.bcf",
#	 output:
#		 "sim_v1/results/HG00096_chr21_subsample_depth{depth}_gl{gl}_angsd_filtered_gt_discordance_doGQ.tsv"
#	 log:
#		 "sim_v1/results/HG00096_chr21_subsample_depth{depth}_gl{gl}_angsd_filtered_gt_discordance_doGQ.log"
#	 shell:
#		 """
#		 ({GET_GT_DISCORDANCE} -t {input.truthgt} -i {input.callgt} -o {output} -doGQ 8 ) 2> {log}
#		 """

# rule tidy_genotype_discordance_doGQ_vcfgl:
#	 input:
#		 "sim_v1/results/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_{method}_gt_discordance_doGQ.tsv"
#	 output:
#		 "sim_v1/results/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_{method}_gt_discordance_doGQ.tsv.tidy"
#	 shell:
#		 """
#		 awk -v REP={wildcards.rep} -v DEPTH={wildcards.depth} -v GL={wildcards.gl} -v METHOD={wildcards.method} 'BEGIN{{FS="\t";OFS="\t"}}{{print $0,REP,DEPTH,GL,METHOD}}' {input} > {output}
#		 """

# rule tidy_genotype_discordance_doGQ_subsamples:
#	 input:
#		 "sim_v1/results/HG00096_chr21_subsample_depth{depth}_gl{gl}_{method}_gt_discordance_doGQ.tsv"
#	 output:
#		 "sim_v1/results/HG00096_chr21_subsample_depth{depth}_gl{gl}_{method}_gt_discordance_doGQ.tsv.tidy"
#	 shell:
#		 """
#		 awk -v REP={wildcards.rep} -v DEPTH={wildcards.depth} -v GL={wildcards.gl} -v METHOD={wildcards.method} 'BEGIN{{FS="\t";OFS="\t"}}{{print $0,REP,DEPTH,GL,METHOD}}' {input} > {output}
#		 """

# rule collect_tidy_genotype_discordance_doGQ_vcfgl:
#	 input:
#		 expand(
#			 "sim_v1/results/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl_gt_discordance_doGQ.tsv.tidy",
#			 depth=DEPTHS,
#			 rep=REPS,
#			 gl=GLS,
#		 ),
#	 output:
#		 "sim_v1/results/HG00096_chr21_gt_discordance_doGQ.tsv"
#	 shell:
#		 """
#		 cat {input} > {output}
#		 """

# rule collect_tidy_genotype_discordance_doGQ_subsamples:
#	 input:
#		 expand(
#			 "sim_v1/results/HG00096_chr21_subsample_depth{depth}_gl{gl}_angsd_filtered_gt_discordance_doGQ.tsv.tidy",
#			 depth=DEPTHS_NO100X,
#			 rep=REPS,
#			 gl=GLS,
#		 ),
#	 output:
#		 "sim_v1/results/HG00096_chr21_subsample_gt_discordance_doGQ.tsv"
#	 shell:
#		 """
#		 cat {input} > {output}
#		 """


# rule call_genotypes_naive_caller_fulldata_gls_angsd:
#	 input:
#		 "sim_v1/realdata/bcf/HG00096_chr21_fulldata_gl{gl}_angsd_filtered.bcf",
#	 output:
#		 "sim_v1/genotype_calling/HG00096_chr21_fulldata_gl{gl}_angsd_filtered_gt.bcf",
#	 shell:
#		 """
#		 {BCFTOOLS} +tag2tag -Ob -o {output} {input} -- --GL-to-GT --threshold 1
#		 """


# # make sure genotype discordance for fulldata for both gl methods
# # header=c("Sample","nSitesTotal","nSitesRetained","nSitesCompared","nSitesCallMis","nDiscordantSites","nConcordantSites","nSitesinTrueNotCall","MissingnessRate","DiscordanceRate","ConcordanceRate","T_Hom","T_Het","F_HomToHom","F_HomToHet","F_HetToHom","F_HetToHet","T_Hom_rate","T_Het_rate","F_HomToHom_rate","F_HomToHet_rate","F_HetToHom_rate","F_HetToHet_rate")
# rule assert_genotype_discordance_for_fulldata:
#	 input:
#		 truthgt="sim_v1/data/HG00096_chr21_gt.bcf",
#		 callgt_gl1="sim_v1/genotype_calling/HG00096_chr21_fulldata_gl1_angsd_filtered_gt.bcf",
#		 callgt_gl2="sim_v1/genotype_calling/HG00096_chr21_fulldata_gl2_angsd_filtered_gt.bcf",
#	 output:
#		 callgt_gl1_sites="sim_v1/genotype_calling/HG00096_chr21_fulldata_gl1_angsd_filtered_gt_sites.bed",
#		 callgt_gl2_sites="sim_v1/genotype_calling/HG00096_chr21_fulldata_gl2_angsd_filtered_gt_sites.bed",
#		 truthgt_sites="sim_v1/genotype_calling/HG00096_chr21_fulldata_truthgt_sites.bed",
#		 common_sites="sim_v1/genotype_calling/HG00096_chr21_fulldata_common_sites.bed",
#		 callgt_gl1_common_bcf="sim_v1/genotype_calling/HG00096_chr21_fulldata_gl1_angsd_filtered_gt_common.bcf",
#		 callgt_gl2_common_bcf="sim_v1/genotype_calling/HG00096_chr21_fulldata_gl2_angsd_filtered_gt_common.bcf",
#		 truthgt_common_bcf="sim_v1/data/HG00096_chr21_gt_common.bcf",
#		 callgt_gl1_discordance="sim_v1/results/HG00096_chr21_fulldata_gl1_angsd_filtered_gt_discordance.tsv",
#		 callgt_gl2_discordance="sim_v1/results/HG00096_chr21_fulldata_gl2_angsd_filtered_gt_discordance.tsv",
#	 log:
#		 "sim_v1/results/HG00096_chr21_fulldata_gls_angsd_filtered_gt_discordance.log",
#	 shell:
#		 """
#		 {BCFTOOLS} query -f '%CHROM\\t%POS0\\t%END\\n' {input.callgt_gl1} > {output.callgt_gl1_sites}
#		 {BCFTOOLS} query -f '%CHROM\\t%POS0\\t%END\\n' {input.callgt_gl2} > {output.callgt_gl2_sites}
#		 {BCFTOOLS} query -f '%CHROM\\t%POS0\\t%END\\n' {input.truthgt} > {output.truthgt_sites}
#		 comm --nocheck-order -12 {output.callgt_gl1_sites} {output.callgt_gl2_sites} | comm --nocheck-order -12 - {output.truthgt_sites} > {output.common_sites}
#		 {BCFTOOLS} view -T {output.common_sites} {input.callgt_gl1} -Ob -o {output.callgt_gl1_common_bcf}
#		 {BCFTOOLS} view -T {output.common_sites} {input.callgt_gl2} -Ob -o {output.callgt_gl2_common_bcf}
#		 {BCFTOOLS} view -T {output.common_sites} {input.truthgt} -Ob -o {output.truthgt_common_bcf}
#		 ({GET_GT_DISCORDANCE} -t {output.truthgt_common_bcf} -i {output.callgt_gl1_common_bcf} -o {output.callgt_gl1_discordance}) 2> {log}
#		 ({GET_GT_DISCORDANCE} -t {output.truthgt_common_bcf} -i {output.callgt_gl2_common_bcf} -o {output.callgt_gl2_discordance}) 2>> {log}
#		 """


#   A C G T
# A
# C x
# G x x
# T x x x
# fix the order if needed
# e.g.
# CA -> AC
# AC -> AC (no change)

# for CA, GA, GC, TA, TC, TG
# sed 's/CA/AC/g; s/GA/AG/g; s/GC/CG/g; s/TA/AT/g; s/TC/CT/g; s/TG/GT/g'

# rule get_genotype_discordance_subsamples_gl1_bcftools:
#	 input:
#		 truthgt="sim_v1/genotype_calling/HG00096_chr21_fulldata_gl1_bcftools_filtered_gt.bcf",
#		 callgt="sim_v1/genotype_calling/HG00096_chr21_subsample_depth{depth}_gl1_bcftools_filtered_gt.bcf",
#	 output:
#		 "sim_v1/results/HG00096_chr21_subsample_depth{depth}_gl1_bcftools_filtered_gt_discordance.tsv",
#	 log:
#		 "sim_v1/results/HG00096_chr21_subsample_depth{depth}_gl1_bcftools_filtered_gt_discordance.log",
#	 shell:
#		 """
#		 ({GET_GT_DISCORDANCE} -t {input.truthgt} -i {input.callgt} -o {output}) 2> {log}
#		 """



# rule tidy_genotype_discordance_subsamples_gl1_bcftools:
#	 input:
#		 "sim_v1/results/HG00096_chr21_subsample_depth{depth}_gl1_bcftools_filtered_gt_discordance.tsv",
#	 output:
#		 "sim_v1/results/HG00096_chr21_subsample_depth{depth}_gl1_bcftools_filtered_gt_discordance.tsv.tidy",
#	 params:
#		 method="bcftools",
#	 shell:
#		 """
#		 awk -v REP={wildcards.rep} -v DEPTH={wildcards.depth} -v GL=1 -v METHOD={params.method} 'BEGIN{{FS="\t";OFS="\t"}}{{print $0,REP,DEPTH,GL,METHOD}}' {input} > {output}
#		 """


# rule get_genotype_discordance_vcfgl:
#	 input:
#		 truthgt="sim_v1/genotype_calling/HG00096_chr21_fulldata_gl1_bcftools_filtered_gt.bcf",
#		 callgt="sim_v1/genotype_calling/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl_gt.bcf",
#	 output:
#		 "sim_v1/results/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl_gt_discordance.tsv",
#	 log:
#		 "sim_v1/results/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl_gt_discordance.log",
#	 shell:
#		 """
#		 ({GET_GT_DISCORDANCE} -t {input.truthgt} -i {input.callgt} -o {output}) 2> {log}
#		 """


# rule tidy_genotype_discordance_vcfgl:
#	 input:
#		 "sim_v1/results/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl_gt_discordance.tsv",
#	 output:
#		 "sim_v1/results/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl_gt_discordance.tsv.tidy",
#	 params:
#		 method="vcfgl",
#	 shell:
#		 """
#		 awk -v REP={wildcards.rep} -v DEPTH={wildcards.depth} -v GL={wildcards.gl} -v METHOD={params.method} 'BEGIN{{FS="\t";OFS="\t"}}{{print $0,REP,DEPTH,GL,METHOD}}' {input} > {output}
#		 """


# rule calculate_gls_angsd_fulldata_gls:
#	 input:
		# bam="resources/HG00096.chr21.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.filtered.bam",
#		 sites="sim_v1/includesites.bed",
#	 output:
#		 bcf="sim_v1/realdata/bcf/HG00096_chr21_fulldata_gl{gl}_angsd.bcf",
#		 filtered_bcf="sim_v1/realdata/bcf/HG00096_chr21_fulldata_gl{gl}_angsd_filtered.bcf",
#	 params:
#		 outprefix="sim_v1/realdata/bcf/HG00096_chr21_fulldata_gl{gl}_angsd",
#	 log:
#		 "sim_v1/realdata/bcf/HG00096_chr21_fulldata_gl{gl}_angsd.log",
#	 shell:
#		 """
#		 {ANGSD} -i {input.bam} -gl {wildcards.gl} -domajorminor 1 -domaf 1 -docounts 1 -dogeno 1 -dopost 1 -doBcf 1 -out {params.outprefix} 2> {log}
#		 {BCFTOOLS} view -T {input.sites} {output.bcf} | grep -v -f resources/rmcontigs.txt | {BCFTOOLS} view -Ob -o {output.filtered_bcf} 
#		 """


# rule calculate_gls_angsd_for_subsamples:
#	 input:
#		 ref="resources/human_g1k_v37_decoy.fasta",
#		 bam="sim_v1/subsamples/bam/HG00096_chr21_subsample_depth{depth}.bam",
#		 sites="sim_v1/includesites.bed",
#	 output:
#		 bcf="sim_v1/subsamples/bcf/HG00096_chr21_subsample_depth{depth}_gl{gl}_angsd.bcf",
#		 filtered_bcf="sim_v1/subsamples/bcf/HG00096_chr21_subsample_depth{depth}_gl{gl}_angsd_filtered.bcf",
#	 params:
#		 outprefix="sim_v1/subsamples/bcf/HG00096_chr21_subsample_depth{depth}_gl{gl}_angsd",
#	 log:
#		 "sim_v1/subsamples/bcf/HG00096_chr21_subsample_depth{depth}_gl{gl}_angsd.log",
#	 shell:
#		 """
#		 {ANGSD} -i {input.bam} -gl {wildcards.gl} -domajorminor 1 -domaf 1 -docounts 1 -dopost 1 -doBcf 1 -out {params.outprefix} 2> {log}
#		 {BCFTOOLS} view -T {input.sites} {output.bcf} | grep -v -f resources/rmcontigs.txt | {BCFTOOLS} view -Ob -o {output.filtered_bcf} 
#		 """


# rule call_genotypes_naive_caller_for_subsamples_angsd:
#	 input:
#		 "sim_v1/subsamples/bcf/HG00096_chr21_subsample_depth{depth}_gl{gl}_angsd_filtered.bcf",
#	 output:
#		 "sim_v1/genotype_calling/HG00096_chr21_subsample_depth{depth}_gl{gl}_angsd_filtered_gt.bcf",
#	 shell:
#		 """
#		 {BCFTOOLS} +tag2tag -Ob -o {output} {input} -- --GL-to-GT --threshold 1
#		 """



# rule extract_sites_to_bed_fulldata:
#	 input:
#		 "sim_v1/realdata/bcf/HG00096_chr21_fulldata_gl1_bcftools.bcf",
#	 output:
#		 "sim_v1/includesites.bed",
#	 shell:
#		 """
#		 {BCFTOOLS} query -f '%CHROM\\t%POS0\\t%END\\n' {input} > {output}
#		 """

# rule get_sitesfiles_from_bed:
#	 input:
#		 "sim_v1/includesites.bed",
#	 output:
#		 "sim_v1/includesites.txt",
#	 shell:
#		 """
#		 awk '{{print $1"\t"$2"\t"}}' {input} > {output}
#		 """

# --------------------------------------------------------------------------- #
# FULL DATA PROCESSING
	
# rule calculate_gl1_bcftools_for_fulldata:
#	 input:
#		 ref="resources/human_g1k_v37_decoy.fasta",
#		 bam="resources/HG00096.chr21.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.filtered.bam",
#	 output:
#		 "sim_v1/realdata/bcf/HG00096_chr21_fulldata_gl1_bcftools.bcf",
#	 log:
#		 "sim_v1/realdata/bcf/HG00096_chr21_fulldata_gl1_bcftools.log",
#	 shell:
#		 """
#		 {BCFTOOLS} mpileup --max-depth 500 --skip-indels -Q 0 --no-BAQ -Ob --fasta-ref {input.ref} {input.bam} > {output}
#		 """

# rule calculate_bcf_avg_persite_depth_fulldata:
#	 input:
#		 "sim_v1/realdata/bcf/HG00096_chr21_fulldata_gl1_bcftools.bcf",
#	 output:
#		 "sim_v1/depth/HG00096_chr21_fulldata_gl1_bcftools.bcf_avg_persite_depth",
#	 shell:
#		 """
#		 {BCFTOOLS} query -f '%DP\\n' {input} | awk '{{sum += $1; n++}} END {{print sum/n}}' > {output}
#		 """

# rule call_genotypes_naive_caller_fulldata_gl1_bcftools:
#	 input:
#		 "sim_v1/realdata/bcf/HG00096_chr21_fulldata_gl1_bcftools.bcf",
#	 output:
#		 gl="sim_v1/realdata/bcf/HG00096_chr21_fulldata_gl1_bcftools_gl.bcf",
#		 gt="sim_v1/genotype_calling/HG00096_chr21_fulldata_gl1_bcftools_gt.bcf",
#	 shell:
#		 """
#		 {BCFTOOLS} +tag2tag -Ob -o {output.gl} {input} -- --PL-to-GL 
#		 {BCFTOOLS} +tag2tag -Ob -o {output.gt} {output.gl} -- --GL-to-GT --threshold 1
#		 """

# rule collect_error_metrics:
# 	input:
# 		bam="resources/HG00096.chr21.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.filtered.bam",
# 		vcf="resources/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
# 		ref="resources/human_g1k_v37_decoy.fasta",
# 		refdict="resources/human_g1k_v37_decoy.dict",
# 	params:
# 		outprefix="sim_v1/error_estimation/HG00096_chr21_errorMetrics",
# 	output:
# 		"sim_v1/error_estimation/HG00096_chr21_errorMetrics.error_by_all",
# 	log:
# 		"sim_v1/error_estimation/HG00096_chr21_errorMetrics.log",
# 	shell:
# 		"""
# 		module load openjdk/17.0.8 python gatk/4.4.0.0
# 		gatk CollectSamErrorMetrics -I {input.bam} --OUTPUT {params.outprefix} --VCF {input.vcf} --REFERENCE_SEQUENCE {input.ref} 2> {log}
# 		"""

# rule calculate_avg_error_rate:
# 	input:
# 		"sim_v1/error_estimation/HG00096_chr21_errorMetrics.error_by_all",
# 	output:
# 		"sim_v1/error_estimation/avg_error_rate.txt",
# 	shell:
# 		"""
# 		cat {input} | grep -A2 '^## METRICS' | tail -1 | awk '{{print $1/$4}}' > {output}
# 		"""

# rule v2_get_subsample_bcf_bed:
#	 input:
#		 "sim_v1/subsamples/bcf/HG00096_chr21_subsample_depth{depth}_gl1_bcftools.bcf",
#	 output:
#		 "sim_v1/subsamples/bcf/HG00096_chr21_subsample_depth{depth}_gl1_bcftools.bed",
#	 shell:
#		 """
#		 {BCFTOOLS} query -f '%CHROM\\t%POS0\\t%END\\n' {input} > {output}
#		 """

# rule v2_collect_subsample_error_metrics:
# 	input:
# 		inbcf="sim_v1/subsamples/bcf/HG00096_chr21_subsample_depth{depth}_gl1_bcftools.bcf",
# 		bam="resources/HG00096.chr21.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.filtered.bam",
# 		vcf="resources/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
# 		ref="resources/human_g1k_v37_decoy.fasta",
# 		refdict="resources/human_g1k_v37_decoy.dict",
# 	output:
# 		error_by_all="sim_v1/error_estimation/HG00096_chr21_subsample_depth{depth}_gl1_bcftools_errorMetrics.error_by_all",
# 		bam="sim_v1/subsamples/bam/HG00096_chr21_subsample_depth{depth}_gl1_bcftools.bam",
# 	params:
# 		outprefix="sim_v1/error_estimation/HG00096_chr21_subsample_depth{depth}_gl1_bcftools_filtered_errorMetrics",
# 	log:
# 		"sim_v1/error_estimation/HG00096_chr21_subsample_depth{depth}_gl1_bcftools_filtered_errorMetrics.log",
# 	shell:
# 		"""
# 		(
# 			{BCFTOOLS} query -f '%CHROM\\t%POS0\\t%END\\n' {input.inbcf} > {output.bed};
# 			# filter bam by bed
# 			{SAMTOOLS} view -b -L {output.bed} {input.bam} > {output.bam};
# 			module load openjdk/17.0.8 python gatk/4.4.0.0;
# 			gatk CollectSamErrorMetrics -I {output.bam} --OUTPUT {params.outprefix} --VCF {input.vcf} --REFERENCE_SEQUENCE {input.ref}
# 		)2> {log}
# 		"""

# rule v2_get_subsample_avg_error_rate:
# 	input:
# 		"sim_v1/error_estimation/HG00096_chr21_subsample_depth{depth}_gl1_bcftools_filtered_errorMetrics.error_by_all",
# 	output:
# 		"sim_v1/error_estimation/HG00096_chr21_subsample_depth{depth}_gl1_bcftools_filtered_avg_error_rate.txt",
# 	shell:
# 		"""
# 		cat {input} | grep -A2 '^## METRICS' | tail -1 | awk '{{print $1/$4}}' > {output}
# 		"""

# rule v2_get_depth_from_subsamples_bam:
# 	input:
# 		bam="sim_v1/subsamples/bam/HG00096_chr21_subsample_depth{depth}.bam",
# 		sites="sim_v1/includesites.txt",
# 	output:
# 		"sim_v1/subsamples/bam/HG00096_chr21_subsample_depth{depth}.depth",
# 	shell:
# 		"""
# 		{SAMTOOLS} depth -a {input.bam} | grep -f {input.sites} | awk '{{sum += $3; n++}} END {{print sum/n}}' > {output}
# 		"""

# rule v2_get_subsample_actual_read_depth_from_subsamples_bcf:
# 	input:
# 		"sim_v1/subsamples/bcf/HG00096_chr21_subsample_depth{depth}_gl1_bcftools.bcf",
# 	output:
# 		"sim_v1/subsamples/bcf/HG00096_chr21_subsample_depth{depth}_gl1_bcftools_filtered.actual_avg_depth",
# 	params:
# 		##contig=<ID=21,length=48129895>	 
# 		chromsize=48129895,
# 	shell:
# 		"""
# 		{BCFTOOLS} query -f '%CHROM\\t%POS\\t%DP\\n' {input} | awk -v TOTN={params.chromsize} '{{sum += $3; n++}} END {{print sum/TOTN}}' > {output}
# 		"""

# rule v2_get_subsample_effective_read_depth_from_subsamples_bcf:
# 	input:
# 		"sim_v1/subsamples/bcf/HG00096_chr21_subsample_depth{depth}_gl1_bcftools.bcf",
# 	output:
# 		"sim_v1/subsamples/bcf/HG00096_chr21_subsample_depth{depth}_gl1_bcftools_filtered.effective_avg_depth",
# 	shell:
# 		"""
# 		{BCFTOOLS} query -f '%CHROM\\t%POS\\t%DP\\n' {input} | awk '{{sum += $3; n++}} END {{print sum/n}}' > {output}
# 		"""

# rule v2_subset_truth_bcf_for_vcfgl:
# 	input:
# 		bed="sim_v1/subsamples/bcf/HG00096_chr21_subsample_depth{depth}_gl1_bcftools_filtered.bed",
# 		fullgt_bcf="sim_v1/genotype_calling/HG00096_chr21_fulldata_gl1_bcftools_filtered_gt.bcf",
# 	output:
# 		"sim_v1/subsamples/bcf/HG00096_chr21_fulldata_truth_subset-of-depth{depth}_gl1_bcftools.bcf",
# 	shell:
# 		"""
# 		{BCFTOOLS} view -T {input.bed} {input.fullgt_bcf} -Ob -o {output}
# 		"""
# rule v2_sim_simulate_vcfgl:
# 	input:
# 		bcf="sim_v1/subsamples/bcf/HG00096_chr21_fulldata_truth_subset-of-depth{depth}_gl1_bcftools.bcf",
# 		error_rate="sim_v1/error_estimation/HG00096_chr21_subsample_depth{depth}_gl1_bcftools_filtered_avg_error_rate.txt",
# 		effective_depth="sim_v1/subsamples/bcf/HG00096_chr21_subsample_depth{depth}_gl1_bcftools_filtered.effective_avg_depth",
# 	output:
# 		"sim_v1/simulations/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl.bcf",
# 	params:
# 		effective_depth=lambda wildcards,input: open(input.effective_depth).read().strip(),
# 		error_rate=lambda wildcards,input: open(input.error_rate).read().strip(),
# 		outprefix="sim_v1/simulations/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl",
# 	log:
# 		"sim_v1/simulations/HG00096_chr21_depth{depth}_rep{rep}_gl{gl}_vcfgl.log",
# 	shell:
# 		"""
# 		({VCFGL} -i {input.bcf} -o {params.outprefix} -O b --depth {params.effective_depth} --error-rate {params.error_rate} -addQS 1 -explode 0 -addInfoAD 1 -addFormatAD 1 -addFormatDP 1 -addPL 1 -GL {wildcards.gl} --rm-invar-sites 0 --rm-empty-sites 0 -doUnobserved 1 --source 1 --seed {wildcards.rep}) 2> {log}
# 		"""
