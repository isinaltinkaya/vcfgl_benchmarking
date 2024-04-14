###############################################################################
# Step 3: Genotype calling
#
# isinaltinkaya
###############################################################################
#
# Naming convention:
# - log files:
# output: sim/{simid}/REST_OF_PATH.EXTENSION
# params.prefix: sim/{simid}/REST_OF_PATH
# log: sim/{simid}/logs/REST_OF_PATH.EXTENSION


simulation_id = config["simulation_id"]


MODEL = config["model"]
ERROR_RATE = config["vcfgl_error_rate"]
DEPTH = config["depth"]
CONTIG = config["contig"]
REP = [*range(config["n_reps"])]


POPS = config["pops"]


ANGSD = config["tools"]["angsd"]
BCFTOOLS = config["tools"]["bcftools"]

QSBETAS = ["0_0", "2_5", "2_6", "2_7"]
VCFGLPARS = ["gl1", "gl2", "precise1_gl2", "platform1_gl1", "platform1_gl2"]

BCFTOOLS_MCALL_PRIOR = list(config["bcftools_mcall_prior"].keys())

###############################################################################
# BEGIN RULES


rule all:
	input:
		expand(
		"sim/{simid}/model_{model_id}/contig_{contig}/gcMethod_mcall_all_prior{mcallPriorID}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_mcall_all_prior{mcallPriorID}.bcf",
			simid=simulation_id,
			model_id=MODEL,
			contig=CONTIG,
			rep=REP,
			depth=DEPTH,
			error_rate=ERROR_RATE,
			vcfgl_qsbeta=["00", "25", "26", "27"],
			vcfgl_gl=["1", "2", "2p"],
			vcfgl_platform=[0, 1],
			mcallPriorID=BCFTOOLS_MCALL_PRIOR,
		),
		expand(
		"sim/{simid}/model_{model_id}/contig_{contig}/gcMethod_mcall_perpop_prior{mcallPriorID}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_mcall_perpop_prior{mcallPriorID}.bcf",
			simid=simulation_id,
			model_id=MODEL,
			contig=CONTIG,
			rep=REP,
			depth=DEPTH,
			error_rate=ERROR_RATE,
			vcfgl_qsbeta=["00", "25", "26", "27"],
			vcfgl_gl=["1", "2", "2p"],
			vcfgl_platform=[0, 1],
			mcallPriorID=BCFTOOLS_MCALL_PRIOR,
		),
		expand(
		"sim/{simid}/model_{model_id}/contig_{contig}/gcMethod_hardcall/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_hardcall.bcf",
			simid=simulation_id,
			model_id=MODEL,
			contig=CONTIG,
			rep=REP,
			depth=DEPTH,
			error_rate=ERROR_RATE,
			vcfgl_qsbeta=["00", "25", "26", "27"],
			vcfgl_gl=["1", "2", "2p"],
			vcfgl_platform=[0, 1],
			mcallPriorID=BCFTOOLS_MCALL_PRIOR,
		),



def get_mcall_prior_val_from_id(wildcards):
	try:
		return config["bcftools_mcall_prior"][wildcards.mcallPriorID]
	except KeyError:
		raise ValueError("mcallPriorID not recognized")


rule bcftools_call_genotypes:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}.bcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/gcMethod_mcall_all_prior{mcallPriorID}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_mcall_all_prior{mcallPriorID}.bcf"
	log:
		"sim/{simid}/model_{model_id}/logs/contig_{contig}/gcMethod_mcall_all_prior{mcallPriorID}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_mcall_all_prior{mcallPriorID}.bcf"
	params:
		mcallPriorVal=get_mcall_prior_val_from_id,
	shell:
		"""
		({BCFTOOLS} call \
			--annotate FORMAT/GQ,FORMAT/GP,INFO/PV4 \
			-m \
			--prior {params.mcallPriorVal} \
			-O b \
			-o {output} \
			{input} )2> {log}
		"""


rule bcftools_call_genotypes_perpop_G:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}.bcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/gcMethod_mcall_perpop_prior{mcallPriorID}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_mcall_perpop_prior{mcallPriorID}.bcf"
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/gcMethod_mcall_perpop_prior{mcallPriorID}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_mcall_perpop_prior{mcallPriorID}.bcf"
	params:
		popinds_list=config["pop_inds_list"],
		mcallPriorVal=get_mcall_prior_val_from_id,
	shell:
		"""
		({BCFTOOLS} call \
			-G {params.popinds_list} \
			--annotate FORMAT/GQ,FORMAT/GP,INFO/PV4 \
			-m \
			--prior {params.mcallPriorVal} \
			-O b \
			-o {output} \
			{input} )2> {log}
		"""


rule bcftools_call_genotypes_hardcall:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}.bcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/gcMethod_hardcall/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_hardcall.bcf"
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/gcMethod_hardcall/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-vcfgl_platform{vcfgl_platform}_gl{vcfgl_gl}_qs{vcfgl_qsbeta}-gcMethod_hardcall.bcf"
	shell:
		"""
		( {BCFTOOLS} +tag2tag -Ob -o {output} {input} -- --GL-to-GT --threshold 1 )2> {log}
		"""


latex:

\subsection{Genotype calling}


\subsubsection{Naive genotype calling method}

For the naive genotype calling method, we used the following command:

\cmd{BCFtools +tag2tag -Ob -o {output} {input} -- --GL-to-GT --threshold 1}

This command converts genotype likelihoods to genotype calls, where the genotype call is the genotype with the highest likelihood. The genotype probability threshold is set to $1$ (i.e., no thresholds are applied). This is the simplest possible genotype calling method, and it is used as a baseline for comparison with other genotype calling methods.


\subsubsection{BCFtools multiallelic genotype caller}

We used BCFtools multiallelic caller \citep{Danecek2016-wc,Danecek2021-rh} to call genotypes with the default theta parameter value $0.0011$ (which is described as an approximate typical value for humans \citep{Danecek2016-wc}), and with the theta parameter disabled. We also compared two approaches: (1) genotype calling with all individuals pooled together, and (2) genotype calling per population. Overall, we used the following four commands:

\cmd{bcftools call -m --prior 0.0011}
\cmd{bcftools call -m --prior 0}
\cmd{bcftools call -m --prior 0.0011 -G {popinds_list}}
\cmd{bcftools call -m --prior 0 -G {popinds_list}}



