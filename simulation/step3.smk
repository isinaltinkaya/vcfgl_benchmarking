###############################################################################
# Step 3: SNP calling
#
# isinaltinkaya
###############################################################################
#
# Naming convention:
# - log files:
# output: sim/{simid}/REST_OF_PATH.EXTENSION
# params.prefix: sim/{simid}/REST_OF_PATH
# log: sim/{simid}/logs/REST_OF_PATH.EXTENSION

import stdpopsim, msprime, tskit
import numpy as np

simulation_id = "sim_vcfgl_2309"


configfile: "config/" + simulation_id + ".yaml"


VCFGL = config["tools"]["vcfgl"]

MODEL = config["model"]

ERROR_RATE = config["vcfgl_error_rate"]

VCFGL = config["tools"]["vcfgl"]

DEPTH = config["depth"]
CONTIG = config["contig"]
REP = [*range(config["n_reps"])]

BETA_VARS = config["beta_variance_values_neg_e"]

###############################################################################
# BEGIN RULES


rule all:
    input:
        expand(
            "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/nogt/betavar{beta_variance}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-betavar{beta_variance}-nogt.bcf",
            simid=simulation_id,
            model_id=MODEL,
            contig=CONTIG,
            rep=REP,
            depth=DEPTH,
            error_rate=ERROR_RATE,
            beta_variance=BETA_VARS,
        ),
        expand(
            "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/nogt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-nogt.bcf",
            simid=simulation_id,
            model_id=MODEL,
            contig=CONTIG,
            rep=REP,
            depth=DEPTH,
            error_rate=ERROR_RATE,
        ),
        # expand(
        #     "sim/{simid}/model_{model_id}/contig_{contig}/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-snp.bcf",
        #     simid=simulation_id,
        #     model_id=MODEL,
        #     contig=CONTIG,
        #     rep=REP,
        #     depth=DEPTH,
        #     error_rate=ERROR_RATE,
        # ),
        # expand(
        #     "sim/{simid}/model_{model_id}/contig_{contig}/snp_calling/betavar{beta_variance}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-betavar{beta_variance}-snp.bcf",
        #     simid=simulation_id,
        #     model_id=MODEL,
        #     contig=CONTIG,
        #     rep=REP,
        #     depth=DEPTH,
        #     error_rate=ERROR_RATE,
        #     beta_variance=BETA_VARS,
        # ),


# SIMULATE GENOTYPE LIKELIHOODS


# remove the simulated genotypes from msprime from vcf files
rule remove_gts_from_vcfgl_vcfs:
    input:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}.bcf",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/nogt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-nogt.bcf",
    log:
        "sim/{simid}/logs/model_{model_id}/contig_{contig}/vcfgl/nogt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-nogt.bcf",
    shell:
        """
        ({BCFTOOLS} annotate -x FORMAT/GT {input} -O b -o {output} )2> {log}
        """


rule remove_gts_from_vcfgl_vcfs_qserror:
    input:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/betavar{beta_variance}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-betavar{beta_variance}.bcf",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/nogt/betavar{beta_variance}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-betavar{beta_variance}-nogt.bcf",
    log:
        "sim/{simid}/logs/model_{model_id}/contig_{contig}/vcfgl/nogt/betavar{beta_variance}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-betavar{beta_variance}-nogt.bcf",
    shell:
        """
        ({BCFTOOLS} annotate -x FORMAT/GT {input} -O b -o {output} )2> {log}
        """


# rule bcftools_call_genotypes:
# 	input:
# 		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/nogt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-nogt.bcf",


# 	shell:
# 	"""
# 	"""


# perform SNP calling
rule call_snps_angsd_vcfgl:
    input:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/nogt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-nogt.bcf",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-snp.bcf",
    params:
        snppval=config["SNP_calling_SNP_pval"],
    log:
        "sim/{simid}/logs/model_{model_id}/contig_{contig}/snp_calling/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-snp.bcf",
    shell:
        """
        ( {ANGSD} \
                -vcf-gl {input} \
                -out {params.prefix} \
                -doMajorMinor 1 \
                -doMaf 1 \
                -skipTriallelic 1 \
                -SNP_pval {params.snppval} \
                -doBcf 1 \
                -doGlf 2 \
                ) 2> {log}
        """


rule call_snps_angsd_vcfgl_qserror:
    input:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl/nogt/betavar{beta_variance}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-betavar{beta_variance}-nogt.bcf",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/snp_calling/betavar{beta_variance}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-betavar{beta_variance}-snp.bcf",
    params:
        snppval=config["SNP_calling_SNP_pval"],
    log:
        "sim/{simid}/logs/model_{model_id}/contig_{contig}/snp_calling/betavar{beta_variance}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-betavar{beta_variance}-snp.bcf",
    shell:
        """
        ( {ANGSD} \
                -vcf-gl {input} \
                -out {params.prefix} \
                -doMajorMinor 1 \
                -doMaf 1 \
                -skipTriallelic 1 \
                -SNP_pval {params.snppval} \
                -doBcf 1 \
                -doGlf 2 \
                ) 2> {log}
        """
