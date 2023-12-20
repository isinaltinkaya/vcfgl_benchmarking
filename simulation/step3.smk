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


simulation_id = "sim_vcfgl_2312"


configfile: "config/" + simulation_id + ".yaml"


MODEL = config["model"]
ERROR_RATE = config["vcfgl_error_rate"]
DEPTH = config["depth"]
CONTIG = config["contig"]
REP = [*range(config["n_reps"])]

# BETA_VARS = config["beta_variance_values_neg_e"]

ANGSD = config["tools"]["angsd"]
BCFTOOLS = config["tools"]["bcftools"]


###############################################################################
# BEGIN RULES

POPS = config["pops"]

GLMODELS = [1, 2]

rule all:
    input:
        expand(
            "sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_perpop_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snp.bcf",
            simid=simulation_id,
            model_id=MODEL,
            contig=CONTIG,
            rep=REP,
            depth=DEPTH,
            error_rate=ERROR_RATE,
            qserr=["0_0", "2_5", "2_6", "2_7"],
			glModel=GLMODELS,
        ),
        expand(
            "sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snp.bcf",
            simid=simulation_id,
            model_id=MODEL,
            contig=CONTIG,
            rep=REP,
            depth=DEPTH,
            error_rate=ERROR_RATE,
            qserr=["0_0", "2_5", "2_6", "2_7"],
			glModel=GLMODELS,
        ),


# ###############################################################################


rule bcftools_call_snp_and_genotypes:
    input:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qserr}.bcf",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snp.bcf",
    log:
        "sim/{simid}/logs/model_{model_id}/contig_{contig}/genotype_calling_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snp.bcf",
    params:
        prior=config["bcftools_call_m_prior"],
    shell:
        """
        ({BCFTOOLS} call --annotate FORMAT/GQ,FORMAT/GP,INFO/PV4 -mv --prior {params.prior} -O b -o {output} {input} )2> {log}
        """


###############################################################################


rule split_pops_for_bcftools_call:
    input:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qserr}.bcf",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qserr}_pop{population}.bcf",
    log:
        "sim/{simid}/logs/model_{model_id}/contig_{contig}/vcfgl_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qserr}_pop{population}.bcf",
    params:
        pop_indlist=lambda wildcards: config["perpop_inds_list"][wildcards.population],
    shell:
        """
        ({BCFTOOLS} view -S {params.pop_indlist} {input} -O b -o {output} ) 2> {log}
        """


rule bcftools_call_snp_and_genotypes_perPop:
    input:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}_qs{qserr}_pop{population}.bcf",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_perpop_gl{glModel}/perpop/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}_pop{population}-snp.bcf",
    log:
        "sim/{simid}/logs/model_{model_id}/contig_{contig}/genotype_calling_perpop_gl{glModel}/perpop/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}_pop{population}-snp.bcf",
    params:
        prior=config["bcftools_call_m_prior"],
    shell:
        """
        ({BCFTOOLS} call --annotate FORMAT/GQ,FORMAT/GP,INFO/PV4 -mv --prior {params.prior} -O b -o {output} {input} )2> {log}
        """


# 231209 annotate -x: hacky solution to bcftools bug
# reported in github.com/samtools/bcftools issue #2055
rule rmgl_and_index_perPop_bcf_for_merge_snp:
    input:
        "sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_perpop_gl{glModel}/perpop/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}_pop{population}-snp.bcf",
    output:
        bcf="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_perpop_gl{glModel}/perpop/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}_pop{population}-snp-rmgl.bcf",
        bcfcsi="sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_perpop_gl{glModel}/perpop/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}_pop{population}-snp-rmgl.bcf.csi",
    log:
        "sim/{simid}/logs/model_{model_id}/contig_{contig}/genotype_calling_perpop_gl{glModel}/perpop/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}_pop{population}-snp-rmgl.bcf.csi",
    shell:
        """
        ({BCFTOOLS} annotate -x "FORMAT/GL" {input} -O b > {output.bcf}; {BCFTOOLS} index {output.bcf} )2> {log}
        """


rule merge_perpop_genotype_calls_snps:
    input:
        bcf=expand(
            "sim/{{simid}}/model_{{model_id}}/contig_{{contig}}/genotype_calling_perpop_gl{{glModel}}/perpop/{{simid}}-{{model_id}}-{{contig}}-rep{{rep}}-d{{depth}}-e{{error_rate}}-qs{{qserr}}_pop{population}-snp-rmgl.bcf",
            population=POPS,
        ),
        csi=expand(
            "sim/{{simid}}/model_{{model_id}}/contig_{{contig}}/genotype_calling_perpop_gl{{glModel}}/perpop/{{simid}}-{{model_id}}-{{contig}}-rep{{rep}}-d{{depth}}-e{{error_rate}}-qs{{qserr}}_pop{population}-snp-rmgl.bcf.csi",
            population=POPS,
        ),
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/genotype_calling_perpop_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snp.bcf",
    log:
        "sim/{simid}/logs/model_{model_id}/contig_{contig}/genotype_calling_perpop_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snp.bcf",
    shell:
        """
        ( {BCFTOOLS} merge {input.bcf} -O b -o {output} )2> {log}
        """
