###############################################################################
# Step 5: Collect and plot genotype discordance metrics
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

simulation_id = "sim_vcfgl_2312"


configfile: "config/" + simulation_id + ".yaml"


VCFGL = config["tools"]["vcfgl"]
MODEL = config["model"]
ERROR_RATE = config["vcfgl_error_rate"]


DEPTH = config["depth"]
CONTIG = config["contig"]
REP = [*range(config["n_reps"])]

BETA_VARS = config["beta_variance_values_neg_e"]


BCFTOOLS = config["tools"]["bcftools"]
PICARDJAR = config["tools"]["picardjar"]


GTCHECK_ERR = config["bcftools_gtcheck_error"]

ploidy = 2
DEF_POPS = config["def_pops"]

haplo_list = []
indv_names = []

for i, (key, value) in enumerate(DEF_POPS.items()):
    haplo_list.append(value * ploidy)
    for ind in range(value):
        # Using PLINK-like format: <Family-ID>_<Individual-ID>
        # to store <Population-ID>_<Individual-ID>
        indv_names.append(f"pop{key}_ind{str(ind+1)}")


###############################################################################
# BEGIN RULES


GC_METHODS = ["genotype_calling", "genotype_calling_perpop"]
GQ_FILTERS = ["-minGQ20", ""]

# GQ_FILTERS = ["-minGQ20", ""]
MINGQ = 20

GLMODELS = [1, 2]


rule all:
    input:
        expand(
            "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance_qs_doGQ3_gl{glModel}/{simid}-{model_id}-{gc_method}.tsv",
            simid=simulation_id,
            model_id=MODEL,
            gc_method=GC_METHODS,
            glModel=GLMODELS,
        ),
        expand(
            "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qserr}-snp_doGQ3_qs_tidy.tsv",
            simid=simulation_id,
            model_id=MODEL,
            contig=CONTIG,
            rep=REP,
            depth=DEPTH,
            error_rate=ERROR_RATE,
            gc_method=GC_METHODS,
            qserr=["0_0", "2_5", "2_6", "2_7"],
            glModel=GLMODELS,
        ),


rule tidy_get_genotype_discordance_doGQ3:
    input:
        "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsval}_{betavar}-snp_doGQ3_qs.tsv",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/gc_evaluation/genotype_discordance/{gc_method}_gl{glModel}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsval}_{betavar}-snp_doGQ3_qs_tidy.tsv",
    shell:
        """
        awk -v REP={wildcards.rep} -v DEPTH={wildcards.depth} -v QSVAL={wildcards.qsval} -v BETAVAR={wildcards.betavar} 'BEGIN{{FS="\t";OFS="\t"}}{{print $0,REP,DEPTH,QSVAL,BETAVAR}}' {input} > {output}
        """


rule collect_get_genotype_discordance_qs_doGQ3:
    input:
        expand(
            "sim/{{simid}}/model_{{model_id}}/contig_{contig}/gc_evaluation/genotype_discordance/{{gc_method}}_gl{{glModel}}/{{simid}}-{{model_id}}-{contig}-rep{rep}-d{depth}-e{error_rate}-qs{qsbeta}-snp_doGQ3_qs_tidy.tsv",
            contig=CONTIG,
            depth=DEPTH,
            rep=REP,
            error_rate=ERROR_RATE,
            betavar=BETA_VARS,
            qsbeta=["0_0", "2_5", "2_6", "2_7"],
        ),
    output:
        "sim/{simid}/model_{model_id}/gc_evaluation/genotype_discordance_qs_doGQ3_gl{glModel}/{simid}-{model_id}-{gc_method}.tsv",
    shell:
        """
        cat {input} > {output}
        """
