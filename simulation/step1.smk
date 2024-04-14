###############################################################################
# Step 1: Simulate data for vcfgl benchmarking
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

import os

simulation_id = config["simulation_id"]


MODEL = config["model"]
CONTIG = config["contig"]

species = stdpopsim.get_species("HomSap")

DEF_POPS = config["def_pops"]
ploidy = 2

haplo_list = []
indv_names = []

samples_pops_dict = {}

for i, (key, value) in enumerate(DEF_POPS.items()):
    haplo_list.append(value * ploidy)
    for ind in range(value):
        # Using PLINK-like format: <Family-ID>_<Individual-ID>
        # to store <Population-ID>_<Individual-ID>
        indv_names.append(f"pop{key}_ind{str(ind+1)}")
        # add to samples_pops_dict
        samples_pops_dict[indv_names[-1]] = key


###############################################################################
# BEGIN RULES


rule all:
    input:
        expand(
            "sim/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}.trees",
            simid=simulation_id,
            model_id=MODEL,
            contig=CONTIG,
        ),
        expand(
            "sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}.vcf",
            simid=simulation_id,
            model_id=MODEL,
            contig=CONTIG,
        ),
        expand(
            "sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}.vcf.isValid",
            simid=simulation_id,
            model_id=MODEL,
            contig=CONTIG,
        ),
        config["pop_inds_list"],


# ________________________________________________________________________________
# BEGIN SIMULATION


rule simulation:
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}.trees",
    params:
        logfile="sim/{simid}/logs/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}.trees",
    log:
        "sim/{simid}/logs/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}.trees",
    run:
        model = species.get_demographic_model(wildcards.model_id)
        contig = species.get_contig(wildcards.contig, mutation_rate=model.mutation_rate)
        samples = DEF_POPS
        engine = stdpopsim.get_engine("msprime")
        seedval = np.random.randint(2**20) + 42
        with open(params.logfile, "w") as logout:
            print(
                "Contig index: "
                + str(int(str(wildcards.contig)[3:5]))
                + ", chr: "
                + str(wildcards.contig),
                file=logout,
            )
            print("Seed value: " + str(seedval), file=logout)
            ts = engine.simulate(
                model, contig, samples, seed=seedval, msprime_model="dtwf"
            )
        with open(output[0], "w") as tsout:
            ts.dump(tsout)


## Convert tree sequence to VCF file
rule tree_to_vcf:
    input:
        "sim/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}.trees",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}.vcf",
    run:
        ts = tskit.load(input[0])
        with open(output[0], "w") as vcfout:
            ts.write_vcf(
                output=vcfout,
                contig_id=str(wildcards.contig),
                individual_names=indv_names,
            )


BCFTOOLS = config["tools"]["bcftools"]


# validate:
# 1) make sure the first position is not 0
# 2) make sure no position is repeated
# 3) make sure the positions are sorted in ascending order
# 0) get the positions from the vcf file
rule validate_vcf:
    input:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}.vcf",
    output:
        "sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}.vcf.isValid",
    shell:
        """
        # 0) get the positions from the vcf file
        {BCFTOOLS} query -f "%CHROM %POS\\n" {input} > {input}.pos
        # 1) make sure the first position is not 0
        first_pos=$(head -n 1 {input}.pos | cut -d " " -f 2)
        if [ $first_pos -eq 0 ]; then
            echo "First position is 0"
            exit 1
        fi
        filelen=$(wc -l < {input}.pos)
        sort -k2 -g {input}.pos > {input}.pos.sorted
        uniq_filelen=$(wc -l < {input}.pos.sorted)
        if [ $filelen -ne $uniq_filelen ]; then
            echo "Position is repeated"
            exit 1
        fi
        # 3) make sure the positions are sorted in ascending order
        if ! cmp -s {input}.pos {input}.pos.sorted; then
            echo "Positions are not sorted in ascending order"
            exit 1
        fi
        rm {input}.pos {input}.pos.sorted
        touch {output}
        """


rule print_pop_inds_list:
    output:
        config["pop_inds_list"],
    run:
        with open(output[0], "w") as pop_inds_list:
            for indv in indv_names:
                print(f"{indv}\t{samples_pops_dict[indv]}", file=pop_inds_list)
