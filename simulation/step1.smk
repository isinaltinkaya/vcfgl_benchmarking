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

simulation_id= "sim_vcfgl_2309"
configfile: "config/"+simulation_id+".yaml"


MODEL=config['model']

CONTIG=config['contig']
REP = [*range(config['n_reps'])]

species=stdpopsim.get_species("HomSap")
recombination_map_id="HapMapII_GRCh37"

DEF_POPS=config['def_pops']


ploidy=2

haplo_list=[]
indv_names=[]

for i, (key, value) in enumerate(DEF_POPS.items()):
	haplo_list.append(value*ploidy)
	for ind in range(value):
		# Using PLINK-like format: <Family-ID>_<Individual-ID>
		# to store <Population-ID>_<Individual-ID>
		indv_names.append(f"pop{key}_ind{str(ind+1)}")

# print(indv_names)
# print(haplo_list)
# print(len(indv_names))
# print(len(haplo_list))


###############################################################################
# BEGIN RULES

rule all:
	input:
		expand("sim/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}.trees",
			simid=simulation_id,
			model_id=MODEL,
			contig=CONTIG,
			rep=REP
			),
		expand("sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf",
			simid=simulation_id,
			model_id=MODEL,
			contig=CONTIG,
			rep=REP
			),



# ________________________________________________________________________________
# BEGIN SIMULATION



rule simulation:
	output: 
		"sim/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}.trees",
	params:
		seedfile="sim/{simid}/model_{model_id}/contig_{contig}/trees/.seed.{simid}-{model_id}-{contig}-rep{rep}.trees",
	run:
		model=species.get_demographic_model(wildcards.model_id)
		contig = species.get_contig(wildcards.contig)
		samples = model.get_samples(*haplo_list)
		engine = stdpopsim.get_engine("msprime")
		seedval=np.random.randint(2**20) + (100*int(wildcards.rep))  + int(str(wildcards.contig)[3:5]) # chr22 -> 22
		with open(params.seedfile,"w") as seedout:
			print("Contig index: "+str(int(str(wildcards.contig)[3:5]))+", chr: "+str(wildcards.contig)+", rep:"+str(wildcards.rep),file=seedout)
			print("Seed value: "+str(seedval), file=seedout)
		ts = engine.simulate(model,contig, samples, seed=seedval, msprime_model="dtwf")
		with open(output[0],"w") as tsout:
			ts.dump(tsout)


## Convert tree sequence to VCF file
# Using legacy format to avoid multiple instances of sites
rule tree_to_vcf:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}.trees",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf"
	run:
		ts=tskit.load(input[0])
		with open(output[0],"w") as vcfout:
			ts.write_vcf(output=vcfout,
					ploidy=ploidy,
					contig_id=str(wildcards.contig),
					position_transform="legacy",
					individual_names=indv_names)




