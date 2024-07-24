import os
import polars as pl

## SAVE RELEVANT PROCESSED CONFIG FILES

# save methods, traitGroups, biosampleGroups, etc.
rule save_config_params:
	params:
		methods = config["methods"],
		biosampleGroups = config["biosampleGroups"],
		traitGroups = config["traitGroups"],
		numPoPSGenes = config["numPoPSGenes"],
		numPredGenes = config["numPredGenes"],
		scratchDir = SCRATCH_DIR
	output:
		outfile =  os.path.join(RESULTS_DIR, "reference_configs", "config_params.tsv")
	shell:
		"""
		set +o pipefail;
		
		echo "METHODS:  {params.methods}" > {output.outfile}
		echo "BIOSAMPLE_GROUPS:  {params.biosampleGroups}" >> {output.outfile}
		echo "TRAIT_GROUPS:  {params.traitGroups}" >> {output.outfile}
		echo "NUM_POPS_GENES:  {params.numPoPSGenes}" >> {output.outfile}
		echo "NUM_PRED_GENES: {params.numPredGenes}" >> {output.outfile}
		echo "SCRATCH_DIR:  {params.scratchDir}" >> {output.outfile}

		"""

rule save_comparisons_table:
	input:
		comparisons_file = config["comparisonsTable"]
	output:
		outfile =  os.path.join(RESULTS_DIR, "reference_configs", "comparisons.tsv")
	shell:
		"""
		set +o pipefail;
		
		cat {input.comparisons_file} > {output.outfile}
		"""

rule save_trait_and_biosample_groups:
	params:
		traitGroups = config["traitGroups"],
		biosampleGroups = config["biosampleGroups"]
	output:
		outTraits = os.path.join(RESULTS_DIR, "reference_configs", "trait_groups.tsv"),
		outBiosamples = os.path.join(RESULTS_DIR, "reference_configs", "biosample_groups.tsv")
	run:
		trait_num = [len(traits) for traits in params.traitGroups.values()]
		max_len = max(trait_num)
		trait_dict = {group : traits + [None] * (max_len-len(traits)) for group,traits in params.traitGroups.items()}
		pl.from_dict(trait_dict).write_csv(output.outTraits, separator="\t")
		
		bs_num = [len(bs) for bs in params.biosampleGroups.values()]
		max_len = max(bs_num)
		bs_dict = {group : bs + [None] * (max_len-len(bs)) for group,bs in params.biosampleGroups.items()}
		pl.from_dict(bs_dict).write_csv(output.outBiosamples, separator="\t")

rule save_methods_config:
	params:
		methods = methods_config
	output:
		outMethods = os.path.join(RESULTS_DIR, "reference_configs", "methods_config.json")
	run:
		params.methods.write_json(output.outMethods)

