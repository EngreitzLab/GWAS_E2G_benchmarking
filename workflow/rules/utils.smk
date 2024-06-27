import os
import polars as pl

MAX_MEM_MB = 250 * 1000  # 250GB

def determine_mem_mb(wildcards, input, attempt, min_gb=8):
	# Memory resource calculator for snakemake rules
	input_size_mb = input.size_mb
	if ".gz" in str(input):
		input_size_mb *= 8  # assume gz compressesed the file <= 8x
	attempt_multiplier = 2 ** (attempt - 1)  # Double memory for each retry
	mem_to_use_mb = attempt_multiplier *  max(4 * input_size_mb, min_gb * 1000)
	return min(mem_to_use_mb, MAX_MEM_MB)

def flatten(list_of_lists):
	return [x for xs in list_of_lists for x in xs]

def add_biosamples_and_files_to_config(methods_config):
	biosampleKeys = methods_config['sampleKey']
	methods = methods_config['method']

	samples = []
	files= []
	groups = []
	for i in range(len(methods)):
		key = pl.read_csv(biosampleKeys[i], separator='\t')
		sampleList = key["biosample"].to_list()
		fileList =  key["predictionFile"].to_list()
		samples.append(sampleList)
		files.append(fileList)

		# check biosample groups
		groupsPresent = ["ALL"]
		for thisGroup, thisBiosamples in config["biosampleGroups"].items():
			if all(b in sampleList for b in thisBiosamples):
				groupsPresent.append(thisGroup)
		groups.append(groupsPresent)

	# create new columns using pl.Series to avoid nested lists
	methods_config = methods_config.with_columns([
		pl.Series('biosamples', samples),
		pl.Series('predFiles', files),
		pl.Series('biosampleGroups', groups)
	])

	return methods_config

# self-explanatory; (does not handle multiple matches, since some values are lists)
def get_col2_from_col1(df, col1, col1_val, col2):
	col2_ser = df.filter(pl.col(col1) == col1_val)[col2]
	return col2_ser.to_list()[0]


def get_pred_file(method, biosample):
	this_biosamples = get_col2_from_col1(methods_config, "method", method, "biosamples")
	this_files = get_col2_from_col1(methods_config, "method", method, "predFiles")

	return this_files[this_biosamples.index(biosample)]

def get_biosamples_for_group(method, group):
	if group=="ALL":
		x = get_col2_from_col1(methods_config, "method", method, "biosamples")
		return get_col2_from_col1(methods_config, "method", method, "biosamples")
	else:
		x = config["biosampleGroups"][group]
		return config["biosampleGroups"][group]

