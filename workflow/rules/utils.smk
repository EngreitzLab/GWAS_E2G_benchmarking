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

def validate_biosample_groups():
	key = pl.read_csv(config["predictionsTable"], separator="\t").drop_nulls("biosample")
	key_biosamples = key["biosample"].to_list()

	for thisGroup, thisBiosamples in config["biosampleGroups"].items():
		for b in thisBiosamples:
			if b not in key_biosamples:
				raise Exception(f"Biosample {b} in group {thisGroup} is not defined in predictions config.")

def validate_trait_groups():
	key_traits = variant_key["trait"].to_list()

	for thisGroup, thisTraits in config["traitGroups"].items():
		for t in thisTraits:
			if t not in key_traits:
				raise Exception(f"Trait {t} in group {thisGroup} is not definted in variant key.")

# add columns to methods config from sample key and config (all lists): biosamples (all individ biosamples for this method), predFiles (corresponding prediction files), biosampleGroups (which groups from config are fully represented)
def add_biosamples_and_files_to_config_old(methods_config):
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

def add_biosamples_and_files_to_config(methods_config):
	key = pl.read_csv(config["predictionsTable"], separator="\t").drop_nulls("biosample")
	methods_config = methods_config.filter(pl.col("method").is_in(config["methods"]))

	samples = []
	files= []
	groups = []
	for i in range(methods_config.height):
		this_method = methods_config["method"].item(i)
		this_key = key.select(["biosample", this_method]).drop_nulls()

		sampleList = this_key["biosample"].to_list()
		fileList = this_key[this_method].to_list()
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

def get_single_biosamples(method, name):
	this_biosamples = get_col2_from_col1(methods_config, "method", method, "biosamples")
	if name in this_biosamples: # if it's a single biosample...
		return name
	elif name=="ALL":
		return get_col2_from_col1(methods_config, "method", method, "biosamples")
	else: # it's a group that's not ALL
		return config["biosampleGroups"][name]

# get list of unique (indiividual) biosamples represented in comparisons provided AND by method... 
# also return the corresponding individual biosamples and associated traits (with repeats)
def get_comparison_biosamples_per_method(comparisons, method):
	this_biosamples = get_col2_from_col1(methods_config, "method", method, "biosamples")
	this_groups = get_col2_from_col1(methods_config, "method", method, "biosampleGroups")

	biosamples_all = comparisons['biosample'].to_list()
	traits_all = comparisons['trait'].to_list()
	groups = [g for g in biosamples_all if g in this_groups]
	trait_groups = [t for t,g in zip(traits_all, biosamples_all) if g in this_groups] 
	singles = [b for b in biosamples_all if b in this_biosamples]
	trait_singles = [t for t,b in zip(traits_all, biosamples_all) if b in this_biosamples]

	groups_exp = flatten([get_single_biosamples(method, g) for g in groups]) # expand group to individual biosamples
	trait_groups_exp = flatten([[t]*len(get_single_biosamples(method, g)) for t,g in zip(trait_groups, groups)]) # propogate traits
	
	all_represented_biosamples = singles + groups_exp
	all_represented_traits = trait_singles + trait_groups_exp
	uniq_list = list(set(all_represented_biosamples))
	return [uniq_list, all_represented_biosamples, all_represented_traits]

# get actual biosample names that are represented per method (from the table)
def get_biosample_names_per_method(comparisons, method):
	this_biosamples = get_col2_from_col1(methods_config, "method", method, "biosamples")
	this_groups = get_col2_from_col1(methods_config, "method", method, "biosampleGroups")
	this_all = this_biosamples + this_groups

	biosamples_all  = comparisons["biosample"].to_list()
	names = [n for n in biosamples_all if n in this_all]

	return names