import os
import pandas as pd

MAX_MEM_MB = 250 * 1000  # 250GB

def determine_mem_mb(wildcards, input, attempt, min_gb=8):
	# Memory resource calculator for snakemake rules
	input_size_mb = input.size_mb
	if ".gz" in str(input):
		input_size_mb *= 8  # assume gz compressesed the file <= 8x
	attempt_multiplier = 2 ** (attempt - 1)  # Double memory for each retry
	mem_to_use_mb = attempt_multiplier *  max(4 * input_size_mb, min_gb * 1000)
	return min(mem_to_use_mb, MAX_MEM_MB)

def add_biosamples_and_files_to_config():
	biosampleKeys = methods_config['sampleKey']
	methods = methods_config['method']

	samples = []
	files = []
	for i in range(len(methods)):
		if biosampleKeys.iloc[i]=="None":
			samples.append("")
			files.append("")
		else:
			key = pd.read_csv(biosampleKeys.iloc[i], sep='\t')
			sampleList = key["biosample"].tolist()
			fileList =  key["predictionFile"].tolist()
			fileDict = dict(zip(sampleList, fileList))
			samples.append(sampleList)
			files.append(fileDict)

	methods_config['biosamples'] = samples # add samples to config
	methods_config['predFiles'] = files 

def get_biosamples_from_tissue(method, tissue):
	sample_key = pd.read_csv(methods_config.loc[method, "sampleKey"], sep="\t")
	if tissue=="ALL":
		biosamples = sample_key['biosample'].tolist()
	else:
		sample_key = sample_key.dropna(subset='tissue') # drop rows with no tissue mapping
		biosamples = []
		for index, row in sample_key.iterrows():
			tissues_this_row = [x.strip() for x in row["tissue"].split(',')] # if tissue is matched with this biosample, add to list
			if tissue in tissues_this_row:
				biosamples.append(row["biosample"])
	return biosamples

def flatten(list_of_lists):
	return [x for xs in list_of_lists for x in xs]