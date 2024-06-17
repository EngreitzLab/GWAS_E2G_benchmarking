import click
import numpy as np
import pandas as pd
from pybedtools import BedTool

def make_annotations(df_bim, predFile):
	iter_bim = [['chr'+str(x1), x2, x2, 1] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
	bimbed = BedTool(iter_bim)
	annotbed = bimbed.intersect(predFile, wb=True)
	bp = [x.start for x in annotbed]
	score = [float(x.fields[7]) for x in annotbed]
	df_int = pd.DataFrame({'BP': bp, 'ANNOT':score})
	df_annot = pd.merge(df_bim, df_int, how='left', on='BP')
	df_annot.fillna(0, inplace=True)
	temp = df_annot[['ANNOT']].astype(float)
	df_annot = pd.concat([df_bim.iloc[:,[0,3,1,2]], temp], axis = 1)
	
	return df_annot

@click.command()
@click.option("--pred_file", required=True)
@click.option("--bim_file", required=True)
@click.option("--output_file", required=True)

def main(pred_file, bim_file, output_file):
	df_bim = pd.read_csv(bim_file, delim_whitespace=True, usecols=[0,1,2,3], names = ['CHR','SNP','CM','BP'])
	df_annot = make_annotations(df_bim, pred_file)
	df_annot.to_csv(output_file, compression="gzip", sep="\t", index=False)

if __name__ == "__main__":
	main()
