import os, argparse, time
import numpy as np
import pandas as pd
from pytorch_tabnet.tab_model import TabNetRegressor
#import torch

parser = argparse.ArgumentParser(description='usage:File path and type of RNAseq')
parser.add_argument('--file_path','-f',type=str, required=True, help="File path of RNAseq: XXX.txt")
parser.add_argument('--norwayname','-n',type=str, required=True,help="File type of RNAseq: tpm or count")
args = parser.parse_args()

def Enrich(file_path, norwayname):
    run="Rscript Requirements/Calculate_enrichment_score.R " + file_path + " " + norwayname
    print(run)
    os.system(run)


def AgePredict():
    reg = TabNetRegressor()
    reg.load_model('Requirements/tabnet_model_0.3.zip')

    features = pd.read_csv('enrich_score.csv',index_col=0)
    x = np.array(features)
    y = reg.predict(np.array(x)).squeeze(1)

    features['Age_pre'] = y
    features.to_csv('outs.csv')


if __name__ == '__main__': 
    start_total_time = time.time() 
    
    Enrich(args.file_path, args.norwayname)
    AgePredict()
    os.remove('enrich_score.csv')
    
    end_total_time = time.time()
    elapsed_total_time = end_total_time - start_total_time
    print(f'Total execution time: {elapsed_total_time:.2f} seconds.')
    print('Finish')