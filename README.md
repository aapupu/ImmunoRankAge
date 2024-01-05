# ImmunoRankAge
A machine learning model for predicting immune age by rank-enrichment algorithm. 
This model, operating in both CPU and GPU modes, quickly predicts over 200 samples in less than 10 seconds.

![image](https://github.com/aapupu/ImmunoRankAge/blob/main/img/ImmunoRankAge.jpg)


## Requirements 

- R == 4.1.0
- data.table == 1.14.6
- matrixStats == 0.63.0
- python == 3.7.3
- pytorch == 1.14.0
- numpy == 1.21.6
- pandas == 1.3.5
- pytorch-tabnet == 4.0
  
**Note: If cuda is available, GPU modes will be run automatically.**
  

Tutorial
-------
Git clone
```bash
git clone git://github.com/aapupu/ImmunoRankAge.git
```

Run in ImmunoRankAge folder
```bash
cd ImmunoRankAge
```

Run ImmunoRankAge
```bash
python ImmunoRankAge_main.py --file_path /path/to/RNAseq.txt --norwayname tpm/count
```

## Output 

- outs.csv: predicted immune age of sample.
- enrich_score.csv: enrichment score of input feature by rank-enrichment algorithm.


Citation
-------
**Decoding Human Immune Aging Dynamics: Insights from PBMC Transcriptome and the ImmunoRankAge Model**


Contacts
-------
kyzy850520@163.com
