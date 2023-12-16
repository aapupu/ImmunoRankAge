from scalex.function import SCALEX
import  numpy as np
import pandas as pd
from pytorch_tabnet.pretraining import TabNetPretrainer
import torch
from pytorch_tabnet.tab_model import TabNetRegressor

# scRNA-seq data integrate
adata = SCALEX(data_list=file_list, # file list of scRNAseq
               batch_categories=batch, # batch information of scRNAseq
              outdir='/path/to/output',
              n_top_features=2000,
              min_cells=3, min_features=0, batch_size=256)

# multicollinearity filter
def set_diagonal_to_zero(df):
    arr = df.values
    np.fill_diagonal(arr, 0)
    df_modified = pd.DataFrame(arr, columns=df.columns, index=df.index)
    return df_modified

def multicollinearity_filter(df, target_cor_df, threshold=0.7):
    cor_matrix = df.corr()
    cor_matrix = set_diagonal_to_zero(cor_matrix)
    feature_list = list(cor_matrix.columns)
    selected_features = []
    for i in feature_list:
        if i in selected_features:
            continue
        feature_col = cor_matrix.loc[:,[i]]
        for j in feature_list:
            if j in selected_features:
                continue
            elif feature_col.loc[j,i] <= threshold:
                continue
            elif feature_col.loc[j,i] > threshold and np.abs(target_cor_df.loc[i,"Correlation"]) > np.abs(target_cor_df.loc[j,"Correlation"]):
                selected_features.append(j)
            elif feature_col.loc[j,i] > threshold and np.abs(target_cor_df.loc[i,"Correlation"]) <= np.abs(target_cor_df.loc[j,"Correlation"]):
                selected_features.append(i)
                break
    return selected_features

selected_features = multicollinearity_filter(GO_es_train, 
                                             GO_correlation_df, 
                                             threshold=0.7)

GO_es_train = GO_es_train.drop(columns=selected_features)

# model train
unsupervised_model = TabNetPretrainer(
    cat_emb_dim=3,
    optimizer_fn=torch.optim.Adam,
    optimizer_params=dict(lr=2e-2),
    mask_type='entmax', 
    n_shared_decoder=1,
    n_indep_decoder=1, 
    verbose=5
)

unsupervised_model.fit(
    X_train=x1,
    eval_set=[x2],
    max_epochs=100 , patience=20,
    batch_size=32, virtual_batch_size=128,
    num_workers=0,
    drop_last=False,
    pretraining_ratio=0.5,
) 

reg = TabNetRegressor()

reg.fit(
    X_train=x1, y_train=np.expand_dims(y1,axis=1),
    eval_set=[(x2, np.expand_dims(y2,axis=1))],
    max_epochs=200,
    patience=20,
    batch_size=8, 
    virtual_batch_size=64,
    num_workers=0,
    drop_last=True,
    pin_memory=False,
    from_unsupervised=unsupervised_model
)
