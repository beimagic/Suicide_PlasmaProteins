
import numpy as np
import pandas as pd
import re
import os
import random
import glob
from tqdm import tqdm
from itertools import product
from joblib import Parallel, delayed
from lightgbm import LGBMClassifier
from sklearn.metrics import roc_auc_score

def get_nb_f(mydf):
    auc_lst = mydf.AUC_mean.tolist()
    i = 0
    while((auc_lst[i]<auc_lst[i+1])|(auc_lst[i]<auc_lst[i+2])|(auc_lst[i]<auc_lst[i+3])):
        i+=1
    return i+1

def select_params_combo(my_dict, nb_items, my_seed):
    combo_list = [dict(zip(my_dict.keys(), v)) for v in product(*my_dict.values())]
    random.seed(my_seed)
    return random.sample(combo_list, nb_items)

def params_iter(mydf, my_f_lst, fold_id_lst, my_seed, my_params):
    auc_cv_lst = []
    my_params0 = my_params.copy()
    for fold_id in fold_id_lst:
        train_idx = mydf['FoldSplit'].index[mydf['FoldSplit'] != fold_id]
        test_idx = mydf['FoldSplit'].index[mydf['FoldSplit'] == fold_id]
        X_train, y_train = mydf.iloc[train_idx][my_f_lst], mydf.iloc[train_idx].target_y
        X_test, y_test = mydf.iloc[test_idx][my_f_lst], mydf.iloc[test_idx].target_y
        my_lgb = LGBMClassifier(objective='binary', metric='auc', is_unbalance=True, seed=my_seed)
        try:
            my_lgb.set_params(**my_params0)
            my_lgb.fit(X_train, y_train)
            y_pred_prob = my_lgb.predict_proba(X_test)[:, 1]
            auc_cv_lst.append(roc_auc_score(y_test, y_pred_prob))
        except:
            pass
    my_params0['AUC_cv_MEAN'] = np.round(np.mean(auc_cv_lst), 5)
    return my_params0

my_params = {'n_estimators': 500,
             'max_depth': 15,
             'num_leaves': 10,
             'subsample': 0.7,
             'learning_rate': 0.01,
             'colsample_bytree': 0.7}

params_dict = {'n_estimators': [100, 200, 300, 400, 500],
               'max_depth': np.linspace(5, 30, 6).astype('int32').tolist(),
               'num_leaves': np.linspace(5, 30, 6).astype('int32').tolist(),
               'subsample': np.linspace(0.6, 1, 9).tolist(),
               'learning_rate': [0.1, 0.05, 0.01, 0.001],
               'colsample_bytree': np.linspace(0.6, 1, 9).tolist()}

nb_params = 99
my_seed = 2024
nb_cpus = 10
type = 'Longitudinal'
imp_type = 'TotalGain'
fold_id_lst = [i for i in range(10)]

candidate_params_lst = select_params_combo(params_dict, nb_params, my_seed)
candidate_params_lst = [my_params] + candidate_params_lst

dpath = '/Volumes/JasonWork/Projects/OtherProjects/ZhangBei/SuicideProtein1005/'
output_file = dpath + 'Results/'+type+'_HoldOut/s3.0_ParameterTuning_'+imp_type+'_pro.csv'

pro_auc_df = pd.read_csv(dpath + 'Results/'+type+'_HoldOut/s20_SFS_'+imp_type+'.csv')
nb_f = get_nb_f(pro_auc_df)
pro_f_lst = pro_auc_df.Pro_code.tolist()[:nb_f]
pro_df = pd.read_csv('/Volumes/JasonWork/Projects/ProDisAtlas/Data/ProteinData/ProteinData.csv', usecols = ['eid'] + pro_f_lst)

target_df = pd.read_csv(dpath + 'Data/suicide_'+type+'.csv', usecols = ['eid', 'target_y', 'TrainTestSplit', 'FoldSplit'])
target_df = target_df.loc[target_df.TrainTestSplit == 0]
mydf = pd.merge(target_df, pro_df, how = 'inner', on = ['eid'])

my_params_lst = Parallel(n_jobs=nb_cpus)(delayed(params_iter)(mydf, pro_f_lst, fold_id_lst, my_seed, my_params) for my_params in candidate_params_lst)
params_df = pd.DataFrame(my_params_lst)
params_df.sort_values(by='AUC_cv_MEAN', ascending=False, inplace=True)
params_df.to_csv(output_file, index=False)

