
import numpy as np
import pandas as pd
from Utility.Training_Utilities import *
from lightgbm import LGBMClassifier
from sklearn.metrics import roc_auc_score
from sklearn.calibration import CalibratedClassifierCV

def get_nb_f(mydf):
    auc_lst = mydf.AUC_mean.tolist()
    i = 0
    while((auc_lst[i]<auc_lst[i+1])|(auc_lst[i]<auc_lst[i+2])|(auc_lst[i]<auc_lst[i+3])):
        i+=1
    return i+1

type = 'Longitudinal'
imp_type = 'TotalGain'

dpath = '/Volumes/JasonWork/Projects/OtherProjects/ZhangBei/SuicideProtein1005/'
outputfile = dpath + 'Results/'+type+'_HoldOut/s4_Predictions.csv'

param_df_pro = pd.read_csv(dpath + 'Results/'+type+'_HoldOut/s3.0_ParameterTuning_'+imp_type+'_pro.csv')
my_params_pro = dict(param_df_pro.iloc[0, :6])
my_params_pro['n_estimators'] = int(my_params_pro['n_estimators'])
my_params_pro['max_depth'] = int(my_params_pro['max_depth'])
my_params_pro['num_leaves'] = int(my_params_pro['num_leaves'])

param_df_pro_cov = pd.read_csv(dpath + 'Results/'+type+'_HoldOut/s3.1_ParameterTuning_'+imp_type+'_pro_cov.csv')
my_params_pro_cov = dict(param_df_pro_cov.iloc[0, :6])
my_params_pro_cov['n_estimators'] = int(my_params_pro_cov['n_estimators'])
my_params_pro_cov['max_depth'] = int(my_params_pro_cov['max_depth'])
my_params_pro_cov['num_leaves'] = int(my_params_pro_cov['num_leaves'])


pro_auc_df = pd.read_csv(dpath + 'Results/'+type+'_HoldOut/s20_SFS_'+imp_type+'.csv')
nb_f = get_nb_f(pro_auc_df)
pro_f_lst = pro_auc_df.Pro_code.tolist()[:nb_f]
pro_df = pd.read_csv('/Volumes/JasonWork/Projects/ProDisAtlas/Data/ProteinData/ProteinData.csv', usecols = ['eid'] + pro_f_lst)
target_df = pd.read_csv(dpath + 'Data/suicide_'+type+'.csv', usecols = ['eid', 'target_y', 'Age', 'Sex', 'TrainTestSplit'])
tgtdf_train, tgtdf_test = target_df.copy(), target_df.copy()
tgtdf_train = tgtdf_train.loc[tgtdf_train.TrainTestSplit == 0]
tgtdf_test = tgtdf_test.loc[tgtdf_test.TrainTestSplit == 1]
traindf = pd.merge(tgtdf_train, pro_df, how = 'inner', on = ['eid'])
testdf = pd.merge(tgtdf_test, pro_df, how = 'inner', on = ['eid'])
my_f_lst = ['Age', 'Sex'] + pro_f_lst

my_lgb1 = LGBMClassifier(objective='binary', metric='auc', is_unbalance=True, n_jobs=4, verbosity=-1, seed=2020)
my_lgb1.set_params(**my_params_pro)
my_lgb1.fit(traindf[pro_f_lst], traindf.target_y)
testdf['y_pred_pro'] = my_lgb1.predict_proba(testdf[pro_f_lst])[:, 1].tolist()

my_lgb2 = LGBMClassifier(objective='binary', metric='auc', is_unbalance=True, n_jobs=4, verbosity=-1, seed=2020)
my_lgb2.set_params(**my_params_pro_cov)
my_lgb2.fit(traindf[my_f_lst], traindf.target_y)
testdf['y_pred_pro_cov'] = my_lgb2.predict_proba(testdf[my_f_lst])[:, 1].tolist()

print(roc_auc_score(testdf.target_y, testdf.y_pred_pro))
print(roc_auc_score(testdf.target_y, testdf.y_pred_pro_cov))

outdf = testdf[['eid', 'target_y', 'TrainTestSplit', 'y_pred_pro', 'y_pred_pro_cov']]
outdf.to_csv(outputfile, index = False)


