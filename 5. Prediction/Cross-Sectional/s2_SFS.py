
import numpy as np
import pandas as pd
from Utility.Training_Utilities import *
from Utility.DelongTest import delong_roc_test
from lightgbm import LGBMClassifier
import scipy.stats as stats
pd.options.mode.chained_assignment = None  # default='warn'

type = 'CrossSectional'
imp_type = 'TotalGain'

dpath = '/Volumes/JasonWork/Projects/OtherProjects/ZhangBei/SuicideProtein1005/'
outfile = dpath + 'Results/'+type+'_HoldOut/s20_SFS_'+imp_type+'.csv'

pro_f_df = pd.read_csv(dpath + 'Results/'+type+'_HoldOut/s1_PROImportance.csv')
pro_f_df.sort_values(by = imp_type+'_cv', ascending=False, inplace = True)
pro_f_lst = pro_f_df.Pro_code.tolist()[:50]
pro_df = pd.read_csv('/Volumes/JasonWork/Projects/ProDisAtlas/Data/ProteinData/ProteinData.csv', usecols = ['eid'] + pro_f_lst)
target_df = pd.read_csv(dpath + 'Data/suicide_'+type+'.csv', usecols = ['eid', 'target_y', 'TrainTestSplit', 'FoldSplit'])
target_df = target_df.loc[target_df.TrainTestSplit == 0]
mydf = pd.merge(target_df, pro_df, how = 'inner', on = ['eid'])

fold_id_lst = [i for i in range(10)]

y_test_full = np.zeros(shape = [1,1])
for fold_id in fold_id_lst:
    test_idx = mydf['FoldSplit'].index[mydf['FoldSplit'] == fold_id]
    y_test_full = np.concatenate([y_test_full, np.expand_dims(mydf.iloc[test_idx].target_y, -1)])

y_pred_full_prev = y_test_full
tmp_f, AUC_cv_lst= [], []

for f in pro_f_lst:
    tmp_f.append(f)
    my_X = mydf[tmp_f]
    AUC_cv = []
    y_pred_full = np.zeros(shape = [1,1])
    for fold_id in fold_id_lst:
        train_idx = mydf['FoldSplit'].index[mydf['FoldSplit'] != fold_id]
        test_idx = mydf['FoldSplit'].index[mydf['FoldSplit'] == fold_id]
        X_train, X_test = mydf.iloc[train_idx][tmp_f], mydf.iloc[test_idx][tmp_f]
        y_train, y_test = mydf.iloc[train_idx].target_y, mydf.iloc[test_idx].target_y
        my_lgb = LGBMClassifier(objective='binary', metric='auc', is_unbalance=True,  n_jobs=4, verbosity=-1, seed=2020)
        my_lgb.fit(X_train, y_train)
        y_pred_prob = my_lgb.predict_proba(X_test)[:, 1]
        AUC_cv.append(np.round(roc_auc_score(y_test, y_pred_prob), 3))
        y_pred_full = np.concatenate([y_pred_full, np.expand_dims(y_pred_prob, -1)])
    log10_p = delong_roc_test(y_test_full[:,0], y_pred_full_prev[:,0], y_pred_full[:,0])
    AUC_full = roc_auc_score(y_test_full[:,0], y_pred_full[:,0])
    y_pred_full_prev = y_pred_full
    tmp_out = np.array([np.round(np.mean(AUC_cv), 5), np.round(np.std(AUC_cv), 5), np.round(AUC_full, 5), 10**log10_p[0][0]] + AUC_cv)
    AUC_cv_lst.append(tmp_out)
    print((f, np.mean(AUC_cv), 10**log10_p[0][0]))

AUC_df = pd.DataFrame(AUC_cv_lst, columns = ['AUC_mean', 'AUC_std', 'AUC_full', 'p_delong'] + ['AUC_' + str(i) for i in range(10)])

AUC_df = pd.concat((pd.DataFrame({'Pro_code':tmp_f}), AUC_df), axis = 1)
AUC_df.to_csv(outfile, index = False)

print('finished')



