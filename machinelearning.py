from numpy import random, array, zeros, empty
import os
import math
import sys
from sklearn.model_selection import KFold
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
#from xgboost import XGBClassifier
from scipy.stats import spearmanr
import joblib
from joblib import dump, load


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn.metrics as metrics

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.model_selection import RandomizedSearchCV
import pandas_profiling
from sklearn import metrics
from sklearn.metrics import RocCurveDisplay
from sklearn.datasets import load_iris
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import PrecisionRecallDisplay
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import GridSearchCV
import xgboost as xgb

import random
from hyperopt import STATUS_OK, Trials, fmin, hp, tpe

input_file = sys.argv[1]
input_folder = sys.argv[2]

print (input_file)
data = pd.read_csv(input_file ,sep="\t",index_col=0)

X = data.drop("class", axis=1)
y = data["class"]

scaler = StandardScaler() 
X_scaled = scaler.fit_transform(X)

logisticRegr = load('./cali/randomforest/firstLogisticRegr.joblib')
mlp = load ('./cali/randomforest/firstMLP.joblib')
second_layer = load('./cali/randomforest/second_layer.joblib')

ntrain = X_scaled.shape[0]
second_input = zeros((ntrain, 2))

predY  = logisticRegr.predict_proba(X_scaled)
second_input[:, 0] = predY[:,1]
predY  = mlp.predict_proba(X_scaled)
second_input[:, 1] = predY[:,1]

Pred_input = second_layer.predict(second_input)
name = list(data.index)
out_file = input_folder + '/filter_final_Variants_mc.vcf'
with open(out_file, 'w') as f:  
    for index, value in enumerate(name):
        elements = value.split("_")
        f.write(elements[1])
        f.write('\t')
        f.write(elements[2])
        f.write('\t')
        f.write('PASS') if Pred_input[index] == 1 else f.write('Low_machine_learning_filter')
        f.write('\n')

