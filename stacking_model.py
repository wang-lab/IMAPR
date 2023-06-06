import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn.metrics as metrics
import sys

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
from sklearn.svm import SVC
import random
from hyperopt import STATUS_OK, Trials, fmin, hp, tpe
from scipy.stats import uniform
import joblib
from joblib import dump, load
from sklearn.model_selection import StratifiedKFold, cross_val_score
pd.options.mode.chained_assignment = None
from sklearn.metrics import confusion_matrix
from numpy import random, array, zeros, empty
from sklearn.preprocessing import PolynomialFeatures

inputFile = sys.argv[1]
outputFile = sys.argv[2]

data = pd.read_csv(inputFile,sep="\t",index_col=0)
columns_to_normalize = ['tumorDepth', 'alterTumorReads', 'mapQuality', 'mapPosition', 'artifactsNormal', 'genoNormal', 'populationAF', 'tumorLog', 'SQ_alt', 'SQ_ref', 'mpileup_alt']
normalizedtrainInput = data
scaler = StandardScaler()
scaler.fit(normalizedtrainInput[columns_to_normalize])
normalized_data  = scaler.transform(normalizedtrainInput[columns_to_normalize])
normalizedtrainInput.loc[:, columns_to_normalize] = normalized_data

Input = pd.read_csv(inputFile,sep="\t",index_col=0)

best_random = joblib.load('./lib/randomForest.joblib')
best_xgb = joblib.load('./lib/XGboost.joblib')
best_mlp = joblib.load('./lib/mlp.joblib')
second = joblib.load('./lib/second.joblib')

nInput = Input.shape[0]
secondInput = zeros((nInput, 3))

predInput  = best_random.predict_proba(Input)
secondInput[:, 0] = predInput[:,1]
predInput  = best_xgb.predict_proba(Input)
secondInput[:, 1] = predInput[:,1]
predInput  = best_mlp.predict_proba(normalizedtrainInput)
secondInput[:, 2] = predInput[:,1]


poly = PolynomialFeatures(degree=3, include_bias=False)
polySecondInput = poly.fit_transform(secondInput)
predictLabel = second.predict(polySecondInput)

nameList = list(Input.index)
index = 0
with open(outputFile, "w") as f:
    for variantName in nameList:
        parts = variantName.split("_")
        string = '\t'.join(parts) + '\t' + str(predictLabel[index])
        print(string, file=f)
        index += 1