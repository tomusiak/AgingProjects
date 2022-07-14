import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn import linear_model
from tensorflow import keras
from tensorflow.keras import Sequential
from tensorflow.keras.layers import Dense
import copy
import sklearn
from sklearn.linear_model import LinearRegression
import math
import os
from sklearn.model_selection import RepeatedKFold
from sklearn.model_selection import cross_val_score
from numpy import absolute
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import ElasticNetCV

def main():
    os.chdir("../../Data/ClockConstruction")
    os.chdir("Data/ClockConstruction")
    ml_sample_table = pd.read_csv('healthy_sample_table.csv',index_col=0)
    test_set = pd.read_csv("test_set.csv",index_col=0)
    training_set = pd.read_csv("training_set.csv",index_col=0)

    cv = RepeatedKFold(n_splits=10, n_repeats=3, random_state=1)
    # define model
    ratios = np.arange(0, 1, 0.1)
    alphas = [ 1e-4, 1e-3, 1e-2, 1e-1, 0.0, 1.0, 10.0]
    model = ElasticNetCV(l1_ratio=ratios, alphas=alphas, cv=cv, n_jobs=-1)
    training_set_data = training_set.values
    X, y = training_set_data[:, :-1], training_set_data[:, -1]
    X = sklearn.preprocessing.scale(X)
    scores = cross_val_score(model, X, y, scoring='neg_mean_absolute_error', cv=cv, n_jobs=-1)
    scores = absolute(scores)
    test_set_data = test_set.values
    X_test, y_test = test_set_data[:, :-1], test_set_data[:, -1]
    model.fit(X,y)
    yhat = model.predict(X_test)
    print("nothing")

def regressionNNPrediction(x_train, y_train, x_test, y_test):
    prediction = regressionNN(x_train, y_train, x_test)
    x_test = pd.DataFrame(x_test)
    y_test = pd.DataFrame(y_test)
    prediction = pd.DataFrame(prediction)
    exportPredictions(x_test,y_test,prediction,"regression_neuralnetwork")
    y_test['prediction'] = prediction
    return y_test

def regressionNN(x_train, y_train, x_test):
    feature_count = len(x_train.columns)
    row_count = len(x_train.index)
    model = Sequential()
    model.add(Dense(row_count, input_dim=feature_count, activation= "relu"))
    model.add(Dense(round(row_count/4), activation= "relu"))
    model.add(Dense(1))
    model.compile(loss= "mae" , optimizer="adam", metrics=["mae"])
    model.fit(x_train, y_train, epochs=5000, verbose=0)
    pred = model.predict(x_test)
    return pred

def cleanData(data):
    list_of_columns = list(data.columns.values)
    list_of_columns.remove("age")
    list_of_columns.append("age")
    data = data[list_of_columns]
    y = data.iloc[:,-1]
    x = data.iloc[:, :-1]
    scaler = preprocessing.StandardScaler().fit(x)
    x = pd.DataFrame(scaler.fit_transform(x.values), columns=x.columns, index=x.index)
    return scaler,x, y

main()
