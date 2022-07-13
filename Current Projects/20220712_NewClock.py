import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import preprocessing
from tensorflow import keras
from tensorflow.keras import Sequential
from tensorflow.keras.layers import Dense
import copy
from sklearn.linear_model import LinearRegression
import math

def main():
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
