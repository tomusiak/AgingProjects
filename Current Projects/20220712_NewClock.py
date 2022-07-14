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
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_squared_error
import sklearn.metrics

def main():
    os.chdir("../../Data/ClockConstruction")
    ml_sample_table = pd.read_csv('healthy_sample_table.csv',index_col=0)
    test_set = pd.read_csv("test_set.csv",index_col=0)
    training_set = pd.read_csv("training_set.csv",index_col=0)
    validation_set = pd.read_csv("validation_set.csv",index_col=0)

    training_set_data = training_set.values
    x_training, y_training = training_set_data[:, :-1], training_set_data[:, -1]
    x_training = sklearn.preprocessing.scale(x_training)

    test_set_data = test_set.values
    x_test, y_test = test_set_data[:, :-1], test_set_data[:, -1]
    x_test = sklearn.preprocessing.scale(x_test)

    validation_set_data = validation_set.values
    x_validation=validation_set_data[:,:]
    x_validation = sklearn.preprocessing.scale(validation_set_data)

    model = ElasticNetCV()
    model.fit(x_training,y_training)
    predicted_training_set_y = model.predict(x_training)
    r2_score_training = sklearn.metrics.r2_score(y_training,predicted_training_set_y)

    predicted_test_set_y = model.predict(x_test)
    r2_score_test = sklearn.metrics.r2_score(predicted_test_set_y,y_test)
    np.savetxt("training_set_predictions.csv",predicted_test_set_y,delimiter=",")

    predicted_validation_set_y = model.predict(x_validation)
    np.savetxt("validation_set_predictions.csv",predicted_validation_set_y,delimeter=",")

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
