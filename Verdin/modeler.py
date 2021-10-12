import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import preprocessing
from keras.models import Sequential
from keras.layers import Dense
import copy

directory = "data/"
data_file = "FinalCombined.csv"
iterations = 150;
alpha = 0.01;

def addIntercept(x):
    row_count = len(x.index)
    initialized_rows = [1] * row_count
    x.insert(0, "Intercept", initialized_rows, True)
    x = pd.DataFrame(x)
    return x

def dropFeatures(data, featureList):
    data = data.drop(columns=featureList)
    return data

def cleanData(data):
    list_of_columns = list(data.columns.values)
    list_of_columns.remove("Predictions Diff")
    list_of_columns.append("Predictions Diff")
    data = data[list_of_columns]
    y = data.iloc[:,-1]
    x = data.iloc[:, :-1]
    scaler = preprocessing.StandardScaler().fit(x)
    x = pd.DataFrame(scaler.fit_transform(x.values), columns=x.columns, index=x.index)
    return scaler,x, y

def splitData(x,y,proportion):
    row_count = len(x.index)
    data = x
    data['y'] = y
    shuffledData= data.sample(frac=1)
    split = round(row_count*proportion)
    testData = shuffledData.head(split)
    trainingData = shuffledData.tail(row_count-split)
    y_train = pd.DataFrame(trainingData.iloc[:,-1])
    x_train = pd.DataFrame(trainingData.iloc[:, :-1])
    y_test = pd.DataFrame(testData.iloc[:,-1])
    x_test = pd.DataFrame(testData.iloc[:,:-1])
    y_train = y_train.reset_index(drop=True)
    y_test = y_test.reset_index(drop=True)
    x_test = x_test.reset_index(drop=True)
    x_train = x_train.reset_index(drop=True)
    return x_train, y_train, x_test, y_test

def main():
    data = pd.read_csv(directory + data_file)
    scaler, x, y = cleanData(data)
    x_train, y_train, x_test, y_test = splitData(x,y,proportion = .20)
    gradientDescentPrediction(x_train,y_train,x_test,y_test)
    print("Gradient Descent Done")
    normalEquationPrediction(x_train,y_train,x_test,y_test)
    print("Normal Equation Done")
    regressionNNPrediction(x_train,y_train, x_test, y_test)
    print("Regression NN done")
    cat_y_train = categoricalTransform(y_train, cutoff=.1)
    cat_y_test = categoricalTransform(y_test, cutoff=.1)
    categoricalNNPrediction(x_train,cat_y_train,x_test,cat_y_test)
    print("Categorical NN done")
    
def gradientDescentPrediction(x_train,y_train, x_test, y_test):
    gradient_x_train = copy.deepcopy(x_train)
    gradient_y_train = copy.deepcopy(y_train)
    gradient_x_test = copy.deepcopy(x_test)
    gradient_y_test = copy.deepcopy(y_test)
    gradient_x_train = addIntercept(gradient_x_train)
    gradient_x_test = addIntercept(gradient_x_test)
    gradient_theta = gradientDescent(gradient_x_train,gradient_y_train,alpha,iterations)
    gradient_prediction = np.dot(gradient_x_test,np.transpose(gradient_theta))
    gradient_prediction = pd.DataFrame(gradient_prediction)
    exportPredictions(gradient_x_test, gradient_y_test, gradient_prediction, 'gradient')

def gradientDescent(x, y, alpha, iterations):
    gradient_x = copy.deepcopy(x)
    gradient_y = copy.deepcopy(y)
    gradient_theta = np.zeros((1,len(gradient_x.columns)))
    gradient_m = len(gradient_y)
    gradient_cost_history = np.zeros(iterations)
    gradient_y = np.matrix(gradient_y)
    for i in range (iterations):
        gradient_predictions = np.dot(gradient_x,np.transpose(gradient_theta))
        gradient_theta = gradient_theta - (1/gradient_m) * alpha * (np.transpose((gradient_predictions-gradient_y)).dot(gradient_x))
        gradient_cost_history[i] = computeCost(gradient_predictions, gradient_y, gradient_m)
    gradient_dataframe = pd.DataFrame(gradient_cost_history)
    #gradient_dataframe.plot()
    #plt.show()
    gradient_theta = pd.DataFrame(gradient_theta)
    return gradient_theta

def computeCost(predictions, y, m):
    cost = (1/(2*m)) * np.sum(np.square(predictions-y))
    return cost

def normalEquationPrediction(x_train, y_train, x_test, y_test):
    normal_x_train = copy.deepcopy(x_train)
    normal_y_train = copy.deepcopy(y_train)
    normal_x_test = copy.deepcopy(x_test)
    normal_y_test = copy.deepcopy(y_test)
    normal_x_train = addIntercept(normal_x_train)
    normal_x_test = addIntercept(normal_x_test)
    normal_theta = np.zeros((1,len(normal_x_train.columns)))
    normal_theta = normalEqn(normal_x_train,normal_y_train)
    normal_prediction = np.dot(normal_x_test,normal_theta)
    exportPredictions(normal_x_test,normal_y_test,normal_prediction,"normaleq")
    
def normalEqn(x,y):
    x_transpose = np.transpose(x)
    x_transpose_dot_x = x_transpose.dot(x)
    temp_1 = np.linalg.pinv(x_transpose_dot_x)
    temp_2= x_transpose.dot(y)
    theta = temp_1.dot(temp_2)
    return theta

def regressionNNPrediction(x_train, y_train, x_test, y_test):
    prediction = regressionNN(x_train, y_train, x_test)
    exportPredictions(x_test,y_test,prediction,"regression_neuralnetwork")

def regressionNN(x_train, y_train, x_test):
    feature_count = len(x_train.columns)
    row_count = len(x_train.index)
    model = Sequential()
    model.add(Dense(row_count, input_dim=feature_count, activation= "relu"))
    model.add(Dense(round(row_count/4), activation= "relu"))
    model.add(Dense(1))
    model.compile(loss= "mae" , optimizer="adam", metrics=["mae"])
    model.fit(x_train, y_train, epochs=500, verbose=1)
    pred = model.predict(x_test)
    return pred

def categoricalNNPrediction(x_train,y_train,x_test,y_test):
    prediction = categoricalNN(x_train,y_train,x_test)
    exportPredictions(x_test,y_test,prediction,"categorical_neuralnetwork")
    
def categoricalNN(x_train, y_train, x_test):
    feature_count = len(x_train.columns)
    row_count = len(x_train.index)
    model = Sequential()
    model.add(Dense(row_count, input_dim=feature_count, activation="relu"))
    model.add(Dense(round(row_count/4), activation= "relu"))
    model.add(Dense(1, activation = "sigmoid"))
    model.compile(loss= "binary_crossentropy" , optimizer="adam", metrics=["accuracy"])
    model.fit(x_train, y_train, epochs=200, verbose = 0)
    pred = model.predict(x_test)
    pred = (pred > .80).astype(int)
    return pred
    
def categoricalTransform(y, cutoff):
    newY = (y > .10).astype(int)
    return newY

def exportPredictions(x_test, y_test, prediction, prediction_type):
    data = copy.deepcopy(x_test)
    data['y'] = y_test
    data['predicted'] = prediction
    data.to_csv(directory + prediction_type + "Prediction.csv")
        
main()