# -*- coding: utf-8 -*-
"""
Spyder Editor

"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from plot_3d_chart import plot_3d_chart


#Import scikit-learn metrics module for accuracy calculation
from sklearn import metrics
# Import train_test_split function
from sklearn.model_selection import train_test_split
#Import svm model
from sklearn import svm
from sklearn.svm import SVC

data = pd.read_csv("D:/Datasets/pocketDescriptorsML.csv")

X = data.drop('pocketQuality', axis=1)
y = data['pocketQuality']
# Split dataset into training set and test set
X_train, X_test, y_train, y_test = train_test_split(data.data, data.target, test_size=0.3,random_state=109) # 70% training and 30% test

#Create a svm Classifier
clf = svm.SVC(kernel='linear') # Linear Kernel

#Train the model using the training sets
clf.fit(X_train, y_train)
svclassifier = SVC(kernel='rbf') #linear, poly (degree = n), sigmoid
svclassifier.fit(X_train, y_train)

#Predict the response for test dataset
y_pred = clf.predict(X_test)

# Model Accuracy: how often is the classifier correct?
print("Accuracy:",metrics.accuracy_score(y_test, y_pred))

# Model Precision: what percentage of positive tuples are labeled as such?
print("Precision:",metrics.precision_score(y_test, y_pred))

# Model Recall: what percentage of positive tuples are labelled as such?
print("Recall:",metrics.recall_score(y_test, y_pred))


plot_3d_chart.plot()

''' Instructions for 3d Axis plot

# Dataset
df=pd.DataFrame({'X': range(1,101), 'Y': np.random.randn(100)*15+range(1,101), 'Z': (np.random.randn(100)*15+range(1,101))*2 })
 
# plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(df['X'], df['Y'], df['Z'], c='skyblue', s=60)
ax.view_init(30, 185)
plt.show()

'''
