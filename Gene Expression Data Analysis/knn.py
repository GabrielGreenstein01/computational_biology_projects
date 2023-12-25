import sys

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt


"""
KNN class computes k-nearest neighbors for a given dataset using leave-one-out cross validation
"""

class KNN(object):

	def __init__(self):

		self.expression_file = ""
		self.sample_file = ""

		self.X_df = pd.DataFrame()
		self.y_df = pd.DataFrame()

		self.X = []
		self.y = []

		self.X_train = []
		self.y_train = []

		self.X_test = []
		self.y_test = []

		self.k = 5
		self.fn = 0.5

	def load_data(self, expfile, sampfile):
		"""
		Takes the paths to the expression and sample file, reads them in, and stores within KNN class.
		"""

		self.expression_file = expfile
		self.sample_file = sampfile

		# Read in X and y
		self.X_df = pd.read_csv(expfile, sep='\t')
		self.y_df = pd.read_csv(sampfile, sep='\t', header = None)

		# Drop gene symbol column
		self.X_df = self.X_df.drop(columns=['SYMBOL'])

		# For each column in X, find corresponding classification value in y
		# Represent each column as one data point and to append X
		# Get corresponding classification value from y and append to y 
		for col in self.X_df.columns:
			self.X.append(list(self.X_df[col]))
			self.y.append(self.y_df.loc[self.y_df.iloc[:,0] == col][1].values[0])

	def split_data(self,train,test):
		"""
		Train and test are arrays containing the indices of data points (i.e. columns of X) for which we want to train and test our model on
		"""

		# Reassign X,y train and test arrays to empty
		self.X_train = []
		self.y_train = []

		self.X_test = []
		self.y_test = []
		
		# Split data into train set
		for index in train:
			self.X_train.append(self.X[index])
			self.y_train.append(self.y[index])

		# Split data into test set
		for index in test:
			self.X_test.append(self.X[index])
			self.y_test.append(self.y[index])

	def get_distance(self, pointA, pointB):
		"""
		Takes in two points and returns the Euclidean distance between them
		"""

		# Make sure points have the same dimensions
		assert len(pointA) == len(pointB)

		dim = len(pointA)
		dist = 0

		# Compute Euclidean distance between the points
		for i in range(0,dim):
			dist += (pointA[i] - pointB[i])**2

		dist = dist**(1/2)

		return dist

	def get_neighbors(self, point):
		"""
		Take in a point and get the k-nearest neighbors
		Note that in implementation, the input "point" will come from the test data
		and the points from which we look for the k-nearest neighbors will come from the train data
		"""
		distances = [] 

		# Get the distance between the input "point" and each point in train
		for i in range(0,len(self.X_train)):

			dist = self.get_distance(point,self.X_train[i])

			distances.append(dist)

		# Sort the distances from smallest to largest and store the index corresponding to that distance
		indices = np.argsort(distances)

		# Return top k indices
		return indices[:self.k]

	def predict_test(self, k, fn):
		"""
		Returns the classification of the test data (trained on the train data) 
		"""

		y_predictions = []

		# Predict each point in test dataset
		for i in range(0,len(self.X_test)):

			# Get the i-th point in test
			x = self.X_test[i]

			# Get its k-nearest neighbors
			x_neighbors_indices = self.get_neighbors(x)

			# Get the classifications of each of the k-nearest neighbors
			self.y_train = np.array(self.y_train)
			neighbors_y_val = self.y_train[x_neighbors_indices]

			# Average the classifications of each of the k-nearest neighbors
			neighbors_y_avg = np.mean(neighbors_y_val)

			# Predict healthy or sick based on fn value
			if neighbors_y_avg > self.fn:
				y_predictions.append(1)
			else:
				y_predictions.append(0)

		return y_predictions

	def generate_LOOCV(self):
		"""
		Generate test and train data per the leave-one-out cross validation
		"""

		LOOCV = []

		# Generate test and train data
		for i in range(0,len(self.X)):

			test_indices = []
			train_indices = []

			for j in range(0,len(self.X)):

				if i == j:
					test_indices.append(j)
				else:
					train_indices.append(j)

			LOOCV.append([train_indices,test_indices])

		return LOOCV

	def get_assignments(self, k, fn):
		"""
		Returns the class assignments for all samples for given values of and as a list of integer 0s and 1s.
		"""

		# Update k, fn values
		self.k = k
		self.fn = fn

		# Assign train,test by generate_LOOCV
		LOOCV = self.generate_LOOCV()

		LOOCV_predictions = []

		# Run prediction on test data
		for train,test in LOOCV:
			self.split_data(train, test)
			LOOCV_predictions.append(self.predict_test(self.k,self.fn)[0])

		predictions_y = []

		# Reorder predictions in terms of rows of y
		for col in list(self.y_df[0]):
			index = self.X_df.columns.get_loc(col)
			predictions_y.append(LOOCV_predictions[index])

		LOOCV_predictions = predictions_y

		return LOOCV_predictions

	def calc_metrics(self, k, fn):
		"""
		Returns a list of float values [sensitivity,specificity] of a KNN classifier using the given values of k and fn

		Sick = 1, Healthy = 0
		True positive: Sick people correctly identified as sick
		False positive: Healthy people incorrectly identified as sick
		True negative: Healthy people correctly identified as healthy
		False negative: Sick people incorrectly identified as healthy

		Sensitivity = TP / (TP + FN)
		Specificity = TN / (TN + FP)
		"""

		# Get classifications from LOOCV
		predictions = self.get_assignments(k,fn)

		predictions_X = []

		# Reorder predictions in the order of the columns of X
		for col in self.X_df.columns:
			index = self.y_df.loc[self.y_df.iloc[:,0] == col].index.values[0]
			predictions_X.append(predictions[index])

		# Update predictions variable
		predictions = predictions_X

		TP = 0
		FP = 0
		TN = 0
		FN = 0

		# Calculate TP, FP, TN, FN
		for i in range(0,len(predictions)):
			if predictions[i] == self.y[i] and self.y[i] == 1:
				TP += 1
			elif predictions[i] == self.y[i] and predictions[i] == 0:
				TN += 1
			elif predictions[i] != self.y[i] and predictions[i] == 0:
				FN += 1
			elif predictions[i] != self.y[i] and predictions[i] == 1:
				FP += 1

		# Calculate specificity and sensitivity
		sensitivity = TP/(TP + FN)
		specificity = TN/(TN + FP)

		return [sensitivity, specificity]


	def plot_ROC(self):

		sensitivities = []
		specificities = []

		k = 3

		fn_vals = [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1]

		for i in fn_vals:
			sensitivity, specificity = self.calc_metrics(k,i)
			sensitivities.append(sensitivity)
			specificities.append(specificity)

		for i in range(0,len(specificities)):
			specificities[i] = 1 - specificities[i]

		print(sensitivities)
		print(specificities)

		plt.style.use('ggplot')
		plt.rcParams["figure.figsize"] = (16,9)

		TPR  = sensitivities
		FPR = specificities
		#you're integrating from right to left. This flips the sign of the result
		auc = -1 * np.trapz(TPR, FPR)

		plt.plot(FPR, TPR, linestyle='--', marker='o', color='darkorange', lw = 2, label='ROC curve', clip_on=False)
		plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
		plt.xlim([0.0, 1.0])
		plt.ylim([0.0, 1.0])
		plt.xlabel('1 - Specificity')
		plt.ylabel('sensitivity')
		plt.title('ROC curve, AUC = %.2f'%auc)
		plt.legend(loc="lower right")
		plt.savefig('AUC_example.png')
		plt.show()



def main():
    """
    Function is called when the file is called in the terminal
    """
    # check that the file is being properly used
    if (len(sys.argv) !=3):
        print("Please specify an input file and an output file as args.")
        return
        
    # input variables
    expression_file = sys.argv[1]
    sample_file = sys.argv[2]

    knn = KNN()
    knn.load_data(expression_file,sample_file)
    knn.plot_ROC()
    # for k in range(1,7):
    # 	print(k)
    # 	knn.calc_metrics(k,0.5)
		# print(knn.calc_metrics(k,0.5))
    # print(knn.generate_LOOCV())
    # print(knn.get_assignments(5,0.7))


if __name__=="__main__":
    main()