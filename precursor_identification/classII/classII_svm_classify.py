#! /usr/bin/env python
#
# Copyright (C) 2015 Jonathan I. Tietz
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2015 Christopher J. Schwalen
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2015 Parth S. Patel
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2015 Douglas A. Mitchell
# University of Illinois
# Department of Chemistry
#
# License: GNU Affero General Public License v3 or later
# Complete license availabel in the accompanying LICENSE.txt.
# or <http://www.gnu.org/licenses/>.
#
# This file is part of RODEO.
#
# RODEO is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RODEO is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
'''
   SVM classification script
   
   Required input:
     a training file CSV
     a set of model data CSV (if fitting)
     
   Output:
     a CSV list of identifiers and classifications
     
   Note that all options must be hard-coded (the script does not take command-line arguments).
   
   This script assumes you have already optimized kernel parameters. Use the included optimization script if not.
   
   RECOMMENDATION:
     Input CSV should ideally have its primary key as column 0, its classification as column 1, and feature as columns [2,...,end]
     For the fitting CSV, there will be no classification; leave it blank if you want (it'll be ignored upon import)
'''

import numpy as np
import csv

from sklearn import svm
from sklearn import preprocessing
from sklearn import model_selection
from sklearn import metrics

# CONFIGURATION OPTIONS
''' change these as desired '''

input_training_file = 'classII_training_set_matrix_20190426.csv'         # the CSV containing the training set
input_fitting_file = 'classII-peps-params-1640-feats.csv'          # the CSV containing the data to be fitted
output_filename = 'classII-peps-classified.csv'     # output filename; this will be a CSV with the first column being the primary key and the second being the classification

primary_key_column = 0;            # the column of the CSV that contains the primary key (identifier) for each record
classification_column = 1;         # the column of the CSV that contains the classification for each record
csv_has_header = True;             # set to true if CSVs have header rows to ignore upon import

# default kernel values
kernel_option = 'rbf'

C_option = 8207.89938876178
gamma_option = 1E-6
class_weight_option = 'balanced'

def parse_CSV_to_dataset(csv_filename, dataset_type):
    '''Parse an input CSV into a data set
    
       Inputs:
            csv_filename            name of CSV file to be parsed
            dataset_type            either 'training' or 'fitting'
    '''
    dataset = []
    with open(csv_filename, 'rb') as csvfile:
      csv_read = csv.reader(csvfile, delimiter=',', quotechar='"')
      if csv_has_header == True:
        next(csv_read,None)
        
      if classification_column < primary_key_column:
        class_first = 1
        pk_first = 0
      else:
        pk_first = 1
        class_first = 0 
        
      for row in csv_read:
        temp_entry = []
        temp_entry.append(row.pop(primary_key_column))
        temp_entry.append(row.pop(classification_column - pk_first))
        for c in row:
          temp_entry.append(float(c))
        # remove all unclassified features if training
        if dataset_type == 'training':
          if int(temp_entry[1]) == 1 or int(temp_entry[1]) == 0:
            temp_entry[1] = int(temp_entry[1])
            dataset.append(temp_entry)
        if dataset_type == 'fitting':
          dataset.append(temp_entry)
    return dataset
    
def write_to_csv(list_of_primary_keys, list_of_classifications, output_file):
    with open(output_filename, 'wb') as csvfile:
      csv_write = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
      for i in range(len(list_of_primary_keys)):
        temp_row = [list_of_primary_keys[i],list_of_classifications[i]]
        csv_write.writerow(temp_row)
    return
    
def main():

    # parse data
    
    print "Importing training and fitting data ..."
    training_data_unrefined = parse_CSV_to_dataset(input_training_file, 'training')
    fitting_data_unrefined = parse_CSV_to_dataset(input_fitting_file, 'fitting')
    
    primary_key_list = []
    fitting_data_just_features = []
    training_data_just_features = []
    training_data_classifications = []
   
    for entry in training_data_unrefined:
        training_data_classifications.append(entry.pop(1))
        entry.pop(0)
        training_data_just_features.append(entry)  
    for entry in fitting_data_unrefined:
        primary_key_list.append(entry.pop(0))
        entry.pop(0)
        fitting_data_just_features.append(entry)
    print "  (done)"
    print "Initiating learning and fitting"
    
    # Scaling -- this ensures standardization of model and target data
    training_data_refined = preprocessing.scale(training_data_just_features)
    scaler = preprocessing.StandardScaler().fit(training_data_refined)
    training_data_refined = scaler.transform(training_data_just_features)
    fitting_data_refined = scaler.transform(fitting_data_just_features)

    # This creates the classifier and learns using the model data
    clf = svm.SVC(kernel=kernel_option,class_weight=class_weight_option,C=C_option,gamma=gamma_option)
    clf.fit(training_data_refined, training_data_classifications)

    # This classifies the input data as a list
    classification_list = clf.predict(fitting_data_refined)
    
    # Output results to file
    write_to_csv(primary_key_list, classification_list, output_filename)
    print " ... Done"
    
    return

if __name__ == '__main__':
  main()
