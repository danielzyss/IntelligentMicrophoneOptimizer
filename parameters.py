import numpy as np
import matplotlib.pyplot as plt
import seaborn
from scipy.spatial.distance import euclidean
import uuid
import pickle
from scipy.spatial import Delaunay
import itertools
from sklearn.neighbors import NearestNeighbors
import tqdm
from sklearn.model_selection import train_test_split
import statsmodels.api as sm
import os
import subprocess
import sys
import sklearn.metrics

#CONTROL:
optimize = True #Load DNA or Run new Optimization process
upload_dna = '1010000001000011110100010101100011000' #DNA sequence to upload

#GENETIC PARAMETERS
n_gen = 10 #number of generation
popsize = 15 #size of population

len_dna = 37 #length of DNA sequence

max_row = 11 #maximum number of rows
max_col = 11 #maximum number of columns
min_row = 4 #minimum number of rows
min_col = 4 #maximum number of columns

n_elit = 3 # size of elit
crossover_size = 8 #size of crossover
p_mut1 = 0.2 #probability of Mutation 1
p_mut2 = 0.02 #probability of Mutation 2
max_children_gen_it = 10 # Maximum number of iteration allowed for children generation
n_seg = 2 #number of segment of mutation 1
l_seg = 3 #sizew of segement of mutation 1

score_attempt = 5 #number of iteration for score

#DNA PARAMETERS
gene_cuts = [[0, 4], [4, 8], [8, 10], [10, 12],
             16, [13, 16], [16, 19], [19, 22],
             [22, 25], [25, 28], [28, 31], [31, 34],
             [34, 37]] # disposition of generes in DNA sequence

gene_names = ['row', 'col', 'xlength', 'ylength', 'del', 'overfeed',
              'overin', 'massgen', 'feedgen', 'ingen', 'stiffgen',
              'dampgen', 'sig'] # gene name

maxminval = [[max_row, min_row], [max_col, min_col], 1, 1] #limit value for genes

#RUNNING PARAMETERS
new = True # True: generate a new membrane, False: use membrane called from 'membrane_id'
membrane_id = 'f8528bf6-b02a-4706-a750-4dafd4e45dbc.pkl' # name of membrane to use if new = False
cplusplus = True # Run in fast mode using C++ binaries
graph = False #Plot graph in real time, only when cplusplus is False (2d only)

#LEARNING PARAMETERS:
train_batch_size = 10 #Size of learning batch
test_batch_size = 5 #Size of testing batch
nb_classes = 2 #Number of classes to learn/test from

#NETWORK GENERATION PARAMETERS
net_dim = 2 #dimension of the membrane (set to 2 or 3)
g = 0.0 #gravitational constant
dt = 0.0025 #time-step
quadratic_spring = True #True: quadratic springs, False: Linear springs

#NETWORK ASSESSMENT PARAMETERS
washout_criteria = 0.1 #Criteria for Washout Assessment standard deviation #ADAPT TO MEMBRANE LENGTH
washout_max_permitted_time = 3 #Maximum time allowed for Membrane to Washout




