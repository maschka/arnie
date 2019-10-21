import numpy as np
import argparse, sys
from arnie.structure_prediction.mea_utils import *

class MEA:
    def __init__(self, bpps, gamma = 1.0, debug=False):
        self.debug = debug
        self.bpps = bpps
        self.N=self.bpps.shape[0]
        self.gamma = gamma
        self.W = np.zeros([self.N,self.N])
        self.MEA_bp_list = []
        self.structure = ['.']*self.N
        self.MEA_bp_matrix = np.zeros([self.N,self.N])
        self.tb = np.zeros([self.N,self.N])
        self.min_hp_length=3
        self.evaluated = False

        self.run_MEA()
        
    def fill_W(self, i, j):
        options = [self.W[i+1, j], self.W[i, j-1], (self.gamma+1)*self.bpps[i,j] + self.W[i+1, j-1] - 1,\
                np.max([self.W[i,k] + self.W[k+1, j] for k in range(i+1,j)])]
        self.W[i,j] = np.max(options) 
        self.tb[i,j] = np.argmax(options) #0: 5' pass, 1: 3' pass, 2: bp, 3: multiloop
        
    def run_MEA(self):
        # fill weight matrix
        for length in range(self.min_hp_length, self.N):
            for i in range(self.N-length):
                j = i + length
                self.fill_W(i,j)
                
        self.traceback(0,self.N-1)
        
        for x in self.MEA_bp_list:
            self.MEA_bp_matrix[x[0],x[1]]=1
            self.structure[x[0]]='('
            self.structure[x[1]]=')'
        
        self.structure = ''.join(self.structure)
        if not self.evaluated: self.evaluated = True
        
    def traceback(self, i, j):
        if j <= i:
            return
        elif self.tb[i,j] == 0: #5' neighbor
            if self.debug: print(i,j, "5'")
            self.traceback(i+1,j)
        elif self.tb[i,j] == 1: #3' neighbor
            if self.debug: print(i,j, "3'")
            self.traceback(i,j-1)
        elif self.tb[i,j] == 2: # base pair
            if self.debug: print(i,j,'bp')
            self.MEA_bp_list.append((i,j))
            self.traceback(i+1,j-1)
        else: #multiloop
            for k in range(i+1,j):
                if self.W[i,j] == self.W[i, k] + self.W[k+1,j]:
                    if self.debug: print(i,j,"multiloop, k=",k)
                    self.traceback(i,k)
                    self.traceback(k+1,j)
                    break

    def score_expected(self):
        '''Compute expected values of TP, FP, etc from predicted MEA structure.

         Returns: 
         pseudoexpected SEN, PPV, MCC, F-score'''

        if not self.evaluated: self.run_MEA()

        pred_m = self.MEA_bp_matrix[np.triu_indices(self.N)]
        probs = self.bpps[np.triu_indices(self.N)]

        TP = np.sum(np.multiply(pred_m, probs)) + 1e-6
        TN = 0.5*self.N*self.N-1 - np.sum(pred_m) - np.sum(probs) + TP + 1e-6
        FP = np.sum(np.multiply(pred_m, 1-probs)) + 1e-6
        FN = np.sum(np.multiply(1-pred_m, probs)) + 1e-6

        a,b = np.triu_indices(self.N)
        cFP = 1e-6
        for i in range(len(pred_m)):
            if np.sum(self.MEA_bp_matrix,axis=0)[a[i]] + np.sum(self.MEA_bp_matrix,axis=0)[b[i]]==0:
               cFP += np.multiply(pred_m[i], 1-probs[i])

        sen = TP/(TP + FN)
        ppv = TP/(TP + FP - cFP)
        mcc = (TP*TN - (FP - cFP)*FN)/np.sqrt((TP + FP - cFP)*(TP + FN)*(TN + FP - cFP)*(TN + FN))
        fscore = 2*TP/(2*TP + FP - cFP + FN)

        return [sen, ppv, mcc, fscore]
