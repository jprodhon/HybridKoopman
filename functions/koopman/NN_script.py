import numpy as np
from random import sample
from math import floor
import tensorflow as tf
import sys

from koopmanlib.dictionary import PsiNN
from koopmanlib.solver import KoopmanDLSolver

N = int(sys.argv[1])
n_funct = int(sys.argv[2])
n_lay = 3
n_neur_perlay = 100
layers = [n_neur_perlay]*n_lay

def extract_data(N):
    """ Data extraction from txt files calculated with Matlab """
    Z = np.loadtxt('functions/dataset/saves/Z.txt',delimiter=',')
    Zt = np.loadtxt('functions/dataset/saves/Zt.txt',delimiter=',')
    eval_raw = np.loadtxt('functions/dataset/saves/eval.txt',delimiter=',')
    m = Z.shape[0]
    n_eval = int(eval_raw.shape[0]/m)
    eval_set = np.zeros((n_eval,Z.shape[0],eval_raw.shape[1]))
    for i in range(n_eval):
        eval_set[i,:,:] = eval_raw[i*m:(i+1)*m,:]
    Zs = np.zeros((1+N,Z.shape[0],Z.shape[1]))
    Zts = np.zeros((1+N,Z.shape[0],Z.shape[1]))
    eval_sets = np.zeros((1+N,n_eval,Z.shape[0],eval_raw.shape[1]))
    Zs[0,:,:] = Z
    Zts[0,:,:] = Zt
    eval_sets[0,:,:,:] = eval_set
    for k in range(N):
        Z = np.loadtxt('functions/dataset/saves/Z_red'+str(k+1)+'.txt',delimiter=',')
        Zt = np.loadtxt('functions/dataset/saves/Zt_red'+str(k+1)+'.txt',delimiter=',')
        eval_raw = np.loadtxt('functions/dataset/saves/eval_red'+str(k+1)+'.txt',delimiter=',')
        eval_set = np.zeros((n_eval,Z.shape[0],eval_raw.shape[1]))
        for i in range(n_eval):
            eval_set[i,:,:] = eval_raw[i*m:(i+1)*m,:]
        Zs[1+k,:,:] = Z
        Zts[1+k,:,:] = Zt
        eval_sets[1+k,:,:,:] = eval_set    
    return Zs, Zts, eval_sets

def export_data(solver,Z,Zt,eval_set,k):
    """Export the data for Matlab"""
    NN_obs = np.transpose(solver.dic.dicNN(np.transpose(Z)).numpy())
    NN_obs_t = np.transpose(solver.dic.dicNN(np.transpose(Zt)).numpy())
    eval_NN = np.zeros((n_funct,n_eval))
    for i in range(n_eval):
        traj = np.squeeze(eval_set[i,:,0]).reshape((1,m))
        psi0 = np.transpose(solver.dic.dicNN(traj).numpy())
        eval_NN[:,i] = np.ravel(psi0)
    np.savetxt('functions/dataset/saves/NN_obs'+str(k)+'.txt',NN_obs,delimiter=',')
    np.savetxt('functions/dataset/saves/NN_obs_t'+str(k)+'.txt',NN_obs_t,delimiter=',')
    np.savetxt('functions/dataset/saves/eval_NN'+str(k)+'.txt',eval_NN,delimiter=',')
    
def model(layers,n_funct,data_train,data_valid,target_dim,batch_size):
    """Generates and train the neural network for EDMD-DL"""
    basis_function = PsiNN(layer_sizes = layers, n_psi_train = n_funct)
    solver = KoopmanDLSolver(dic = basis_function,
                             target_dim = target_dim,
                             reg = 0.1)
    solver.build(data_train = data_train,
                 data_valid = data_valid,
                 epochs = 1000,
                 batch_size = batch_size,
                 lr = 1e-4,
                 log_interval = 20,
                 lr_decay_factor = 0.8)
    return solver

Zs, Zts, eval_sets = extract_data(N)
a, n_eval, m, nt = eval_sets.shape
for k in range(1+N):
    Z = np.squeeze(Zs[k,:,:])
    Zt = np.squeeze(Zts[k,:,:])
    ind = floor(0.8*Z.shape[1])
    eval_set = np.squeeze(eval_sets[k,:,:,:])
    data_train = [tf.transpose(Z[:,:ind]), tf.transpose(Zt[:,:ind])]
    data_valid = [tf.transpose(Z[:,ind:]), tf.transpose(Zt[:,ind:])]
    solver = model(layers,n_funct,data_train,data_valid,m,Z.shape[1])
    export_data(solver,Z,Zt,eval_set,k)