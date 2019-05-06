
from robustness import *
from struct_formula import *
from load_data import *
#from bayesian_nn import *
import matplotlib.pyplot as plt
from scipy.io import loadmat
import edward as ed
import numpy as np
import tensorflow as tf
from sklearn.model_selection import train_test_split
from edward.models import Normal
from sklearn.gaussian_process.kernels import RBF
from numpy.linalg import inv, pinv

def generate_train(data_name, Num):
    signal, time2, name, label = load_data(data_name)
    label = np.squeeze(label)
    X =[]
    y = []
    for _ in range(Num):
        vector, formula = Init_state(name, signal, time2, 3)
        tree = formula.get_tree()
        rewards = reward(tree, name, signal,time2 )
        rob = [a*b for a, b in zip(rewards, label)]
        X.append(vector)
        y.append(min(rob))
    header = "X_train, y_train \n"
    data = np.column_stack((X, y))
    dat_name = 'nn_'+data_name[0:-3]+'dat'
    np.savetxt(dat_name, data, header=header)
    return X , y

def generate_vector(data_name, Num):
    signal, time2, name, label = load_data(data_name)
    X =[]
    for _ in range(Num):
        vector, formula = Init_state(name,signal, time2,3)
        X.append(vector)
    header = "X_sample\n"
    dat_name = 'Xs_'+data_name[0:-3]+'dat'
    np.savetxt(dat_name, X, header=header)
    return X


def robust_function(vector,data_name, pool):
    signal, time1, name, label = load_data(data_name)
    state, formulas = Init_state(name, signal, time1,3)
    tree = formulas.vector_tree(vector)
    if pool:
        rewards = poolreward(tree,name,signal,time1)
    else:
        rewards= reward(tree,name,signal,time1)

    robust = [a*b for a, b in zip(rewards, label)]
    return min(robust)

def find_a_candidate(X_train, Xs,ys, kappa):
    var = calculate_variance(X_train,Xs)
    ucb = [m+kappa*v for m, v in zip(ys,var)]
    idx = np.argsort(ucb)
    x  = Xs[idx[-10:]].tolist()
    return x

def is_invertible(a):
    return a.shape[0] == a.shape[1] and np.linalg.matrix_rank(a) == a.shape[0]
    
def calculate_variance(X,Xs):
    kernel = RBF()
    kernel.__init__()
    sigma = 0.01
    ks = kernel.diag(Xs)
    Kt = kernel.__call__(X,X)+sigma**2
    if is_invertible(Kt):
        Kt = inv(Kt)
    else:
        Kt = pinv(Kt)
    K = kernel.__call__(Xs,X)
    Ktt = np.matmul(K, Kt)
    Ktt = np.matmul(Ktt,K.transpose())
    Ktt = Ktt.diagonal()
    var = ks - Ktt.transpose()
    return var


def del_all_flags(FLAGS):
    flags_dict = FLAGS._flags()
    keys_list = [keys for keys in flags_dict]
    for keys in keys_list:
        FLAGS.__delattr__(keys)



def bayesian_optimization(X_train,y_train,Xs):
    del_all_flags(tf.flags.FLAGS)
    N_sample = len(X_train)
    N_feature = len(X_train[0])
    tf.flags.DEFINE_integer("N", default=N_sample, help="Number of data points.")
    tf.flags.DEFINE_integer("D", default=N_feature, help="Number of features.")
    tf.flags.DEFINE_integer("H", default=50, help="Number of hidden node.")
    FLAGS = tf.flags.FLAGS

    def neural_network(X):
        h = tf.tanh(tf.matmul(X, W_0) + b_0)
        h = tf.tanh(tf.matmul(h, W_1) + b_1)
        h = tf.matmul(h, W_2) + b_2
        return tf.reshape(h, [-1])


  # MODEL
    with tf.name_scope("model"):
        W_0 = Normal(loc=tf.zeros([FLAGS.D, FLAGS.H]), scale=tf.ones([FLAGS.D, FLAGS.H]),
                 name="W_0")
        W_1 = Normal(loc=tf.zeros([FLAGS.H, FLAGS.H]), scale=tf.ones([FLAGS.H, FLAGS.H]), name="W_1")
        W_2 = Normal(loc=tf.zeros([FLAGS.H, 1]), scale=tf.ones([FLAGS.H, 1]), name="W_2")
        b_0 = Normal(loc=tf.zeros(FLAGS.H), scale=tf.ones(FLAGS.H), name="b_0")
        b_1 = Normal(loc=tf.zeros(FLAGS.H), scale=tf.ones(FLAGS.H), name="b_1")
        b_2 = Normal(loc=tf.zeros(1), scale=tf.ones(1), name="b_2")

    X = tf.placeholder(tf.float32, [FLAGS.N, FLAGS.D], name="X")
    y = Normal(loc=neural_network(X), scale=0.1 * tf.ones(FLAGS.N), name="y")

  # INFERENCE
    with tf.variable_scope("posterior", reuse=tf.AUTO_REUSE):
        with tf.variable_scope("qW_0", reuse=tf.AUTO_REUSE):
            loc = tf.get_variable("loc", [FLAGS.D, FLAGS.H])
            scale = tf.nn.softplus(tf.get_variable("scale", [FLAGS.D, FLAGS.H]))
            qW_0 = Normal(loc=loc, scale=scale)
        with tf.variable_scope("qW_1", reuse=tf.AUTO_REUSE):
            loc = tf.get_variable("loc", [FLAGS.H, FLAGS.H])
            scale = tf.nn.softplus(tf.get_variable("scale", [FLAGS.H, FLAGS.H]))
            qW_1 = Normal(loc=loc, scale=scale)
        with tf.variable_scope("qW_2", reuse=tf.AUTO_REUSE):
            loc = tf.get_variable("loc", [FLAGS.H, 1])
            scale = tf.nn.softplus(tf.get_variable("scale", [FLAGS.H, 1]))
            qW_2 = Normal(loc=loc, scale=scale)
        with tf.variable_scope("qb_0", reuse=tf.AUTO_REUSE):
            loc = tf.get_variable("loc", [FLAGS.H])
            scale = tf.nn.softplus(tf.get_variable("scale", [FLAGS.H]))
            qb_0 = Normal(loc=loc, scale=scale)
        with tf.variable_scope("qb_1", reuse=tf.AUTO_REUSE):
            loc = tf.get_variable("loc", [FLAGS.H])
            scale = tf.nn.softplus(tf.get_variable("scale", [FLAGS.H]))
            qb_1 = Normal(loc=loc, scale=scale)
        with tf.variable_scope("qb_2", reuse=tf.AUTO_REUSE):
            loc = tf.get_variable("loc", [1])
            scale = tf.nn.softplus(tf.get_variable("scale", [1]))
            qb_2 = Normal(loc=loc, scale=scale)

    inference = ed.KLqp({W_0: qW_0, b_0: qb_0,
                       W_1: qW_1, b_1: qb_1,
                       W_2: qW_2, b_2: qb_2}, data={X: X_train, y: y_train})
    
    print(qW_0)

              

    inference.run(n_iter=10)
    
    y_post = ed.copy(y, {
        W_0: qW_0, b_0: qb_0,
        W_1: qW_1, b_1: qb_1,
        W_2: qW_2, b_2: qb_2,})

    sess = ed.get_session()
    predictions=np.array([])
    for index in range(len(Xs)//FLAGS.N):
        xs = Xs[FLAGS.N*index:(index+1)*FLAGS.N]
        predictions = np.append(predictions, sess.run(y_post, feed_dict={X: xs}))


    return predictions




def run_program():
    data_name ='train_norm.mat'
    T = 200
    X, y = generate_train(data_name,200)
    dat_name = 'nn_'+data_name[0:-3]+'dat'
    X, y  = load_py_data(dat_name)
    y  = np.squeeze(y)
    Xs   = generate_vector(data_name,10000)
    dat_name = 'Xs_'+data_name[0:-3]+'dat'
    Xs  = load_tree_vector(dat_name)
    X = X.tolist()
    res =[max(y)]
    ed.set_seed(42)
    for idx  in range(T):
        ys = bayesian_optimization(X,y,Xs)
        x_new = find_a_candidate(X,Xs,ys,2)
        for ii in range(len(x_new)):
            xx = x_new[ii]
            X.append(xx)
            y_new = robust_function(xx,data_name, False)
            y = np.concatenate((y,y_new))

        res.append(max(y))
        print('Iteration', idx)

    plt.figure()
    plt.plot(res)
    plt.show()




 



if __name__ == '__main__':
	run_program()





