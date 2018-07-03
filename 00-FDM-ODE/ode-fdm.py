
# coding: utf-8

# $$ -U''(x) = f(x) = 9sin(3x)-6x $$
# two-order center finite difference method
# 
# $$ \frac{U_{i-1}-2U_i + U_{i+1}}{h^2} = f_i $$

# In[1]:


import numpy as np
import matplotlib.pyplot as plt


# In[2]:


def Thomas(La, Mb, Uc, b):
    """
    Arguments:
        La -- [lower item for tri-diagonal matrix]
        Mb -- [mian item for tri-diagonal matrix]
        Uc -- [upper item for tri-diagonal matrix]
        b -- [AX = b, where A is the tri-diagonal matrix]
    """
    n = len(Mb)
    Uc[0] = Uc[0] / Mb[0]
    for i in range(1, n-1):
        Uc[i] = Uc[i] / (Mb[i] - La[i - 1] * Uc[i - 1])
    b[0] = b[0] / Mb[0]
    for i in range(1, n):
        b[i] = (b[i] - La[i-1]*b[i-1]) / (Mb[i] - La[i-1] * Uc[i-1])
    ls = list(range(n-1))[::-1]
#     print(b)
    for i in ls:
        b[i] = b[i] - Uc[i]*b[i+1]
    return b

def Exact(X):
    return np.sin(3*X) + np.power(X,3)

def f(X):
    return 9.0*np.sin(3.*X) - 6.*(X)

def FDMode(n):
    Xa = 0
    Xb = np.pi
    U0 = 0
    Upi = np.pi**3
    h = np.array((Xb-Xa)/n)
    h2 = h**2
    N = n+1 # all point
    U_numerical = np.zeros((N,1))
    U_numerical[0] = U0
    U_numerical[-1] = Upi
    X = np.linspace(Xa,Xb,N)
    La = np.zeros((n-2,))
    Mb = np.zeros((n-1,))
    Uc = np.zeros((n-2,))
    La[:] = -1.0/h2
    Mb[:] = 2.0/h2
    Uc[:] = -1.0/h2
    b = f(X)[1:-1]

    b[0] = b[0] + U0
    b[-1] = b[-1] + Upi/h2

    res = Thomas(La,Mb,Uc,b)
    exact = Exact(X)
    # errors 
    norm_max = np.max(np.abs(exact[1:-1]-res))
    norm_2 = np.sqrt(h*np.sum(np.abs(exact[1:-1]-res)**2))
#     print(n,norm_max, norm_2)
    return norm_max, norm_2
    
# FDMode(30)


# In[6]:


nlist = list(range(10,101,10))
# nlist.remove(30)
NN = len(nlist)
Norm2 = np.zeros((NN,))
NormMax = np.zeros((NN,))
for i in range(NN):
    norm_max, norm_2 = FDMode(nlist[i])
    Norm2[i] = norm_2
    NormMax[i] = norm_max

print("norm_2", Norm2)
print("norm_max", NormMax)
fig, ax = plt.subplots()
ax.loglog(nlist,Norm2, basex=np.e, basey=np.e)
ax.loglog(nlist,NormMax, basex=np.e, basey=np.e)

ax.set_xlabel('n', fontsize=10)
ax.set_ylabel('Error', fontsize=10)
ax.legend(["2 norm", "max norm"])
plt.show()

