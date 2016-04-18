#Importing what I need
import networkx as nx
from scipy.sparse import *
from scipy import *
import numpy as np
import pycpx 
import scipy.io as sio
#Read in the files
options=np.genfromtxt('leks_options.txt',names=True,dtype=('f8','f8','f8','a21','a21','a51'))
mipgap=float(options['mipgap'])
number_solutions=float(options['number_solutions'])
timelimit=float(options['timelimit'])
measurements=np.loadtxt(str(options['measurements_file']),dtype='S')
inputs=np.loadtxt(str(options['inputs_file']),dtype='S')
net=np.loadtxt(str(options['network_file']),dtype='S')

#Defining the reactants and the products
reactants=net[:,0].astype('str')
products=net[:,2].astype('str')
species=set(reactants)|set(products)
species=list(species)

nsp=len(species)
nr=net.shape[0]
nexp=measurements.shape[0]-1
sigma=np.reshape(net[:,1].astype('int'),(nr,nexp))
R=lil_matrix((nr,nsp))
P=lil_matrix((nr,nsp))
#Make the measurements matrix: maybe look for ways to make this  part nicer
#M=lil_matrix((nsp,nexp))
#for i in range(nexp):
#        for j in range(nsp):
#                for k in range(measurements[0].size):
#                        if measurements[0][k]==species[j]:M[j,i]=measurements[i+1][k]

#Make the R and P matrix
#for i in range(nr):
#     for j in range(nsp):
#             if reactants[i]==species[j]:R[i,j]=1
#	     if products[i]==species[j]:P[i,j]=1


#Make the input matrix 
Inputs=np.zeros((nsp,nexp))
for i in range(nexp):
	for j in range(nsp):
		for k in range(inputs[0].size):
			if species[j]==inputs[0][k]:Inputs[j,i]=inputs[i+1][k]


P=sio.mmread('testP.mtx')
R=sio.mmread('test.mtx')
M=sio.mmread('M.mtx')

I=abs(Inputs) #Helper matrix for constraints 11/12
USR=np.zeros((nsp,1)) #Helper matrix for constraint 13 
USR[np.where(P.sum(axis=0).T==0)]=1
#Define the model
m = pycpx.CPlexModel()
#Define the variables
B=m.new((nsp,nexp),lb=-1,ub=1,vtype=int,name='B')
X=m.new((nsp,nexp),lb=-1,ub=1,vtype=int,name='X')
Xp=m.new((nsp,nexp),lb=0,ub=1,vtype=bool,name='Xp')
Xm=m.new((nsp,nexp),lb=0,ub=1,vtype=bool,name='Xm')
Up=m.new((nr,nexp),lb=0,ub=1,vtype=bool,name='Up')
Um=m.new((nr,nexp),lb=0,ub=1,vtype=bool,name='Um')
ABS=m.new((nsp,nexp),lb=0,ub=2,vtype=int,name='ABS')
dR=m.new((nr,nexp),lb=0,ub=100,vtype=int,name='dR')
dP=m.new((nr,nexp),lb=0,ub=100,vtype=int,name='dP')
#Define constraints
for j in range(X.shape[1]):
	for i in range(R.shape[0]):
		rowR=R.getrow(i).todense()
		rowP=P.getrow(i).todense()
		m.constrain(Up[i,j].A-X[:,j]*sigma[i]*rowR>=0)
		m.constrain(Um[i,j].A+X[:,j]*sigma[i]*rowR>=0)
		m.constrain(Up[i,j]+Um[i,j]-1<=0)
		m.constrain(-Up[i,j]+Um[i,j]+X[:,j]*sigma[i]*rowR>=0)
		m.constrain(Up[i,j]-Um[i,j]-X[:,j]*sigma[i]*rowR>=0)
 		m.constrain(dR[i,j].A*np.dot(rowR,(1-I[:,j]))-dP[i,j].A*np.dot(rowP,(1-I[:,j]))-101+100*Up[i,j]<=0)
           	m.constrain(dR[i,j].A*np.dot(rowR,(1-I[:,j]))-dP[i,j].A*np.dot(rowP,(1-I[:,j]))-101+100*Um[i,j]<=0)

for j in range(X.shape[1]):	
	for k in range(X.shape[0]):
		rowM=M.getrow(k).todense()
		colP=P.getcol(k).todense()
           	m.constrain(Up[:,j].T*colP-Xp[k,j]>=0)
 		m.constrain(Um[:,j].T*colP-Xm[k,j]>=0)
		m.constrain(X[k,j]-Xp[k,j]+Xm[k,j]-B[k,j]==0)
		m.constrain(X[k,j]*I[k,j]-Inputs[k,j]==0)
		m.constrain(B[k,j]*(1-I[k,j])==0)
		m.constrain(X[k,j]*USR[k,j]-B[k,j].A*USR[k,j]==0)
		m.constrain(abs(rowM[j]-X[k,j])-ABS[k,j]<=0)
              

#minimize the objective function with default parameters: algorithm used: auto
m.minimize(ABS.sum())

#Saving 
names=np.array(species)
names.shape=(nsp,nexp)
np.savetxt('Xs.txt',np.concatenate((names,m[X].A),axis=1),delimiter=" ", fmt="%s")
np.savetxt('Bs.txt',np.concatenate((names,m[B].A),axis=1),delimiter=" ", fmt="%s")
np.savetxt('Xps.txt',np.concatenate((names,m[Xm].A),axis=1),delimiter=" ", fmt="%s")
np.savetxt('Xms.txt',np.concatenate((names,m[Xp].A),axis=1),delimiter=" ", fmt="%s")
names=np.arange(nr)
names.shape=(nr,nexp)
np.savetxt('Ums.txt',np.concatenate((names,m[Um].A),axis=1),delimiter=" ", fmt="%s")
np.savetxt('Ups.txt',np.concatenate((names,m[Up].A),axis=1),delimiter=" ", fmt="%s")



