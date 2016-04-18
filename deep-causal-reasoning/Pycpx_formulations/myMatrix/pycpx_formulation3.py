#Importing what I need
import networkx as nx
from scipy import stats
import numpy as np
import pycpx 
from scipy.sparse import *
import scipy.io as sio
#Read in the files
options=np.genfromtxt('leks_options.txt',names=True,dtype=('f8','f8','f8','a21','a21','a21'))
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
R=np.zeros((nr,nsp))
P=np.zeros((nr,nsp))
#Make the measurements matrix: maybe look for ways to make this  part nicer
M=np.zeros((nsp,1))
for i in range(nexp):
        for j in range(nsp):
                for k in range(measurements[0].size):
                        if measurements[0][k]==species[j]:M[j,i]=measurements[i+1][k]

#Make the R and P matrix
for i in range(nr):
     for j in range(nsp):
             if reactants[i]==species[j]:R[i,j]=1
	     if products[i]==species[j]:P[i,j]=1


#Make the input matrix 
Inputs=np.zeros((nsp,nexp),dtype=float)
for i in range(nexp):
	for j in range(nsp):
		for k in range(inputs[0].size):
			if species[j]==inputs[0][k]:Inputs[j,i]=inputs[i+1][k]


I=abs(Inputs) #Helper matrix for constraints 11/12
USR=np.zeros((nsp,1)) #Helper matrix for constraint 13 
USR[np.where(np.sum(P,axis=0)==0)]=1
#Define the model
m = pycpx.CPlexModel()
#Define the variables
#Regarding U+ U- : 1 to 5
r=np.transpose(np.nonzero(R))
p=np.transpose(np.nonzero(P))
inp=np.transpose(np.nonzero(Inputs))
usr=np.transpose(np.nonzero(USR))
ninp=np.transpose(np.nonzero(1-I))

m = pycpx.CPlexModel()
B=m.new((nsp,nexp),lb=-1,ub=1,vtype=int,name='B')
X=m.new((nsp,nexp),lb=-1,ub=1,vtype=int,name='X')
Xp=m.new((nsp,nexp),lb=0,ub=1,vtype=bool,name='Xp')
Xm=m.new((nsp,nexp),lb=0,ub=1,vtype=bool,name='Xm')
Up=m.new((nr,nexp),lb=0,ub=1,vtype=bool,name='Up')
Um=m.new((nr,nexp),lb=0,ub=1,vtype=bool,name='Um')
ABS=m.new((nsp,nexp),lb=0,ub=2,vtype=int,name='ABS')
dR=m.new((nr,nexp),lb=0,ub=100,vtype=int,name='dR')
dP=m.new((nr,nexp),lb=0,ub=100,vtype=int,name='dP')

for j in range(X.shape[1]):
	for i in range(r.shape[0]):
		Nr=int(r[i][0])
		Rs=int(r[i][1])
		Ps=int(p[i][1])
		m.constrain(Up[Nr,j]-X[Rs,j]*sigma[Nr]*R[Nr,Rs]>=0)
		m.constrain(Um[Nr,j]+X[Rs,j]*sigma[Nr]*R[Nr,Rs]>=0)
		m.constrain(Up[Nr,j]+Um[Nr,j]-1<=0)
		m.constrain(-Up[Nr,j]+Um[Nr,j]+X[Rs,j]*sigma[Nr]*R[Nr,Rs]>=0)
		m.constrain(Up[Nr,j]-Um[Nr,j]-X[Rs,j]*sigma[Nr]*R[Nr,Rs]>=0)
 		m.constrain(dR[Nr,j]*np.dot(R[Nr,Rs],(1-I[Rs,j]))-dP[Nr,j]*np.dot(P[Nr,Ps],(1-I[Ps,j]))-101+100*Up[Nr,j]<=0)
           	m.constrain(dR[Nr,j]*np.dot(R[Nr,Rs],(1-I[Rs,j]))-dP[Nr,j]*np.dot(P[Nr,Ps],(1-I[Ps,j]))-101+100*Um[Nr,j]<=0)	
	for i in range(X.shape[0]):
                #m.constrain(Up[:,j].A*P[:,i]-Xp[i,j]>=0)
              	m.constrain((Up[:,j].A*P[:,i]).sum()-Xp[i,j]>=0)
		m.constrain((Um[:,j].A*P[:,i]).sum()-Xm[i,j]>=0)
	m.constrain(X[:,j]-Xp[:,j]+Xm[:,j]-B[:,j]==0)
	for i in range(inp.shape[0]):
		Ns=int(inp[i][0])
		m.constrain(X[Ns,j]-Inputs[Ns,j]==0)
	for i in range(usr.shape[0]):
		Ns=int(usr[i][0])
		m.constrain(X[Ns,j]-B[Ns,j]==0)
	for i in range(ninp.shape[0]):
		Ns=int(ninp[i][0])
		m.constrain(B[Ns,j]==0)  		

m.constrain(abs(M-X)-ABS<=0)
              
#minimize the objective function with default parameters: algorithm used: auto
m.minimize(ABS.sum())

#Saving 
names=np.array(species)
names.shape=(nsp,nexp)
np.savetxt('Xs.txt',np.concatenate((names,m[X].A),axis=1),delimiter=" ", fmt="%s")
np.savetxt('Bs.txt',np.concatenate((names,m[B].A),axis=1),delimiter=" ", fmt="%s")
np.savetxt('Xps.txt',np.concatenate((names,m[Xm].A),axis=1),delimiter=" ", fmt="%s")
np.savetxt('Xms.txt',np.concatenate((names,m[Xp].A),axis=1),delimiter=" ", fmt="%s")
np.savetxt('Ds.txt',np.concatenate((names,np.dot(m[dP].A.T,P).T),axis=1),delimiter=" ", fmt="%s")
names=np.arange(nr)
names.shape=(nr,nexp)
np.savetxt('Ums.txt',np.concatenate((names,m[Um].A),axis=1),delimiter=" ", fmt="%s")
np.savetxt('Ups.txt',np.concatenate((names,m[Up].A),axis=1),delimiter=" ", fmt="%s")



