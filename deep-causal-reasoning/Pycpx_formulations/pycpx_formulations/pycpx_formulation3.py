#Importing what I need
import networkx as nx
from scipy.sparse import *
from scipy import *
import numpy as np
import pycpx 
import scipy.io as sio
#Read in the files
options=np.genfromtxt('leks_options.txt',names=True,dtype=('f8','f8','f8','a200','a200','a200'))
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


#Make the input matrix 
Inputs=np.zeros((nsp,nexp))
for i in range(nexp):
	for j in range(nsp):
		for k in range(inputs[0].size):
			if species[j]==inputs[0][k]:Inputs[j,i]=inputs[i+1][k]

#Loading pre-made matrices which are sparse
P=sio.mmread('testP.mtx')    #Products matrix
R=sio.mmread('test.mtx')     #Reactants matrix
M=sio.mmread('M.mtx')        #Measurements matrix

#converting to list on lists format for sparse matrices
P=lil_matrix(P)
R=lil_matrix(R)
M=lil_matrix(M)

I=abs(Inputs) #Matrix indicating if a node is input or not
USR=np.zeros((nsp,1)) #Matrix indicating if a node has an UpStream Reaction or not
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

#Helper variables indicating non-zero elements of the array
r=np.transpose(np.nonzero(R))
p=np.transpose(np.nonzero(P))
inp=np.transpose(np.nonzero(Inputs))
usr=np.transpose(np.nonzero(USR))
ninp=np.transpose(np.nonzero(1-I))
Pa=P.A

#Define constraints
print('Reading in now...only using non zero elements except in 8/9')
for j in range(X.shape[1]):
	for i in range(r.shape[0]):#reactions
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
	for i in range(X.shape[0]):#species (The memory runs out in this constraint)
        exp1=(Up[:,j].A*Pa[:,i]).sum()
        exp2=(Um[:,j].A*Pa[:,i]).sum()
        m.constrain(exp1-Xp[i,j]>=0)
        m.constrain(exp2-Xm[i,j]>=0)
	m.constrain(X[:,j]-Xp[:,j]+Xm[:,j]-B[:,j]==0)
	for i in range(inp.shape[0]):#species
		Ns=int(inp[i][0])
		m.constrain(X[Ns,j]-Inputs[Ns,j]==0)
	for i in range(usr.shape[0]):#species
		Ns=int(usr[i][0])
		m.constrain(X[Ns,j]-B[Ns,j]==0)
	for i in range(ninp.shape[0]):#species
		Ns=int(ninp[i][0])
		m.constrain(B[Ns,j]==0)  		
	m.constrain(abs(M[:,j].todense()-X[:,j])-ABS[:,j]<=0)

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

time=m.getSolverTIme()
print(time)
