import networkx as nx
from scipy.sparse import *
import numpy as np
import scipy.io as sio
from sympy import *
execfile('sympybasedfns.py')

## OPTIONS for running/ Measurements and Network (all loaded as arrays)
options=np.genfromtxt('leks_options.txt',names=True,dtype=('f8','f8','f8','a200','a200','a200'))
mipgap=float(options['mipgap'])
number_solutions=float(options['number_solutions'])
timelimit=float(options['timelimit'])
measurements=np.loadtxt(str(options['measurements_file']),dtype='S')
inputs=np.loadtxt(str(options['inputs_file']),dtype='S')
net=np.loadtxt(str(options['network_file']),dtype='S')

#Defining the reactants and the products and input matrix
reactants=net[:,0].astype('str')
products=net[:,2].astype('str')
#species=set(reactants)|set(products)
species=set()
for i in range(len(net)):
    species.add(net[i][0])
    species.add(net[i][2])


species=list(species)
#inputs[1][np.where(inputs[1]=='NaN')]='0'
nsp=len(species)
nr=net.shape[0]
nexp=measurements.shape[0]-1
sigma=net[:,1].astype('int')

# Reading already stored P R M M_ind and Inputs matrices (These are made from the net,species...)
P=sio.mmread('Pf.mtx')
R=sio.mmread('Rf.mtx')
M=sio.mmread('Mf.mtx')
Inputs=sio.mmread('Inputsf.mtx')
M_ind=sio.mmread('Mindf.mtx')
P=lil_matrix(P)
R=lil_matrix(R)

#Making the Helper Matrices
I=Inputs.copy()
I[I!=0]=1
USR=np.zeros((nsp,1))
USR[np.where(P.sum(axis=0).T==0)]=1
I_nan=np.isnan(Inputs).astype('int')
I_nan=1-I_nan
I_nan=I_nan*I
notI=(1-I)
P_nI=lil_matrix((nr,nr))
PnI=P*notI
P_nI.setdiag(PnI)
R_d=P_nI*R
P_d=P_nI*P
PR=sigma.copy()
NR=sigma.copy()
PR[PR==-1]=0
NR[NR==1]=0

sigmaR=lil_matrix((nr,nr))
sigmaP=lil_matrix((nr,nr))
sigmaN=lil_matrix((nr,nr))
sigmaP.setdiag(PR)
sigmaN.setdiag(NR)
sigmaR.setdiag(sigma)
sigmaP=sigmaP*R
sigmaN=sigmaN*R
sigmaR=sigmaR*R
Pt=P.T

PR=PR.reshape((nr,1))
NR=NR.reshape((nr,1))
# recording time
import time as t

# Convert to sympy
t0=t.time()
sigmaN,sigmaP,NR,PR,notI,P,P_d,R_d,I_nan,USR,R,I,Inputs,M_ind,M,sigmaR,Pt,P_nI =[sci2sym(i) for i in [sigmaN,sigmaP,NR,PR,notI,P,P_d,R_d,I_nan,USR,R,I,Inputs,M_ind,M,sigmaR,Pt,P_nI]]
t01=t.time()

t1=t.time()
model=tocplex()
X=model.Makevars('x',species,nexp,1,-1,'I',0)
B=model.Makevars('B',species,nexp,1,-1,'I',0)
ABS=model.Makevars('abs',species,nexp,2,0,'I',100)
Xp=model.Makevars('x+',species,nexp,1,0,'B',20)
Xm=model.Makevars('x-',species,nexp,1,0,'B',20)
Up=model.Makevars('u+',nr,nexp,1,0,'B',0,name=False)
Um=model.Makevars('u-',nr,nexp,1,0,'B',0,name=False)
D=model.Makevars('dist',species,nexp,100,0,'C',0)
## the elements will be the dictionary variables (Call This before making any constraints)
t2=t.time()

## convert it to the format CPLEX can handle (unlist everything into a dict)
## Constraints
t3=t.time()                                                                   #C1                                                                   #C2
#model.Makecons(Up-sigmaR*X,'G',0)
#model.Makecons(Um+sigmaR*X,'G',0)
#model.Makecons(Um+Up,'L',1)                                                                             #C3
#model.Makecons(-Up+Um+sigmaR*X,'L',0)                                                                   #C4
#model.Makecons(Up-Um-sigmaR*X,'L',0)                                                                    #C5
model.Makecons(multiply_elem(Up,PR)-sigmaP*X,'G',0)
##model.Makecons(-multiply_elem(Up,NR)-sigmaN*X,'G',0)
##model.Makecons(multiply_elem(Um,PR)+sigmaP*X,'G',0)
##model.Makecons(-multiply_elem(Um,NR)+sigmaN*X,'G',0)
##model.Makecons(Um+Up,'L',1) 
##model.Makecons(-multiply_elem(Up,PR)+multiply_elem(Um,PR)+sigmaP*X,'L',0)
##model.Makecons(multiply_elem(Up,NR)-multiply_elem(Um,NR)+sigmaN*X,'L',0)
##model.Makecons(multiply_elem(Up,PR)-multiply_elem(Um,PR)-sigmaP*X,'L',0)
##model.Makecons(-multiply_elem(Up,NR)+multiply_elem(Um,NR)-sigmaN*X,'L',0)    
##model.Makecons(-Pt*Up+Xp,'L',0)                                                                         #C6
##model.Makecons(-Pt*Um+Xm,'L',0)                                                                         #C7
##model.Makecons(X-Xp+Xm-B,'E',0)                                                                         #C8
##model.Makecons(multiply_elem(ABS+X,M_ind),'G',multiply_elem(M,M_ind))                                   #C12
##model.Makecons(multiply_elem(ABS-X,M_ind),'G',multiply_elem(-M,M_ind))                                  #C13
##model.Makecons(P_d*D-R_d*D-101*P_nI*Up,'G',-100)                                                         #C14
##model.Makecons(P_d*D-R_d*D-101*P_nI*Um,'G',-100)                                                         #C15
##model.Makecons(multiply_elem(X,I_nan),'E',multiply_elem(Inputs,I_nan))                                  #C9
##model.Makecons(multiply_elem(X-B,USR),'E',0)                                                            #C10
##model.Makecons(multiply_elem(B,notI),'E',0)                                                             #C11
t4=t.time()


#Making cplex object
import cplex

m=cplex.Cplex()
#m.variables.add(obj=model.obj,ub=model.ub,lb=model.lb,types=model.vtypes,names=model.names)
#m.linear_constraints.add(lin_expr=model.cons_rows,senses=model.cons_senses,rhs=model.cons_rhs)
m.read("main_sympy2.lp")
# Set Model Parameters
m.parameters.threads.set(1)
m.parameters.mip.limits.populate.set(int(number_solutions))
m.parameters.mip.pool.absgap.set(0.0)
m.parameters.mip.pool.intensity=1
m.parameters.mip.pool.replace=0
m.parameters.mip.tolerances.mipgap.set(float(mipgap))
m.parameters.mip.strategy.search=1
#m.parameters.timelimit.set(float(timelimit))
m.parameters.randomseed=100
#m.write("main_sympy2.lp")
#m.parameters.preprocessing.presolve.set(0)
## Solve Problem
m.solve()
m.populate_solution_pool()

# parse the solution
nsol=m.solution.pool.get_num()

#get values
sol=[ m.solution.pool.get_values(i) for i in range(nsol)]
xs=getsols(sol,X)
Bs=getsols(sol,B)
ups=getsols(sol,Up)
ums=getsols(sol,Um)

#get names
xnames=X.tolist()
Bnames=B.tolist()
upnames=Up.tolist()
umnames=Um.tolist()

m.solution.pool.write("my_sol_matrix2.sol")
# calculate average xs and us

avx=[ float(sum( [ xs[i][j] for i in range(len(xs)) ]))/len(xs) for j in range(len(xs[0]))  ]
avup=[ float(sum( [ ups[i][j] for i in range(len(ups)) ]))/len(ups) for j in range(len(ups[0]))  ]
avum=[ float(sum( [ ums[i][j] for i in range(len(ums)) ]))/len(ums) for j in range(len(ums[0]))  ]


f=open("xs_matrix2.txt", "w")
[ f.write("%s\t%f\n" % (xnames[i], avx[i])) for i in range(len(xnames)) if avx[i]!=0.0 ]
f.close()

active_species=set([ species[i] for i in range(len(xnames)) if avx[i]!=0.0 ])
up_species=set([ species[i] for i in range(len(xnames))  if avx[i]>0.0 ])
down_species=set([ species[i] for i in range(len(xnames))  if avx[i]<0.0 ])

#Make this simpler
f=open("us_matrix2.txt", "w")
[ f.write("%s->%s\t%f\n" % (net[i][0], net[i][2], avup[i])) for i in range(len(avup)) if net[i][2] in up_species and avup[i]>0 ]
[ f.write("%s->%s\t-%f\n" % (net[i][0], net[i][2], avum[i])) for i in range(len(avum)) if net[i][2] in down_species and avum[i]>0 ]
f.close()


#create graphviz file

f=open("network_out_matrix2.dot", "w")
f.write("digraph {\n\n")
[ f.write("\"%s\"->\"%s\"\t[penwidth=%f, color=%s]\n" % (net[i][0], net[i][2], avup[i], "black")) for i in range(len(avup)) if net[i][2] in up_species and avup[i]>0 ]
[ f.write("\"%s\"->\"%s\"\t[penwidth=%f, color=%s]\n" % (net[i][0], net[i][2], avum[i], "red")) for i in range(len(avum)) if net[i][2] in down_species and avum[i]>0 ]
f.write("\n\n}")
f.close()


#write times to a file
f=open("time_used2.txt","w")
f.write("to convert scipy to sympy sparse matrices=%f \n" %(t01-t0) )
f.write("to make variable=%f \n" %(t2-t1) )
f.write("to make constraints=%f " %(t4-t3) )
f.close()

