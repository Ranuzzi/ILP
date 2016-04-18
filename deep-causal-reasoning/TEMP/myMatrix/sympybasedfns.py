##CLASS tocplex() helps to make calling constraints to cplex easier: uses scympy matrices
def make_symbolic(a,n, m,name=True):
    if(name==False):
        rows = []
        for i in xrange(n):
            col = []
            for j in xrange(m):
                col.append(Symbol('%s_{%d,%d}' % (a,i,j)))
            rows.append(col)
        return SparseMatrix(rows)
    else:
        rows = []
        for i in n:
            col = []
            for j in xrange(m):
                col.append(Symbol('%s_{%s,%d}' % (a,i,j)))
            rows.append(col)
    return SparseMatrix(rows)



class tocplex:
    def __init__(self):
        self.data = []
        self.ub=[]
        self.names=[]
        self.lb=[]
        self.vtypes=[]
        self.obj=[]
        self.cons_rows=[]
        self.cons_senses=[]
        self.cons_rhs=[]
        self.NAMEdict=dict()
        self.cons_len=[]
        self.helpers=[]
    def Makevars(self,a,n, m,UB,LB,Vtype,OBJ,name=True):
        MAT=make_symbolic(a,n, m,name)
        mat=range(len(MAT))
        Names=MAT.tolist()
        self.ub.extend([UB for j in mat])
        self.lb.extend([LB for j in mat])
        self.vtypes.extend([Vtype for j in mat])
        self.obj.extend([OBJ for j in mat])
        self.names.extend([str(i) for i in MAT])
        for i in mat:self.NAMEdict[MAT[i]]=len(self.NAMEdict)
        return(MAT)
        #def Helpers(self,objs):
        #ans=[sci2sym(i) for i in objs]
        #self.helpers.extend(ans)
        #return(ans)
    def Makecons(self,expr,sens,rhs):
        rows=[[[self.NAMEdict[i] for i in J.as_coefficients_dict().keys()],[int(k) for k in J.as_coefficients_dict().values()]] for J in expr if J!=0 ]
        sens=[sens]*len(rows)
        if(isinstance(rhs,int)): rhs=[rhs]*len(rows)
        else: rhs=[int(rhs[i]) for i in range(len(rhs)) if expr[i]!=0 ]
        self.cons_rows.extend(rows)
        self.cons_senses.extend(sens)
        self.cons_rhs.extend(rhs)
        self.cons_len.extend([len(rows)])

def multiply_elem(X,I,makeint=True):
    if(isinstance(I,MutableSparseMatrix)==False):I=sci2sym(I)
    if(isinstance(X,MutableSparseMatrix)==False):X=sci2sym(X)
    expr=SparseMatrix(X.shape[0],X.shape[1],0)
    if(makeint==True):
        for i in range(I.shape[0]):
            for j in range(X.shape[1]):
                expr[i,j]=X[i,j]*int(I[i,j])
    else:
        for i in range(I.shape[0]):
            for j in range(X.shape[1]):
                expr[i,j]=X[i,j]*int(I[i,j])
    return(expr)


def sqMat(index,value,n,m):
    if(index.shape[1]==1):index=np.tile(index,2)
    if(isinstance(value,int)):smat=SparseMatrix(n,m,{(index[i][0],index[i][1]):value for i in range(index.shape[0])})
    else:smat=SparseMatrix(n, m, {(index[i][0],index[i][1]):value[index[i][0],index[i][1]] for i in range(index.shape[0])})
    return(smat)


def sci2sym(sci):
    if(isinstance(I,MutableSparseMatrix)==False):
        n=sci.shape[0]
        m=sci.shape[1]
        index=np.transpose(np.nonzero(sci))
        smat=sqMat(index,sci,n,m)
    else: print('Already a Sympy Sparse Matrix')
    return(smat)


def getsols(sol,S):
    nsol=len(sol)
    xs=[[ sol[i][model.NAMEdict[k]] for k in S ] for i in range(nsol)]
    return(xs)


