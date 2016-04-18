# setting objective weights

obj_x<-1.0
obj_B<-3.0
obj_fx<-0.0


# add variables

variables<-list()
variables[["ids"]]<-array(NA, dim=c(nexp*nspc,1))
variables[["names"]]<-character()
variables[["upper_bound"]]<-numeric()
variables[["lower_bound"]]<-numeric()
variables[["vtype"]]<-character()
variables[["obj"]]<-numeric()


for (i in 1:nexp) {
  
  # add the x_j^k
  variables[["names"]]<-c(variables[["names"]], apply(array(spc$V1), 1, function(x) sprintf("x_{%d,%d}",x,i)))
  variables$upper_bound<-c(variables$upper_bound, rep_len(1, length.out=nspc))
  variables$lower_bound<-c(variables$lower_bound, rep_len(0, length.out=nspc))
  variables$vtype<-c(variables$vtype, rep_len("B", length.out=nspc))
  variables$obj<-c(variables$obj, rep_len(obj_x, length.out=nspc))
  
  #add the abs_j^k
  variables[["names"]]<-c(variables[["names"]], apply(array(spc$V1), 1, function(x) sprintf("abs_{%d,%d}",x,i)))
  variables$upper_bound<-c(variables$upper_bound, rep_len(1, length.out=nspc))
  variables$lower_bound<-c(variables$lower_bound, rep_len(0, length.out=nspc))
  variables$vtype<-c(variables$vtype, rep_len("B", length.out=nspc))
  variables$obj<-c(variables$obj, rep_len(0, length.out=nspc))
  
  #add the dist_j^k
  variables[["names"]]<-c(variables[["names"]], apply(array(spc$V1), 1, function(x) sprintf("dist_{%d,%d}",x,i)))
  variables$upper_bound<-c(variables$upper_bound, rep_len(50, length.out=nspc))
  variables$lower_bound<-c(variables$lower_bound, rep_len(0, length.out=nspc))
  variables$vtype<-c(variables$vtype, rep_len("C", length.out=nspc))
  variables$obj<-c(variables$obj, rep_len(0, length.out=nspc))
  
  #add the B_j^k
  variables[["names"]]<-c(variables[["names"]], apply(array(spc$V1), 1, function(x) sprintf("B_{%d,%d}",x,i)))
  variables$upper_bound<-c(variables$upper_bound, rep_len(1, length.out=nspc))
  variables$lower_bound<-c(variables$lower_bound, rep_len(0, length.out=nspc))
  variables$vtype<-c(variables$vtype, rep_len("B", length.out=nspc))
  variables$obj<-c(variables$obj, rep_len(obj_B, length.out=nspc))
  
  #add the z_i^k
  variables[["names"]]<-c(variables[["names"]], apply(array(1:length(net$V1)), 1, function(x) sprintf("z_{%d,%d}",x,i)))
  variables$upper_bound<-c(variables$upper_bound, rep_len(1, length.out=nreactions))
  variables$lower_bound<-c(variables$lower_bound, rep_len(0, length.out=nreactions))
  variables$vtype<-c(variables$vtype, rep_len("B", length.out=nreactions))
  variables$obj<-c(variables$obj, rep_len(0, length.out=nreactions))
  
  #add the fx_j^k
  variables[["names"]]<-c(variables[["names"]], apply(array(spc$V1), 1, function(x) sprintf("fx_{%d,%d}",x,i)))
  variables$upper_bound<-c(variables$upper_bound, rep_len(nsignals, length.out=nspc))
  variables$lower_bound<-c(variables$lower_bound, rep_len(0, length.out=nspc))
  variables$vtype<-c(variables$vtype, rep_len("I", length.out=nspc))
  variables$obj<-c(variables$obj, rep_len(0, length.out=nspc))
  
  #add the fz_i^k
  variables[["names"]]<-c(variables[["names"]], apply(array(1:length(net$V1)), 1, function(x) sprintf("fz_{%d,%d}",x,i)))
  variables$upper_bound<-c(variables$upper_bound, rep_len(nsignals, length.out=nreactions))
  variables$lower_bound<-c(variables$lower_bound, rep_len(0, length.out=nreactions))
  variables$vtype<-c(variables$vtype, rep_len("I", length.out=nreactions))
  variables$obj<-c(variables$obj, rep_len(0, length.out=nreactions))
  
}

#add the y_i
variables[["names"]]<-c(variables[["names"]], apply(array(1:length(net$V1)), 1, function(x) sprintf("y_{%d}",x)))
variables$upper_bound<-c(variables$upper_bound, rep_len(1, length.out=nreactions))
variables$lower_bound<-c(variables$lower_bound, rep_len(0, length.out=nreactions))
variables$vtype<-c(variables$vtype, rep_len("B", length.out=nreactions))
variables$obj<-c(variables$obj, rep_len(0, length.out=nreactions))

# creating hashtable for variable names 

for (i in 1:nexp) {
  begin_index<-(i-1)*(5*nspc+2*nreactions)
  for (j in 1:nspc) {
    variables$ids[[variables$names[begin_index+j]]]<-begin_index+j
  }
}


# define structures for Rcplex

# nvals<-nexp*(5*nspc+2*nreactions)+nreactions
bvec<-array(0, dim=c(200000,1))           # righthand side
sensevec<-array(0, dim=c(200000,1))       # sense of the constraints

# add constraints

counter<-0
nconstraints<-0
nzrow<-array(data=NA, dim=c(1000000,1))  # row indices for A matrix (constraints matrix)
nzcol<-array(data=NA, dim=c(1000000,1))  # col indices for A matrix
ncoef<-array(data=NA, dim=c(1000000,1))  # coefficients for A matrix 
for (i in 1:nexp) {
  
  begin_index<-(i-1)*(5*nspc+2*nreactions)
  
  #   z = R + y - 1
  for (j in 1:nreactions) {
    nconstraints<-nconstraints+1
    bvec[nconstraints]<--1
    sensevec[nconstraints]<-"E"
    nzrow[(counter+1):(counter+3)]<-rep_len(nconstraints, length.out=3)
    nzcol[(counter+1):(counter+3)]<-c(begin_index+net$V1[j], begin_index+4*nspc+j, nexp*(5*nspc+2*nreactions)+j)
    ncoef[(counter+1):(counter+3)]<-c(-1,1,-1)
    counter<-counter+3
  }
  
  
  # d_P \ge d_R - 50 + 51*z 
  
  for (j in 1:nreactions) {
    nconstraints<-nconstraints+1
    bvec[nconstraints]<--50
    sensevec[nconstraints]<-"G"
    nzrow[(counter+1):(counter+3)]<-rep_len(nconstraints, length.out=3)
    nzcol[(counter+1):(counter+3)]<-c(begin_index+2*nspc+net$V3[j], begin_index+2*nspc+net$V1[j], begin_index+4*nspc+j)
    ncoef[(counter+1):(counter+3)]<-c(1, -1, -51)
    counter<-counter+3
  }
    
  # x \leq d
  
  for (j in 1:nspc) {
    nconstraints<-nconstraints+1
    bvec[nconstraints]<-0
    sensevec[nconstraints]<-"L"
    nzrow[(counter+1):(counter+2)]<-rep_len(nconstraints, length.out=2)
    nzcol[(counter+1):(counter+2)]<-c(begin_index+j, begin_index+2*nspc+j)
    ncoef[(counter+1):(counter+2)]<-c(1, -1)
    counter<-counter+2
  }
  
  
  # d \leq 50
  for (j in 1:nspc) {
    nconstraints<-nconstraints+1
    bvec[nconstraints]<-50
    sensevec[nconstraints]<-"L"
    nzrow[(counter+1)]<-rep_len(nconstraints, length.out=1)
    nzcol[(counter+1)]<-begin_index+2*nspc+j
    ncoef[(counter+1)]<-1
    counter<-counter+1
  }
  
  
  # P \leq 1 - z
  
  for (j in 1:nreactions) {
    if ((net$V2[j]=="-1") | (net$V2[j]=="-2")) {
      nconstraints<-nconstraints+1
      bvec[nconstraints]<-1
      sensevec[nconstraints]<-"L"
      nzrow[(counter+1):(counter+2)]<-rep_len(nconstraints, length.out=2)
      nzcol[(counter+1):(counter+2)]<-c(variables$ids[[sprintf("x_{%d,%d}", net$V3[j], i)]], begin_index+4*nspc+j)
      ncoef[(counter+1):(counter+2)]<-c(1, 1)
      counter<-counter+2
    }  
  }
  
  # x \leq sum(z) + B
  
  for (j in 1:nspc) {
    upact<-upactivations[[j]]
    nupact<-length(upact)
    if (nupact >= 1) {
      nconstraints<-nconstraints+1
      bvec[nconstraints]<-0
      sensevec[nconstraints]<-"L"
      nzrow[(counter+1):(counter+(nupact+2))]<-rep_len(nconstraints, length.out=(nupact+2))
      nzcol[(counter+1):(counter+(nupact+2))]<-c(begin_index+j, apply(array(1:nupact), 1, function(x) {begin_index+4*nspc+upact[[x]]}), begin_index+3*nspc+j)
      ncoef[(counter+1):(counter+(nupact+2))]<-c(1, rep_len(-1, length.out=(nupact)), -1)
      counter<-counter+(nupact+2)
    }
  }
  
  # x \geq R - sum(zI)

  for (j in 1:nspc) {
    upact<-upactivations[[j]]
    nupact<-length(upact)
    if (nupact >= 1) {
      upinh<-upinhibitions[[j]]
      nupinh<-length(upinh)
      if (nupinh>=1) {
        for (k in 1:nupact) {
          nconstraints<-nconstraints+1
          bvec[nconstraints]<-0
          sensevec[nconstraints]<-"G"
          nzrow[(counter+1):(counter+(nupinh+2))]<-rep_len(nconstraints, length.out=(nupinh+2))
          nzcol[(counter+1):(counter+(nupinh+2))]<-c(begin_index+j, apply(array(1:nupinh), 1, function(x) {begin_index+4*nspc+upinh[[x]]}),  begin_index+4*nspc+upact[[k]])
          ncoef[(counter+1):(counter+(nupinh+2))]<-c(1, rep_len(1, length.out=(nupact)), -1)
          counter<-counter+nupinh+2
        }        
      } else {
        for (k in 1:nupact) {
          nconstraints<-nconstraints+1
          bvec[nconstraints]<-0
          sensevec[nconstraints]<-"G"
          nzrow[(counter+1):(counter+2)]<-rep_len(nconstraints, length.out=(2))
          nzcol[(counter+1):(counter+2)]<-c(begin_index+j, begin_index+4*nspc+upact[[k]])
          ncoef[(counter+1):(counter+2)]<-c(1, -1)
          counter<-counter+2
        }        
      }
    }
  }
  
  # sum(B) \leq 100
  
  nconstraints<-nconstraints+1
  bvec[nconstraints]<-100
  sensevec[nconstraints]<-"L"
  nzrow[(counter+1):(counter+nspc)]<-rep_len(nconstraints, length.out=(nspc))
  nzcol[(counter+1):(counter+nspc)]<-apply(array(1:nspc), 1, function(x) {begin_index+3*nspc+x})
  ncoef[(counter+1):(counter+nspc)]<-rep_len(1, length.out=(nspc))
  counter<-counter+nspc
  
  
  # x = B for all terminal input nodes
  
  terminal<-apply(array(1:nspc), 1, function(x) { nupact<-length(upactivations[[x]]) ; 
                                                  if (nupact==0) { 1 } else { 0 } })
  tinp<-apply(array(which(terminal==1)), 1, function(x) { if (is.na(inputs_hash[x])==FALSE) { x } else { 0 }} )
  tinp_pos<-tinp[tinp > 0]
  if (length(tinp_pos)>=1) {
    for (j in 1:length(tinp_pos)) {
      nconstraints<-nconstraints+1
      bvec[nconstraints]<-0
      sensevec[nconstraints]<-"E"
      nzrow[(counter+1):(counter+2)]<-rep_len(nconstraints, length.out=(2))
      nzcol[(counter+1):(counter+2)]<-c(begin_index+tinp_pos[j], begin_index+3*nspc+tinp_pos[j])
      ncoef[(counter+1):(counter+2)]<-c(1,-1)
      counter<-counter+2
    }  
  }
  
  
  # x = 0 for all terminal non input nodes
  
  terminal<-apply(array(1:nspc), 1, function(x) { nupact<-length(upactivations[[x]]) ; 
                                                  if (nupact==0) { 1 } else { 0 } })
  tinp<-apply(array(which(terminal==1)), 1, function(x) { if (is.na(inputs_hash[x])==TRUE) { x } else { 0 }} )
  tinp_pos<-tinp[tinp > 0]
  if (length(tinp_pos)>=1) {
    for (j in 1:length(tinp_pos)) {
      nconstraints<-nconstraints+1
      bvec[nconstraints]<-0
      sensevec[nconstraints]<-"E"
      nzrow[(counter+1)]<-nconstraints
      nzcol[(counter+1)]<-begin_index+tinp_pos[j]
      ncoef[(counter+1)]<-1
      counter<-counter+1
    }
  }
    
  #  abs_{j,k} + x_{j,k} \geq x{j,k,m}  
  
  nmeas<-length(measurements[1,])
  for (j in 1:nmeas) {
    nconstraints<-nconstraints+1
    bvec[nconstraints]<- measurements[i+1,j]
    sensevec[nconstraints]<-"G"
    nzrow[(counter+1):(counter+2)]<-rep_len(nconstraints, length.out=(2))
    nzcol[(counter+1):(counter+2)]<-c(begin_index+1*nspc+measurements[1,j], begin_index+measurements[1,j])
    ncoef[(counter+1):(counter+2)]<-c(1, 1)
    counter<-counter+2
  }
  
  
  #  abs_{j,k} - x_{j,k} \geq -x{j,k,m}
  
  nmeas<-length(measurements[1,])
  for (j in 1:nmeas) {
    nconstraints<-nconstraints+1
    bvec[nconstraints]<- -1*measurements[i+1,j]
    sensevec[nconstraints]<-"G"
    nzrow[(counter+1):(counter+2)]<-rep_len(nconstraints, length.out=(2))
    nzcol[(counter+1):(counter+2)]<-c(begin_index+1*nspc+measurements[1,j], begin_index+measurements[1,j])
    ncoef[(counter+1):(counter+2)]<-c(1, -1)
    counter<-counter+2
  }
  
  
  #  sum(abs_{j,k}) = 0 
  
  nmeas<-length(measurements[1,])
  nconstraints<-nconstraints+1
  bvec[nconstraints]<-0
  sensevec[nconstraints]<-"E"
  nzrow[(counter+1):(counter+nmeas)]<-rep_len(nconstraints, length.out=(nmeas))
  nzcol[(counter+1):(counter+nmeas)]<-apply(array(1:nmeas), 1, function(x) {begin_index+1*nspc+measurements[1,x]})
  ncoef[(counter+1):(counter+nmeas)]<-rep_len(1, length.out=nmeas)
  counter<-counter+nmeas
  
  
  # boundary conditions at the genes 
  # fz = 1
  
  nmeas<-length(measurements[1,])
  for (j in 1:nmeas) {
    upact<-upactivations[[measurements[1,j]]]
    nupact<-length(upact)
    if (nupact>=1) {
      nconstraints<-nconstraints+1
      bvec[nconstraints]<-measurements[i+1,j]
      sensevec[nconstraints]<-"E"
      nzrow[(counter+1):(counter+nupact)]<-rep_len(nconstraints, length.out=(nupact))
      nzcol[(counter+1):(counter+nupact)]<-apply(array(1:nupact), 1, function(x) {begin_index+5*nspc+nreactions+upact[[x]]})
      ncoef[(counter+1):(counter+nupact)]<-rep_len(1, length.out=nupact)
      counter<-counter+nupact
    }
  }
  
  # fx = sum(z_out)
  
  for (j in 1:nspc) {
    dact<-dactivations[[j]]
    ndact<-length(dact)
    if (ndact >= 1) {
      nconstraints<-nconstraints+1
      bvec[nconstraints]<-0
      sensevec[nconstraints]<-"E"
      nzrow[(counter+1):(counter+ndact+1)]<-rep_len(nconstraints, length.out=(ndact+1))
      nzcol[(counter+1):(counter+ndact+1)]<-c(apply(array(1:ndact), 1, function(x) {begin_index+5*nspc+nreactions+dact[[x]]}), begin_index+4*nspc+nreactions+j)
      ncoef[(counter+1):(counter+ndact+1)]<-c(rep_len(-1, length.out=ndact), 1)
      counter<-counter+ndact+1
    }
  }
  
  # sum(Zin) = sum(Zout)
  
  interm<-apply(array(1:nspc), 1, function(x) { nupact<-length(upactivations[[x]]) ; 
                                                if (nupact==0) { 0 } else { 1 } })
  interm2<-apply(array(which(interm==1)), 1, function(x) { ndact<-length(dactivations[[x]]) ; 
                                                           if (ndact==0) { 0 } else { x } })
  interm_pos<-interm2[interm2>0]
  for (j in 1:length(interm_pos)) {
    if (interm_pos[j]!=0) {
      nconstraints<-nconstraints+1
      bvec[nconstraints]<-0
      sensevec[nconstraints]<-"E"
      upact<-upactivations[[interm_pos[j]]]
      nupact<-length(upact)
      dact<-dactivations[[interm_pos[j]]]
      ndact<-length(dact)
      nzrow[(counter+1):(counter+ndact+nupact)]<-rep_len(nconstraints, length.out=(ndact+nupact))
      nzcol[(counter+1):(counter+ndact+nupact)]<-c(apply(array(1:nupact), 1, function(x) { begin_index+5*nspc+nreactions+upact[[x]] }), apply(array(1:ndact), 1, function(x) { begin_index+5*nspc+nreactions+dact[[x]] })   )
      ncoef[(counter+1):(counter+ndact+nupact)]<-c(rep_len(1, length.out=nupact), rep_len(-1, length.out=ndact))
      counter<-counter+ndact+nupact
    }
  }
  
  
  # fx <= x * 1000
  
  for (j in 1:nspc) {
    nconstraints<-nconstraints+1
    bvec[nconstraints]<-0
    sensevec[nconstraints]<-"L"
    nzrow[(counter+1):(counter+2)]<-rep_len(nconstraints, length.out=(2))
    nzcol[(counter+1):(counter+2)]<-c(begin_index+j, begin_index+4*nspc+nreactions+j)
    ncoef[(counter+1):(counter+2)]<-c(-1000, 1)
    counter<-counter+2
  }
  
    
  # fz <= z * 1000
  for (j in 1:nreactions)  {
    nconstraints<-nconstraints+1
    bvec[nconstraints]<-0
    sensevec[nconstraints]<-"L"
    nzrow[(counter+1):(counter+2)]<-rep_len(nconstraints, length.out=(2))
    nzcol[(counter+1):(counter+2)]<-c(begin_index+4*nspc+j, begin_index+5*nspc+nreactions+j)
    ncoef[(counter+1):(counter+2)]<-c(-1000, 1)
    counter<-counter+2
  }
  
  
  # z <= fz 
  
  for (j in 1:nreactions)  {
    nconstraints<-nconstraints+1
    bvec[nconstraints]<-0
    sensevec[nconstraints]<-"L"
    nzrow[(counter+1):(counter+2)]<-rep_len(nconstraints, length.out=(2))
    nzcol[(counter+1):(counter+2)]<-c(begin_index+4*nspc+j, begin_index+5*nspc+nreactions+j)
    ncoef[(counter+1):(counter+2)]<-c(1, -1)
    counter<-counter+2
  }  
  
}


# y_fixed = 1

fixed_reactions<-array(which((net$V2==2) | (net$V2==-2)))
nfixed<-length(fixed_reactions)
if (nfixed >= 1) {
  for (j in 1:nfixed)  {
    nconstraints<-nconstraints+1
    bvec[nconstraints]<-1
    sensevec[nconstraints]<-"E"
    nzrow[(counter+1)]<-rep_len(nconstraints, length.out=(1))
    nzcol[(counter+1)]<-c(begin_index+4*nspc+j, (nexp)*(5*nspc+2*nreactions)+fixed_reactions[x])
    ncoef[(counter+1)]<-1
    counter<-counter+1
  }    
}

# put everything into A matrix

A<-sparseMatrix(i=nzrow[1:counter], j=nzcol[1:counter], x=ncoef[1:counter])

# # print constraints
# 
# constraints_start<-98366
# constraints_end<-98366
# 
# lpfile<-character(length=constraints_end-constraints_start+1)
# for (i in constraints_start:constraints_end) {
#   nz<-which(A[i,]!=0) 
#   cons<-character(length=(length(nz)+2))
#   cons<-c(apply(array(1:length(nz)), 1, function(x) {
#     sprintf("%d %s ", A[i, nz[x]], variables$names[nz[x]])
#   }), "=", sprintf("%f", bvec[i]))
#   lpfile[i]<-paste(cons, sep="", collapse="")
# }


# solve problem

sol<-Rcplex(cvec=variables$obj, Amat=A, bvec=bvec[1:nconstraints], Qmat=NULL, lb = variables$lower_bound, ub = variables$upper_bound, control = list(epgap=mipgap, tilim=timelimit, solnpoolgap=0.0, solnpoolintensity=4),
         objsense = "min", sense = sensevec[1:nconstraints], vtype = variables$vtype, n = populatelim)

# parse solution and break into vectors

nsol<-length(sol)
if (nsol>=1) {
  xs<-array(NA, dim=c(nexp*nspc, nsol))
  Bs<-array(NA, dim=c(nexp*nspc, nsol))
  zs<-array(NA, dim=c(nexp*nreactions, nsol))
  fxs<-array(NA, dim=c(nexp*nspc, nsol))
  fzs<-array(NA, dim=c(nexp*nreactions, nsol))
  ys<-array(NA, dim=c(nreactions, nsol))
  for (j in 1:nsol) {
    for (i in 1:nexp) {
      begin_index<-(i-1)*(5*nspc+2*nreactions)
      xs[((i-1)*nspc+1):(i*nspc), j]<-apply(array(1:nspc), 1, function(x) sol[[j]]$xopt[begin_index+x])
      zs[((i-1)*nreactions+1):(i*nreactions), j]<-apply(array(1:nreactions), 1, function(x) sol[[j]]$xopt[begin_index+4*nspc+x])
      fxs[((i-1)*nspc+1):(i*nspc), j]<-apply(array(1:nspc), 1, function(x) sol[[j]]$xopt[begin_index+4*nspc+nreactions+x])
      fzs[((i-1)*nreactions+1):(i*nreactions), j]<-apply(array(1:nreactions), 1, function(x) sol[[j]]$xopt[begin_index+5*nspc+nreactions+x])
      Bs[((i-1)*nspc+1):(i*nspc), j]<-apply(array(1:nspc), 1, function(x) sol[[j]]$xopt[begin_index+3*nspc+x])
    }  
    ys[1:nreactions, j]<-apply(array(1:nreactions), 1, function(x) sol[[j]]$xopt[nexp*(5*nspc+2*nreactions)+x])
  }
}


