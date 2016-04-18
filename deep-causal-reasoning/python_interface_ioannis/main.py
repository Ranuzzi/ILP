# -*- coding: utf-8 -*-
"""
Created on Thu May 08 16:34:58 2014

@author: Teo Sakel
"""

__author__ = 'Teo Sakel (teo.sakel@gmail.com)'
__version__= '0.3'

# Modules Used
import cplex
from math import isnan
from os import path as op
from random import shuffle
from operator import itemgetter, add

# Define Variable nomenclature
_x = lambda spc, k: "x_{{{},{}}}".format(spc, k)
_B = lambda spc, k: "B_{{{},{}}}".format(spc, k)
_d = lambda spc, k: "d_{{{},{}}}".format(spc, k)
_y = lambda i: "y_{}".format(i)
_z = lambda i, k: "z_{{{},{}}}".format(i, k)
_F = lambda i, k: "F_{{{},{}}}".format(i, k)
_Fx = lambda spc, k: "Fx_{{{},{}}}".format(spc, k)


# Core Function
def write_DCR_model(filename, reactions, N_up, species, x_m, B_0, write_model=True,
                    x_weights=None, B_weights=None, longest_cascade=50, max_perturbations=100):
    """
    A function that creates a Cplex model and writes it without solving it.

    """

    # Setup Variables

    I = len(reactions)
    N = len(species)
    K = len(x_m[x_m.keys()[0]])
    # The maximum allowed flow is equal to the maximum sum of the flows that spring from all the sources in across all samples
    Fmax = max([sum(x) for x in zip(*[[0 if isnan(v) else v for v in val] for val in x_m.values()])])
    # Parse the Reaction list to find the reactions associated of every node.
    neighbors = {spc: [[], [], []] for spc in species}
    for i, r in enumerate(reactions[:N_up]):
        src, tgt = r
        neighbors[src][0].append(i) # src activates i
        neighbors[tgt][1].append(i) # tgt is activated by i
    for i, r in enumerate(reactions[N_up:]):
        src, tgt = r
        neighbors[tgt][2].append(N_up + i) # tgt is inhibitied by i

    # Define Model
    m = cplex.Cplex()

    # Add Variables

    # y_i : is reaction (i) present in the final Network?
    vnames = map(_y, range(I))
    numVar = I
    coef = [0]*numVar
    lwb = [0]*numVar
    upb = [1]*numVar
    vtypes = "B"*numVar
    m.variables.add(obj=coef, lb=lwb, ub=upb, types=vtypes, names=vnames)

    # z_{i,k}: is reactions (i) active in experiment/sample (k)?
    vnames = [_z(i, k) for i in range(I) for k in range(K)]
    numVar = I*K
    coef = [0]*numVar
    lwb = [0]*numVar
    upb = [1]*numVar
    vtypes = "B"*numVar
    m.variables.add(obj=coef, lb=lwb, ub=upb, types=vtypes, names=vnames)

    # F_{i,k}: flow through reaction (i) in experiment/sample (k)
    # (Flows are only defined for possitive reactions)
    vnames = [_F(i, k) for i in range(N_up) for k in range(K)]
    numVar = len(vnames)
    coef = [0]*numVar
    lwb = [0]*numVar
    upb = [Fmax]*numVar
    vtypes = "C"*numVar
    m.variables.add(obj=coef, lb=lwb, ub=upb, types=vtypes, names=vnames)

#    # Fx_{j,k}: flow through node (j) in experiment/sample (k)
#    # Can be computed later but added for extra functionality (in case we would like to maximize the flow through a node)
#    vnames = [_Fx(spc, k) for spc in species for k in range(K)]
#    numVar = len(vnames)
#    coef = [0]*numVar
#    lwb = [0]*numVar
#    upb = [Fmax]*numVar
#    vtypes = "C"*numVar
#    m.variables.add(obj=coef, lb=lwb, ub=upb, types=vtypes, names=vnames)

    # x_{j,k}: is node (j) active in experiment/sample (k)?
    vnames = [_x(spc, k) for spc in species for k in range(K)]
    numVar = N*K
    if type(x_weights) is list:
        coef = x_weights 
    elif isinstance(x_weights,(float, int, long)):
        coef = [x_weights]*numVar
    else: 
        coef = [1]*numVar
    lwb = [0]*numVar
    upb = [1]*numVar
    vtypes = "C"*numVar
    m.variables.add(obj=coef, lb=lwb, ub=upb, types=vtypes, names=vnames)

    # B_{j,k}: is node (j) perturbed in experiment/sample (k)
    vnames = [_B(spc, k) for spc in B_0 for k in range(K)]
    numVar = len(vnames)
    if type(B_weights) is list:
        coef = B_weights 
    elif isinstance(B_weights,(float, int, long)):
        coef = [B_weights]*numVar
    else: 
        coef = [3]*numVar
    lwb = [0]*numVar
    upb = [1]*numVar
    vtypes = "B"*numVar
    m.variables.add(obj=coef, lb=lwb, ub=upb, types=vtypes, names=vnames)

    # d_{j,k}: distance of node (j) in experiment/sample (k) (from a perturbed one)
    vnames = [_d(spc, k) for spc in species for k in range(K)]
    numVar = N*K
    coef = [0]*numVar
    lwb = [0]*numVar
    upb = [longest_cascade]*numVar
    vtypes = "C"*numVar
    m.variables.add(obj=coef, lb=lwb, ub=upb, types=vtypes, names=vnames)

    # Add Constrains

    # x_{j,k} = m_{j,k}: Fit the data perfectly
    lhs = [[[_x(spc, k)], [1]] for spc in x_m for k in range(K)]
    # if a value is missing the use a trivial constrain (good for reshuffling the constrains latter)
    sns = ''.join(['G' if isnan(x_m[spc][k]) else 'E' for spc in x_m for k in range(K)])
    rhs = [0 if isnan(x_m[spc][k]) else x_m[spc][k] for spc in x_m for k in range(K)]
    row_names = ["Fit_{%s,%d}" % (spc, k) for spc in x_m for k in range(K)]
    m.linear_constraints.add(lin_expr=lhs, senses=sns, rhs=rhs, names=row_names)

    # z_{i,k} <= y_i: A reaction cannot be active is it's not present
    lhs = [[[_z(i, k), _y(i)], [1, -1]] for i in range(I) for k in range(K)]
    sns = 'L'*len(lhs)
    rhs = [0]*len(lhs)
    row_names = ["Feasible_{%d,%d}" % (i, k) for i in range(I) for k in range(K)]
    m.linear_constraints.add(lin_expr=lhs ,senses=sns, rhs=rhs, names=row_names)

    # z_{i,k} <= x{src,k}: A reaction cannot be active if the source node is not active
    lhs = [[[_z(i, k), _x(r[0], k)], [1, -1]] 
            for i, r in enumerate(reactions) for k in range(K)]
    sns = 'L'*len(lhs)
    rhs = [0]*len(lhs)
    row_names = ["Reactant_%s_{%d,%d}" % (r[0], i, k)
                 for i, r in enumerate(reactions) for k in range(K)]
    m.linear_constraints.add(lin_expr=lhs ,senses=sns, rhs=rhs, names=row_names)

    # z_{i,k} >= x_{j,k}+y_i-1: A reaction is active if it's feasible and the reactant is active
    lhs = [[[_z(i, k), _x(r[0], k), _y(i)], [1, -1, -1]] 
            for i, r in enumerate(reactions) for k in range(K)]
    sns = 'G'*len(lhs)
    rhs = [-1]*len(lhs)
    row_names = ["ActiveZ_{%d,%d}" % (i, k) for i in range(I) for k in range(K)]
    m.linear_constraints.add(lin_expr=lhs ,senses=sns, rhs=rhs, names=row_names)
    
    # y_i <= sum_k(z_{i,k}): A reaction that is never active should be removed from the network.
    lhs = [[[_y(i)] + [_z(i, k) for k in range(K)], [1] + [-1]*K] for i in range(I)]
    sns = 'L' * I
    rhs = [0] * I
    row_names = ["Useful_%d" % i for i in range(I)]
    m.linear_constraints.add(lin_expr=lhs ,senses=sns, rhs=rhs, names=row_names)

    # x_{j,k} <= 1-z_{i,k} for all j in P_i if sign(i)=-1 (inhibitory reaction)
    lhs = [[[_x(r[1], k), _z(N_up+i, k)], [1, 1]] 
            for i, r in enumerate(reactions[N_up:]) for k in range(K)]
    sns = 'L' * len(lhs)
    rhs = [1] * len(lhs)
    row_names = ["Inhibit_$s_{%d,%d}" % (r[1], N_up+i, k) 
                    for i, r in enumerate(reactions[N_up:]) for k in range(K)]
    m.linear_constraints.add(lin_expr=lhs, senses=sns, rhs=rhs, names=row_names)

    # x_{j,k} >= z_{i,k} - sum_{j \in H_i}(z_{i,k}): if reactions (i) is active then its product must be active as well, unless it's inhibited by an other reaction
    lhs = [[[_x(r[1], k), _z(i, k)] + [_z(j, k) for j in neighbors[r[1]][2]],
            [1, -1] + [1]*len(neighbors[r[1]][2])]
            for i, r in enumerate(reactions[:N_up]) for k in range(K)]
    sns = 'G'*len(lhs)
    rhs = [0]*len(lhs)
    row_names = ["Produce_%s_{%d,%d}" % (r[1],i,k)
                 for i, r in enumerate(reactions[:N_up]) for k in range(K)]
    m.linear_constraints.add(lin_expr=lhs, senses=sns, rhs=rhs, names=row_names)
    # del inhibitors

    # x_{j,k} <= sum(z_{i,k})+B_{j,k}: a node must remain inactive unless it's activated by a reaction or it's pertrurbed
    lhs = [[[_x(spc, k), _B(spc, k)] + [_z(i, k) for i in neighbors[spc][1]],
             [1, -1] + [-1]*len(neighbors[spc][1])] for spc in species]
    sns = 'L' * len(lhs)
    rhs = [0] * len(lhs)
    row_names = ["Silence_{%s,%d}" %(spc,k) for spc in species for k in range(K)]
    m.linear_constraints.add(lin_expr=lhs, senses=sns, rhs=rhs, names=row_names)

    # d_{tgt,k} >= d_{src,k}+1+(z_{i,k}-1)*M: Remove Cicles and long cascades in the final graph
    lhs = [[[_d(r[1], k), _d(r[0], k), _z(i, k)], [1, -1, -longest_cascade]] 
             for i,r in enumerate(reactions) for k in range(K)]
    sns = 'G'*len(lhs)
    rhs = [1-longest_cascade]*len(lhs)
    row_names = ["Distance_{%d,%d}" %(i,k) for i in range(I) for k in range(K)]
    m.linear_constraints.add(lin_expr=lhs, senses=sns, rhs=rhs, names=row_names)

    # Flow conservation at every node: sum(F_in) = sum(F_out)
    lhs = [[[_F(i, k) for i in n[0]] + [_F(i, k) for i in n[1]],
             [1]*len(n[0]) + [-1]*len(n[1])]
             for n in neighbors.values() if n[0] and n[1] 
             for k in range(K)]
    sns = 'E' * len(lhs)
    rhs = [0] * len(lhs)
    row_names = ["Flux_{%s,%d}" % (spc, k) for spc, n in neighbors.iteritems()
                 if n[0] and n[1] for k in range(K)]
    m.linear_constraints.add(lin_expr=lhs, senses=sns, rhs=rhs, names=row_names)
    
#    # Node Flow calculation: Fx(spc, k) = sum(F_out)
#    lhs = [[[_Fx(spc, k)] + [_F(i, k) for i in neighbors[spc][0]], 
#             [1] + [-1]*len(neighbors[spc][0])] for spc in species for k in range(K)]
#    sns = 'E' * len(lhs)
#    rhs = [0] * len(lhs)
#    row_names = ["Flow_{%s,%d}" % (spc, k) for spc in species for k in range(K)]
#    m.linear_constraints.add(lin_expr=lhs, senses=sns, rhs=rhs, names=row_names)

    # Flow sources
    lhs = [[[_F(i, k) for i in neighbors[spc][1]], [1]*len(neighbors[spc][1])]
             for spc in x_m for k in range(K)]
    sns = ''.join(['G' if isnan(x_m[spc][k]) else 'E' for spc in x_m for k in range(K)])
    rhs = [0 if isnan(x_m[spc][k]) else x_m[spc][k] for spc in x_m for k in range(K)]
    row_names = ["Source_{%s,%d}" % (spc, k) for spc in x_m for k in range(K)]
    
    lhs += [[[_F(i, k) for i in neighbors[spc][1]],
              [1]*len(neighbors[spc][1])] for spc in species 
              if not neighbors[spc][0] and spc not in x_m for k in range(K)]
    sns += 'E' * (len(lhs)-len(sns))             
    rhs += [0] * (len(lhs)-len(rhs))
    row_names += ["Silence_{%s,%d}" % (spc, k) 
                  for spc in species if not neighbors[spc][0] and spc not in x_m 
                  for k in range(K)]
    m.linear_constraints.add(lin_expr=lhs, senses=sns, rhs=rhs, names=row_names)

    # Active Flow
    lhs = [[[_F(i, k), _z(i, k)], [1, -Fmax]] 
                for i, r in enumerate(reactions[:N_up]) for k in range(K)]
    sns = 'L'*len(lhs)
    rhs = [0]*len(lhs)
    row_names = ["ActiveF_{%d,%d}" % (i, k) 
                    for i, r in enumerate(reactions[:N_up]) for k in range(K)]
    m.linear_constraints.add(lin_expr=lhs, senses=sns, rhs=rhs, names=row_names)

    # Imposed Stimuli
    lhs = [[[_B(spc, k)], [1]] for spc, b0 in B_0.iteritems() for k, val in b0]
    sns = 'E'*len(lhs)
    rhs = [val for spc, b0 in B_0.iteritems() for k, val in b0]
    row_names = ["Impose_{%s,%d}" %(spc,k) for spc, b0 in B_0.iteritems() for k, val in b0]
    m.linear_constraints.add(lin_expr=lhs, senses=sns, rhs=rhs, names=row_names)
    #del B_0

    # sum(B) <= max_perturbations: limit the possible number of stimuli
    lhs = [[[_B(spc, k) for spc in B_0], [1]*N]  for k in range(K)]
    sns = 'L'*K
    rhs = [max_perturbations]*K
    row_names = ["Max Petrubations_%d" % k for k in range(K)]
    m.linear_constraints.add(lin_expr=lhs, senses=sns, rhs=rhs, names=row_names)

    if write_model: m.write(filename)

    return m

# Auxiliary Functions
def read_PKN(filename):
    """
    A function that reads a Simple Interaction Format file (.sif) and parses it into a list of tuples.
    Each tuple consists of 2 elements (reactant, product). 
    The first tuples of the list represent up-regulating reactions while the following down-regulating.
    The function also return the number of up-regulating reactions (that effectively splits the list)
    and a list of all the species/nodes contained in the network.   
    """
    assert filename[-4:].lower() == '.sif', "Currently parsing only SIF format. Sorry for the inconvenience!"
    activations = []
    inhibitions = []
    N_unknown = 0
    species = set()
    with open(filename, 'r') as fid: edges = fid.read().splitlines()
    for edge in edges:
        src, sign, tgt = edge.split("\t")
        sign = int(sign)
        if sign == 1:
            species.update([src, tgt])
            activations.append((src, tgt))
        elif sign == -1: # Giannis is also including -2
            species.update([src, tgt])
            inhibitions.append((src, tgt))
        else:
            N_unknown += 1
    if N_unknown > 0:
        print "There are %d reactions of unknown type in the reaction list that will be ignored." % N_unknown
    return activations + inhibitions, len(activations), list(species)


def read_measurements(filename, samples_row=True):
    """
    A function that reads the experimental measurements from a space delimited file.
    Measurements for each sample can be stored either in rows or in columns.
    In either case, the first element of the list must be the name of the measured node while the rest must the measured values.
    All signals must have the same number of measurements and NaNs should be used to fill in the missing values.
    The function returns a dictionary of with signals as its keys and lists containing the measured values as its values.
    """
    with open(filename, 'r') as fid:
        lines = fid.read().splitlines()
    if samples_row:
        Signals = lines[0].split()
        x_m = {spc: [] for spc in Signals}
        for line in lines[1:]:
            measurement = map(float, line.split())
            for i, spc in enumerate(Signals):
                x_m[spc].append(measurement[i])
    else:
        x_m = {}
        for line in lines:
            measurement = line.split()
            x_m[measurement[0]] = map(float, measurement[1:])
    # Check if all list have the same number of samples
    K = len(x_m.values()[0])
    error_msg = "Erroneous file format of %s.\nAll signals must have the same number of measurements.\n Use NaN to fill in the missing values." % filename
    assert all([len(m) == K for m in x_m.itervalues()]), error_msg

    return x_m
# Maybe I could use another format to read sparse data.


def read_stimuli(filename, samples_row=True):
    """
    A function that reads the stimulations imposed on each sample.
    Stimulation can be arranged either by rows or by columns.
    In either case, the first element of the list must the ne name of the stimulus and the rest must be a list of tuples.
    Each tuple contains 2 elements describing in which samples was the stimulus present or absent.
    The first element of each tuple is the number of the sample while the second is either 1 or 0 (present of absent)
    """
    with open(filename, 'r') as fid:
        lines = fid.read().splitlines()
    if samples_row:
        Stimuli = lines[0].split()
        B_0 = {spc: [] for spc in Stimuli}
        for k,line in enumerate(lines[1:]):
            imposed = map(float, line.split())
            for j, spc in enumerate(Stimuli):
                if not isnan(imposed[j]):
                    B_0[spc].append((k,imposed[j]))
    else:
        B_0 = {}
        for line in lines:
            imposed = line.split()
            B_0[imposed[0]] = [(k, val) for k, val in enumerate(map(float, imposed[1:])) if not isnan(val)]
    return B_0


def compute_node_flow(flow, reactions, species):
    """
    Computes the flow passing through each node.
    Inputs:
        - a list of flows passing through each reaction
        - a list of reactions (tuples)
        - a list of species
    Outputs:
        - a dictionary with the species as keys and the flow as values
    """
    flow_x = {spc: 0 for spc in species}
    for i, r in enumerate(reactions):
        flow_x[r[0]] += flow[i]
    return flow_x


def average_solutions(m, indices):
    """
    Computes the average value of a variable across multiple solutions.
    Inputs:
        - the cplex model m
        - a list containing the indices of the variables to be averaged. 
            Each element of the list must be a list containing the indices of the grouped variables
    Outputs:
        - a list of values corresponding to the indices specified
    """
    Var = [[0]*len(ind) for ind in indices]
    Nsol = m.solution.pool.get_num()
    for i, index in enumerate(indices):
        for s in range(Nsol):
            Var[i] = map(add, Var[i], m.solution.pool.get_values(s, index))
        Var[i] = [v / Nsol for v in Var[i]]

    return Var


def randomize_labels(m, labels, K, Var_ind, N_shuf, pool_size=10):
    """
    A function that takes as input a DCR CPLEX model, randomizes the labels to calculate an average random solution 
    and returns the average value of every Variable in the randomized datasets.
    Inputs:
        - a CPLEX model formalizing the DCR formulation
        - the list of labels (signals) that will be shuffled
        - the number of samples
        - the indices of the Variables that we wish to average. 
            Format: it's a list every element of which is a sublist containing the indices of a variable grouped by label
        - the number of reshuffling before returning the values
    Outputs:
        - a list of lists, where its sublist contains the average value of each variable (1 per label)
    """
    
    # Find the indices of the constrains that are affected by the Randomization grouped by the label
    # to make it more generic, maybe the list of indices of the constrains should be used as inputs
    Fit = [m.linear_constraints.get_indices(["Fit_{%s,%d}" % (l, k) for k in range(K)]) for l in labels]
    Flow = [m.linear_constraints.get_indices(["Source_{%s,%d}" % (l, k) for k in range(K)]) for l in labels]
    # Find the sense and the RHS of the affected Constraints grouped
    sns = [m.linear_constraints.get_senses(indx) for indx in Fit]
    rhs = [m.linear_constraints.get_rhs(indx) for indx in Fit]
    # Spread out the indices to make just 1 call per shuffling
    indices = [tmp for sublist in Fit for tmp in sublist] + [tmp for sublist in Flow for tmp in sublist]
    order = range(len(labels))
    Var = [[[0]*len(k) for k in ind] for ind in Var_ind]
    m.parameters.mip.limits.populate.set(pool_size) # too much??
    m.set_results_stream(None)
    m.set_log_stream(None)
    for n in range(N_shuf):
        shuffle(order)
        # Spread out the senses and the RHSs but keep the order within its label
        sns_rand = [tmp for i_rand in order for tmp in sns[i_rand]]
        rhs_rand = [tmp for i_rand in order for tmp in rhs[i_rand]]
        # the 2 constraints have the same senses and RHSs
        m.linear_constraints.set_senses(zip(indices, sns_rand*2))
        m.linear_constraints.set_rhs(zip(indices, rhs_rand*2))
        # Calculate new solution
        m.populate_solution_pool()
        for i, ind in enumerate(Var_ind):
            New = average_solutions(m, ind)
            Var[i] = [map(add, tmp, New[j]) for j, tmp in enumerate(Var[i])]
    Var = [[[n/N_shuf for n in group] for group in var] for var in Var]
    return Var


def write_dot(filename, edgelist, sources, sinks):
    """
    
    """
    with open(filename, 'w') as fid:
        fid.write("digraph {\n\n")

        fid.write("// Reactions (Edges)\n")
        for edge in edgelist:
            fid.write("\"%s\"->\"%s\"\t[penwidth=%.3f]\n" % edge)

        fid.write("\nStimuli (Sources)\n")
        for n in sources:
            fid.write("\"%s\"\t[style=filled, shape=oval, color=green]\n" % n)

        fid.write("\nSignals (Sinks)\n")
        for n in sinks:
            fid.write("\"%s\"\t[style=filled, shape=oval, color=blue]\n" % n)

        fid.write("\n}")


# Stand-Alone Functionality
def main(name, sifile, measurements, stimuli, MIP_options=None, sol=True,
         x_weights=None, B_weights=None, randomize=None, penwidth=(0.1, 5), out='', dot=True):


    # Read Inputs    
    reactions, N_up, species = read_PKN(sifile)
    x_m = read_measurements(measurements)
    B_0 = read_stimuli(stimuli)
    model_name = op.splitext(name)[0]

    # Print Problem Description
    print
    print "Initial Network (PKN):"
    total = len(reactions)
    print "Number of binary-reactions:\t%d (%d activations + %d inhibitions)" % (total, N_up, total-N_up)
    print "Number of species:\t%d" % len(species)
    print

    K = len(x_m.values()[0])
    print "Experimental Dataset:"
    print "Optimization will run using %d measured nodes over %d samples" % (len(x_m), K)
    print "Considering %d possible stimuli" %len(B_0)

    # Create Model
    print "Constructing CPLEX Model..."
    m = write_DCR_model(name, reactions, N_up, species, x_m, B_0, write_model=True, 
                        x_weights=x_weights, B_weights=B_weights)
    m.set_problem_name(op.basename(model_name))

    # Set Model Parameters
    if MIP_options: m.parameters.read_file(MIP_options)

    # Solve Problem
    m.populate_solution_pool()

    # Write Solution to files

    # Standard .sol file
    print
    print "Producing output files:"
    if sol:
        print "Writing .sol..."
        m.solution.pool.write(model_name + '_pool.sol')
    
    # Average Flow through a reaction across all solutions
    print "Calculating mean flow through each reaction..."
    flow_ind = [m.variables.get_indices([_F(i, k) for k in range(K)]) for i in range(N_up)]
    flow_z = average_solutions(m, flow_ind)
    print "Writing %s_flow_reactions.txt..." % model_name
    flow = map(lambda x: sum(x) / len(x), flow_z)
    edgelist = [(r[0], r[1], flow[i]) for i, r in enumerate(reactions[:N_up])]
    with open(model_name + '_flow_reactions.txt','w') as fid:
        for edge in sorted(edgelist, key=itemgetter(2), reverse=True):
            fid.write("%s\t%s\t%f\n" % edge)

    # Average Flow through a node across all solutions
    print "Calculating mean flow through each node..."
    flow_x = compute_node_flow(flow, reactions[:N_up], species)
    print "Writing %s_flow_nodes.txt..." % model_name
    with open(model_name + '_flow_nodes.txt','w') as fid:
        for f in sorted(flow_x.items(), key=itemgetter(1), reverse=True):
            fid.write("%s\t%f\n" % f)

    # Print Variables requested by the user
    if 'x' in out:
        filename = model_name + '_x.txt'
        print "Writing 'x' values in %s..." % filename
        ind = [m.variables.get_indices([_x(spc, k) for k in range(K)]) for spc in species]
        val = average_solutions(m, ind)
        with open(filename, 'w') as fid:
            for j, spc in enumerate(species):
                fid.write(spc+':\t' + '\t'.join(map(str, val[j])) + '\n')
    if 'd' in out:
        filename = model_name + '_d.txt'
        print "Writing 'd' values in %s..." % filename
        ind = [m.variables.get_indices([_d(spc, k) for k in range(K)]) for spc in species]
        val = average_solutions(m, ind)
        with open(filename, 'w') as fid:
            for j, spc in enumerate(species):
                fid.write(spc+':\t' + '\t'.join(map(str, val[j])) + '\n')
    if 'b' in out:
        filename = model_name + '_B.txt'
        print "Writing 'B' values in %s..." % filename
        ind = [m.variables.get_indices([_B(spc, k) for k in range(K)]) for spc in species]
        val = average_solutions(m, ind)
        with open(filename, 'w') as fid:
            for j, spc in enumerate(species):
                fid.write(spc+':\t' + '\t'.join(map(str, val[j])) + '\n')
    if 'z' in out:
        filename = model_name + '_z.txt'
        print "Writing 'z' values in %s..." % filename
        ind = [m.variables.get_indices([_z(i, k) for k in range(K)]) for i in range(total)]
        val = average_solutions(m, ind)
        with open(filename, 'w') as fid:
            for i, r in enumerate(reactions):
                fid.write("%s->%s:\t" % r)
                fid.write("\t".join(map(str, val[i])))
                fid.write('\n')
    if 'f' in out:
        filename = model_name + '_Fz.txt' 
        print "Writing 'Fz' values in %s..." % filename
        with open(filename, 'w') as fid:
            for i, r in enumerate(reactions):
                fid.write("%s->%s:\t" % r)
                fid.write("\t".join(map(str, flow_z[i])))
                fid.write('\n')
    if 'y' in out:
        filename = model_name + '_y.txt'
        print "Writing 'y' values in %s..." % filename
        ind = [m.variables.get_indices(_y(i)) for i in range(total)]
        val = average_solutions(m, ind)
        with open(filename, 'w') as fid:
            for i, r in enumerate(reactions):
                fid.write("{}\t{}\t{}\n".format(r[0], r[1], val[i]))

    # Plot Concensus Network
    if dot:
        Nsol = m.solution.pool.get_num()
        Fmin = min(flow, key=lambda x: x if x > 0 else float('inf'))
        Fmax = max(flow)
        Wmin, Wmax = penwidth
        interpolate = lambda f: Wmax - (Fmax-f) * (Wmax-Wmin) / (Fmax-Fmin)
        edgelist = []
        sinks = set()
        sources = set()
        for i, r in enumerate(reactions[:N_up]):
            if flow[i] > 0:
                edgelist.append((r[0], r[1], interpolate(flow[i])))
                if r[0] in x_m:
                    sinks.add(r[0])
                if r[1] in B_0 and r[1] not in sources: # the double checking is to avoid the cost of traversing all the solutions
                    B_ind = m.variables.get_indices([_B(r[1], k) for k in range(K)])
                    for s in range(Nsol):
                        if sum(m.solution.pool.get_values(s, B_ind)) > 0:
                            sources.add(r[1])
                            break
        write_dot(model_name + "_concensus_net.dot", edgelist, sources, sinks)

    # Find Flow in Random Network
    if randomize:
        N_shuf, rand_pool = randomize
        flow_rand = randomize_labels(m, x_m.keys(), K, [flow_ind], N_shuf, rand_pool)
        edgelist = [(r[0], r[1], flow_rand[i]) for i, r in enumerate(reactions[:N_up])]
        with open(model_name + "_flow_reactions_Rand.txt", 'w') as fid:
            for edge in sorted(edgelist, key=itemgetter(2), reverse=True):
                fid.write("%s\t%s\t%f\n" % edge)
        
        flow_x_rand = compute_node_flow(flow_rand, reactions[:N_up], species)
        with open(model_name + "_flow_node_Rand.txt", 'w') as fid:
            for f in sorted(flow_x_rand.items(), key=itemgetter(1), reverse=True):
                fid.write("%s\t%f\n" % f)


if __name__ == '__main__':
    
    import argparse
    # Parsing Arguments
    text1 = "A program for the identification of Signaling Networks based on \
    Boolean logic using the Melas Deep Causal Reasoning approach"
    text2 = "For more details see individual functions or read \
    https://github.com/mschubert/deep-causal-reasoning/blob/master/DCR_report.pdf"
    parser = argparse.ArgumentParser(description=text1, epilog=text2)

    # Necessary Arguments
    parser.add_argument("name", help="the name of the model (with a valid suffix mps or lp)")
    parser.add_argument('PKN', help='a .sif file containing a list of reactions')
    parser.add_argument('measurements', help='a tab delimited file containing \
    the names of the measured nodes (1st row) and the states they were \
    in each sample.\nEvery node must have the same number of measuremnts (use of NaNs is allowed)')
    parser.add_argument('inputs', help='a tab delimited file containing \
    the names of the possible stimuli (1st row) and the state they were \
    in each sample.\nEvery node must have the same number of measuremnts (use of NaNs is allowed)')

    # Optional Arguments
    parser.add_argument('-opt', '--options', dest='opt',
                        help='a file containing the options of the CPLEX model to be read by parameters.read_file method')
    parser.add_argument('-R', '--Randomize', nargs=2, type=int, dest='rand',
                        help="The program will randomize the labels of the measurements and output the results using the '_Rand' suffix.\
                        \nIf used, the number of random shuffles of the labels and the number of random solution per shuffle must be specified")
    parser.add_argument('-xw', '--XWeights', dest='xw',
                        help='either a single number (uniform weighting) or path to a file containing the weights of each \
                        variable x in each sample in the same order as they appear in the measurements, repeated as many times as the samples (single line).')
    parser.add_argument('-bw', '--BWeights', dest='bw',
                        help='either a single number (uniform weighting) or path to a file containing the weights of each \
                        variable B in each sample in the same order as they appear in the inputs, repeated as many times as the samples (single line).')
    parser.add_argument('-out', default='', help='Specifies with variables would be written in txt file \
    apart from the average flow through each node and reaction.\nThe variables are defined as \
    y->reactions, z->reactions state, f->flow through reactions, x->node state, b->stimulus state, d->node distance')
    parser.add_argument('-dot', action='store_false', help='Prevents the program from exporting a .dot file for graphviz')
    parser.add_argument('-sol', action='store_false', help='Prevents the program from exporting the solution in a .sol file')
    
    args = parser.parse_args()
    if args.xw:
        try:
            xw = int(args.xw)
        except ValueError:
            with open(args.xw, 'r') as fid:
                xw = map(float, fid.read().split())
    else:
        xw = None
    if args.bw:
        try:
            bw = int(args.bw)
        except ValueError:
            with open(args.bw, 'r') as fid:
                bw = map(float, fid.read().split())
    else:
        bw = None

    main(args.name, args.PKN, args.measurements, args.inputs, MIP_options=args.opt, sol=args.sol,
         x_weights=xw, B_weights=bw, randomize=args.rand, out=args.out.lower(), dot=args.dot)