import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import xpress as xp
import pandas as pd
import itertools
import math
from preflibtools.instances import MatchingInstance
import random
from itertools import combinations
import gurobipy as gp
from gurobipy import GRB

EPSILON = 1e-5
M = 100000000

# define the distribution of age
age_distribution = {
    "0-15": 0.0119,
    "16-55": 0.5809,
    "55-64": 0.2756,
    "65+": 0.1316
}

# Generate random age ranges
def simulate_age_group(n,seed):
    age_groups = list(age_distribution.keys())
    probabilities = list(age_distribution.values())
    
    # Use random.choices to generate n age ranges according to the distribution
    random.seed(seed)
    simulated_ages = random.choices(age_groups, probabilities, k=n)
    return simulated_ages

def import_data(file, num = None):
    instance = MatchingInstance()
    wmd_file = file + ".wmd"
    instance.parse_file(wmd_file)

    if num is None or num >= instance.num_alternatives:
        num_of_nodes = instance.num_alternatives
    else:
        num_of_nodes = num

    edges = instance.edges()

    weights = np.zeros((num_of_nodes, num_of_nodes))
    for i,j,w in edges:
        if i <= num_of_nodes and j <= num_of_nodes:
            weights[i-1, j-1] = w
    weights

    dat_file = file + ".dat"
    df = pd.read_csv(dat_file, delimiter=',', header=0)

    s_p = df[df['%Pra'] >= 0.8]['Pair'].tolist()
    sensitive_patients = [sp for sp in s_p if sp <= num_of_nodes]

    seed = file[file.rfind('0') + 1:].rstrip(''.join(filter(lambda x: not x.isdigit(), file)))
    ages = simulate_age_group(num_of_nodes,seed)
    eldest_patients = [index+1 for index, value in enumerate(ages) if value in ['0.15','65+']]

    data_with_s_and_e = {'weights': weights, 'sensitive_patients': sensitive_patients, 'eldest_patients': eldest_patients}
    return data_with_s_and_e


# Build graph
def build_graph(matrix):
    # number of donor-patient pairs
    n = len(matrix[0])
    # name of nodes
    id = np.array(range(1, n + 1))

    # weights of edges
    weight_values = matrix

    # create a directed graph
    G = nx.DiGraph()

    # add node
    for i in range(n):
        G.add_node(id[i])

    # add edges with weights
    for j in range(n):
        for i in range(n):
            if weight_values[i][j] != 0:
                if weight_values[j][i] != 0:
                    G.add_edge(id[i], id[j], weight=weight_values[i][j], direction='both')
                else:
                    G.add_edge(id[i], id[j], weight=weight_values[i][j], direction='single')

    return G

# Show figure of graph
def show_KEP(G, selected_edges=None):
    # seed = 100
    # pos = nx.random_layout(G, seed=seed)
    # pos = nx.spring_layout(G)
    pos = nx.spring_layout(G, seed=42)
    # pos={1:(0,5),2:(10,15),3:(20,5),4:(20,-5),5:(10,-15),6:(0,-5)}

    # draw nodes
    nx.draw_networkx_nodes(G, pos, node_size=1000, node_color='skyblue', alpha=0.6)

    # draw nodes label
    nx.draw_networkx_labels(G, pos, font_size=15)

    # draw edges
    if selected_edges is not None:      # If we have selected edges
        for (u, v, d) in G.edges(data=True):
            if (u, v) in selected_edges:
                if d['direction'] == 'single':
                    nx.draw_networkx_edges(G, pos, edgelist=[(u, v)],
                                        edge_color='red', width=2
                                        )
                elif d['direction'] == 'both':
                    nx.draw_networkx_edges(G, pos, edgelist=[(u, v)],
                                        edge_color='red', width=2,
                                        connectionstyle='arc3, rad = 0.1')
                else:
                    raise ValueError("Unsupported direction type")
            else:
                if d['direction'] == 'single':
                    nx.draw_networkx_edges(G, pos, edgelist=[(u, v)],
                                        edge_color='black', width=2
                                        )
                elif d['direction'] == 'both':
                    nx.draw_networkx_edges(G, pos, edgelist=[(u, v)],
                                        edge_color='black', width=2,
                                        connectionstyle='arc3, rad = 0.1')
                else:
                    raise ValueError("Unsupported direction type")
    else:
        for (u, v, d) in G.edges(data=True):
            if d['direction'] == 'single':
                nx.draw_networkx_edges(G, pos, edgelist=[(u, v)],
                                    edge_color='black', width=2
                                    )
            elif d['direction'] == 'both':
                nx.draw_networkx_edges(G, pos, edgelist=[(u, v)],
                                    edge_color='black', width=2,
                                    connectionstyle='arc3, rad = 0.1')
            else:
                raise ValueError("Unsupported direction type")

    # def compute_label_pos(pos, edge, offset=0.1):
    #     u, v = edge
    #     x1, y1 = pos[u]
    #     x2, y2 = pos[v]
    #     x = x1 * (1 - offset) + x2 * offset
    #     y = y1 * (1 - offset) + y2 * offset
    #     return (x, y)

    # edge_label_pos = {key: compute_label_pos(pos, (key[0], key[1])) for key in edge_labels.keys()}


    # # add edges label
    # nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=15, label_pos=0.3)

    # show graph
    # plt.title("KEP Graph")
    plt.axis('off') 
    plt.show()

    return 0


def find_cycles(w,k=3):
    graph = {}
    num_of_nodes = len(w[0])
    name_of_nodes = range(1,num_of_nodes+1)
    for i,node_i in enumerate(name_of_nodes):
        neighbour = []
        for j,node_j in enumerate(name_of_nodes):
            if w[i,j] == 1:
                neighbour.append(node_j)
        graph[node_i] = neighbour
    
    cycles = []

    def dfs(node, start, path, visited):
        if len(path) > 4:  # cycles exceeding 3 are no longer checked
            return
        if node in visited:
            if node == start and len(path) > 1:
                cycles.append(path[:])
            return
        visited.add(node)
        path.append(node)
        for neighbor in graph.get(node, []):
            dfs(neighbor, start, path, visited)
        visited.remove(node)
        path.pop()

    for start_node in graph.keys():
        dfs(start_node, start_node, [], set())
    
    # Filter out cycles with length greater than 3
    short_cycles = [cycle for cycle in cycles if len(cycle) <= 3]
    
    # Using a collection to store unique cycles
    unique_cycles_set = set()
    unique_cycles = []
    
    for cycle in short_cycles:
        normalized_cycle = tuple(sorted(cycle))
        if normalized_cycle not in unique_cycles_set:
            unique_cycles_set.add(normalized_cycle)
            unique_cycles.append(cycle)

    return unique_cycles
    

def cycles_to_edges(cycles):
    # initial empty set for edges
    edges = set()

    for cycle in cycles:
        cycle_edges = [(cycle[i], cycle[i+1]) for i in range(len(cycle)-1)]
        cycle_edges.append((cycle[-1], cycle[0]))  # Add an edge from the last node to the first node to form a closed loop
        edges.update(cycle_edges)
    return edges


def xpress_solve_cycle(w,K):
    G = build_graph(w)

    # Xpress model

    # index
    number_of_nodes = len(w[0])
    nodes = list(range(number_of_nodes))
    name_of_nodes = range(1, number_of_nodes+1)

    cycles_K = find_cycles(w,K)
    weights_c = [len(cycle) for cycle in cycles_K]
    number_of_cycles_K = len(cycles_K)
    index_of_cycles_K = range(number_of_cycles_K)

    max_cycle_length = K

    # define the problem
    model_2 = xp.problem('KEP_simple_cycle')

    # define variables
    is_selected_c = np.array([xp.var(name='x_{0}'.format(i+1), vartype=xp.binary)
                                    for i in index_of_cycles_K], dtype=xp.npvar).reshape(number_of_cycles_K)

    model_2.addVariable(is_selected_c)

    # constraints

    # each node can exist in at most one cycle
    for node in name_of_nodes:
        model_2.addConstraint(xp.Sum([is_selected_c[i] for i, cycle in enumerate(cycles_K) if node in cycle]) <= 1)

    # objective
    model_2.setObjective(
        xp.Sum(weights_c[i]*is_selected_c[i] for i in index_of_cycles_K), sense = xp.maximize)

    # solve
    model_2.solve()

    model_2.write('KEP_simple_cycle','lp')

    x_df = model_2.getSolution(is_selected_c)

    x_result = [cycle for i, cycle in enumerate(cycles_K) if abs(x_df[i]-1)<=EPSILON]

    return x_result


def xpress_solve(w,K,method = 'cycle'):
    method_lower = method.lower()
    if method_lower == 'arc':
        # return xpress_solve_arc(w,K)
        return 0
    elif method_lower == 'cycle':
        return xpress_solve_cycle(w,K)
    else:
        raise ValueError("Parameter 'method' must be 'arc' or 'cycle'.")


# Combination attempt on probability principle

def find_all_solutions(all_cycles):

    def to_bitmask(c):
        """Convert c to a bit mask"""
        bitmask = 0
        for num in c:
            bitmask |= (1 << num)
        return bitmask

    def find_all_valid_selections(c):
        n = len(c)
        bitmasks = [to_bitmask(ci) for ci in c]
        all_selections = []
        
        def backtrack(start_index, current_bitmask, current_selection):
            if current_selection:
                all_selections.append(list(current_selection))
            
            for i in range(start_index, n):
                if current_bitmask & bitmasks[i] == 0:

                    current_selection.append(c[i])
                    new_bitmask = current_bitmask | bitmasks[i]
                    
                
                    backtrack(i + 1, new_bitmask, current_selection)
                    
                    
                    current_selection.pop()
        
        backtrack(0, 0, [])
        return all_selections


    all_valid_selections = find_all_valid_selections(all_cycles)

    def length_of_selection(selection):
        l = 0
        for cycle in selection:
            l += len(cycle)
        return l

    len_of_selection = [length_of_selection(selection) for selection in all_valid_selections]

    max_len = max(len_of_selection)

    used_selections = []

    for i,selection in enumerate(all_valid_selections):
        if len_of_selection[i] >= 0.5*max_len:
            used_selections.append(selection)

    return used_selections



def num_of_specific_patients_in_solution(feasible_solution,specific_patients=None):
    if specific_patients is None:
        num = sum([len(cycle) for i, cycle in enumerate(feasible_solution)])
    else:
        num = 0
        for c, cycle in enumerate(feasible_solution):
            for v, node in enumerate(cycle):
                if node in specific_patients:
                    num += 1
    
    return num


def is_included(patient, feasible_solution):
    for cycle in feasible_solution:
        if patient in cycle:
            return True
    
    return False


def all_selected_patients(feasible_solutions,patients):
    re = []

    for i,p in enumerate(patients):
        for s, fs in enumerate(feasible_solutions):
            if is_included(p,fs):
                re.append(p)
                break

    return re    

   
def xpress_solve_lottery_one_principle(data,Feasible_solutions,K,type):
        
    # 
    w = data.get('weights')
    sensitive_patients = data.get('sensitive_patients', None)
    
    G = build_graph(w)

    # Xpress model

    # index
    # vertex
    number_of_nodes = len(w[0])
    index_of_nodes = list(range(number_of_nodes))
    name_of_nodes = range(1, number_of_nodes+1)

    # cycles
    # cycles_K = find_cycles(w,K)
    # weights_c = [len(cycle) for cycle in cycles_K]
    # number_of_cycles_K = len(cycles_K)
    # index_of_cycles_K = range(number_of_cycles_K)

    # all feasible solutions
    feasible_solutions = Feasible_solutions
    weights_fs = [num_of_specific_patients_in_solution(fs) for s, fs in enumerate(feasible_solutions)]
    number_of_fs = len(feasible_solutions)
    index_of_fs = range(number_of_fs)

    # all selected patients
    selected_patients = all_selected_patients(feasible_solutions,name_of_nodes)
    number_of_sp = len(selected_patients)
    index_of_sp = range(number_of_sp)

    # IF principle

    # Aristotle’s
    # sensitive_patients = [1,2,3,4]
    weights_fs_sp = [num_of_specific_patients_in_solution(fs,sensitive_patients) for s, fs in enumerate(feasible_solutions)]


    # Parameter
    p = 1

    # define the problem
    model = xp.problem('KEP_lottery_one_principle')

    # define variables
    prop_s = np.array([xp.var(name='prop_{0}'.format(s+1), vartype=xp.continuous, lb=0, ub=1)
                                    for s in index_of_fs], dtype=xp.npvar).reshape(number_of_fs)
  
    u = xp.var(name='u', vartype=xp.continuous)

    model.addVariable(prop_s)
    model.addVariable(u)

    #### constraints

    # probability
    model.addConstraint(xp.Sum(prop_s[s] for s in index_of_fs) == 1)

    # function
    if type == 0:
        # utilitarian
        model.addConstraint(u == xp.Sum([prop_s[s]*weights_fs[s] for s in index_of_fs]))
    else:
        prop_v = np.array([xp.var(name='prop_v_{0}'.format(v+1), vartype=xp.continuous, lb=EPSILON, ub=1)
                                        for v in index_of_sp], dtype=xp.npvar).reshape(number_of_sp)
        model.addVariable(prop_v)

        for i,sp in enumerate(selected_patients):
            model.addConstraint(prop_v[i] == xp.Sum([prop_s[s] for s, fs in enumerate(feasible_solutions) if is_included(sp,fs)]))

        if type == 1:
            # IF
            mean_of_prop_v = xp.var(name='mean_of_prop_v', vartype=xp.continuous)
            d_v = np.array([xp.var(name='d_v_{0}'.format(v+1), vartype=xp.continuous)
                                        for v in index_of_sp], dtype=xp.npvar).reshape(number_of_sp)

            model.addVariable(mean_of_prop_v)
            model.addVariable(d_v)

            model.addConstraint(xp.Sum([prop_v[i] for i, sp in enumerate(selected_patients)]) == number_of_sp*mean_of_prop_v)

            for i,sp in enumerate(selected_patients):
                model.addConstraint(d_v[i] >= prop_v[i] - mean_of_prop_v)
                model.addConstraint(d_v[i] >= mean_of_prop_v - prop_v[i])

            model.addConstraint(xp.log(u) == - (xp.Sum([d_v[i]**p for i, sp in enumerate(selected_patients)])) ** (1/p))
            # model.addConstraint(u == (xp.Sum([xp.abs(prop_v[i] - mean_of_prop_v)**p for i, sp in enumerate(selected_patients)])) ** (1/p))
        elif type == 2:
            # Aristotle’s
            model.addConstraint(u == xp.Sum([prop_s[s]*weights_fs_sp[s] for s, fs in enumerate(feasible_solutions)]))
        elif type == 3:
            # Nash
            model.addConstraint(xp.log(u) == xp.Sum([ xp.log(prop_v[i]) for i, sp in enumerate(selected_patients)]))
        elif type == 4:
            # Rawlsian justice
            r = xp.var(name='r', vartype=xp.continuous)
            model.addVariable(r)

            for i, sp in enumerate(selected_patients):
                if prop_v[i] >= EPSILON:
                    model.addConstraint(r <= prop_v[i])
            
            model.addConstraint(u == r)
        else:
            raise ValueError("Unsupported parameter 'type' (0,1,2,3,4)")

    # objective
    model.setObjective(u, sense = xp.maximize)

    # solve
    model.solve()

    model.write('KEP_lottery_one_principle','lp')

    # x_df = model.getSolution(prop_s)
    # x_df_r = np.round(model.getSolution(prop_s),3)
    # # y_df = model.getSolution(prop_v)

    # optimal_value = model.getObjVal()

    # return optimal_value

    obj_max = model.getObjVal()

    return obj_max


def xpress_solve_lottery_one_fairness_principle(data,Feasible_solutions,K,lamda,type):
            
    # 
    w = data.get('weights')
    sensitive_patients = data.get('sensitive_patients', None)

    G = build_graph(w)

    # Xpress model

    # index
    # vertex
    number_of_nodes = len(w[0])
    index_of_nodes = list(range(number_of_nodes))
    name_of_nodes = range(1, number_of_nodes+1)

    # cycles
    # cycles_K = find_cycles(w,K)
    # weights_c = [len(cycle) for cycle in cycles_K]
    # number_of_cycles_K = len(cycles_K)
    # index_of_cycles_K = range(number_of_cycles_K)

    # all feasible solutions
    feasible_solutions = Feasible_solutions
    weights_fs = [num_of_specific_patients_in_solution(fs) for s, fs in enumerate(feasible_solutions)]
    number_of_fs = len(feasible_solutions)
    index_of_fs = range(number_of_fs)

    # all selected patients
    selected_patients = all_selected_patients(feasible_solutions,name_of_nodes)
    number_of_sp = len(selected_patients)
    index_of_sp = range(number_of_sp)

    # # principle
    number_of_principles = 4
    index_of_principles = range(number_of_principles)

    # Parameter
    # Aristotle’s
    # sensitive_patients = [1,2,3,4]
    weights_fs_sp = [num_of_specific_patients_in_solution(fs,sensitive_patients) for s, fs in enumerate(feasible_solutions)]

    p = 1
    lambda_k=lamda

    # define the problem
    model = xp.problem('KEP_lottery_one_fairness_principle')

    # define variables
    prop_s = np.array([xp.var(name='prop_{0}'.format(s+1), vartype=xp.continuous, lb=EPSILON, ub=1)
                                    for s in index_of_fs], dtype=xp.npvar).reshape(number_of_fs)
    prop_v = np.array([xp.var(name='prop_v_{0}'.format(v+1), vartype=xp.continuous, lb=EPSILON, ub=1)
                                    for v in index_of_sp], dtype=xp.npvar).reshape(number_of_sp)

    u = xp.var(name='u', vartype=xp.continuous)
    f = xp.var(name='f', vartype=xp.continuous)
    h = xp.var(name='H', vartype=xp.continuous)

    model.addVariable(prop_s)
    model.addVariable(prop_v)
    model.addVariable(u)
    model.addVariable(f)
    model.addVariable(h)

    #### constraints

    # probability
    model.addConstraint(xp.Sum(prop_s[s] for s in index_of_fs) == 1)

    # utilitarian
    model.addConstraint(u == xp.Sum([prop_s[s]*weights_fs[s] for s in index_of_fs]))

    # 
    for i,sp in enumerate(selected_patients):
        model.addConstraint(prop_v[i] == xp.Sum([prop_s[s] for s, fs in enumerate(feasible_solutions) if is_included(sp,fs)]))
    
    # f function
    if type == 1:
        # IF
        mean_of_prop_v = xp.var(name='mean_of_prop_v', vartype=xp.continuous)
        d_v = np.array([xp.var(name='d_v_{0}'.format(v+1), vartype=xp.continuous)
                                    for v in index_of_sp], dtype=xp.npvar).reshape(number_of_sp)

        model.addVariable(mean_of_prop_v)
        model.addVariable(d_v)

        model.addConstraint(xp.Sum([prop_v[i] for i, sp in enumerate(selected_patients)]) == number_of_sp*mean_of_prop_v)

        for i,sp in enumerate(selected_patients):
            model.addConstraint(d_v[i] == abs(prop_v[i] - mean_of_prop_v))
            # model.addConstraint(d_v[i] >= prop_v[i] - mean_of_prop_v)
            # model.addConstraint(d_v[i] >= mean_of_prop_v - prop_v[i])

        model.addConstraint(xp.log(f) == - (xp.Sum([d_v[i]**p for i, sp in enumerate(selected_patients)])) ** (1/p))
        # model.addConstraint(u == (xp.Sum([xp.abs(prop_v[i] - mean_of_prop_v)**p for i, sp in enumerate(selected_patients)])) ** (1/p))
    elif type == 2:
        # Aristotle’s
        model.addConstraint(f == xp.Sum([prop_s[s]*weights_fs_sp[s] for s, fs in enumerate(feasible_solutions)]))
    elif type == 3:
        # Nash
        model.addConstraint(xp.log(f) == xp.Sum([ xp.log(prop_v[i]) for i, sp in enumerate(selected_patients)]))
    elif type == 4:
        # Rawlsian justice
        r = xp.var(name='r', vartype=xp.continuous)
        model.addVariable(r)

        for i, sp in enumerate(selected_patients):
            if prop_v[i] >= EPSILON:
                model.addConstraint(r <= prop_v[i])
        
        model.addConstraint(f == r)
    else:
        raise ValueError("Unsupported parameter 'type' (1,2,3,4)")

    # h function
    model.addConstraint(h == u + lambda_k * f)

    # objective
    model.setObjective(h, sense = xp.maximize)

    # solve
    model.solve()

    model.write('KEP_lottery_one_fairness_principle','lp')

    # x_df = np.round(model.getSolution(prop_s),3)
    # y_df = np.round(model.getSolution(prop_v),3)
    # alpha_best = np.round(model.getObjVal(),3)
    # x_df = model.getSolution(prop_s)
    # y_df = model.getSolution(prop_v)
    # alpha_best = model.getObjVal()
                          
    # return feasible_solutions, x_df,y_df,alpha_best

    obj_max = model.getObjVal()

    return obj_max


def xpress_solve_lottery(data,Feasible_solutions,lambda_k,U_max,F_max,H_opt):
           
    # 
    w = data.get('weights')
    sensitive_patients = data.get('sensitive_patients', None)
    G = build_graph(w)

    # Xpress model

    # index
    # vertex
    number_of_nodes = len(w[0])
    index_of_nodes = list(range(number_of_nodes))
    name_of_nodes = range(1, number_of_nodes+1)

    # cycles
    # cycles_K = find_cycles(w,K)
    # weights_c = [len(cycle) for cycle in cycles_K]
    # number_of_cycles_K = len(cycles_K)
    # index_of_cycles_K = range(number_of_cycles_K)

    # all feasible solutions
    feasible_solutions = Feasible_solutions
    weights_fs = [num_of_specific_patients_in_solution(fs) for s, fs in enumerate(feasible_solutions)]
    number_of_fs = len(feasible_solutions)
    index_of_fs = range(number_of_fs)

    # all selected patients
    selected_patients = all_selected_patients(feasible_solutions,name_of_nodes)
    number_of_sp = len(selected_patients)
    index_of_sp = range(number_of_sp)

    # # principle
    number_of_principles = 4
    index_of_principles = range(number_of_principles)

    # Parameter
    p = 1
    u_max = U_max

    # sensitive_patients = [1,2,3,4]
    weights_fs_sp = [num_of_specific_patients_in_solution(fs,sensitive_patients) for s, fs in enumerate(feasible_solutions)]

    f_max = F_max
    
    u_max_divided_by_f_max = [u_max/v if v!= 0 else 0 for i,v in enumerate(f_max)]

    h_opt = H_opt
    # for i in index_of_principles:
    #     lamda = lambda_k*u_max_divided_by_f_max[i]
    #     h_opt.append(xpress_solve_lottery_one_fairness_principle(data,K,lamda,i+1))

    # define the problem
    model = xp.problem('KEP_lottery')

    # define variables
    prop_s = np.array([xp.var(name='prop_{0}'.format(s+1), vartype=xp.continuous, lb=0, ub=1)
                                    for s in index_of_fs], dtype=xp.npvar).reshape(number_of_fs)
    
    prop_v = np.array([xp.var(name='prop_v_{0}'.format(v+1), vartype=xp.continuous, lb=0, ub=1)
                                    for v in index_of_sp], dtype=xp.npvar).reshape(number_of_sp)
    
    u = xp.var(name='u', vartype=xp.continuous)

    mean_of_prop_v = xp.var(name='mean_of_prop_v', vartype=xp.continuous)
    
    f = np.array([xp.var(name='f_{0}'.format(k+1), vartype=xp.continuous)
                                    for k in index_of_principles], dtype=xp.npvar).reshape(number_of_principles)
    
    h = np.array([xp.var(name='H_{0}'.format(k+1), vartype=xp.continuous)
                                    for k in index_of_principles], dtype=xp.npvar).reshape(number_of_principles)
    
    z = xp.var(name='z', vartype=xp.continuous)

    model.addVariable(prop_s)
    model.addVariable(prop_v)
    model.addVariable(u)
    model.addVariable(mean_of_prop_v)
    model.addVariable(f)
    model.addVariable(h)
    model.addVariable(z)

    # Auxiliary variables
    # IF
    d_v = np.array([xp.var(name='d_v_{0}'.format(v+1), vartype=xp.continuous)
                                    for v in index_of_sp], dtype=xp.npvar).reshape(number_of_sp)

    # Aristotle’s

    # Nash

    # Rawlsian justice
    r = xp.var(name='r', vartype=xp.continuous)

    model.addVariable(d_v)
    model.addVariable(r)

    #### constraints

    # probability
    model.addConstraint(xp.Sum(prop_s[s] for s in index_of_fs) == 1)

    # utilitarian
    model.addConstraint(u == xp.Sum([prop_s[s]*weights_fs[s] for s in index_of_fs]))

    # # # IF
    for i,sp in enumerate(selected_patients):
        model.addConstraint(prop_v[i] == xp.Sum([prop_s[s] for s, fs in enumerate(feasible_solutions) if is_included(sp,fs)]))

    model.addConstraint(xp.Sum([prop_v[i] for i, sp in enumerate(selected_patients)]) == number_of_sp*mean_of_prop_v)

    for i,sp in enumerate(selected_patients):
        model.addConstraint(d_v[i] >= prop_v[i] - mean_of_prop_v)
        model.addConstraint(d_v[i] >= mean_of_prop_v - prop_v[i])

    model.addConstraint(xp.log(f[0]) == - (xp.Sum([d_v[i]**p for i, sp in enumerate(selected_patients)])) ** (1/p))
    # model.addConstraint(f[0] == (xp.Sum([xp.abs(prop_v[i] - mean_of_prop_v)**p for i, sp in enumerate(selected_patients)])) ** (1/p))

    # Aristotle’s
    model.addConstraint(f[1] == xp.Sum([prop_s[s]*weights_fs_sp[s] for s, fs in enumerate(feasible_solutions)]))

    # # Nash
    model.addConstraint(xp.log(f[2]) == xp.Sum([ xp.log(prop_v[i]) for i, sp in enumerate(selected_patients)]))
    
    # Rawlsian justice
    for i, sp in enumerate(selected_patients):
        if prop_v[i] >= EPSILON:
            model.addConstraint(r <= prop_v[i])
    
    model.addConstraint(f[3] == r)

    # h function
    for k in index_of_principles:
        model.addConstraint(h[k] == u + lambda_k*u_max_divided_by_f_max[k]*f[k])

    # # z function
    for k in index_of_principles:
        model.addConstraint(z * h_opt[k] <= h[k])
    
    # objective
    model.setObjective(z, sense = xp.maximize)

    # solve
    model.solve()

    model.write('KEP_lottery','lp')

    x_df = np.round(model.getSolution(prop_s),3)
    y_df = np.round(model.getSolution(prop_v),3)

    alpha_best = np.round(model.getObjVal(),3)

    u_with_current_solution = np.dot(x_df,weights_fs)

    pof = (u_max - u_with_current_solution)/u_max

    return x_df,alpha_best,u_with_current_solution,pof
    # return feasible_solutions, x_df,y_df,alpha_best,np.round(u_max,3),np.round(f_max,3),np.round(h_opt,3)


# Determining results for cycle formulas

def number_of_specific_people_in_cycle(cycle,specific_people):
    return sum(1 for c in cycle if c in specific_people)

# Consider only one principle to obtain the maximum function
def xpress_solve_deterministic_one_principle(data,K,type):
    
    # 获取权重矩阵w
    w = data.get('weights')
    sensitive_patients = data.get('sensitive_patients', None)
    eldest_patients = data.get('eldest_patients', None)

    G = build_graph(w)

    # Xpress model

    # index
    number_of_nodes = len(w[0])
    nodes = list(range(number_of_nodes))
    name_of_nodes = range(1, number_of_nodes+1)

    cycles_K = find_cycles(w,K)
    weights_c = [len(cycle) for cycle in cycles_K]
    number_of_cycles_K = len(cycles_K)
    index_of_cycles_K = range(number_of_cycles_K)

    # Parameter
    # sensitive_patients = [1,2,3,4]
    if sensitive_patients is not None:
        weights_s = [number_of_specific_people_in_cycle(cycle,sensitive_patients) for cycle in cycles_K]
    # eldest_patients = [1,6]
    if eldest_patients is not None:
        weights_e = [number_of_specific_people_in_cycle(cycle,eldest_patients) for cycle in cycles_K]

    # define the problem
    model = xp.problem('KEP_deterministic')

    # define variables
    is_selected_c = np.array([xp.var(name='x_{0}'.format(i+1), vartype=xp.binary)
                                    for i in index_of_cycles_K], dtype=xp.npvar).reshape(number_of_cycles_K)
    
    u = xp.var(name='u', vartype=xp.continuous)

    model.addVariable(is_selected_c)
    model.addVariable(u)


    # constraints

    # each node can exist in at most one cycle
    for node in name_of_nodes:
        model.addConstraint(xp.Sum([is_selected_c[i] for i, cycle in enumerate(cycles_K) if node in cycle]) <= 1)

    # f function
    if type == 0:
        # utilitarian function
        model.addConstraint(u == xp.Sum(weights_c[i]*is_selected_c[i] for i in index_of_cycles_K))
    elif type == 1:
        # sensitive patients
        model.addConstraint(u == xp.Sum(weights_s[i]*is_selected_c[i] for i in index_of_cycles_K))
    elif type == 2:
        # eldest patients
        model.addConstraint(u == xp.Sum(weights_e[i]*is_selected_c[i] for i in index_of_cycles_K))
    else:
        raise ValueError("Unsupported parameter 'type' (0,1,2)")

    model.setObjective(u, sense = xp.maximize)

    # solve
    model.solve()

    model.write('model_deterministic_one_principle','lp')

    x_df = model.getSolution(is_selected_c)
    
    x_result = [cycle for i, cycle in enumerate(cycles_K) if abs(x_df[i]-1)<=EPSILON]

    obj_max = model.getObjVal()

    return obj_max

# Considering efficiency + a certain fairness principle, obtain the maximum function
def xpress_solve_deterministic_one_fairness_principle(data,K,lamda,type):
    
    # 获取权重矩阵w
    w = data.get('weights')
    sensitive_patients = data.get('sensitive_patients', None)
    eldest_patients = data.get('eldest_patients', None)

    G = build_graph(w)

    # Xpress model

    # index
    number_of_nodes = len(w[0])
    nodes = list(range(number_of_nodes))
    name_of_nodes = range(1, number_of_nodes+1)

    cycles_K = find_cycles(w,K)
    weights_c = [len(cycle) for cycle in cycles_K]
    number_of_cycles_K = len(cycles_K)
    index_of_cycles_K = range(number_of_cycles_K)

    # Parameter
    # sensitive_patients = [1,2,3,4]
    if sensitive_patients is not None:
        weights_s = [number_of_specific_people_in_cycle(cycle,sensitive_patients) for cycle in cycles_K]
    # eldest_patients = [1,6]
    if eldest_patients is not None:
        weights_e = [number_of_specific_people_in_cycle(cycle,eldest_patients) for cycle in cycles_K]

    # define the problem
    model = xp.problem('KEP_deterministic_one_fairness_principle')

    # define variables
    is_selected_c = np.array([xp.var(name='x_{0}'.format(i+1), vartype=xp.binary)
                                    for i in index_of_cycles_K], dtype=xp.npvar).reshape(number_of_cycles_K)
    
    u = xp.var(name='u', vartype=xp.continuous)
    f = xp.var(name='f', vartype=xp.continuous)
    h = xp.var(name='H', vartype=xp.continuous)

    model.addVariable(is_selected_c)
    model.addVariable(u)
    model.addVariable(f)
    model.addVariable(h)

    # constraints

    # each node can exist in at most one cycle
    for node in name_of_nodes:
        model.addConstraint(xp.Sum([is_selected_c[i] for i, cycle in enumerate(cycles_K) if node in cycle]) <= 1)
    
    # utilitarian function
    model.addConstraint(u == xp.Sum(weights_c[i]*is_selected_c[i] for i in index_of_cycles_K))
    
    # f function
    if type == 1:
        # sensitive patients
        model.addConstraint(f == xp.Sum(weights_s[i]*is_selected_c[i] for i in index_of_cycles_K))
    elif type == 2:
        # eldest patients
        model.addConstraint(f == xp.Sum(weights_e[i]*is_selected_c[i] for i in index_of_cycles_K))
    else:
        raise ValueError("Unsupported parameter 'type' (0,1,2)")

    # h function
    model.addConstraint(h == u + lamda*f)
    
    # objective
    model.setObjective(h, sense = xp.maximize)

    # solve
    model.solve()

    model.write('KEP_deterministic_one_fairness_principle','lp')

    x_df = model.getSolution(is_selected_c)
    
    x_result = [cycle for i, cycle in enumerate(cycles_K) if abs(x_df[i]-1)<=EPSILON]

    obj_max = model.getObjVal()

    return obj_max

# Consider all
def xpress_solve_deterministic(data,Cycles_K,lambda_k,U_max,F_max,H_opt):
    
    # 获取权重矩阵w
    w = data.get('weights')
    sensitive_patients = data.get('sensitive_patients', None)
    eldest_patients = data.get('eldest_patients', None)

    G = build_graph(w)

    # Xpress model

    # index
    number_of_nodes = len(w[0])
    nodes = list(range(number_of_nodes))
    name_of_nodes = range(1, number_of_nodes+1)

    cycles_K = Cycles_K
    weights = [len(cycle) for cycle in cycles_K]
    number_of_cycles_K = len(cycles_K)
    index_of_cycles_K = range(number_of_cycles_K)

    number_of_principles = 2
    index_of_principles = range(number_of_principles)

    # Parameter
    u_max = U_max
    # sensitive_patients = [1,2,3,4]
    weights_s = [number_of_specific_people_in_cycle(cycle,sensitive_patients) for cycle in cycles_K] 
    # eldest_patients = [1,6]
    weights_e = [number_of_specific_people_in_cycle(cycle,eldest_patients) for cycle in cycles_K] 

    f_max = F_max

    u_max_divided_by_f_max = [u_max/v if v!= 0 else 0 for i,v in enumerate(f_max)]

    h_opt = H_opt

    # define the problem
    model = xp.problem('KEP_deterministic')

    # define variables
    is_selected_c = np.array([xp.var(name='x_{0}'.format(i+1), vartype=xp.binary)
                                    for i in index_of_cycles_K], dtype=xp.npvar).reshape(number_of_cycles_K)
    
    u = xp.var(name='u', vartype=xp.continuous)
    
    f = np.array([xp.var(name='f_{0}'.format(k+1), vartype=xp.continuous)
                                    for k in index_of_principles], dtype=xp.npvar).reshape(number_of_principles)
    
    h = np.array([xp.var(name='H_{0}'.format(k+1), vartype=xp.continuous)
                                    for k in index_of_principles], dtype=xp.npvar).reshape(number_of_principles)

    z = xp.var(name='z', vartype=xp.continuous)

    model.addVariable(is_selected_c)
    model.addVariable(u)
    model.addVariable(f)
    model.addVariable(h)
    model.addVariable(z)

    # constraints

    # each node can exist in at most one cycle
    for node in name_of_nodes:
        model.addConstraint(xp.Sum([is_selected_c[i] for i, cycle in enumerate(cycles_K) if node in cycle]) <= 1)
    
    # utilitarian function
    model.addConstraint(u == xp.Sum(weights[i]*is_selected_c[i] for i in index_of_cycles_K))
    
    # f1 function
    model.addConstraint(f[0] == xp.Sum(weights_s[i]*is_selected_c[i] for i in index_of_cycles_K))
    # f2 function
    model.addConstraint(f[1] == xp.Sum(weights_e[i]*is_selected_c[i] for i in index_of_cycles_K))
    
    # h function
    for k in index_of_principles:
        model.addConstraint(h[k] == u + lambda_k*u_max_divided_by_f_max[k]*f[k])
    
    # z function
    for k in index_of_principles:
        model.addConstraint(z * h_opt[k] <= h[k])
    # objective
    model.setObjective(z, sense = xp.maximize)
    
    # solve
    model.solve()

    model.write('model_deterministic','lp')

    x_df = model.getSolution(is_selected_c)
    
    x_result = [cycle for i, cycle in enumerate(cycles_K) if abs(x_df[i]-1)<=EPSILON]

    alpha_best = model.getObjVal()

    u_with_current_solution = np.dot(weights,x_df)

    pof = (u_max - u_with_current_solution) / u_max

    return x_result,alpha_best,u_with_current_solution,pof
    


def lambda_test_deterministic(data,data_name):
    # Parameter
    # lambda_test = [0,1,10,100,1000,6252]
    lambda_test = np.arange(0.0, 4.1, 0.1)
    K = 3
    w = data.get('weights')
    sensitive_patients = data.get('sensitive_patients', None)
    eldest_patients = data.get('eldest_patients', None)

    # return value
    x_result = []
    alpha_best = []
    u_with_s = []
    pof = []

    # find all cycles with K
    cycles_K = find_cycles(w,K)

    # calculate u_max and f_max
    u_max = xpress_solve_deterministic_one_principle(data,K,0)

    f_max = []
    for i in (0,1):
        f_max.append(xpress_solve_deterministic_one_principle(data,K,i+1))

    for l in lambda_test:
        
        # Calculate h_opt
        h_opt = []
        for i in (0,1):
            lamda = l * u_max / f_max[i]
            h_opt.append(xpress_solve_deterministic_one_fairness_principle(data,K,lamda,i+1))


        x, a, u_w_s, p = xpress_solve_deterministic(data,cycles_K,l,u_max,f_max,h_opt)

        

        x_result.append(x)
        alpha_best.append(a)
        u_with_s.append(u_w_s)
        pof.append(p)

    f_n = 'data/deterministic_result-' + data_name + '.npz'
    
    np.savez(f_n, x = lambda_test, y = pof, z = alpha_best)

    # return x_result,alpha_best,u_with_s,u_max,f_max,pof
    # return u_max,f_max,u_max_divided_by_f_max
    return 0


def lambda_test_lottery(data,data_name):
    # Parameter
    # lambda_test = [0,1,10,100,1000,6252]
    lambda_test = np.arange(1,4.1,111)
    K = 3
    w = data.get('weights')
    sensitive_patients = data.get('sensitive_patients', None)
    # eldest_patients = data.get('eldest_patients', None)

    # return value
    x_result = []
    alpha_best = []
    u_with_s = []
    pof = []
    
    # find all cycles with K
    cycles_K = find_cycles(w,K)

    # find all feasible solution
    feasible_solutions = find_all_solutions(cycles_K)

    # calculate u_max and f_max
    # u_max = xpress_solve_deterministic_one_principle(data,K,0)
    u_max = xpress_solve_lottery_one_principle(data,feasible_solutions,K,0)

    f_max = []
    for i in range(0,4):
        f_max.append(xpress_solve_lottery_one_principle(data,feasible_solutions,K,i+1))

    for l in lambda_test:
        
        # Calculate h_opt
        h_opt = []
        for i in range(0,4):
            lamda = l * u_max / f_max[i]
            h_opt.append(xpress_solve_lottery_one_fairness_principle(data,feasible_solutions,K,lamda,i+1))

        x, a, u_w_s, p = xpress_solve_lottery(data,feasible_solutions,l,u_max,f_max,h_opt)

        x_result.append(x)
        alpha_best.append(a)
        u_with_s.append(u_w_s)
        pof.append(p)
    
    f_n = 'data/lottery_result-' + data_name + '.npz'
    
    np.savez(f_n, x = lambda_test, y = pof, z = np.round(alpha_best,3))

    return 0


def main_deterministic():
    # Calculate result with different lambda in some instances 

    list_index = list(range(31, 40)) + list(range(71, 81))

    for i in list_index:
        file = 'data/kidney/00036-000000' + str(i)

        data_with_s_and_e = import_data(file)
        data_name = file[-14:]

        lambda_test_deterministic(data_with_s_and_e,data_name)
    
    # draw plot
    lambda_test = np.arange(0.0, 4.1, 0.1)
    # lambda_test = [0,1,10,100,1000,6252]
    alpha_best = []
    pof = []
    for i in list_index:
    # for i in list(range(80, 81)):
        file = 'data/deterministic_result-00036-000000' + str(i) + '.npz'
        data = np.load(file)
        alpha_best.append(data['z'])
        pof.append(data['y'])

    
    # 32-38,71,73,80 pof-lambda
    plt.figure(figsize=(10, 6))
    list1 = [1,2,3,4,5,6,7,9,11,18]

    for i in list1:
        plt.plot(lambda_test, pof[i], label=list_index[i])

    plt.xlabel(r'$\lambda$')
    plt.ylabel('PoF')
    plt.legend()
    plt.grid(True)
    plt.show()

    # 33,35-38,71,73 alpha_best-lambda
    plt.figure(figsize=(10, 6))
    list1 = [2,4,5,6,7,9,11]

    for i in list1:
        plt.plot(lambda_test, alpha_best[i], label=list_index[i])

    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$\alpha_{best}$')
    plt.legend()
    plt.grid(True)
    plt.show()
    

    # 31,39,72,74-79 pof-lambda
    plt.figure(figsize=(10, 6))
    list2 = [0,8,10,12,13,14,15,16,17]

    for i in list2:
        plt.plot(lambda_test, pof[i], label=list_index[i])

    plt.xlabel(r'$\lambda$')
    plt.ylabel('PoF')
    plt.legend()
    plt.grid(True)
    plt.show()
    

def main_lottery():
    # Calculate result with different lambda in some instances 

    list_index = ['02', '09', '10','37','72','74','77']

    for i in list_index:
        file = 'data/kidney/00036-000000' + str(i)

        data_with_s_and_e = import_data(file, 25)
        data_name = file[-14:]

        # lambda_test_lottery(data_with_s_and_e,data_name)

    
    # draw plot
    lambda_test = np.arange(0.1, 4.1, 0.1)
    # lambda_test = [0,1,10,100,1000,6252]
    alpha_best = []
    pof = []
    for i in list_index:
        file = 'data/lottery_result-00036-000000' + str(i) + '.npz'
        data = np.load(file)
        alpha_best.append(data['z'])
        pof.append(data['y'])
    

    for i,index in enumerate(list_index):
        file = 'data/kidney/00036-000000' + str(index)

        data_with_s_and_e = import_data(file, 25)
        data_name = file[-14:]

        weight_values = data_with_s_and_e.get('weights')
        all_cycles_K = find_cycles(weight_values)
        feasible_solutions = find_all_solutions(all_cycles_K)


        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

        # pof-lambda
        mask = alpha_best[i]>=-1
        x = lambda_test[mask]
        y = pof[i][mask]

        ax1.plot(x, y, label=index)
        ax1.set_xlabel(r'$\lambda$')
        ax1.set_ylabel('Pof')


        count = [0]*len(all_cycles_K)
        for j,solution in enumerate(feasible_solutions):
            for k,cycle in enumerate(solution):
                index = all_cycles_K.index(cycle)
                count[index] += 1

        ax2.bar(range(len(all_cycles_K)), count)
        ax2.set_xlabel('Cycle')
        ax2.set_ylabel('Frequency')

        plt.tight_layout()  
        plt.show()
    

    
    # alpha_best-lambda
    plt.figure(figsize=(10, 6))

    for i,j in enumerate(list_index):
        mask = alpha_best[i]>=-1
        x = lambda_test[mask]
        z = alpha_best[i][mask]
        plt.plot(x, z, label=j)

    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$\alpha_{best}$')
    plt.legend()
    plt.grid(True)
    plt.show()

 

if __name__ == '__main__':
    # main_deterministic()
    main_lottery()


