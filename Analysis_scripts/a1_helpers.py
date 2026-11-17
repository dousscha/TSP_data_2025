"""Helper functions used by A1_DataAnalyse_cleaned_up.ipynb."""

import numpy as np
import pandas as pd
import networkx as nx

def get_time_taken_sec(t):

    import time
    from datetime import datetime
        
    # Truncate the microseconds to 6 digits
    t = t[:-1]
    
    # Parse using datetime
    try :
        parsed_time = datetime.strptime(t, '%H:%M:%S.%f')
        # Calculate the total number of seconds
        seconds = parsed_time.hour * 3600 + parsed_time.minute * 60 + parsed_time.second + parsed_time.microsecond / 1e6
    except :
        print('error',t)
        parsed_time = datetime.strptime(t, '%H:%M:%S')
        seconds = parsed_time.hour * 3600 + parsed_time.minute * 60 + parsed_time.second 
    
    
    
    return seconds

def extract_strategies(row):
    '''from survey get first, second and third prefered strategy'''
    try:
        list_of_strategies = [s.strip() for s in row.split(',')]
        first_strat    = list_of_strategies[0] if len(list_of_strategies) > 0 else None
        second_strat   = list_of_strategies[1] if len(list_of_strategies) > 1 else None
        third_strat    = list_of_strategies[2] if len(list_of_strategies) > 2 else None

    except Exception as e:

        print(f"Error processing row, invalid answer: {row} -> {e}")
        return None, None, None  # Ensure three values are always returned 

    return first_strat, second_strat, third_strat


def get_percentages(grp):
    total = grp.sum()    
    percentage = grp / total * 100  # Actual percentage calculation
    return pd.DataFrame({'percentage': percentage})


def parse_values(x):
    item_ = x.split('-')[1]
    if item_ in ['10', '9']:
       return 'L'
    elif item_ in ['2', '3']:
       return 'R'
    else:
       return 'other'


def directed_to_undirected_with_summed_weights(G):
    
    UG = nx.Graph()  # Create an empty undirected graph

    for u, v, d in G.edges(data=True):

        weight = d.get("weight", 1)  # Default weight is 1 if not provided
        
        if UG.has_edge(u, v):

            UG[u][v]["weight"] += weight  # Sum weights if edge already exists
        else:
            UG.add_edge(u, v, weight=weight)

    return UG


def graph_primer(all_nodes, group, ax, sortby='trial_nb',bi = True,sizing=10): 
    
    '''set the parameters needed to build a weighted and directed networkX graph to describe route following
    all_nodes = list of nodes name e.g. ['A','B','C']
    data = dataframe column with 'trials_nb' and one column with 'sequence'
    ax = ax of figure
    '''
    
    index_df = [(p1, p2) for p1 in all_nodes for p2 in all_nodes]
    
    #create a bi directed graph
    #if bi : DG = nx.DiGraph()

    DG = nx.DiGraph()

    #create empty sequence for all consecutive routes
    list_pairs             = pd.DataFrame(0,index=index_df,columns=['weight'])

    list_pairs.index       = pd.MultiIndex.from_tuples(list_pairs.index)

    list_pairs.index.names = ['F1','F2']

    all_seq=list()
    
    for bout_id,bout in group.groupby(sortby):
 
        #get sequence of the bout
        raw_seq = bout.Node_ID.values 
        seq = [item for item in raw_seq if item in all_nodes]

        #put all in one
        all_seq.append(seq)

        # get the travelled segments between visited nodes out of the sequence of nodes'visits from that bout
        for f1,f2 in zip(seq[:-1], seq[1:]):
            
            #bee flew f1 to f2
            tup = (f1,f2)
            #add weight to seg f1-> f2
            list_pairs.loc[tup].weight += 1
            
    #flatten the list to get all nodes        
    all_seq = [item for sublist in all_seq for item in sublist]
    
    #list_pairs= list_pairs[list_pairs['weight'] != 0]       
    #print(all_seq)

    for edge in list_pairs.index:

        DG.add_edge(edge[0], edge[1], weight = list_pairs.loc[edge[0],edge[1]].weight)
        
      
    if bi : edges, weights = zip(*nx.get_edge_attributes(DG,'weight').items()) 

    else : 

        DG = directed_to_undirected_with_summed_weights(DG)
        edges,weights = zip(*nx.get_edge_attributes(DG,'weight').items()) 

        edges, weights = zip(*nx.get_edge_attributes(DG,'weight').items()) 


     
    list_pairs_nozero= list_pairs[list_pairs['weight'] != 0]       
        
        
    G_nozero = nx.DiGraph()

    for edge in list_pairs_nozero.index:

        G_nozero.add_edge(edge[0],edge[1],weight=list_pairs_nozero.loc[edge[0],edge[1]].weight)
        
    if bi : edges_nozero,weights_nozero = zip(*nx.get_edge_attributes(G_nozero,'weight').items()) 

    else : 

        G_nozero = directed_to_undirected_with_summed_weights(G_nozero)
        edges_nozero,weights_nozero = zip(*nx.get_edge_attributes(G_nozero,'weight').items()) 

    
    node_sizes = [all_seq.count(i)*sizing for i in all_nodes]
    
    return DG,edges,weights,node_sizes,edges_nozero,G_nozero,weights_nozero


def create_dic(data_array,condition,scaling):

    chest_pos = data_array.loc[(data_array.cond == condition) & (data_array.scale == scaling)].set_index('node').drop(columns=['scale','cond','optimal'])
    pos = chest_pos.transpose().to_dict('list')

    return pos,chest_pos.index


def percentage_scaler(X):
    return [ x*100/np.sum(X) for x in  X]


def get_transition_df(G):

    data = [(u, v, d.get('weight', 1)) for u, v, d in G.edges(data=True)]
    return pd.DataFrame(data, columns=['From', 'To', 'Weight'])


def Euclidean_Dist(df1, df2, cols=['x','y']):
    
    return np.linalg.norm(df1[cols].values - df2[cols].values,
                   axis=1)
