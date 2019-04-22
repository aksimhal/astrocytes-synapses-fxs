"""FXS Data Analysis"""

import os
import copy
import scipy
import scipy.stats
import numpy as np
import pandas as pd
from skimage import exposure
from skimage import color
from at_synapse_detection import dataAccess as da

class mouse_region:
    """This class encompasses a region within a mouse"""

    def __init__(self, data_list, num_queries, z_spans, section_name, layer_name):
        """

        Parameters
        ---------------
        data_list : list
        num_queries : int
        z_spans : list
        section_name : str
        layer_name : str

        """
        self.section_name = section_name
        self.layer_name = layer_name

        z_spans = list(range(0, len(z_spans)))

        all_queries = {}

        for n in range(0, num_queries):
            key = 'q' + str(n)

            query_data = {}

            start_z = 0
            for z in range(0, len(z_spans)):
                query_data[z + 1] = data_list[n + z * num_queries]

            all_queries[key] = query_data

        self.all_queries = all_queries

    def get_small_synapse(self, query_id):
        """return the number of 'small synapses', 1 slice - 2 slice"""

        return self.all_queries[query_id][1] - self.all_queries[query_id][2]

    def get_medium_synapse(self, query_id):
        """return the number of 'medium synapses', 2 slice - 3 slice"""

        return self.all_queries[query_id][2] - self.all_queries[query_id][3]


class fxs_mouse:
    """
    Contains the information associated with each mouse.  
    Making this a class because the mouse images are formatted 
    differently and this will hopefully make things less confusing.

    Properties 
    --------------

    """

    def __init__(self, name):
        """
        """
        self.name = name

    def create_mouse_fn(self, name, fn, mouse_type, layer_order, num_queries, query_names, region_names, layer_names, z_spans):
        """
        Parameters
        ---------------
        fn : str - excel sheet filepath 
        layer_order : str - 'forward' means F000 is layer 1. 'backward' means F000 is layer 4
        num_queries : int - number of queries 
        z_spans : list[int] - [1, 2, 3]

        """
        self.name = name
        self.layer_order = layer_order
        self.num_queries = num_queries
        self.z_spans = z_spans
        self.query_names = query_names
        self.region_names = region_names
        self.type = mouse_type

        region_data = []
        df = pd.read_excel(fn)
        r_step = num_queries * len(z_spans)

        if layer_order == 'forward':
            ind_list = [0, 1 * r_step, 2 * r_step, 3 * r_step]
        else:
            ind_list = [3 * r_step, 2 * r_step, 1 * r_step, 0]

        rn_ind = 0
        for region_ind in ind_list:
            print(region_ind, region_ind + r_step)
            data_list = df.iloc[region_ind:region_ind + r_step, 0].values
            region = mouse_region(
                data_list, num_queries, z_spans, region_names[rn_ind], layer_names[rn_ind])
            region_data.append(region)
            rn_ind = rn_ind + 1

        self.region_data = region_data


    def create_mouse_df(self, name, df, mouse_type, layer_order, num_queries, query_names, region_names, layer_names, z_spans):
        """
        Parameters
        ---------------
        """

        self.name = name
        self.layer_order = layer_order
        self.num_queries = num_queries
        self.z_spans = z_spans
        self.query_names = query_names
        self.region_names = region_names
        self.type = mouse_type

        region_data = []
        #df = pd.read_excel(fn)
        r_step = num_queries * len(z_spans)

        if layer_order == 'forward':
            ind_list = [0, 1 * r_step, 2 * r_step, 3 * r_step]
        else:
            ind_list = [3 * r_step, 2 * r_step, 1 * r_step, 0]

        rn_ind = 0
        for region_ind in ind_list:
            print(region_ind, region_ind + r_step)
            data_list = df.iloc[region_ind:region_ind + r_step, 0].values
            region = mouse_region(
                data_list, num_queries, z_spans, region_names[rn_ind], layer_names[rn_ind])
            region_data.append(region)
            rn_ind = rn_ind + 1

        self.region_data = region_data

    def average_layers(self, z_span, query_id=None):

        if query_id == None:
            # get list of keys
            keylist = self.region_data[0].all_queries.keys()
            average_dict = {}

            if z_span == 'small':
                for key in keylist:
                    foo = 0
                    for n in range(0, 4):
                        foo = foo + self.region_data[n].get_small_synapse(key)
                    average_dict[key] = np.mean(foo)
            else:
                for key in keylist:
                    foo = 0
                    for n in range(0, 4):
                        foo = foo + \
                            self.region_data[n].all_queries[key][z_span]
                    average_dict[key] = np.mean(foo)

            return average_dict

        else:
            if z_span == 'small':
                sum_density = 0
                for n in range(0, 4):
                    sum_density = sum_density + \
                        self.region_data[n].get_small_synapse(query_id)
                average_density = sum_density / 4
            else:
                sum_density = 0
                for n in range(0, 4):
                    sum_density = sum_density + \
                        self.region_data[n].all_queries[query_id][z_span]
                average_density = sum_density / 4

            return average_density


def create_df(mouse_2ss, row_labels):
    """This takes a fxs_mouse and formats it into yi format """
    region_names = ['F000', 'F001', 'F002', 'F003']
    size_labels = ['small', 'one', 'two', 'three']
    column_labels = []
    
    for sizelabel in size_labels:
        for regionname in region_names:
            columnlabel = sizelabel + '-' + regionname
            column_labels.append(columnlabel)

    df = pd.DataFrame(np.nan, index=row_labels, columns=column_labels)
    for query_itr, query_id in enumerate(row_labels):

        size_itr = 0
        q = query_id.lower()

        df.iloc[query_itr, size_itr * 4 +
                0] = mouse_2ss.region_data[0].get_small_synapse(q)
        df.iloc[query_itr, size_itr * 4 +
                1] = mouse_2ss.region_data[1].get_small_synapse(q)
        df.iloc[query_itr, size_itr * 4 +
                2] = mouse_2ss.region_data[2].get_small_synapse(q)
        df.iloc[query_itr, size_itr * 4 +
                3] = mouse_2ss.region_data[3].get_small_synapse(q)

        for size_itr in range(1, 4):
            df.iloc[query_itr, size_itr * 4 +
                    0] = mouse_2ss.region_data[0].all_queries[q][size_itr]
            df.iloc[query_itr, size_itr * 4 +
                    1] = mouse_2ss.region_data[1].all_queries[q][size_itr]
            df.iloc[query_itr, size_itr * 4 +
                    2] = mouse_2ss.region_data[2].all_queries[q][size_itr]
            df.iloc[query_itr, size_itr * 4 +
                    3] = mouse_2ss.region_data[3].all_queries[q][size_itr]
    
    df.name = mouse_2ss.name + '-' + mouse_2ss.type
    
    return df

def create_ratio_df(mouse_2ss, row_labels):
    """This takes a fxs_mouse and formats it into yi format """
    region_names = ['F000', 'F001', 'F002', 'F003']
    size_labels = ['small', 'one', 'two', 'three']
    column_labels = []

    for sizelabel in size_labels:
        for regionname in region_names:
            columnlabel = sizelabel + '-' + regionname
            column_labels.append(columnlabel)

    df = pd.DataFrame(np.nan, index=row_labels, columns=column_labels)
    for query_itr, query_id in enumerate(row_labels):

        q = query_id.lower()

        for size_itr in range(0, 4):
            df.iloc[query_itr, size_itr * 4 +
                    0] = mouse_2ss.region_data[0].all_queries[q][size_itr]
            df.iloc[query_itr, size_itr * 4 +
                    1] = mouse_2ss.region_data[1].all_queries[q][size_itr]
            df.iloc[query_itr, size_itr * 4 +
                    2] = mouse_2ss.region_data[2].all_queries[q][size_itr]
            df.iloc[query_itr, size_itr * 4 +
                    3] = mouse_2ss.region_data[3].all_queries[q][size_itr]
    
    df.name = mouse_2ss.name + '-' + mouse_2ss.type
    
    return df


def create_editedquery_df(mouse_2ss): 
    """This takes a fxs_mouse and formats it into yi format.
        It also ensures query uniqueness-ish """
    
    
    row_labels = ['Q0', 'Q1', 'Q2', 'Q3']
    region_names = ['F000', 'F001', 'F002', 'F003']
    size_labels = ['small', 'one', 'two', 'three']
    column_labels = []
    
    for sizelabel in size_labels: 
        for regionname in region_names: 
            columnlabel = sizelabel+'-'+regionname
            column_labels.append(columnlabel)
    
    df = pd.DataFrame(np.nan, index=row_labels, columns=column_labels)
    
    ####### QUERY 0 : PSD SYNAPSIN ONLY 
    n=0
    size_ind = 0    

    for layer_ind in range(0, 4): 
        df.iloc[n, size_ind*4+layer_ind] = mouse_2ss.region_data[layer_ind].get_small_synapse('q0') \
            - mouse_2ss.region_data[layer_ind].get_small_synapse('q1') \
            - mouse_2ss.region_data[layer_ind].get_small_synapse('q2') \
            + 2*mouse_2ss.region_data[layer_ind].get_small_synapse('q3')
    
    

    for size_ind in range(1, 4): 
        for layer_ind in range(0, 4): 
            df.iloc[n, size_ind*4+layer_ind] = mouse_2ss.region_data[layer_ind].all_queries['q0'][size_ind] \
                - mouse_2ss.region_data[layer_ind].all_queries['q1'][size_ind] \
                - mouse_2ss.region_data[layer_ind].all_queries['q2'][size_ind] \
                + 2*mouse_2ss.region_data[layer_ind].all_queries['q3'][size_ind]

        
    
    ######### QUERY 1 : PSD SYNAPSIN V1 ONLY 
    n=1
    size_ind = 0        
    for layer_ind in range(0, 4): 
        df.iloc[n, size_ind*4+layer_ind] = mouse_2ss.region_data[layer_ind].get_small_synapse('q1') \
            - mouse_2ss.region_data[layer_ind].get_small_synapse('q3')
    

    for size_ind in range(1, 4): 
        for layer_ind in range(0, 4): 
            df.iloc[n, size_ind*4+layer_ind] = mouse_2ss.region_data[layer_ind].all_queries['q1'][size_ind] \
                - mouse_2ss.region_data[layer_ind].all_queries['q3'][size_ind]
        
        
    
    ######### QUERY 2 : PSD SYNAPSIN V2 ONLY 

    n=2
    size_ind = 0        
    for layer_ind in range(0, 4): 
        df.iloc[n, size_ind*4+layer_ind] = mouse_2ss.region_data[layer_ind].get_small_synapse('q2') \
            - mouse_2ss.region_data[layer_ind].get_small_synapse('q3')
    

    for size_ind in range(1, 4): 
        for layer_ind in range(0, 4): 
            df.iloc[n, size_ind*4+layer_ind] = mouse_2ss.region_data[layer_ind].all_queries['q2'][size_ind] \
                - mouse_2ss.region_data[layer_ind].all_queries['q3'][size_ind]
        
    
    
    ######### QUERY 3 : PSD SYNAPSIN VGLUT1 VGLUT2 
    n=3
    size_ind = 0    
    for layer_ind in range(0, 4): 
        df.iloc[n, size_ind*4+layer_ind] = mouse_2ss.region_data[layer_ind].get_small_synapse('q3')
    

    for size_ind in range(1, 4): 
        for layer_ind in range(0, 4): 
            df.iloc[n, size_ind*4+layer_ind] = mouse_2ss.region_data[layer_ind].all_queries['q3'][size_ind]
        
    return df
    


def update_mouse(mouse, fn, layer_order): 
    """update mouse data structure"""
    q = 'q3'
    df = pd.read_excel(fn)

    r_step = 3

    if layer_order == 'forward':
        ind_list = [0, 1 * r_step, 2 * r_step, 3 * r_step]
    else:
        ind_list = [3 * r_step, 2 * r_step, 1 * r_step, 0]

    rn_ind = 0
    for n, region_ind in enumerate(ind_list):
        
        data_list = df.iloc[region_ind:region_ind + r_step, 0].values
        mouse.region_data[n].all_queries[q][1] = data_list[0]
        mouse.region_data[n].all_queries[q][2] = data_list[1]
        mouse.region_data[n].all_queries[q][3] = data_list[2]
        
    return mouse

def add_query_to_mouse(mouse, fn, num_queries_in_fn, q_keylist, layer_order): 
    """Add a query to the mouse data structure. 4 layers and 3 sizes are assumed. 
    
    Parameters
    ---------------
    mouse : mouse object 
    fn : str
    q_key : list of str
    layer_order : str
    
    Return
    ---------------
    mouse : mouse object 
    """

    df = pd.read_excel(fn)

    r_step = 3*num_queries_in_fn

    if layer_order == 'forward':
        ind_list = [0, 1 * r_step, 2 * r_step, 3 * r_step]
    else:
        ind_list = [3 * r_step, 2 * r_step, 1 * r_step, 0]


    for n, region_ind in enumerate(ind_list):
        data_list = df.iloc[region_ind:region_ind + r_step, 0].values
        itr = 0 

        # Add keys to the dictionary
        for q_key in q_keylist: 
            mouse.region_data[n].all_queries[q_key] = {} 

        for size_key in range(1, 4): 
            for q_n, q_key in enumerate(q_keylist): 
                mouse.region_data[n].all_queries[q_key][size_key] = data_list[itr]
                itr = itr + 1
            
    return mouse


def average_mice(mouse_list, mouse_name, row_labels): 
    """Average multiple mice
    
    Parameters
    -------------
    mouse_list
    syanpse_size
    
    Returns
    -------------
    average_mouse
    """
    region_names = ['F000', 'F001', 'F002', 'F003']
    z_spans = [0, 1, 2, 3]

    sizedict = {'0': -1, '1': -1, '2': -1, '3': -1}
    regiondict = {'F000': sizedict.copy(), 'F001': sizedict.copy(), 'F002': sizedict.copy(), 'F003': sizedict.copy()}
    average_mouse = {'Q0': copy.deepcopy(regiondict), 'Q1': copy.deepcopy(regiondict), 'Q2': copy.deepcopy(regiondict), 'Q3': copy.deepcopy(regiondict), 
                        'Q4': copy.deepcopy(regiondict), 'Q5': copy.deepcopy(regiondict), 'Q6': copy.deepcopy(regiondict)}

    num_of_mice = len(mouse_list)
    
    for qID in row_labels: 

        for r_n, regionID in enumerate(region_names): 
            for sizeID in z_spans:
                average_density = 0 
                for mouse in mouse_list: 
                    #print(mouse.name)
                    #print(mouse.region_data[r_n].all_queries.keys())

                    average_density = average_density + mouse.region_data[r_n].all_queries[qID.lower()][sizeID]
                
                average_density = average_density/num_of_mice
                average_mouse[qID][regionID][str(sizeID)] = average_density
   
    
    return average_mouse
    

def average_mouse_to_df(average_mouse, row_labels, sizeID, mouse_name): 
    """Convert nested dictionary to dataframe
    Parameters
    -------------
    average_mouse
    
    Return
    -------------
    """
    region_names = ['F000', 'F001', 'F002', 'F003']
    column_labels = []
    sizeID = str(sizeID)
    for regionname in region_names:
        columnlabel = sizeID + '-' + regionname
        column_labels.append(columnlabel)
    df = pd.DataFrame(np.nan, index=row_labels, columns=column_labels)
    
    for q_n, qID in enumerate(row_labels):
        for r_n, regionID in enumerate(region_names):
            df.iloc[q_n, r_n] = average_mouse[qID][regionID][sizeID]
            
    df.name = mouse_name
            
    return df



    
    
def get_p_value(queryID, regionNum, sizeNum, ko_mouse_list, wt_mouse_list, pairkey): 
    """
    
    """
    ko_values = [] 
    for mouse in ko_mouse_list: 
        val = mouse.region_data[regionNum].all_queries[queryID.lower()][sizeNum]
        ko_values.append(val)
    
    wt_values = [] 
    for mouse in wt_mouse_list: 
        val = mouse.region_data[regionNum].all_queries[queryID.lower()][sizeNum]
        wt_values.append(val)

    if pairkey == 'paired':
        val = scipy.stats.ttest_rel(wt_values, ko_values)
    elif pairkey =='unpaired': 
        val = scipy.stats.ttest_ind(wt_values, ko_values)

    pval = val.pvalue

    return pval
        
def create_pval_df(query_list, sizeNum, ko_mouse_list, wt_mouse_list, df_name, pairkey): 
    """
    Parameters
    ---------------
    query_list : list of strings 
    sizeNum : int
    ko_mouse_list : list of mouse objects 
    
    Returns
    ---------------
    df: dataframe with p values 
    """
    
    region_names = ['F000', 'F001', 'F002', 'F003']
    column_labels = []
    sizeID = str(sizeNum)
    for regionname in region_names:
        columnlabel = sizeID + '-' + regionname
        column_labels.append(columnlabel)
        
    df = pd.DataFrame(np.nan, index=query_list, columns=column_labels)
    
    for q_n, queryID in enumerate(query_list): 
        for regionNum, regionID in enumerate(region_names): 
            pval = get_p_value(queryID, regionNum, sizeNum, ko_mouse_list, wt_mouse_list, pairkey)
            df.iloc[q_n, regionNum] = pval
    df.name = df_name
    return df
            

def write_dfs_to_excel(df_list, sheets, file_name, spaces=1):
    """
    Write multiple dataframes to an excel file

    Parameters
    --------------
    df_list : list of dataframes
    sheets : str - name of sheet in excel
    file_name : str
    spaces : int - number of rows to skip, default = 1
    """

    writer = pd.ExcelWriter(file_name, engine='xlsxwriter')
    row = 1
    col = 0
    for dataframe in df_list:
        dataframe.to_excel(writer, sheet_name=sheets, startrow=row, startcol=col)
        
        if hasattr(dataframe, 'name'):
            print(dataframe.name)    
            worksheet = writer.sheets[sheets]
            worksheet.write_string(row-1, col, dataframe.name)
        row = row + len(dataframe.index) + spaces + 2
        
    
    writer.save()



def compute_small_synapses(mouse): 
        """
        compute small synapses
        """
        for region in mouse.region_data:
            for key in region.all_queries.keys(): 
                region.all_queries[key][0] = region.get_small_synapse(key)
        
        return mouse 


def compute_medium_synapses(mouse): 
        """
        compute medium synapses
        """
        for region in mouse.region_data:
            for key in region.all_queries.keys(): 
                region.all_queries[key][2] = region.get_medium_synapse(key)
        
        return mouse 



def divide_mouse(mouse1, mouse2, query_list): 
    """mouse2/mouse1"""
    
    result = copy.deepcopy(mouse1)
    
    
    for r_n in range(0, 4): 
        for qID in query_list: 
            for size_n in range(0, 4):
                a = mouse1.region_data[r_n].all_queries[qID.lower()][size_n]
                b = mouse2.region_data[r_n].all_queries[qID.lower()][size_n]
                c = b/a
                result.region_data[r_n].all_queries[qID.lower()][size_n] = c
                #print(c)
    return result

def divide_mouse_by_query(mouse1, blob_mouse2, blob_queryID): 
    """mouse2/mouse1
    mouse1 and mouse2 have different number of queries 
    """
    
    result = copy.deepcopy(mouse1)
    
    query_list = mouse1.region_data[0].all_queries.keys()
    
    for r_n in range(0, 4): 
        for qID in query_list: 
            for size_n in range(0, 4):
                a = mouse1.region_data[r_n].all_queries[qID][size_n]
                b = blob_mouse2.region_data[r_n].all_queries[blob_queryID][size_n] #this holds the query constant
                c = b/a
                result.region_data[r_n].all_queries[qID][size_n] = c
                #print(c)
    return result

def divide_mouse_by_query2(mouse1, blob_mouse2, blob_queryID): 
    """mouse2/mouse1
    mouse1 and mouse2 have different number of queries 
    """
    
    result = copy.deepcopy(mouse1)
    
    query_list = mouse1.region_data[0].all_queries.keys()
    
    for r_n in range(0, 4): 
        for qID in query_list: 
            for size_n in range(0, 4):
                a = mouse1.region_data[r_n].all_queries[qID][size_n]
                b = blob_mouse2.region_data[r_n].all_queries[blob_queryID][size_n] #this holds the query constant
                c = a/b
                result.region_data[r_n].all_queries[qID][size_n] = c
                #print(c)
    return result

def average_mouse_layers(mouse, query_list): 
    """Average multiple mice
    
    Parameters
    -------------
    mouse_list
    syanpse_size
    
    Returns
    -------------
    average_mouse
    """
    region_names = ['F000', 'F001', 'F002', 'F003']
    z_spans = [0, 1, 2, 3]

    sizedict = {'0': -1, '1': -1, '2': -1, '3': -1}
    regiondict = {'average': sizedict.copy()}
    
    average_mouse = {'Q0': copy.deepcopy(regiondict), 'Q1': copy.deepcopy(regiondict), 'Q2': copy.deepcopy(regiondict), 'Q3': copy.deepcopy(regiondict), 
                        'Q4': copy.deepcopy(regiondict), 'Q5': copy.deepcopy(regiondict), 'Q6': copy.deepcopy(regiondict)}

    
    for qID in query_list: 
        average_mouse[qID]['average_L23'] = {} 
        average_mouse[qID]['average_L1'] = {}
        average_mouse[qID]['average_L4'] = {} 
        for sizeID in z_spans: 
            layer2_3_avg = (mouse.region_data[1].all_queries[qID.lower()][sizeID] + mouse.region_data[2].all_queries[qID.lower()][sizeID])/2
            layer1 = mouse.region_data[0].all_queries[qID.lower()][sizeID]
            layer4 = mouse.region_data[3].all_queries[qID.lower()][sizeID]
            
            avg_layer = (layer1 + layer2_3_avg + layer4)/3 
            
            average_mouse[qID]['average'][str(sizeID)] = avg_layer
            average_mouse[qID]['average_L23'][str(sizeID)] = layer2_3_avg
            average_mouse[qID]['average_L1'][str(sizeID)] = layer1
            average_mouse[qID]['average_L4'][str(sizeID)] = layer4

    
    return average_mouse


def average_layer_mice(mouse_list, mouse_name, query_list): 
    """Average multiple mice of the same type
    
    Parameters
    -------------
    mouse_list
    syanpse_size
    
    Returns
    -------------
    average_mouse
    """
    
    z_spans = [0, 1, 2, 3]
    sizedict = {'0': -1, '1': -1, '2': -1, '3': -1}
    regiondict = {'average': sizedict.copy(), 'average_std_error': sizedict.copy()}
    average_mouse = {'Q0': copy.deepcopy(regiondict), 'Q1': copy.deepcopy(regiondict), 'Q2': copy.deepcopy(regiondict), 'Q3': copy.deepcopy(regiondict), 
                        'Q4': copy.deepcopy(regiondict), 'Q5': copy.deepcopy(regiondict), 'Q6': copy.deepcopy(regiondict)}

    num_of_mice = len(mouse_list)
    # average_mouse = {} 
    key_list = ['average_L1', 'average_L4', 'average_L23', 'average']

    
    for qID in query_list: 
        for key in key_list: 
            average_mouse[qID][key] = {} 
            average_mouse[qID][key+"_std_error"] = {} 
            average_mouse[qID][key+"_data"] = {} 
            for sizeID in z_spans:
                average_density = 0 
                data_list = [] 
                for mouse in mouse_list: 
                    data_list.append(mouse[qID][key][str(sizeID)])

                average_density = np.mean(data_list)
                avg_std = np.std(data_list)
                average_mouse[qID][key][str(sizeID)] = average_density
                average_mouse[qID][key+"_std_error"][str(sizeID)] = avg_std/np.sqrt(len(data_list))
                average_mouse[qID][key+"_data"][str(sizeID)] = data_list
    return average_mouse


def average_layer_mice_to_df(average_mouse, row_labels, sizeID, mouse_name): 
    """Convert nested dictionary to dataframe
    Parameters
    -------------
    average_mouse
    
    Return
    -------------
    """
    
    column_labels = ['layers averaged']
    sizeID = str(sizeID)

    df = pd.DataFrame(np.nan, index=row_labels, columns=column_labels)
    
    for q_n, qID in enumerate(row_labels):
        df.iloc[q_n, 0] = average_mouse[qID]['average'][sizeID]
            
    df.name = mouse_name
            
    return df


def create_layer_avg_pval_df(query_list, layer_key, sizeNum, ko_mouse_list, wt_mouse_list, df_name, pairkey): 
    """
    Parameters
    ---------------
    query_list : list of strings 
    sizeNum : int
    ko_mouse_list : list of mouse objects 
    
    Returns
    ---------------
    df: dataframe with p values 
    """
    
    column_labels = []
    sizeID = str(sizeNum)
    columnlabel = sizeID + '-' + 'average'
    column_labels.append(columnlabel)
        
    df = pd.DataFrame(np.nan, index=query_list, columns=column_labels)
    
    for q_n, queryID in enumerate(query_list): 
        pval = get_average_layer_p_value(queryID, sizeNum, layer_key, ko_mouse_list, wt_mouse_list, pairkey)
        df.iloc[q_n, 0] = pval
    df.name = df_name
    return df
            

def get_average_layer_p_value(queryID, sizeNum, layer_key, ko_mouse_list, wt_mouse_list, pairkey): 
    """
    
    """
    ko_values = [] 
    for mouse in ko_mouse_list: 
        val = mouse[queryID][layer_key][str(sizeNum)]
        ko_values.append(val)
    
    wt_values = [] 
    for mouse in wt_mouse_list: 
        val = mouse[queryID][layer_key][str(sizeNum)]
        wt_values.append(val)

    if pairkey == 'paired':
        val = scipy.stats.ttest_rel(wt_values, ko_values)
    elif pairkey =='unpaired': 
        val = scipy.stats.ttest_ind(wt_values, ko_values)

    pval = val.pvalue

    return pval

def create_layer_avg_pval_dict(ko_mouse_list, wt_mouse_list, query_list, pairkey): 
    """
    Parameters
    ---------------
    sizeNum : int
    ko_mouse_list : list of mouse objects 
    
    Returns
    ---------------
    """
    key_list = ['average_L1', 'average_L4', 'average_L23', 'average']
    
    pval_dict = {} 
    zspans = [0, 1, 2, 3]
    for qkey in query_list: 
        pval_dict[qkey] = {} 
        for layer_key in key_list: 
            pval_dict[qkey][layer_key] = {} 
            for zkey in zspans: 
                zkey = str(zkey)
                pval = get_average_layer_p_value(qkey, zkey, layer_key, ko_mouse_list, wt_mouse_list, pairkey)
                pval_dict[qkey][layer_key][zkey] = pval

    return pval_dict    


def get_yticks(y_data): 
    yticks24 = np.array([0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4])

    yticks20 = np.array([0, 0.4, 0.8, 1.2, 1.6, 2.0])
    yticks16 = np.array([0, 0.4, 0.8, 1.2, 1.6])
    yticks14 = np.array([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4])
    yticks12 = np.array([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
    yticks10 = np.array([0, 0.2, 0.4, 0.6, 0.8, 1])
    yticks08 = np.array([0, 0.2, 0.4, 0.6, 0.8])
    yticks05 = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    yticks04 = np.array([0, 0.1, 0.2, 0.3, 0.4])
    yticks025 = np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25])
    yticks02 = np.array([0, 0.05, 0.1, 0.15, 0.2])
    yticks0125 = np.array([0, 0.025, 0.05, 0.075, 0.1, 0.125])
    yticks003 = np.array([0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03])
    yticks0025 = np.array([0, 0.005, 0.01, 0.015, 0.02, 0.025])

    if y_data>=1.6: 
        y = 2.2
        yticks = yticks24
    elif y_data >= 1.4 and y_data<1.6: 
        y = 1.8
        yticks = yticks20
    elif y_data>= 1.2 and y_data<1.4: 
        y = 1.4
        yticks = yticks16
    elif y_data>= 1 and y_data<1.2: 
        y = 1.3
        yticks = yticks14
    elif y_data>=0.8 and y_data<1: 
        y = 1.1
        yticks = yticks12
    elif y_data >= 0.6 and y_data<0.8: 
        y = 0.9
        yticks = yticks10
    elif y_data >=0.4 and y_data<0.6: 
        y = 0.7 
        yticks = yticks08
    elif y_data >=0.3 and y_data<0.4: 
        y = 0.45
        yticks = yticks05
    elif y_data >=0.2 and y_data < 0.3: 
        y = 0.35
        yticks = yticks04
    elif y_data >=0.15 and y_data <0.2: 
        y = 0.225
        yticks = yticks025
    elif y_data >=0.1 and y_data < 0.15: 
        y =0.175
        yticks = yticks02
    elif y_data >= 0.025 and y_data < 0.1: 
        y = 0.1125
        yticks = yticks0125
    elif y_data >= 0.02 and y_data < 0.025: 
        y = 0.0275
        yticks = yticks003
    elif y_data > 0 and y_data < 0.02: 
        y = 0.0225
        yticks = yticks0025

    
    return y, yticks



def colorize(image, hue, saturation=1,v=1):
    ### Add color of the given hue to an RGB greyscale image.
    hsv = color.rgb2hsv(image)
    hsv[:, :, 2] *= v
    hsv[:, :, 1] = saturation
    hsv[:, :, 0] = hue
    return color.hsv2rgb(hsv)

def setup_slice(fn, z): 
    vol = da.imreadtiff(fn)
    img = vol[:, :, z]
    p2, p98 = np.percentile(img, (5, 98)) #contrast stretching
    img_rescale = exposure.rescale_intensity(img, in_range=(p2, p98))
    img_rgb = color.gray2rgb(img_rescale)
    
    return img_rgb

def update_vglut2_queries(mouse, fn, layer_order): 
    """update mouse data structure"""
    q2 = 'q2'
    q3 = 'q3'
    
    df = pd.read_excel(fn)


    if layer_order == 'forward':
        q2_L1 = 0
        q3_L1 = 1
        q2_L2 = 2
        q3_L2 = 3
        q2_L3 = 4
        q3_L3 = 5
        q2_L4 = 6
        q3_L4 = 7
        
    else:
        q2_L1 = 6
        q3_L1 = 7
        q2_L2 = 4
        q3_L2 = 5
        q2_L3 = 2
        q3_L3 = 3
        q2_L4 = 0
        q3_L4 = 1

    data_list = df.iloc[0:8, 0].values
    mouse.region_data[0].all_queries[q2][1] = data_list[q2_L1]
    mouse.region_data[0].all_queries[q3][1] = data_list[q3_L1]
    mouse.region_data[1].all_queries[q2][1] = data_list[q2_L2]
    mouse.region_data[1].all_queries[q3][1] = data_list[q3_L2]
    mouse.region_data[2].all_queries[q2][1] = data_list[q2_L3]
    mouse.region_data[2].all_queries[q3][1] = data_list[q3_L3]
    mouse.region_data[3].all_queries[q2][1] = data_list[q2_L4]
    mouse.region_data[3].all_queries[q3][1] = data_list[q3_L4]
    
        
    return mouse