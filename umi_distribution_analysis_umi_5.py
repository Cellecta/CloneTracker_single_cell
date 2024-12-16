#!/usr/bin/env python
# coding: utf-8

# In[332]:


import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import re


# In[398]:


data = pd.read_csv('barcode_assignment_umi_2.xls', sep='\t')


# In[399]:


data = data.fillna(0)


# In[400]:


pattern = pd.read_csv('tag.csv', header=None)
bc14 = pattern[0][:100]
bc30 = pattern[0][100:]


# In[401]:


bc14_name = pattern[1][:100]
bc30_name = pattern[1][100:]


# In[402]:


bc14_dict = dict(zip(bc14, bc14_name))
bc30_dict = dict(zip(bc30, bc30_name))


# In[484]:


def bc_type(bc14_str, bc30_str):
    if bc14_str == 0:
        bc_type = "bc14_only"
    if bc30_str == 0:
        bc_type = "bc30_only"
    if bc14_str != 0 and bc30_str != 0:
        bc_type = 'bc14_bc30'
    return bc_type


def barcode_name(bc14_str, bc30_str, bc14_dict, bc30_dict, R2_seq):
    if bc14_str != 0 and bc30_str != 0:
        barcode = bc14_dict.get(bc14_str)+str("_")+bc30_dict.get(bc30_str)
    if bc14_str == 0:
        barcode = R2_seq[25:39]+str("_")+bc30_dict.get(bc30_str)
    if bc30_str == 0:
        barcode = bc14_dict.get(bc14_str)+str("_")+R2_seq[43:73]
    return barcode

def final_assgined_barcode_func(barcode_list, umi_count_list):
    if max(umi_count_list) <5:
        final_assigned_barcode = 'NA'
    if (max(umi_count_list) >= 5) and (max(umi_count_list) >= (sum(umi_count_list)/2)) and len(umi_count_list)>1:
        final_assigned_barocode = barcode_list[0] 
    if len(umi_count_list) > 1 and (umi_count_list[1]/umi_count_list[0]) >= 0.7 and umi_count_list[0]>= 5:
        final_assigned_barcode = 'multi_barcode'
    if len(umi_count_list) > 1 and (umi_count_list[1]/umi_count_list[0]) < 0.7 and umi_count_list[0]>= 5:
        final_assigned_barcode = barcode_list[0]
    if len(umi_count_list)==1 and umi_count_list[0]>=5:
        final_assigned_barcode = barcode_list[0]
    return final_assigned_barcode

def final_assigned_type(final_assigned_barcode_str):
    pattern = re.compile(r'^bc14-(\d+)_bc30-(\d+)')
    if final_assigned_barcode_str == "NA":
        barcode_type = "Undetermined due to low UMI"
    elif final_assigned_barcode_str == "multi_barcode":
        barcode_type = 'multi_barcode'
    elif pattern.match(final_assigned_barcode_str):
        barcode_type = 'One Barcode'
    else:
        barcode_type = 'One Barcode with mutant'
    return barcode_type


def final_barcode_distribution(bc_type_list):
    # Example Counter
    my_counter = Counter(bc_type_list)
    # Extract labels and values from the Counter
    labels = list(my_counter.keys())
    sizes = list(my_counter.values())

    barcode_stat = pd.DataFrame(dict(my_counter), index=[0])
    barcode_stat.to_csv('barcode_stat.xls', sep='\t', index=False)

    plt.figure(figsize=(12, 6))

    # Create a pie chart
    plt.pie(sizes, labels=['']*len(labels), autopct='', startangle=90)

    # Add a title
    plt.title('Pie Chart of Categories')
    
    plt.gca().set_aspect('equal')  # Equal aspect ratio ensures that pie is drawn as a circle
    plt.gca().legend(labels = [f'{label}: {count}' for label, count in zip(labels, sizes)], loc='lower left', bbox_to_anchor=(1, 0.5), fontsize=10)  # Adjust font size and label position


    # Display the pie chart
    plt.savefig('barcode_assignment.png')
    
def umi_distribution_pie(umi_list):
     # Example Counter
    data = umi_list
    
    # Generate bins with integer boundaries
    bins = np.arange(np.floor(np.min(umi_list)), np.ceil(np.max(umi_list)) + 21, 20)

    # Create histogram
    hist, bin_edges = np.histogram(data, bins=bins)
    
    # Convert histogram to percentages
    percentages = (hist / sum(hist)) * 100


    # Create labels for each category
    labels = [f'{int(bin_edges[i])}-{int(bin_edges[i+1])}' for i in range(len(bin_edges)-1)]
    
    # Create a pie chart without showing labels
    plt.figure(figsize=(12, 6))
    patches, texts, autotexts = plt.pie(hist, labels=None, startangle=90, counterclock=False, autopct='', pctdistance=1.5)


    # Ensure the pie chart is drawn as a circle
    plt.axis('equal')
    
    # Create a legend manually with bin information
    legend_labels = [f'{label}: {hist[i]}' for i, label in enumerate(labels)]
    plt.legend(legend_labels, title='UMI bins', loc='center left', bbox_to_anchor=(0.8, 0.5), fontsize=10)

    # Display the pie chart
    plt.savefig('umi_distribution.png')


# In[417]:


data['bc_type'] = data.apply(lambda x: bc_type(x['bc14'], x['bc30']), axis=1)


# In[418]:


data_sorted = data.sort_values(by=['cell', 'umi'], ascending=[False, False])


# In[419]:


data_sorted['barcode_name'] = data_sorted.apply(lambda x: barcode_name(x['bc14'], x['bc30'], bc14_dict, bc30_dict, x['R2 sequnece']), axis=1)


# In[420]:


data_umi = data_sorted.groupby('cell')['umi'].apply(list).reset_index(name = 'umi')
data_barcode = data_sorted.groupby('cell')['barcode_name'].apply(list).reset_index(name = 'barcode')
data_final = data_umi.merge(data_barcode, on='cell', how='left')


# In[485]:


data_final['final_assigned_barcode'] = data_final.apply(lambda x: final_assgined_barcode_func(x['barcode'], x['umi']), axis=1)


# In[353]:


data_final['umi_count'] = data_final['umi'].apply(lambda x: x[0])


# In[354]:


data_final['barcode_type'] = data_final['final_assigned_barcode'].apply(lambda x: final_assigned_type(x))


# In[388]:


final_barcode_distribution(data_final['barcode_type'])


# In[393]:


umi_distribution_pie(data_final['umi_count'])


# In[395]:


data_final = data_final[['cell', 'final_assigned_barcode', 'umi_count', 'barcode_type']]


# In[397]:


data_final.to_csv('barcode_assignment_summary.xls', sep='\t', index=False)


# In[ ]:




