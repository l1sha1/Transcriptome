# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 23:26:56 2022

@author: Timofei
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.font_manager import FontProperties

file = pd.read_csv('C:\Data\gene_used_norm.txt') #менять путь чтения и имя чтения
new_file = []

for string in range(len(file['name'])):
    
    new_string = []
    for simb in range(len(file['name'][string])):
        if file['name'][string][simb] == '$' or file['name'][string][simb] == '[' or file['name'][string][simb] == ']' or file['name'][string][simb] == ' ':
            continue        
        else:
            new_string.append(file['name'][string][simb])
            #new_string.extend()
            
    print('до ind', new_string)
    
    ind = 0
    if new_string[0] != 'C':
        for i in range(len(new_string)):
            if new_string[i] == '"':
                ind = i+1
                break
        if string+1 != len(file['name']):
            if 'Cluster' in file['name'][string+1]:
                temp = ','.join((''.join(new_string[ind:-1])).split('""'))
                new_file.append(temp)
            else:
                temp = ','.join((''.join(new_string[ind:-1])).split('""')) + ','
                new_file.append(temp)

        elif string+1 == len(file['name']):
            temp = ','.join((''.join(new_string[ind:-1])).split('""'))
            new_file.append(temp)
    else:
        print('после', new_file)
        temp = ''.join(new_string)
        #temp = ''.join(''.join(new_string)).split('""')
        new_file.append(temp)    
    #print(string)



saved_file = open('C:\Data\\New_gene_used_norm.txt', 'w') #менять путь записи и имя записи

for string in (new_file):
    if 'Cluster' in string:
        saved_file.write('\n')
