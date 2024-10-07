# -*- coding: utf-8 -*-
"""
Created the 2024.10.07
@author: Joé Petrazoller, LEM3, Metz, France : joe.petrazoller@univ-lorraine.fr
"""

#################################
#Tensor name example : Pij_11_22 : means : "the elastic dipole P_22 component after a strain in the 11 direction
#################################
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats 
from matplotlib.ticker import ScalarFormatter

########################    To display number in indexes in plot  ###########################
dict1 = {'0':'\u2070',
         '1':'\u00b9',
         '2':'\u00b2',
         '3':'\u00b3',
         '4':'\u2074',
         '5':'\u2075',
         '6':'\u2076',
         '7':'\u2077',
         '8':'\u2078',
         '9':'\u2079',
         '+':'\u207A',
         '-':'\u207B',
         '=':'\u207C',
         '(':'\u207D',
         ')':'\u207E',
         'n':'\u2084',}
dict2 = {'0':'\u2080',
         '1':'\u2081',
         '2':'\u2082',
         '3':'\u2083',
         '4':'\u2084',
         '5':'\u2085',
         '6':'\u2086',
         '7':'\u2087',
         '8':'\u2088',
         '9':'\u2089',
         '+':'\u208A',
         '-':'\u208B',
         '=':'\u208C',
         '(':'\u208D',
         ')':'\u208E',
         'a':'\u2090',
         'e':'\u2091',
         'o':'\u2092',
         'x':'\u2093',
         'h':'\u2095',
         'k':'\u2096',
         'l':'\u2097',
         'm':'\u2098',
         'n':'\u2099',
         'p':'\u209A',
         's':'\u209B',
         't':'\u209C'}
def subs(base,x):
    z = '{}'.format(dict2.get(x))
    return base + z

def plot_results(y,y_name,x):
    x=x*100 #because strains are displayed in %
    plt.figure()
    plt.plot(x,y,'o',color='blue')
    (a, b, R2, _, uA ) = stats.linregress (x, y)
    plt.plot(x, a*x+b,'blue')
    plt.gca().yaxis.set_major_formatter(ScalarFormatter())
    plt.gca().ticklabel_format(style='plain', axis='y')
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.2)    
    x_label_tempo=y_name[4:6]
    if x_label_tempo=='11':
        x_label1='1'
        x_label2='1'
    elif x_label_tempo=='22':
        x_label1='2'
        x_label2='2'        
    elif x_label_tempo=='33':
        x_label1='3'
        x_label2='3'
    elif x_label_tempo=='23':
        x_label1='2'
        x_label2='3'
    elif x_label_tempo=='13':
        x_label1='1'
        x_label2='3'
    elif x_label_tempo=='12':
        x_label1='1'
        x_label2='2'
    plt.xlabel(subs('ε',x_label1)+subs('',x_label2)+' (%)',fontsize=20)
    plt.ylabel(subs('P',y_name[-2:-1])+subs('',y_name[-1])+' (eV)',fontsize=20)
    plt.gca().get_yaxis().get_offset_text().set_visible(False)
    plt.gca().get_xaxis().get_offset_text().set_visible(False)
    #Create a "Figure" folder with all the P_ij (elastic dipole) in terms of epsilon_ij (strains)
    if not os.path.isdir(figures_dir):
        os.makedirs(figures_dir)
    plt.text(max(x)-0.05, min(y),'R² = '+"%.3f" % R2,fontsize='20')
    plt.savefig(figures_dir+y_name+".png", dpi=200) #This line to tab or untab 
    plt.close()
    
##################################################
#############  Computation starting   ############
##################################################
global E
global nu
WORKDIR = glob.glob(os. getcwd())[0] #Set Working directory on 
figures_dir=os.path.join(WORKDIR, 'Figures/')
ymin_list=[]
#List all .txt files results from previous LAMMPS computation
file_raw=os.listdir()
file_txt=[]
for i in range(0, len(file_raw)):
    if file_raw[i][-4:]=='.txt' and file_raw[i][0]!='1':
        file_txt.append(file_raw[i])
alpha_list=[]

for k in [1,2]: #positive and negative strains
    if k==1:
        file_txt2=file_txt[0:int(len(file_txt)/2)] #negative values
    if k==2:
        file_txt2=file_txt[int(len(file_txt)/2):int(len(file_txt))] #positive values
    alpha=np.array([[0.0,0.0,0.0,0.0,0.0,0.0],
          [0.0,0.0,0.0,0.0,0.0,0.0],
          [0.0,0.0,0.0,0.0,0.0,0.0],
          [0.0,0.0,0.0,0.0,0.0,0.0],
          [0.0,0.0,0.0,0.0,0.0,0.0],
          [0.0,0.0,0.0,0.0,0.0,0.0]])
    name_list=[]
    Pij_list=[]
################  Assignin all P_ij values in lists  ###############
    Pij_11_11=[]
    Pij_11_22=[]
    Pij_11_33=[]
    Pij_11_23=[]
    Pij_11_13=[]
    Pij_11_12=[]
    Pij_22_11=[]
    Pij_22_22=[]
    Pij_22_33=[]
    Pij_22_23=[]
    Pij_22_13=[]
    Pij_22_12=[]
    Pij_33_11=[]
    Pij_33_22=[]
    Pij_33_33=[]
    Pij_33_23=[]
    Pij_33_13=[]
    Pij_33_12=[]
    Pij_23_11=[]
    Pij_23_22=[]
    Pij_23_33=[]
    Pij_23_23=[]
    Pij_23_13=[]
    Pij_23_12=[]
    Pij_13_11=[]
    Pij_13_22=[]
    Pij_13_33=[]
    Pij_13_23=[]
    Pij_13_13=[]
    Pij_13_12=[]
    Pij_12_11=[]
    Pij_12_22=[]
    Pij_12_33=[]
    Pij_12_23=[]
    Pij_12_13=[]
    Pij_12_12=[]
    
    for j in range(0,len(file_txt2)): #going through all negative (k=1) or positive (k=2)  strains
        filename=file_txt2[j]
        with open(filename,'r') as file:
            line=file.readlines()
        P11=float(line[0].split()[1])
        P22=float(line[1].split()[1])
        P33=float(line[2].split()[1])
        P23=float(line[3].split()[1])
        P13=float(line[4].split()[1])
        P12=float(line[5].split()[1])
        Pij=np.array([[P11,P12,P13],
              [P12,P22,P23],
              [P13,P23,P33]])
        #Initilisation of the strains values
        if filename[-5]=='x':
            name_list.append(float(filename[:-6][4:]))
        if filename[-6]=='_' and filename[-5]=='x' and filename[-4]=='.':
            Pij_11_11.append(float(Pij[0][0]))
            Pij_11_22.append(float(Pij[1][1]))
            Pij_11_33.append(float(Pij[2][2]))
            Pij_11_23.append(float(Pij[1][2]))
            Pij_11_13.append(float(Pij[0][2]))
            Pij_11_12.append(float(Pij[0][1]))
        if filename[-6]=='_' and filename[-5]=='y' and filename[-4]=='.':
            Pij_22_11.append(float(Pij[0][0]))
            Pij_22_22.append(float(Pij[1][1]))
            Pij_22_33.append(float(Pij[2][2]))
            Pij_22_23.append(float(Pij[1][2]))
            Pij_22_13.append(float(Pij[0][2]))
            Pij_22_12.append(float(Pij[0][1]))
        if filename[-6]=='_' and filename[-5]=='z' and filename[-4]=='.':
            Pij_33_11.append(float(Pij[0][0]))
            Pij_33_22.append(float(Pij[1][1]))
            Pij_33_33.append(float(Pij[2][2]))
            Pij_33_23.append(float(Pij[1][2]))
            Pij_33_13.append(float(Pij[0][2]))
            Pij_33_12.append(float(Pij[0][1]))
        if filename[-6]=='y' and filename[-5]=='z':
            Pij_23_11.append(float(Pij[0][0]))
            Pij_23_22.append(float(Pij[1][1]))
            Pij_23_33.append(float(Pij[2][2]))
            Pij_23_23.append(float(Pij[1][2]))
            Pij_23_13.append(float(Pij[0][2]))
            Pij_23_12.append(float(Pij[0][1]))
        if filename[-6]=='x' and filename[-5]=='z':
            Pij_13_11.append(float(Pij[0][0]))
            Pij_13_22.append(float(Pij[1][1]))
            Pij_13_33.append(float(Pij[2][2]))
            Pij_13_23.append(float(Pij[1][2]))
            Pij_13_13.append(float(Pij[0][2]))
            Pij_13_12.append(float(Pij[0][1]))
        if filename[-6]=='x' and filename[-5]=='y':
            Pij_12_11.append(float(Pij[0][0]))
            Pij_12_22.append(float(Pij[1][1]))
            Pij_12_33.append(float(Pij[2][2]))
            Pij_12_23.append(float(Pij[1][2]))
            Pij_12_13.append(float(Pij[0][2]))
            Pij_12_12.append(float(Pij[0][1]))
    x=np.array(name_list)
    for i in range(0,len(name_list)):
        if k==1: 
            x[i]=x[i]
        if k==2: 
            x[i]=x[i]
            
    list_pij=[]
    list_pij.append(Pij_11_11)
    list_pij.append(Pij_11_22)
    list_pij.append(Pij_11_33)
    list_pij.append(Pij_11_23)
    list_pij.append(Pij_11_13)
    list_pij.append(Pij_11_12)
    list_pij.append(Pij_22_11)
    list_pij.append(Pij_22_22)
    list_pij.append(Pij_22_33)
    list_pij.append(Pij_22_23)
    list_pij.append(Pij_22_13)
    list_pij.append(Pij_22_12)
    list_pij.append(Pij_33_11)
    list_pij.append(Pij_33_22)
    list_pij.append(Pij_33_33)
    list_pij.append(Pij_33_23)
    list_pij.append(Pij_33_13)
    list_pij.append(Pij_33_12)
    list_pij.append(Pij_23_11)
    list_pij.append(Pij_23_22)
    list_pij.append(Pij_23_33)
    list_pij.append(Pij_23_23)
    list_pij.append(Pij_23_13)
    list_pij.append(Pij_23_12)
    list_pij.append(Pij_13_11)
    list_pij.append(Pij_13_22)
    list_pij.append(Pij_13_33)
    list_pij.append(Pij_13_23)
    list_pij.append(Pij_13_13)
    list_pij.append(Pij_13_12)
    list_pij.append(Pij_12_11)
    list_pij.append(Pij_12_22)
    list_pij.append(Pij_12_33)
    list_pij.append(Pij_12_23)
    list_pij.append(Pij_12_13)
    list_pij.append(Pij_12_12)
    
    #Concatenate into a big list Pij
    list_pij_name=[]
    list_pij_name.append('11_11')
    list_pij_name.append('11_22')
    list_pij_name.append('11_33')
    list_pij_name.append('11_23')
    list_pij_name.append('11_13')
    list_pij_name.append('11_12')
    list_pij_name.append('22_11')
    list_pij_name.append('22_22')
    list_pij_name.append('22_33')
    list_pij_name.append('22_23')
    list_pij_name.append('22_13')
    list_pij_name.append('22_12')
    list_pij_name.append('33_11')
    list_pij_name.append('33_22')
    list_pij_name.append('33_33')
    list_pij_name.append('33_23')
    list_pij_name.append('33_13')
    list_pij_name.append('33_12')
    list_pij_name.append('23_11')
    list_pij_name.append('23_22')
    list_pij_name.append('23_33')
    list_pij_name.append('23_23')
    list_pij_name.append('23_13')
    list_pij_name.append('23_12')
    list_pij_name.append('13_11')
    list_pij_name.append('13_22')
    list_pij_name.append('13_33')
    list_pij_name.append('13_23')
    list_pij_name.append('13_13')
    list_pij_name.append('13_12')
    list_pij_name.append('12_11')
    list_pij_name.append('12_22')
    list_pij_name.append('12_33')
    list_pij_name.append('12_23')
    list_pij_name.append('12_13')
    list_pij_name.append('12_12')
    
############    Linear regression to calculate the alpha values #############
    for i in range(0,len(list_pij)):
        y=np.array(list_pij[i])
        (a, b, R2, _, uA ) = stats.linregress (x, y)
        if i<6:
            alpha[i][0]=a
        if i<12 and i>5:
            alpha[i-6][1]=a
        if i<18 and i>11:
            alpha[i-12][2]=a
        if i<24 and i>17:
            alpha[i-18][3]=a
        if i<30 and i>23:
            alpha[i-24][4]=a
        if i<36 and i>29:
            alpha[i-30][5]=a        
            
        #Plot
        if k==1:
            plot_results(list_pij[i], 'neg_'+list_pij_name[i],x)
        if k==2:
            plot_results(list_pij[i], 'pos_'+list_pij_name[i],x)
            
    alpha_list.append(alpha)
        
alpha=0.5*(alpha_list[0]+alpha_list[1])

# For visualization, all alpha tensor values closer than 0.0001 eV are set manually to 0
for i in range(0,6):
    for j in range(0,6):
        if abs(alpha[i][j])<0.0001:
            alpha[i][j]=0

# Write results in a .txt file. Values with more digits after comma can be found within the  "alpha" tensor
filename='1.polarizability_tensor.txt'
with open(filename,'w') as file:
    file.write(str("%.1f" % alpha[0][0]))
    file.write('     ')
    file.write(str("%.1f" % alpha[0][1]))
    file.write('     ')
    file.write(str("%.1f" % alpha[0][2]))
    file.write('     ')
    file.write(str("%.1f" % alpha[0][3]))
    file.write('     ')
    file.write(str("%.1f" % alpha[0][4]))
    file.write('     ')
    file.write(str("%.1f" % alpha[0][5]))
    file.write('\n')
    file.write(str("%.1f" % alpha[1][0]))
    file.write('     ')
    file.write(str("%.1f" % alpha[1][1]))
    file.write('     ')
    file.write(str("%.1f" % alpha[1][2]))
    file.write('     ')
    file.write(str("%.1f" % alpha[1][3]))
    file.write('     ')
    file.write(str("%.1f" % alpha[1][4]))
    file.write('     ')
    file.write(str("%.1f" % alpha[1][5]))
    file.write('\n')
    file.write(str("%.1f" % alpha[2][0]))
    file.write('     ')
    file.write(str("%.1f" % alpha[2][1]))
    file.write('     ')
    file.write(str("%.1f" % alpha[2][2]))
    file.write('     ')
    file.write(str("%.1f" % alpha[2][3]))
    file.write('     ')
    file.write(str("%.1f" % alpha[2][4]))
    file.write('     ')
    file.write(str("%.1f" % alpha[2][5]))
    file.write('\n')
    file.write(str("%.1f" % alpha[3][0]))
    file.write('     ')
    file.write(str("%.1f" % alpha[3][1]))
    file.write('     ')
    file.write(str("%.1f" % alpha[3][2]))
    file.write('     ')
    file.write(str("%.1f" % alpha[3][3]))
    file.write('     ')
    file.write(str("%.1f" % alpha[3][4]))
    file.write('     ')
    file.write(str("%.1f" % alpha[3][5]))
    file.write('\n')
    file.write(str("%.1f" % alpha[4][0]))
    file.write('     ')
    file.write(str("%.1f" % alpha[4][1]))
    file.write('     ')
    file.write(str("%.1f" % alpha[4][2]))
    file.write('     ')
    file.write(str("%.1f" % alpha[4][3]))
    file.write('     ')
    file.write(str("%.1f" % alpha[4][4]))
    file.write('     ')
    file.write(str("%.1f" % alpha[4][5]))
    file.write('\n')
    file.write(str("%.1f" % alpha[5][0]))
    file.write('     ')
    file.write(str("%.1f" % alpha[5][1]))
    file.write('     ')
    file.write(str("%.1f" % alpha[5][2]))
    file.write('     ')
    file.write(str("%.1f" % alpha[5][3]))
    file.write('     ')
    file.write(str("%.1f" % alpha[5][4]))
    file.write('     ')
    file.write(str("%.1f" % alpha[5][5])) 
file.close()