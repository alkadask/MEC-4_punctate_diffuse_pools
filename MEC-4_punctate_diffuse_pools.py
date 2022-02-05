# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 05:01:01 2021

@author: alkad
"""

#import basic libraries for plotting, data structures and signal processing
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas
from scipy.signal import find_peaks, peak_widths
from scipy.ndimage import gaussian_filter, minimum_filter
import imageio
import os
import fnmatch
import datetime
import seaborn as sns
import shutil

#%%
fpath = 'G:/My Drive/ECM manuscript/github codes/MEC-4_punctate_diffuse_pools/sample_data/input_files/' #filepath where the data is
imgfiles = fnmatch.filter(os.listdir(fpath), '*.tif')

dfpath = 'G:/My Drive/ECM manuscript/github codes/MEC-4_punctate_diffuse_pools/sample_data/output_files/'                
toa = str(datetime.datetime.today()).split()
today = toa[0]
now = toa[1]
timestamp = today.replace('-','')+'-'+now.replace(':','')[:6]

def height_cutoff(nf):
    avnoise = np.mean(nf[nf < np.percentile(nf, 50)])
    stdnoise = np.std(nf[nf < np.percentile(nf, 50)])
    height = avnoise + 5*stdnoise
    return(height)

def peakfinder(nf,height):
    peaks_pos=find_peaks(nf, height=height, prominence=0.5*height)[0]
    pw=peak_widths(nf, peaks_pos, rel_height=0.5)[0]
    pmi=[nf[i] for i in peaks_pos]
    
    #only keep wide peaks if they are also bright
    indices=[]
    for i in np.arange(0, len(peaks_pos)):
        if pw[i]>8 and pw[i]<16 and pmi[i]<1.5*height: indices = indices+[i]
        elif pw[i]>=16 and pmi[i]<3*height: indices = indices+[i]
    peaks_pos=np.delete(peaks_pos,indices)
    pw=np.delete(pw, indices)
    pmi=np.delete(pmi,indices)
    pd=peaks_pos*mu_per_px
    
    #compute nonpeakf
    nonpeakf=nf
    pw2=peak_widths(nf, peaks_pos, rel_height=0.75)
    for i in np.arange(0, len(peaks_pos)):
        indices=np.arange(int(pw2[2][i])-1,int(pw2[3][i])+2)
        np.put(nonpeakf, indices, 0)
        
    return(peaks_pos, pd, pw, pmi, nonpeakf)

#%%
code_version = os.path.basename(__file__)                  #stores the filename of the current script in f
mu_per_px = 0.126     #pixels to microns conversion factor

#specify columns of the pandas dataframe and excel sheets
cols_Data =     ['Date','Strain', 'Allele', 'ImageID', 'Distance', 'Normalized distance', 'Neurite intensity', 'Diffuse intensity', 'Puncta intensity']
cols_Analysis = ['Date','Strain', 'Allele', 'ImageID', 'Segment', 'Total fluorescence', 'Total diffuse fluorescence', 'Total puncta fluorescence', 'Average fluorescence', 'Average diffuse fluorescence', 'Average puncta fluorescence']

#initialize Pandas DataFrames
df_Data = pandas.DataFrame()
df_Analysis = pandas.DataFrame()


strain_key=pandas.DataFrame({('GN753', 'WT', 0),
                             ('GN923', 'nid-1(cg119)', 0),
                             ('GN935','mec-9(u437)', 0),
                             ('GN922','mec-1(e1738)', 0),
                             ('GN941','WT', 0),
                             ('GN1004','mec-1(e1526)', 0),
                             ('GN1062','mec-1(pg154)', 0),
                             ('inj64', 'mec-1(pg164)', 0)
                             }, columns=['Strain name','Allele', 'n'])

#%%    

for x in imgfiles:                            #create loop for number of images in folder
    img = imageio.imread(fpath+x)[:,:,1]    #import image and store it in a list of lists

    imsize = np.shape(img)                  #calculate image size
    dist = np.arange(0,imsize[1])*mu_per_px
    normdist = np.arange(0,imsize[1])/imsize[1]
    
    #extract info from filename
    date = x.split('_')[0]
    strain = x.split('_')[1].split('-')[0]
    row_index=strain_key[(strain_key['Strain name']==strain)].index[0]
    allele = strain_key.loc[row_index,'Allele']
    count = strain_key.loc[row_index,'n'] + 1
    strain_key.at[row_index,'n']=count

    n = img[6:14, 0:]                       #extract rows to use for neurite
    bg = np.concatenate((img[0:4, :], img[16:,:]))   #extract rows to use for background
    rawf = np.mean(n, axis=0)               #calculate average raw neurite fluorescence
    bgf = gaussian_filter(np.mean(bg, axis=0), sigma=10)             #calculate average background fluorescence
    nf = rawf - bgf                         #calculate background subtracted neurite fluorescence
    for i in range(0,len(nf)): 
        if nf[i]<0: nf[i]=0

    bsf = gaussian_filter(minimum_filter(nf, 20),5)
    
    nf2 = nf-bsf
    
    height = 5
    peaks_pos, pd, pw, pmi, nonpeakf = peakfinder(nf2, height)
    
    df = bsf+nonpeakf
    pf = nf-df
    
    
    plt.figure(1, figsize=(20,5))
    plt.title(x)
    plt.axis([0, max(dist), 0, 100])
    plt.plot(dist,nf,'y-')
    plt.plot(dist,df,'b-')
    plt.savefig(dfpath+timestamp+'_'+x+'_traces.png')
    plt.show()
    plt.close()

    all_data1 = pandas.DataFrame({'Date':[date]*imsize[1], 'Strain':[strain]*imsize[1], 'Allele':[allele]*imsize[1], 'ImageID':[x]*imsize[1], 'Distance':dist, 'Normalized distance':normdist, 'Neurite intensity':nf, 'Diffuse intensity':df, 'Puncta intensity':pf}, columns=cols_Data)
    df_Data=df_Data.append(all_data1)

    frame = pandas.DataFrame([[date, strain, allele, x, 'Full', np.sum(nf), np.sum(df), np.sum(pf), np.mean(nf), np.mean(df), np.mean(pf)]], columns=cols_Analysis)
    df_Analysis = df_Analysis.append(frame)
    
#%%
wb = pandas.ExcelWriter(dfpath+timestamp+'_Analysis.xlsx', engine='xlsxwriter')
df_Analysis.to_excel(wb, sheet_name='Analysis')
strain_key.to_excel(wb, sheet_name='Strain summary')
wb.save()

df_Data.to_pickle(dfpath+timestamp+'_Data.pkl')
df_Analysis.to_pickle(dfpath+timestamp+'_Analysis.pkl')
