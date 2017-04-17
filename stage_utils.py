#reindex_stage_ids.py
import re
import os
import glob
import json

def reindex_json_stage_ids():
    stage_meta = []
    #gs =  glob.glob(os.path.dirname(os.path.abspath(__file__))+'/stages/*.json')
    gs =  glob.glob('./stages/*.json')   
    for g in gs:
        if not g.endswith('stage_wrapper.json'):
            stage_meta += [g]
            
    for i in range(0,len(stage_meta)):
        with open(stage_meta[i],'rb') as stage:
            lines = stage.readlines()
            
        with open(stage_meta[i],'wb') as stage:
            print('processing %s'%stage_meta[i])
            for j in range(0,len(lines)):
                m = re.search('stage_id',lines[j])
                if m is not None and m.group()=='stage_id':
                    lines[j] = lines[j].replace(lines[j].split('"stage_id":')[1].split(',')[0],str(i))
                #print(lines[j])
            stage.writelines(lines)

def get_stage_meta():
    stage_meta = {} 
    gs =  glob.glob(os.path.dirname(os.path.abspath(__file__))+'/stages/*.json')
    for stage in gs:
        if not stage.endswith('stage_wrapper.json'):
            with open(stage) as stage_json:
                try:
                    data = json.load(stage_json)
                    k = data.pop('stage_id')
                    stage_meta[k] = data
                except Exception:
                    print('error reading stage meta data from %s'%stage)
    return stage_meta
    
def get_stage_name_id(stage_meta):
    return {stage_meta[k]['name']:'_S'+str(k) for k in stage_meta}
            
#given a list of n strings
#return the longest common string to the left
def get_common_string_left(L):
    S = ''
    if len(L)>1:
        j,n = 0,min([len(i) for i in L])
        for i in range(n):
            if not all([L[0][i]==m[i] for m in L[1:]]): break
            else:                                       j+=1
        S = L[0][0:j]
    if len(L)==1: S = L
    return S

#parses the target string, checks that the files exist and maps each file
#to the opropriate stage_id converted file that the SVE will use internally
#usually this will involve copying and renaming files into a ref directory
#sids = su.get_stage_name_id(stage_meta)
def get_target_map(target_str):
    T,sids = {},get_stage_name_id(get_stage_meta())
    targets = target_str.split(';')
    for t in targets:
        k,s,R,S = t.split(':')[0],t.split(':')[1].split(','),[],[]
        for i in s:
            R += glob.glob(i)
        for r in R:
            if os.path.exists(r): S += [r]
            else: print('%s was not located, please check your path!'%r)
        if k in sids: T[sids[k]] = R
    #check each file for existance
    return T
    
#given a full reference_path to a SVE reference directory
#and a set of stage
def map_stage_names_targets(target_str,ref_fa_path):
    copy_map = {}
    refbase = ref_fa_path.rsplit('/')[-1].replace('.fasta','').replace('.fa','')
    ref_dir = '/'.join(ref_fa_path.rsplit('/')[:-1])+'/'
    T = get_target_map(target_str)
    #new_path = '/'.join(ref_fa_path.rsplit('/')[0:-1])+'/'+refbase+sids['breakseq']
    for t in T:
        for i in range(len(T[t])):
            file_ext = '.'+'.'.join(T[t][i].rsplit('/')[-1].rsplit('.')[1:])
            copy_map[T[t][i]] = ref_dir+refbase+t+file_ext
    return copy_map
    
    
    
    