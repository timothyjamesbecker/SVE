#!/usr/bin/env python
import argparse
import time
import sys
import os
import numpy as np
import multiprocessing as mp
import subprocess32 as subprocess

des = """1000 genomes ftp download management script"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-l', '--ftp_list',type=str, help='path to the tsv ftp download list file\t[None]')
parser.add_argument('-p', '--pop_list',type=str, help='path to the tsv population list file\t[None]')
parser.add_argument('-o', '--out_dir',type=str, help='outputdirectory to save ...bam/ into\t[None]')
parser.add_argument('-c', '--connections',type=int,help='number of connections to use at once\t[1]')
parser.add_argument('-s', '--num_samples',type=int, help='total number of samples to download\t[1]')
args = parser.parse_args()

if args.ftp_list is not None:
    sample_ftp_path = args.ftp_list
else:
    raise IOError
if args.pop_list is not None:
    sample_list_path = args.pop_list
else:
    raise IOError
if args.out_dir is not None:
    out_dir = args.out_dir
    if not os.path.exists(out_dir): os.makedirs(out_dir)
else:
    raise IOError
if args.connections is not None:
    cpus = args.connections
else:
    cpus = 1
if args.num_samples is not None:
    num_samples = args.num_samples
else:
    num_samples = 1

#wget from the ftp a specified sample, unless it is already in the download.check file
def wget(base_url,log_path,sample):
        print('starting sample %s'%sample)
        output,err = '',''
        #[1]unmapped index
        url = base_url+'/%s/alignment/%s.unmapped*.bam.bai'%(sample,sample)
        print(url)
        command = ['cd','/'.join(log_path.rsplit('/')[0:-1])+'/','&&','wget','-c',url]
        try:
            output += subprocess.check_output(' '.join(command), stderr=subprocess.STDOUT, shell=True)
        except Exception:
            err += '\t'.join(command)+'\n'
            pass
        
        #[2]unmapped bam
        url = base_url+'/%s/alignment/%s.unmapped*.bam'%(sample,sample)
        print(url)
        command = ['cd','/'.join(log_path.rsplit('/')[0:-1])+'/','&&','wget','-c',url]
        try:
            output += subprocess.check_output(' '.join(command), stderr=subprocess.STDOUT, shell=True)
        except Exception:
            err += '\t'.join(command)+'\n'
            pass
        
        url = base_url+'/%s/alignment/%s.mapped*.bam.bai'%(sample,sample)
        print(url)
        command = ['cd','/'.join(log_path.rsplit('/')[0:-1])+'/','&&','wget','-c',url]
        #[3]mapped index
        try:
            output += subprocess.check_output(' '.join(command), stderr=subprocess.STDOUT, shell=True)
        except Exception:
            err += '\t'.join(command)+'\n'
            pass
        
        url = base_url+'/%s/alignment/%s.mapped*.bam'%(sample,sample)
        print(url)
        command = ['cd','/'.join(log_path.rsplit('/')[0:-1])+'/','&&','wget','-c',url]
        #[4]mapped bam
        try:
            output += subprocess.check_output(' '.join(command), stderr=subprocess.STDOUT, shell=True)
        except Exception:
            err += '\t'.join(command)+'\n'
            pass
        
            
        print('output:\n'+output)
        #[3a]execute the command here----------------------------------------------------
        #[3b]do a os directory/data check or a ls type command to check the size
        #[3b]of the produced data to ensure that everything is fine...        
        if err == '':
            with open(log_path+'_'+sample,'w') as f:
                f.write(output)
            return output
        else:
            with open(log_path+'_err_'+sample,'w') as f:
                f.write(err)
            return 'error on sample %s'%sample

results = []
def collect_results(result):
    results.append(result)  

if __name__ == '__main__':
    base_url = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/'
    log_path = out_dir
    
    print('using sample_list: %s'%sample_list_path)
    print('using sample ftp: %s'%sample_ftp_path)
    print('writing logs to %s'%log_path)
    print('using %s cpus and wget commands'%cpus)
    
    P,N = {},{}
    with open(sample_list_path,'r') as f:
        sample_list = f.readlines()
    sample_list = [s.replace('\n','').split('\t') for s in sample_list]
    with open(sample_ftp_path, 'r') as f:
        sample_ftp = f.readlines()
    sample_ftp = [s.replace('\n','') for s in sample_ftp]

    for sample in sample_list:
        if sample[0] in sample_ftp:
            if P.has_key(sample[1]):
                P[sample[1]] += [sample[0]]
            else:
                P[sample[1]]  = [sample[0]]
        else:
            if N.has_key(sample[1]):
                N[sample[1]] += [sample[0]]
            else:
                N[sample[1]]  = [sample[0]]

    #------------------------------------------------------
    pops = list(np.random.choice(P.keys(),num_samples,replace=True))
    pick_list = list(np.random.choice(list(set([y for k in P for y in P[k]])),num_samples,replace=False))

    #start || wget calls
    # p1 = mp.Pool(processes=cpus)
    # for sample in pick_list: #for each sample download both mapped and unmapped patterns
    #     p1.apply_async(wget, args=(base_url,log_path,sample), callback=collect_results)
    #     time.sleep(1)
    # p1.close()
    # p1.join()
    #
    # L = []
    # for i in results:
    #     if not i.startswith('error on sample'): L += [i]
    #     else: print(i)
    # print('%s samples were successfully downloaded'%len(L))
