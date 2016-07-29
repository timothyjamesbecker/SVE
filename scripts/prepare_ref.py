#!/usr/bin/env python
#prepare_ref.py v 0.5, 5/23/2015
#Timothy Becker, UCONN/SOE/CSE Phd Candidate
#preprocessing of reference sequences as well as random sequence generation
#[EX] python ~/path/prepare_ref.py -r -o ~/node_data/ -L 5E4 -C 5
import argparse 
import os
import sys
import socket
import subprocess32 as subprocess
import time
import random
import hashlib
relink = os.path.dirname(os.path.abspath(__file__))+'/../'
sys.path.append(relink) #go up one in the modules
import read_utils as sr
import svedb
import stage

def path(path):
    return os.path.abspath(path)+'/'
des = """
Prepares a random fasta reference or existing fasta 
reference file for downstream SV analysis.
This includes chrom copying, indexing, masking and 
all associated one-time processecies that need to 
be done before || per sample or other || means
can be started in indepencence."""
parser = argparse.ArgumentParser(description=des)
group = parser.add_mutually_exclusive_group()
group.add_argument('-g', '--gen_rnd', action='store_true',help='generate a rgXXX.fa file using -L and -C')
group.add_argument('-f', '--file', type=str, help='fasta reference file')
parser.add_argument('-o', '--out_dir', type=str, help='output directory to store resulting files')
parser.add_argument('-l','--len',type=float, help='random reference total length')
parser.add_argument('-c','--chr',type=float,help='random reference number of chroms')
parser.add_argument('-d', '--database', type=str, help='UCONN,JAX')
parser.add_argument('-e','--erase_db',action='store_true',help='wipe the SVEDB')
args = parser.parse_args()

directory = path('~/refs')    #users base home folder as default plus refs folder
if args.out_dir is not None:  #optional reroute
    directory = args.out_dir
if not os.path.exists(directory): os.makedirs(directory)
if args.gen_rnd: #validate flag and params
    #do the random generation with length and chrom number params
    start = time.time()
    refbase = 'rg'+str(int(random.uniform(0,1E3))).zfill(3)
    ref_fa_path = directory+refbase+'.fa'
    if args.len>float(1E3) and args.chr>2 and args.len>float(1E2)*args.chr: #baseline limits
        L = int(args.len)                                            #overall length
        num_chroms   = int(args.chr)                                 #2+chr (X,Y for ploidy)
    else:
        L = int(args.len)                                            #overall length
        num_chroms   = 3                                             #2+chr (X,Y for ploidy)
    chrom_names  = ['chr'+str(i) for i in range(1,num_chroms-1)] #chrom_num,chrom_len
    chrom_names += ['chrX','chrY']                               #can simulate sex/ploidy here...
    print('building random reference of length %s'%L)
    print('and geometric-like distribution of %s chroms'%len(chrom_names))
    
    #start generation and print run time metrics  
    chroms = gv.gen_chrom_lens(L,chrom_names) #{'chr1':200,'chr2':100, ... ,'chrn':10}
    print(chroms)
    seqs = gv.gen_random_seqs(chroms,method='fast')
    seqs = sorted(seqs,key=lambda x: len(x.seq),reverse=True)
    print('random reference written: %s'%sr.write_fasta(seqs,ref_fa_path)) #this is in directory
    seqs = sr.read_fasta(ref_fa_path,dictionary=True)
    chr_lens = [len(seqs[k]) for k in seqs]
    stop = time.time()
    print('real cpu time: %s seconds'%(stop-start))
elif args.file is not None: #using existing ref, now just process
    ref_fa_path = args.file
    print('using ref_path = '+ref_fa_path)
    out = ''
    try:
        out = subprocess.check_output(' '.join(['cp',ref_fa_path,directory]),shell=True)
    except Exception as E:
        print('I/O Copy Error')
        
    if ref_fa_path[-3:]=='.gz':
        print('.gz extension detected...')
        gzref = directory+ref_fa_path.rsplit('/')[-1]
        refbase = ref_fa_path.rsplit('/')[-1].split('.fa')[0]
        ref_fa_path = directory+refbase+'.fa'
        try:
            print('uncompressing '+gzref)
            out = subprocess.check_output(' '.join(['gzip','-d',gzref]),shell=True)
            print('renaming %s to %s'%(gzref[:-3],ref_fa_path))
            out = subprocess.check_output(' '.join(['mv',gzref[:-3],ref_fa_path]),shell=True)
        except Exception as E:
            print('I/O Copy Error')
    else:
        refbase = ref_fa_path.rsplit('/')[-1].split('.fa')[0]
        ref_fa_path = directory+refbase+'.fa'
    print('reading reference sequence...')
    seqs = sr.read_fasta(ref_fa_path,dictionary=True)     #finally read the fasta with read_utils
    chr_lens = [len(seqs[k]) for k in seqs]
else:
    raise IOError #return and Error
    
print('reference validated, moving on to preprossesing...')
#start DB and stages for processing now

dbc = {'srv':'','db':'','uid':'','pwd':''}
if args.database is not None:
    dbc['db']  = 'sve'
    dbc['uid'] = 'sv_calibrator'
    dbc['pwd'] = 'sv_calibrator'
    if args.database.upper() == 'UCONN':
        dbc['srv'] = 'arc-gis.ad.engr.uconn.edu'
    elif args.database.upper() == 'JAX':
        dbc['srv'] = 'ldg-jgm003.jax.org'

else:
    print('invalid database configuration')
    raise KeyError
    
with svedb.SVEDB(dbc['srv'], dbc['db'], dbc['uid'], dbc['pwd']) as dbo:
    ks = sorted(seqs.keys())
    host = socket.gethostname()    
    #start db entriesDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
    if args.erase_db:
        print('reseting SVEDB :::::::::::::::::::::')
        dbo.new() #conditional reset
    dbo.embed_schema() #embed the db schema into the dbo object
    #get table dict
    tables = dbo.select_tables()
    for t in tables['names']: print t #good you can get table names
    if args.gen_rnd: #simulation saves the exact rnd seqs used
        rsb = dbo.obj_to_blob(seqs)
    else: #well known seqs get a hash of the sorted concatenated contigs
        hsh = hashlib.sha512(''.join([seqs[i].seq for i in ks]))
        rsb = dbo.obj_to_blob(hsh.digest())
    #get next non-assigned ref_id
    ref_id = 0 #run_id == get next int form DB? 
#    name,ref_len,seq_names,seq_lens,seqs='',url=''
    dbo.new_ref(name=refbase,ref_len=sum(chr_lens),
                seq_names=','.join(seqs.keys()),
                seq_lens=','.join([str(l) for l in chr_lens]),
                seqs=rsb,url='')
    ref_id = dbo.get_max_key('refs')
    print('using registered ref_id=%s'%str(ref_id))
    print('with contig list:')
    print(seqs.keys())
    print(''.join(['-' for i in range(60)]))
    
    dbo.new_run('illumina',host,ref_id)
    run_id = dbo.get_max_key('runs')
    #do this once for each ref ... RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR    
    #split by chrom here and index for random access...
    ks = sorted(seqs.keys())
    ss = [seqs[i] for i in ks] #convert the dict to seq list...
    print('writting fasta by chrom...')
    chrom_names = sr.write_fasta_by_chrom(ss,directory,'') #shared to the directory for some callers
    ss_len = len(ss)
    tt = [x.name for x in ss] #get names
    for x in tt:
        #index for alignment
        st = stage.Stage('bwa_index',dbc)
        outs = st.run(run_id,{'.fa':[directory+x+'.fa']})
        print(outs)
        st = stage.Stage('picard_dict',dbc)
        outs = st.run(run_id,{'.fa':[directory+x+'.fa']})
        print(outs)
    seqs,ss = {},[] #clear out mem?
    #indexes for alignment
    st = stage.Stage('samtools_fasta_index',dbc)
    out = st.run(run_id, {'.fa':[ref_fa_path]})
    st = stage.Stage('bwa_index',dbc)
    outs = st.run(run_id,{'.fa':[ref_fa_path]})
    st = stage.Stage('picard_dict',dbc)
    outs = st.run(run_id,{'.fa':[ref_fa_path]})

    #genomestrip index and genome masking
    st = stage.Stage('genome_strip_prepare_ref',dbc)
    outs = st.run(run_id, {'.fa':[ref_fa_path]})
    print(outs)

    #add a fa2bit stage for use with the hydra to vcf converter...
    st = stage.Stage('fa_to_2bit',dbc)
    outs = st.run(run_id, {'.fa':[ref_fa_path]})
    print(outs)

    #mrfast/variation hunter fasta ref indexing
    st = stage.Stage('mrfast_index')
    outs = st.run(run_id,{'.fa':[ref_fa_path]})
    print(outs)
    
    
    #do this once for each ref ... RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR