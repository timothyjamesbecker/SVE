#!/usr/bin/env python

#auto.py --ref_path|--ref_file --fqs --out_dir --sample --database
#prep and variant calling steps___________________________________________________________________________________
#(A) --pr_cpus --ref_path (optional if indexes are not found)
#(B) --pb_cpus --pb_threads --pb_mem | uses the bwa_split alogorithm
#(C) --vp_cpus --vp_threads --vp_mem | default ordering with a few || steps
#(D) [1A: breakdancer, 1B: breakseq, 1C: cnmops, 2A:cnvnator, 2B:hydra, 4 lumpy, 5 delly, 6 GS, 7 GATK (optional)]
#merging steps----------------------------------------------------------------------------------------------------
#(E) --fsv_cpus --in_dir --out_dir --model
#(F) --vp_cpus --vp_threads --vp_mem [8 tigra-ext ith fsv input for --target]
#(G) --fsv_cpus --in_dir --contig_dir --model --cross_map_chain --out_dir
#----------------------------------------------------------------------------------------------------------------- 
import argparse
import os
import sys
import glob
import socket
import time
import subprocess32 as subprocess
relink = os.path.dirname(os.path.abspath(__file__))+'/../'
sys.path.append(relink) #go up one in the modules
import read_utils as ru
import stage_utils as su
import svedb

def path(path):
    return os.path.abspath(path)[:-1]

#[1]parse command arguments
des = """
Autonomous Structural Variation Engine:
Given a .fa reference file and a pair: NA12878_1.fq.gz,NA12878_2.fq.gz, 
produce a FusorSV VCF file all_samples_merge.vcf with comprehensive merged SV calls using a fusion model.
#[USAGE] auto.py --ref --fqs --out_dir
#----------------------------------------------------------------------------------------------------------------- """
parser = argparse.ArgumentParser(description=des,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-o', '--out_dir', type=str, help='output directory to store resulting files\t[None]')
parser.add_argument('-r', '--ref', type=str, help='fasta reference path (if indexes are not found, run (A))\t[None]')
parser.add_argument('-d','--database',type=str, help='database configuration file\t[SVE/data]')
fqs_help = """
fq comma-sep file path list\t[None]
[EX PE] --fqs ~/data/sample1_FWD.fq,~/data/sample1_REV.fq"""
parser.add_argument('-f', '--fqs',type=str, help=fqs_help)
parser.add_argument('-b', '--bam',type=str, help='BAM format file for use without alignment or for realignmnet\t[None]')
parser.add_argument('-s', '--sample',type=str, help='Sample identifier\t[infered from BAM|fastq names trimmed by . and -]')
target_help = """
maps a stage name to a comma seperated list of filenames with a colon : (wildcards will work for multiple files)
multiple mappings are sperated by a semi colon ; This is a nessesary step when using breakseq2 and some other callers 
that rely on specific libraries/files with a new reference. The SVE image has hg19 libraries for breakseq2 included.
[EX] -t delly:/data/exclude_list.bed;breakseq:/data/breakseq2_bplib_20150129*\n\t[None]
"""
parser.add_argument('-t', '--target',type=str, help=target_help)
parser.add_argument('-m', '--model',type=str, help='data fusion model\t[FusorSV/data/models/g1k_v37_decoy.P3.pickle]')
parser.add_argument('-P', '--cpus',type=str, help='number of cpus for alignment and sorting, ect\t[1]')
parser.add_argument('-T', '--threads',type=str, help='number of threads per CPU\t[4]')
parser.add_argument('-M', '--mem',type=str, help='ram in GB units to use for processing per cpu/thread unit\t[4]')
parser.add_argument('-a', '--algorithm',type=str, help='alignment/realignment algorithm pipeline\t[speed_seq]')
parser.add_argument('-R', '--test_run',type=str, help='generates paired end FASTQ files and tests instalation\t[place holder]')
args = parser.parse_args()

#read the database configuration file
dbc = {'srv':'','db':'','uid':'','pwd':''}
if args.database is not None:
    with open(args.database, 'r') as f:
        params = f.read().split('\n') #newline seperated configuration file
    try:
        dbc['srv']  = params[0].split('srv=')[-1]
        dbc['db'] = params[0].split('srv=')[-1]
        dbc['uid'] = params[0].split('srv=')[-1]
        dbc['pwd'] = params[0].split('srv=')[-1]
    except Exception:
        print('invalid database configuration')
        print('running the SVE without the SVEDB')
        pass
else:
    with open(os.path.dirname(os.path.abspath(__file__))+'/../data/svedb.config', 'r') as f:
        params = f.read().split('\n') #newline seperated configuration file
    try:
        dbc['srv'] = params[0].split('srv=')[-1]
        dbc['db']  = params[1].split('db=')[-1]
        dbc['uid'] = params[2].split('uid=')[-1]
        dbc['pwd'] = params[3].split('pwd=')[-1]
        schema = {}
        with svedb.SVEDB(dbc['srv'], dbc['db'], dbc['uid'], dbc['pwd']) as dbo:
            dbo.embed_schema()   #check the schema for a valid db
            schema = dbo.schema
        if len(schema)<1:
            print('dbc:%s' % [c + '=' + dbc[c] for c in dbc])
            print('invalid database configuration')
            print('running the SVE without the SVEDB')
        else:
            print('dbc:%s'%[c+'='+dbc[c] for c in dbc])
            print('valid database configuration found')
            print('running the SVE with the SVEDB')
    except Exception:
        print('invalid database configuration')
        print('running the SVE without the SVEDB')
        pass

host = socket.gethostname()
directory = path('~/'+host+'/') #users base home folder as default plus hostname
scripts_path = os.path.dirname(os.path.abspath(__file__))+'/'
if args.out_dir is not None:    #optional reroute
    directory = args.out_dir
    if directory[-1]=='/': directory = directory[:-1]
if not os.path.exists(directory): os.makedirs(directory)
#if the reference has a reference folder set up and you pass its path, it will skip prepare_ref.py
if args.ref is not None and os.path.exists(args.ref):
    ref_fa_path = args.ref
else:
    print('no ref pattern found: args=%s'%args.ref)
    raise IOError
    ref_fa_path = ''
refbase = ref_fa_path.rsplit('/')[-1].split('.fa')[0]
if args.fqs is not None and all([os.path.exists(f) for f in args.fqs.split(',')]):
    reads = args.fqs.split(',') #CSL returns a list of 1+    
else:
    print('no fqs pattern found: arg=%s'%args.fqs)
    reads = []
if not all([os.path.exists(r) for r in reads]):
    print('fastq files not found!')
if args.bam is not None and os.path.exists(args.bam):
    bam = args.bam
else:
    print('no bam pattern found: arg=%s'%args.bam)
    bam = ''
if bam=='' and len(reads)<1:
    print('fastq files or a BAM are needed to proceed!')
    raise IOError    
#check that either a bam is found:with a matching header to reference
#optionally employ a realignment pipeline using samtools fastq to bwa mem workflow

if args.sample is not None:
    SM = args.sample
    print('using SM = %s'%SM)
else:
    if len(reads)>0:
        SM = reads[0].rsplit('/')[-1].rsplit('-')[0].rsplit('.')[0][:-1]
        print('using SM = %s'%SM)
if args.target is not None:
    target_str = args.target #target mapping string, need to parse it
else:
    target_str = 'breakseq:'+os.path.dirname(os.path.abspath(__file__))+'/../data/breakseq2*.fna'
if args.cpus is not None:
    cpus = int(args.cpus)
else:
    cpus = 4
if args.threads is not None:
    threads = args.threads
else:
    threads = 6
if args.mem is not None:
    mem = int(args.mem)
else:
    mem = 12
if args.algorithm is not None:
    algorithm = args.algorithm
else:
    algorithm = 'speed_seq'
    
#(A) check for the reference preparations in that directory
#[1] read ref as dict
seqs = set(ru.get_fasta_seq_names_lens(ref_fa_path).keys())
#[2] glob the reference directory to see if files are present
ref_prep = True
for f in ['.fa','.dict','.amb','.ann','.bwt','.pac','.sa','.svmask.fasta']:
    g = set([g.rsplit('/')[-1].replace(f,'').replace('.fa','').rsplit('_S15_')[-1] for g in \
             glob.glob('/'.join((ref_fa_path.rsplit('/')[:-1]))+'/*'+f)])
    if len(seqs.difference(g))>0:
        print('reference files did not match the expectation!\n%s'%seqs.difference(g))
        ref_prep = False
#check for the target_str now
output = ''
target_name_map = su.map_stage_names_targets(target_str,ref_fa_path)
print('reference target file mapping: %s\n'%target_name_map)
for t in target_name_map:
    if not os.path.exists(target_name_map[t]):
        try:
            file_copy = ['cp',t,target_name_map[t]]
            output += subprocess.check_output(' '.join(file_copy),shell=True)
        except Exception as E:
                print(E)
                ref_prep = False

#[3] if not run prepare_ref.py (easiest way to include breakseq bp-lib is a target map that is input)
r_start = time.time()
if ref_prep: #ref_directory is the ref_path inside a valid ref_folder
    ref_directory = '/'.join(ref_fa_path.rsplit('/')[:-1])+'/'
    print('reference has already been prepared, skipping...')
else:
    output = ''
    ref_directory = '/'.join(ref_fa_path.rsplit('/')[:-1])+'/'+\
                         ref_fa_path.rsplit('/')[-1].rsplit('.')[0]+'/'
    print('reference not prepared:\nbuilding reference directory at location: %s'%ref_directory)
    prepare_ref = [scripts_path+'prepare_ref.py','-r',ref_fa_path, #add target map for breakseq
                   '-o',ref_directory,'-P',str(cpus)]
    try:
        output += subprocess.check_output(' '.join(prepare_ref),shell=True)
        ref_fa_path = ref_directory+ref_fa_path.rsplit('/')[-1] #update new location
    except Exception as E:
        print(E)
    print(output)
r_stop = time.time()
if os.path.exists(ref_fa_path):
    print('prepare_ref has completed in %s minutes'%round((r_stop-r_start)/60.0,2))
else:
    print('prepare_ref failure shuting down execution!')
    raise IOError

    #(B) prepare bam files using the fastest method
b_start = time.time()
if len(reads)>0 and len(glob.glob(directory+'/bam/*%s*.bam'%SM))<1:
    output = ''
    bam_directory = directory+'/bam/'
    print('BAM file not aligned:\nstarting alignment of FASTQ files at location: %s'%bam_directory)
    prepare_bam = [scripts_path+'prepare_bam.py','-r',ref_fa_path,'-f',','.join(reads),'-o',bam_directory,'-s',SM,
                   '-P',str(cpus),'-T',str(threads),'-M',str(mem),'-a',algorithm]
    try:
        output += subprocess.check_output(' '.join(prepare_bam),shell=True)
        bam_path = glob.glob(bam_directory+'%s.bam'%SM)[0]
    except Exception as E:
        print(E)
    print(output)
else:
    print('bam files were found, skipping alignment step')
    bam_path = glob.glob(directory+'/bam/*%s*.bam'%SM)[0]
b_stop = time.time()
if os.path.exists(ref_fa_path):
    print('prepare_bam has completed in %s minutes'%round((b_stop-b_start)/60.0,2))
else:
    print('prepare_bam failure shuting down execution!')
    raise IOError
    
#(C) start as many variant callers as you can (in series ...)
v_start = time.time()
output = ''
vcf_directory = directory+'/vcf/'
if not os.path.exists(vcf_directory): os.mkdir(vcf_directory)
vcf_sample = vcf_directory+'/%s/'%SM
print('starting variant calling')
variant_processor = [scripts_path+'variant_processor.py','-r',ref_fa_path,'-b',bam_path,
                     '-o',vcf_sample,'-s',
                     'breakdancer,breakseq,cnmops,cnvnator,delly,hydra,lumpy,genome_strip,gatk_haplo']
try:
    output += subprocess.check_output(' '.join(variant_processor),shell=True)
except Exception as E:
    print(E)
print(output)
v_stop  = time.time()

#(D) Run FusorSV to merge the results
output = ''
print('starting FusorSV processing with a prior model')
fusorsv_out = directory+'/fusorsv_out/'

fusor_sv = [scripts_path+'../../FusorSV/FusorSV.py',
            '-r',ref_fa_path,'-i',vcf_directory,
            '-f',scripts_path+'../../FusorSV/data/models/human_g1k_v37_decoy.P3.pickle',
            '-o',fusorsv_out,'-p',str(cpus),'-M',str(0.5),'-L','DEFAULT']
try :
    output += subprocess.check_output(' '.join(fusor_sv),shell=True)
except Exception as E:
    print(E)
print(output)

#(D) check the number of VCF files and total SV calls per caller
#[1] check FusorSV file folder structure
#[2] issue some statistics per sample found in terms of call inputs
#[3] if targets are specified, use them, otherwise (if new ref is issued, do coordinate map on it)
#[4] check the model to see if it is in a valid ref space
#[5] run the FusorSV analysis
#[6] optional cross_map coordinates and gene_finding tools



  