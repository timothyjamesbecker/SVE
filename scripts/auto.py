#!/usr/bin/env python
import argparse
import os
import sys
import glob
import socket
import time
import subprocess32 as subprocess
import multiprocessing as mp
relink = os.path.dirname(os.path.abspath(__file__))+'/../'
sys.path.append(relink) #go up one in the modules
import read_utils as ru
import stage_utils as su

result_list = []
def collect_results(result):
    result_list.append(result)

def path(path):
    return os.path.abspath(path)[:-1]

def run_variant_processor(variant_processor,stages):
    output = ''
    v_start = time.time()
    try:
        output += subprocess.check_output(' '.join(variant_processor+[stages]),shell=True)
    except Exception as E:
        print(E)
        raise IOError
    v_stop  = time.time()
    return [[output,round((v_stop-v_start)/60.0,2)]] #time in hours

#[1]parse command arguments
des = """Autonomous Structural Variation Engine:
Given a .fa reference file and a pair: NA12878_1.fq.gz,NA12878_2.fq.gz, 
produce a FusorSV VCF file all_samples_merge.vcf with comprehensive merged SV calls using a fusion model.
#[USAGE] auto.py --ref --fqs|bam --out_dir
#----------------------------------------------------------------------------------------------------------------- """
parser = argparse.ArgumentParser(description=des,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-o', '--out_dir', type=str, help='output directory to store resulting files\t[None]')
parser.add_argument('-r', '--ref', type=str, help='fasta reference path (if indexes are not found, run (A))\t[None]')
fqs_help = """fq comma-sep file path list\t[None]
[EX PE] --fqs ~/data/sample1_FWD.fq,~/data/sample1_REV.fq"""
parser.add_argument('-f', '--fqs',type=str, help=fqs_help)
parser.add_argument('-b', '--bam',type=str, help='BAM format file for use without alignment or for realignmnet\t[None]')
parser.add_argument('-s', '--sample',type=str, help='Sample identifier\t[infered from BAM|fastq read 1 trimmed by . and _]')
target_help = """maps a stage name to a comma seperated list of filenames with a colon : (wildcards will work for multiple files)
multiple mappings are sperated by a semi colon ; This is a nessesary step when using breakseq2 and some other callers 
that rely on specific libraries/files with a new reference. The SVE image has hg19 libraries for breakseq2 included.
[EX] -t delly:/data/exclude_list.bed;breakseq:/data/breakseq2_bplib_20150129*\n\t[None]"""
parser.add_argument('-t', '--target',type=str, help=target_help)
parser.add_argument('-m', '--model',type=str, help='data fusion model\t[FusorSV/data/models/g1k_v37_decoy.P3.pickle]')
parser.add_argument('-P', '--cpus',type=str, help='number of cpus for alignment and sorting, ect\t[1]')
parser.add_argument('-T', '--threads',type=str, help='number of threads per CPU\t[4]')
parser.add_argument('-M', '--mem',type=str, help='ram in GB units to use for processing per cpu/thread unit\t[4]')
parser.add_argument('-a', '--algorithm',type=str, help='alignment/realignment algorithm pipeline\t[speed_seq]')
parser.add_argument('-R', '--realign',type=str, help='does realignment of the BAM file to the reference\t[place holder]')
parser.add_argument('-v', '--verbose',action='store_true', help='stderr and stdout from all stages\t[False]')
args = parser.parse_args()

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
    print('no ref pattern found from: args=%s'%args.ref)
    raise IOError
    ref_fa_path = ''
refbase = ref_fa_path.rsplit('/')[-1].split('.fa')[0]
if args.fqs is not None and all([os.path.exists(f) for f in args.fqs.split(',')]):
    reads = args.fqs.split(',') #CSL returns a list of 1+    
else:
    print('no valid fqs pattern found from: arg=%s'%args.fqs)
    print('fastq files should be comma seperated and end with 1.fq for read 1 or 2.fq for read2 (can be .gz compressed)')
    reads = []
if not all([os.path.exists(r) for r in reads]):
    print('fastq files do not exist at the path specified!')
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
    print('SM parameter not specified, stripping off . and _ char with the extension on read 1 or bam input')
    if len(reads)>0:
        SM = reads[0].rsplit('/')[-1].rsplit('.')[0].rsplit('_')[0]
        print('using SM = %s from fastq file name'%SM)
    elif bam!='':
        SM = bam.rsplit('/')[-1].rsplit('.')[0].rsplit('_')[0]
        print('using SM = %s from bam file name'%SM)
    else:    
        print('fastq or bam files were not located!')
        raise IOError
if args.target is not None:
    target_str = args.target #target mapping string, need to parse it
else: #default targets
    target_str = 'breakseq:'+os.path.dirname(os.path.abspath(__file__))+'/../data/breakseq2*.fna'
    #assume genome_strip is in the same directory as referece and ends with .tar.gz
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
    mem = 6
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
    if args.verbose: print(output)
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
    print('BAM files not located:\nstarting alignment of FASTQ files at location: %s'%bam_directory)
    if algorithm=='speed_seq': a_cpus,a_threads,a_mem = 1,cpus*threads,cpus*mem
    else:                      a_cpus,a_threads,a_mem = cpus,threads,mem
    prepare_bam = [scripts_path+'prepare_bam.py','-r',ref_fa_path,'-f',','.join(reads),'-o',bam_directory,'-s',SM,
                   '-P',str(a_cpus),'-T',str(a_threads),'-M',str(a_mem),'-a',algorithm]
    try:
        output += subprocess.check_output(' '.join(prepare_bam),shell=True)
        bam_path = sorted(glob.glob(bam_directory+'%s*.bam'%SM),key=lambda x: len(x))[0] #more permissive
    except Exception as E:
        print(E)
        raise IOError
    if args.verbose: print(output)
else:
    bam_directory = '/'.join(bam.rsplit('/')[:-1])+'/'
    bam_path = bam
    print('BAM files : %s was found, skipping alignment step'%bam)
b_stop = time.time()
if os.path.exists(bam_path):
    print('prepare_bam has completed in %s minutes'%round((b_stop-b_start)/60.0,2))
else:
    print('prepare_bam failure shuting down execution!')
    raise IOError
    
#(C) start as many variant callers as you can (in series ...)
v_start = time.time()
output,i = '',1
vcf_directory = directory+'/vcf/'
if not os.path.exists(vcf_directory): os.mkdir(vcf_directory)
vcf_sample = vcf_directory+'/%s/'%SM
print('starting variant calling')
variant_processor = [scripts_path+'variant_processor.py','-r',ref_fa_path,'-b',bam_path,'-o',vcf_sample,'-s']
#redo this when every you add new callers to create faster workflows
jobs = [['bam_stats','breakdancer','lumpy','breakseq'],
        ['cnmops','cnvnator','genome_strip','hydra'],
        ['gatk_haplo','delly']]
for job in jobs:                                      
    print('dispatching %s of %s jobs'%(i,len(jobs)))   #each group will wait for each other to finish
    p1 = mp.Pool(processes=cpus)                                 #so bam_stats can pass on D L values
    for stage in job:                                         #run one at a time in || and then block
        print('dispatching stage: %s'%stage)
        p1.apply_async(run_variant_processor,
                       args=(variant_processor,stage),
                       callback=collect_results)
        time.sleep(5.0)
    p1.close()
    p1.join()
    print('|| execution pool completed on job %s'%i)
    if args.verbose: print(result_list)
    result_list,i = [],i+1
#(D) Run FusorSV to merge the results
output = ''
print('starting FusorSV processing with a prior model')
fusorsv_out = directory+'/fusorsv_out/'

fusor_sv = [scripts_path+'../../FusorSV/FusorSV.py',
            '-r',ref_fa_path,'-i',vcf_directory,
            '-f',scripts_path+'../../FusorSV/data/models/human_g1k_v37_decoy.P3.pickle',
            '-o',fusorsv_out,'-p',str(cpus),'-L','DEFAULT']
try :
    output += subprocess.check_output(' '.join(fusor_sv),shell=True)
except Exception as E:
    print(E)
    raise IOError
if args.verbose: print(output)

#(E) Run TigraSV with FusorSv output
#get the FusorSV VCF file and set as target
output = ''
targeted_assembly = [scripts_path+'variant_processor.py','-r',ref_fa_path,'-b',bam_path,
                     '-o',vcf_sample,'-s','tigra','-t']
try:
    fusor_vcf = glob.glob(fusorsv_out+'/vcf/*_S*.vcf')[0]
    output += subprocess.check_output(' '.join(targeted_assembly+[fusor_vcf]),shell=True)
except Exception as E:
    print(E)
    raise IOError
if args.verbose: print(output)
##(F) Rerun FusorSV with TigraSV and ctg directory
#try :
#    output += subprocess.check_output(' '.join(fusor_sv),shell=True)
#except Exception as E:
#    print(E)
#    raise IOError
#if args.verbose: print(output)
