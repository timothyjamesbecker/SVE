#!/usr/bin/env python
import argparse
import os
import sys
import csv
import datetime
import subprocess32 as subprocess #to call qsub a bunch of times

des = """
script that auto generates PBS scripts for simple bam splitting and submits them via qsub"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-l', '--ftp_list',type=str, help='path to the csv ftp download list file')
parser.add_argument('-o', '--output_dir',type=str, help='outputdirectory to save ...chrA.bam/...chrA.bam.bai into')
parser.add_argument('-c', '--chroms',type=str, help='comma seperated chrom list matching the bam header')
parser.add_argument('-d', '--data_base',type=str, help='database to log bam_split_x into')
parser.add_argument('-m', '--email_address',type=str,help='cluster email results to this email address')
parser.add_argument('-t', '--time',type=str,help='time to execute one set of m connections: 01:00:00 = 1 hour')
parser.add_argument('-n', '--connections',type=int,help='number of connections to use at once')
args = parser.parse_args()

#add seconds to a time utilizing the local timezone and DST info on the server
def addSecs(t0, secs):
    fulldate = datetime.datetime(100, 1, 1, t0.hour, t0.minute, t0.second)
    fulldate = fulldate + datetime.timedelta(seconds=secs)
    return fulldate.time()

#given the batch number j and time for each batch t
#compute the time to insert in qsub -a in HHMM.SS string format
#t is HH:MM:SS format string, b is a base delay offset for qsub
def get_time_offset(j,t,b=30):
    t0 = datetime.datetime.now().time()
    t = [int(i) for i in t.split(':')] #HH:MM:SS = [1,12,56]
    secs = j*(t[0]*60*60+t[1]*60+t[2])+b #get batch number j multiplied by each batch in seconds
    a = addSecs(t0,secs) #now HHMM.SS for qsube -a
    a = str(a.hour).zfill(2)+str(a.minute).zfill(2)+'.'+str(a.second).zfill(2)
    return a

#sample_name.pbs
#/home/tbecker/software/anaconda/bin/python2.7 /home/tbecker/software/SVE/tests/variant_processor
# -d jax -c 20 -s bam_split -o /data/tbecker/sample_name/

#check ftp file list and make it a dict sample_name:ftp-URL 
if args.ftp_list is not None:
    ftp_list = args.ftp_list
    samples = {}
    with open(ftp_list, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader: #row[0] was just an index, row[1] is sample name, row[2] is link
            samples[row[1]] = row[2]
else:
    print('ftp list not formed correctly')
    raise IOError

if args.output_dir is not None:
    out_dir = args.output_dir
    if not os.path.exists(out_dir): os.makedirs(out_dir)
else:
    print('out put directory not specified')
    raise IOError

if args.email_address is not None:
    email = args.email_address
else:
    print('email not valid')
    raise AttributeError

if args.data_base is not None:
    dbc = args.data_base
else:
    print('no data base specified')

if args.chroms is not None:
    chroms = args.chroms
else:
    print('chroms not listed')
    raise AttributeError

if args.time is not None:
    time = args.time
else:
    print('time for running m connections not specified')
    raise AttributeError

if args.connections is not None:
    connections = args.connections
else:
    print('number of connections at once not specified')
    raise AttributeError
    
#make a temp directory for pbs scripts
pbs_dir = out_dir+'PBS'
if not os.path.exists(pbs_dir): os.makedirs(pbs_dir)
    
#write out the .pbs scripts
PBS = []
python = '/home/tbecker/software/anaconda/bin/python'
variant_processor = '/home/tbecker/software/SVCP/tests/variant_processor.py'
for k in samples:
    sample_pbs = pbs_dir+'/'+k+'.pbs' #name of pbs script for sample k
    PBS += [sample_pbs]
    with open(sample_pbs,'w') as pbs:
        #variant processor string, samples[k] is the ftp.bam file name to use
        cut = [python,variant_processor,'-d',dbc,'-b',samples[k],'-c',chroms,
               '-s','bam_split_simple','-o',out_dir+k+'/']
        pbs.write('#!/bin/bash\n'+' '.join(cut))
#execute qsub with the scripts, getting the jids back (can display these)
output,err = '',{}
n,i,j = len(samples),1,0
for pbs in PBS:
    print('doing %s'%pbs)
    try:
        command = ['qsub','-l','walltime=8:00:00,mem=2gb','-a',get_time_offset(j,time),
                   '-m','e','-M',email,'-o','split.log','-j oe',pbs]
        print(' '.join(command)) #don't run it yet!
        output += subprocess.check_output(' '.join(command),stderr=subprocess.STDOUT,shell=True)
    #catch all errors that arise under normal call behavior
    except subprocess.CalledProcessError as E:
        print('call error: '+E.output)        #what you would see in the term
        err['output'] = E.output
        #the python exception issues (shouldn't have any...
        print('message: '+E.message)          #?? empty
        err['message'] = E.message
        #return codes used for failure....
        print('code: '+str(E.returncode))     #return 1 for a fail in art?
        err['code'] = E.returncode
    except OSError as E:
        print('os error: '+E.strerror)        #what you would see in the term
        err['output'] = E.strerror
        #the python exception issues (shouldn't have any...
        print('message: '+E.message)          #?? empty
        err['message'] = E.message
        #the error num
        print('code: '+str(E.errno))
        err['code'] = E.errno
    i += 1
    if i%connections == 0: j += 1
print('output:\n'+output)
#we will simply use variant_processor -s bam_split here!

#remove/delete the intermediate pbs scripts TO DO...