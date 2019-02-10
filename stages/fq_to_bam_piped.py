import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper
import stage_utils as su

#function for auto-making svedb stage entries and returning the stage_id
class fq_to_bam_piped(stage_wrapper.Stage_Wrapper):
    #path will be where a node should process the data using the in_ext, out_ext
    #stage_id should be pre-registered with db, set to None will require getting
    #a new stage_id from the  db by writing and registering it in the stages table
    def __init__(self,wrapper,dbc,retrieve,upload,params):
        #inheritance of base class stage_wrapper    
        stage_wrapper.Stage_Wrapper.__init__(self,wrapper,dbc,retrieve,upload,params)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        return 0  

    def get_records(self,m,t):
        r = 250000*16
        try:
            r = int((int(m)*250000)/int(t))
        except Exception as E:
            pass
        return r

    #override this function in each wrapper...
    #bwa sampe ref.fa L.sai R.sai L.fq R.fq -f out.sam
    def run(self,run_id,inputs):
        #workflow is to run through the stage correctly and then check for error handles
        #[1b]
        in_names  = {'.fa':inputs['.fa'][0],'.fq':inputs['.fq']}
        out_ext = self.split_out_exts()[0]
        #add tigra.ctg.fa bam by allowing a single '.fa' input key
        if len(in_names['.fq']) == 2: csl = su.get_common_string_left(in_names['.fq'])
        else:                         csl = in_names['.fq'][0]
        stripped_name = self.strip_path(csl)
        if inputs.has_key('out_dir'):
            out_dir = inputs['out_dir'][0]
            out_name = out_dir+stripped_name
            out_name = out_name.rstrip('_')
        else: #untested...
            right = in_names['.fa'][0].rsplit('/')[-1]
            out_dir = in_names['.fq'][0].replace(right,'')
            cascade = self.strip_in_ext(in_names['.fq'][0],'.fq')
            out_name = cascade
            out_name = out_name.rstrip('_')
        if inputs.has_key('SM'):
            SM = inputs['SM'][0]
            print('found SM tag = %s'%SM)
        else:
            SM = stripped_name
            print('missing a SM tag, using %s'%SM)
            
        #[2]build command args
        #:::TO DO::: ALLOW THE USER TO GENERATE AND ENTER RG
        if not os.path.exists(out_dir): os.makedirs(out_dir)
        #if not os.path.exists(out_dir+'/sort/'): os.makedirs(out_dir+'/sort/')
        threads = str(self.get_params()['-t']['value'])
        bwa = self.software_path+'/bwa-master/bwa' #latest release
        samtools = self.software_path+'/samtools-1.3/samtools'
        sambamba = self.software_path+'/sambamba/sambamba'
        sample = stripped_name+'RG'

        #'@RG\tID:H7AGF.2\tLB:Solexa-206008\tPL:illumina\tPU:H7AGFADXX131213.2\tSM:HG00096\tCN:BI'
        RG = r'\t'.join(["'@RG",'ID:'+sample,'LB:'+'Solexa'+sample,'PL:'+inputs['platform_id'][0],
                         'PU:'+sample,'SM:'+SM+"'"])
        bwa_mem = [bwa,'mem','-M','-t',threads,'-R',RG,in_names['.fa']]+in_names['.fq']+['|']
        view = [samtools,'view','-Shb','-','>',out_name+'.bam']
        sort   =  [sambamba,'sort','-l','9','-m',str(self.get_params()['-m']['value'])+'GB',
                   '-t',threads,'-o',out_name+'.sorted.bam','--tmpdir=%s'%out_dir+'/sort', out_name+'.bam']
        mark   =  [sambamba,'markdup','-l','9','-t',threads,'--tmpdir=%s'%out_dir+'/mark',
                   out_name+'.sorted.bam',out_name+'.mark.bam']
        clean = ['rm','-rf',out_dir+'/sort',out_dir+'/mark']

        self.command = bwa_mem+view
        print(self.get_command_str())
        self.db_start(run_id,in_names['.fq'][0])
        #[3a]execute the commands here----------------------------------------------------
        output,err = '',{}
        bam_size = 0
        try:
            if not os.path.exists(out_name+'.bam') or os.path.getsize(out_name+'.bam') <= 4096:
                output += subprocess.check_output(' '.join(bwa_mem+view),stderr=subprocess.STDOUT,shell=True)
            bam_size = os.path.getsize(out_name+'.bam')
        except subprocess.CalledProcessError as E:
            print('call error: '+E.output)
            err['output'] = E.output
            print('message: '+E.message)
            err['message'] = E.message
            print('code: '+str(E.returncode))
            err['code'] = E.returncode
        except OSError as E:
            print('os error: '+E.strerror)
            err['output'] = E.strerror
            print('message: '+E.message)
            err['message'] = E.message
            print('code: '+str(E.errno))
            err['code'] = E.errno
        except Exception as E:
            print(E)
        print('output:\n'+output)

        sorted_size = 0 #sorting step
        try:
            if not os.path.exists(out_name + '.bam.bai') :
                if not os.path.exists(out_name + '.sorted.bam'):
                    output += subprocess.check_output(' '.join(sort), stderr=subprocess.STDOUT, shell=True)
                elif os.path.getsize(out_name + '.sorted.bam') <= bam_size/2.0:
                    output += subprocess.check_output(' '.join(['rm',out_name+'.sorted.bam*']),
                                                      stderr=subprocess.STDOUT, shell=True)
                    output += subprocess.check_output(' '.join(sort), stderr=subprocess.STDOUT, shell=True)
            sorted_size = os.path.getsize(out_name+'.sorted.bam')
        except subprocess.CalledProcessError as E:
            print('call error: ' + E.output)
            err['output'] = E.output
            print('message: ' + E.message)
            err['message'] = E.message
            print('code: ' + str(E.returncode))
            err['code'] = E.returncode
        except OSError as E:
            print('os error: ' + E.strerror)
            err['output'] = E.strerror
            print('message: ' + E.message)
            err['message'] = E.message
            print('code: ' + str(E.errno))
            err['code'] = E.errno
        except Exception as E:
            print(E)
        print('output:\n' + output)

        mark_size = 0
        try:
            if os.path.exists(out_name+'.sorted.bam') and sorted_size > bam_size/2.0:
                output += subprocess.check_output(' '.join(['rm',out_name+'.bam']),
                                                  stderr=subprocess.STDOUT,shell=True)
                output += subprocess.check_output(' '.join(mark),stderr=subprocess.STDOUT,shell=True)
            mark_size = os.path.getsize(out_name+'.mark.bam')
        except subprocess.CalledProcessError as E:
            print('call error: ' + E.output)
            err['output'] = E.output
            print('message: ' + E.message)
            err['message'] = E.message
            print('code: ' + str(E.returncode))
            err['code'] = E.returncode
        except OSError as E:
            print('os error: ' + E.strerror)
            err['output'] = E.strerror
            print('message: ' + E.message)
            err['message'] = E.message
            print('code: ' + str(E.errno))
            err['code'] = E.errno
        except Exception as E:
            print(E)
        print('output:\n' + output)

        final_size = 0
        try:
            if os.path.exists(out_name+'.mark.bam') and mark_size > bam_size/2.0:
                output += subprocess.check_output(' '.join(['rm',out_name+'.sorted.bam*']),
                                                  stderr=subprocess.STDOUT,shell=True)
                output += subprocess.check_output(' '.join(['mv',out_name+'.mark.bam',out_name+'.bam']),
                                                  stderr=subprocess.STDOUT,shell=True)
                output += subprocess.check_output(' '.join(['mv',out_name+'.mark.bam.bai',out_name+'.bam.bai']),
                                                  stderr=subprocess.STDOUT, shell=True)
                output += subprocess.check_output(' '.join(clean), stderr=subprocess.STDOUT,shell=True)
            final_size = os.path.getsize(out_name+'.bam')
        except subprocess.CalledProcessError as E:
            print('call error: ' + E.output)
            err['output'] = E.output
            print('message: ' + E.message)
            err['message'] = E.message
            print('code: ' + str(E.returncode))
            err['code'] = E.returncode
        except OSError as E:
            print('os error: ' + E.strerror)
            err['output'] = E.strerror
            print('message: ' + E.message)
            err['message'] = E.message
            print('code: ' + str(E.errno))
            err['code'] = E.errno
        except Exception as E:
            print(E)
        print('output:\n' + output)

        #[3b]check results--------------------------------------------------
        if err == {}:
            self.db_stop(run_id,{'output':output},'',True)
            results = [out_name+'.bam']
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("sucessfull........")
                return results   #return a list of names
            else:
                print("failure...........")
                return False
        else:
            self.db_stop(run_id,{'output':output},err['message'],False)
            return None