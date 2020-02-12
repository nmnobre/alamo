#!/usr/bin/env python3
import os
import re
from glob import glob
import subprocess
from datetime import datetime
import filecmp
import argparse
import configparser
import sqlite3
import simba
import ansi2html
import random



timestamp = datetime.today().strftime('%Y-%m-%d_%H%M%S')

parser = argparse.ArgumentParser(description='Sift through outputs')
parser.add_argument('inifile', help='Configuration file')
parser.add_argument('--benchmark',action='store_true',default=False,help='Set this run as benchmark for all tests')
args = parser.parse_args()

db = sqlite3.connect('regtest.db')
db.text_factory = str
cur= db.cursor()
types = dict()

config = configparser.ConfigParser()
config.read(args.inifile)
print(config)
print(config.sections())
for s in config.sections():
    simba.updateTable(cur,s,types,verbose=False)
    print(config[s])
    for r in config[s]:
        print(r,config[s][r])

alamo_path = os.path.abspath('.')
alamo_configure_flags = ''
regtest_dir = os.path.abspath('.')
branches = ['']
nprocs_build = 1
benchmark_run = None
if 'main' in config:
    if 'alamo_path'            in config['main']: alamo_path = os.path.abspath(config['main']['alamo_path'])
    if 'alamo_configure_flags' in config['main']: alamo_configure_flags = config['main']['alamo_configure_flags']
    if 'regtest_dir'           in config['main']: regtest_dir = os.path.abspath(config['main']['regtest_dir'])
    if 'branches'              in config['main']: branches = config['main']['branches'].split(' ')
    if 'nprocs_build'          in config['main']: nprocs_build = int(config['main']['nprocs_build'])

simba.updateRegTestTable(cur,verbose=False)

conv = ansi2html.Ansi2HTMLConverter()

for branch in branches:

    returncode = 0
    run_id = timestamp
    build_stdout = ""

    if branch != '':
        run_id += "-" + branch

        ret = subprocess.run(['make','realclean'],cwd=alamo_path,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        build_stdout += ret.stdout.decode()
        if (ret.returncode):
            print("\033[31m"+ret.stdout.decode())
            print(ret.stderr.decode()+"\033[0m")
            raise(Exception("There was an error cleaning"))

        ret = subprocess.run(['git','checkout',branch],cwd=alamo_path,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        build_stdout += ret.stdout.decode()
        if (ret.returncode):
            print("\033[31m"+ret.stdout.decode())
            print(ret.stderr.decode()+"\033[0m")
            raise(Exception("There was an error checking out branch \""+branch+"\":"))

        ret = subprocess.run(['git','pull'],cwd=alamo_path,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        build_stdout += ret.stdout.decode()
        if (ret.returncode):
            print("\033[31m"+ret.stdout.decode())
            print(ret.stderr.decode()+"\033[0m")
            raise(Exception("There was an error pulling"))

        # Configure 2D
        ret = subprocess.run(['./configure','--dim=2']+alamo_configure_flags.split(),cwd=alamo_path,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        build_stdout += ret.stdout.decode()
        print(ret.stdout.decode())
        simba.updateRegTestRun(cur,run_id,ret.returncode,conv.convert(build_stdout))
        db.commit()
        if (ret.returncode): 
            print("Encountered error configuring {} in 2D".format(branch))
            continue
        
        # Compile 2D
        ret = subprocess.run(['make','-j{}'.format(nprocs_build)],cwd=alamo_path,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        build_stdout += ret.stdout.decode()
        print(ret.stdout.decode())
        simba.updateRegTestRun(cur,run_id,ret.returncode,conv.convert(build_stdout))
        db.commit()
        if (ret.returncode): 
            print("Encountered error making {} in 2D".format(branch))
            continue

        # Configure 3D
        ret = subprocess.run(['./configure','--dim=3']+alamo_configure_flags.split(),cwd=alamo_path,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        build_stdout += ret.stdout.decode()
        print(ret.stdout.decode())
        simba.updateRegTestRun(cur,run_id,ret.returncode,conv.convert(build_stdout))
        db.commit()
        if (ret.returncode): 
            print("Encountered error configuring {} in 3D".format(branch))
            continue
        
        # Compile 3D
        ret = subprocess.run(['make','-j{}'.format(nprocs_build)],cwd=alamo_path,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        build_stdout += ret.stdout.decode()
        print(ret.stdout.decode())
        simba.updateRegTestRun(cur,run_id,ret.returncode,conv.convert(build_stdout))
        db.commit()
        if (ret.returncode): 
            print("Encountered error making {} in 3D".format(branch))
            continue

    simba.updateRegTestRun(cur,run_id,0,conv.convert(build_stdout))
    db.commit()

    if 'benchmark_run' in config['main']:
        benchmark_run = config['main']['benchmark_run']
    else:
        benchmark_run = None


    #continue

    for test in config.sections():
        if test in ['main']: continue

        status = simba.Status()

        #
        # Settings from Config File
        # 

        if 'input' in config[test]: input_file = config[test]['input']
        else:                       input_file = "tests/"+test+"/input"

        if 'dim' in config[test]:   dim = int(config[test]['dim'])
        else:                       dim = 3

        if 'nprocs' in config[test]: nprocs = int(config[test]['nprocs'])
        else:                        nprocs = 1



        #test_dir = "tests/" + test['name'] + "/"
        rt_plot_dir = regtest_dir + "/rt-" + run_id + "/" + test + "/"
        print("------------- rt_plot_dir: ", rt_plot_dir)

        if 'benchmark_run' in config[test]:
            benchmark_run = config[test]['benchmark_run']

        if benchmark_run:
            bm_plot_dir = regtest_dir + "/rt-" + benchmark_run + "/" + test + "/"
            print("------------- bm_plot_dir: ", bm_plot_dir)

        subprocess.run(["mkdir", "-p", rt_plot_dir])

        print(rt_plot_dir+"/output/")
        print("Running ...... ")
        ret = subprocess.run(["mpirun", "-np", str(nprocs), "./bin/alamo-{}d-g++".format(dim), input_file, "plot_file={}/output".format(rt_plot_dir)],
                             cwd=alamo_path,
                             stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        status.runcode = ret.returncode
        print(rt_plot_dir+"/output/")
        if not os.path.isdir(rt_plot_dir+"/output/"):
            print("CREATING DIRECTORY!\n")
            print(status.runcode)
            print(os.path.isdir(rt_plot_dir+"/output/"))
            subprocess.run(["mkdir","-p",rt_plot_dir+"/output/"])
            with open(rt_plot_dir+"/output/metadata","w") as f:
                f.write("Simulation_run_time = 0\n")
                f.write("HASH = " + ''.join(random.choice('0123456789') for i in range(20))+'\n')
        with open(rt_plot_dir+"/output/stdout","w") as f: f.writelines(conv.convert(ret.stdout.decode()))
        with open(rt_plot_dir+"/output/stderr","w") as f: f.writelines(conv.convert(ret.stderr.decode()))

        print("[PASS]" if status.runcode == 0 else "[FAIL]")


        #
        # Do a direct file-by-file comparison
        #
        if benchmark_run:
            print("Comparing: [", rt_plot_dir, "] <==> [", bm_plot_dir, "]")
            match = True
            for rt in sorted(glob(rt_plot_dir+"/output/**",recursive=True)):
                if not os.path.isfile(rt): continue
                if os.path.basename(rt) in ["output", "metadata", "diff.html", "stdout", "stderr"]: continue
                bm = rt.replace(rt_plot_dir,bm_plot_dir)
                if not os.path.isfile(bm):
                    print("Error - mismatched files")
                    match = False
                    break
                if not filecmp.cmp(bm,rt):
                    print(bm, " does not match ", rt)
                    match = False
                    break
            if (match) : print("OK - files match")
            else: print("Error - files do not match")
            if (match) : status.compare = "YES"
            else : status.compare = "NO"
        else:
            status.compare = "NONE"

        # Get timing
        rt_metadata_f = open(rt_plot_dir+"/output/metadata","r")
        rt_run_time = float(re.findall("Simulation_run_time = (\S*)","".join(rt_metadata_f.readlines()))[0])
        status.runtime = rt_run_time

        #
        # Get metadata and check timing
        #
            
        if benchmark_run:
            bm_metadata_f = open(bm_plot_dir+"/output/metadata","r").readlines()
            bm_hash       = re.findall("HASH = (\S*)","".join(bm_metadata_f))[0]
            print("Benchmark hash is ", bm_hash)
            bm_run_time   = float(re.findall("Simulation_run_time = (\S*)","".join(bm_metadata_f))[0])
            status.bm_runtime = bm_run_time
        else:
            bm_hash = "NONE"
            status.bm_runtime = 0

        data = simba.parse(rt_plot_dir+'/output')
        types = simba.getTypes(data)
        simba.updateTable(cur,test,types,verbose=False)
        simba.updateRecord(cur,test,data,verbose=False)

        simba.updateRegTestRecord(cur,data['HASH'],run_id,test,status,bm_hash,benchmark_run)
        db.commit()


db.close()
#exit()
    



    
