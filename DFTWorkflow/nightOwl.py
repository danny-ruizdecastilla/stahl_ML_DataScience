import subprocess
import time
from pathlib import Path
#Danny Ruiz de Castilla | 10.08.2025
#Monitors and submits jobs on kestrel when the number of jobs dips past 500 queues 
def slurmJobChecker(user):
    result = subprocess.run(["squeue", "-u", user],capture_output=True,text=True,check=True)
    lines = result.stdout.strip().split("\n")
    return len(lines) - 1 if len(lines) > 1 else 0
def kestrelSubmit(fileName , action):
    subprocess.run([action , fileName])
def nightOwl(mainDir, rootName , action ,user ,  outputDir):#submission framework to kestrel when you have more files to run than slots available

    initialFiles = sorted(list(f'*{rootName}'))
    count = 0 
    for file in initialFiles: #submit the first 500
        kestrelSubmit(file , action)
        initialFiles.remove(file)
        count +=1
        if count == 500:
            break
    #now we initialize the process of indefinite checking 
    while True:
        time.sleep(2)
        if len(initialFiles) == 0:
            break
        else:
            numJobs = slurmJobChecker(user)
            if numJobs < 500:
                diff = 500 - numJobs 
                count = 0 
                while True:
                    if count == diff:
                        break
                    else:
                        file = initialFiles[0]
                        kestrelSubmit[file , action]
                        initialFiles.remove(file)
                        count +=1
            