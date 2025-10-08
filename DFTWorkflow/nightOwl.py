import os
import glob
import shutil
import time
import subprocess
import sys
#Danny Ruiz de Castilla | 10.08.2025
#Monitors slurm job queues while you sleep
def slurmJobChecker(user):
    result = subprocess.run(["squeue", "-u", user], capture_output=True, text=True, check=True)
    lines = result.stdout.strip().split("\n")
    return len(lines) - 1 if len(lines) > 1 else 0

def kestrelSubmit(fileName, action):
    subprocess.run([action, fileName])

def basicTerm(logFile, errorPhrase, termPhrase):
    termCount = 0
    with open(logFile, 'r') as f:
        for line in f:
            if errorPhrase in line:
                # Log the error line to a file
                with open("errorLog.dat", 'a') as log:
                    log.write(f"{logFile}: {line}")
                return True
            elif termPhrase in line:
                termCount += 1
    return termCount == 0

def nightOwl(action, user, outputDir):  
    os.makedirs(outputDir, exist_ok=True)

    initialFiles = sorted(glob.glob("*.log"))
    count = 0
    for file in list(initialFiles):  
        kestrelSubmit(file, action)
        initialFiles.remove(file)
        count += 1
        if count == 500:
            break

    while True:
        time.sleep(2)

        if len(initialFiles) == 0:
            break

        numJobs = slurmJobChecker(user)
        if numJobs < 500:
            diff = 500 - numJobs
            count = 0
            while count < diff and len(initialFiles) > 0:
                file = initialFiles.pop(0)  
                kestrelSubmit(file, action)  
                count += 1
        termsList = sorted(glob.glob("*.chk"))
        if len(termsList) != 0: 
            for file in termsList:
                fileName = os.path.splitext(file)[0]
                logFile = f"{fileName}.log"

                falseTerm = basicTerm(logFile, "Error termination", "Normal termination")
                allFiles = glob.glob(f"{fileName}*")

                if not falseTerm:
                    for f in allFiles:
                        shutil.move(f, os.path.join(outputDir, os.path.basename(f)))
                else:
                    for f in allFiles:
                        if f.endswith(".chk") or f.endswith(".log"):
                            os.remove(f)

if __name__ == "__main__":
    action = sys.argv[1]
    user = sys.argv[2]
    outputDir = sys.argv[3]
    nightOwl(action, user, outputDir)
