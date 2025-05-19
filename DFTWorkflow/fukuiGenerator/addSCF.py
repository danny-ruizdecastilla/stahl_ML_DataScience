import os 
import sys
import glob
import shutil
#Danny Ruiz de Castilla
##takes in any error files and rewrites to add scf=qc to the setting line
def readDat(filename):
    with open(filename, 'r') as file:
        line = file.readline().strip()
        return line.split(',')

def addtoCom(comFile, phrase, lineInd, phraseInd):
    with open(comFile , 'r') as file:
        lines = file.readlines()
    for i , line in enumerate(lines):
        if i == lineInd:
            lineMAST = ""
            allFracts = line.split(" ")
            for j , fract in enumerate(allFracts):
                if j == phraseInd:
                    lineMAST += phrase
                lineMAST += fract + " "
            lines[lineInd] = lineMAST
    return lines
def makeMoves(ionPhrase, errFile):
    errFiles = readDat(errFile)
    for file in errFiles:

        name = file.split(ionPhrase)[0]
        chks = glob.glob("chks/" + name+ "*")
        for chk in chks:
            if "start" in chk:
                shutil.move(chk ,  name + ".chk")
            else:
                os.remove(chk)
        logs = glob.glob("logs/" + name+ "*")
        for log in logs:
            os.remove(log)
        comName = file.split(".")[0] + ".com"
        comLines = addtoCom("finished/" + comName , "scf=qc ", 3, 5)
        with open(comName, 'w') as file:
            for line in comLines:
                file.write(line)
        os.remove("finished/" + comName)
if __name__ == "__main__":
    ion = str(sys.argv[1])
    errFile = str(sys.argv[2])
    makeMoves(ion,errFile)
