import os 
import sys
import glob
import shutil
import re
import pandas as pd 
import numpy as np

parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(parentDir)
from DFTWorkflow.ionComGenerator import locateinLog
from DFTWorkflow.featureMaping import createCSV

#Danny Ruiz de Castilla
#Generates Fukui maps for each molecule 
def extractEnergies(log: str, energyStr):

    # Energy extraction patterns for Gaussian
    patterns = {
        'electronic': r'SCF Done:.*?E\([^)]+\)\s*=\s*([-\d.]+)',
        'gibbs': r'Sum of electronic and thermal Free Energies=\s*([-\d.]+)',
        'enthalpy': r'Sum of electronic and thermal Enthalpies=\s*([-\d.]+)',
        'zpe': r'Sum of electronic and zero-point Energies=\s*([-\d.]+)'
    }
    if energyStr not in patterns:
        raise ValueError(f"Unknown energy type: {energyStr}. "
                        f"Available types: {list(patterns.keys())}")
    
    try:
        with open(log, 'r') as f:
            content = f.read()
        matches = re.findall(patterns[energyStr], content)
        
        if matches:
            return float(matches[-1])
               
    except FileNotFoundError:
        print(f"Error: File {log} not found")
        return None
    except Exception as e:
        print(f"Error reading {log}: {str(e)}")
        return None
def getBoltzmannWeightsGauss(logDirs, temperature, energyStr):
    R = 1.987204e-3
    HARTREE_TO_KCAL = 627.5094740631
    results = []
    
    for log in logDirs:
        if not os.path.exists(log):
            print(f"Warning: File {log} not found. Skipping.")
            continue        
        try:
            energy = extractEnergies(log, energyStr)
            if energy is not None:
                results.append({'logID': str(log.split("/")[-1].split(".")[0]),  'E_Ha': energy})
            else:
                print(f"Warning: Could not extract {energyStr} energy from {log}")     
        except Exception as e:
            print(f"Error processing {log}: {str(e)}")
            continue
    if len(results) == 0:
        raise ValueError("No energies could be extracted from any log files")
    # Convert to DataFrame
    df = pd.DataFrame(results)
    # Convert energies to kcal/mol
    df['E_kCal'] = df['E_Ha'] * HARTREE_TO_KCAL
    
    # Calculate relative energies (relative to lowest energy)
    minE = df['E_kCal'].min()
    df['rE_kCal'] = df['E_kCal'] - minE
    
    # Calculate Boltzmann factors
    df['boltzmannFacts'] = np.exp(-df['rE_kCal'] / (R * temperature))
    
    # Normalize to get weights (sum to 1)
    normFacts = df['boltzmannFacts'].sum()
    df['boltzWeights'] = df['boltzmannFacts'] / normFacts

    df = df.sort_values('E_kCal').reset_index(drop=True)
    
    return df[['logID', 'E_Ha', 'rE_kCal', 'boltzWeights']]


def boltzmannAveraging(weights, featureList): #same type features returns final DFMAST number
    if len(weights) != len(featureList):
        raise ValueError("Number of weights must equal number of features")
    if abs(sum(weights) - 1.0) > 1e-6:
        print(f"Warning: Weights sum to {sum(weights):.6f}, not 1.0")
    
    return sum(w * f for w, f in zip(weights, featureList))

    
def fukuiFunction(substrateDF):
    populations = ['pop_cat', 'pop_an', 'pop_neut']
    missingCols = [col for col in populations if col not in substrateDF.columns]
    
    if missingCols:
        print(f"Error: Missing required columns: {missingCols}")
        print(f"Available columns: {list(substrateDF.columns)}")
        return None
    try:
        substrateDF.insert(len(substrateDF.columns), "f_neg", 0)
        substrateDF.insert(len(substrateDF.columns), "f_pos", 0) 
        substrateDF.insert(len(substrateDF.columns), "f_neut", 0)
        for index, row in substrateDF.iterrows():
            # f- (nucleophilic attack): q(N) - q(N+1) 
            substrateDF.loc[index, "f_neg"] = row["pop_neut"] - row["pop_cat"]
            # f+ (electrophilic attack): q(N-1) - q(N)
            substrateDF.loc[index, "f_pos"] = row["pop_an"] - row["pop_neut"]
            # f0 (radical attack): 0.5 * [q(N-1) - q(N+1)]
            substrateDF.loc[index, "f_neut"] = 0.5 * (float(row["pop_an"]) - float(row["pop_cat"])) 
                
    except Exception as e:
        print(f"Unexpected error in fukuiFunction: {e}")
        return None
    
    print("Fukui function calculation completed successfully")
    return substrateDF
def extractOccupanciesNBO(logFile:str, extract1:str, extract2:str , location:str):
    lowerInd = locateinLog(logFile , extract1, location)
    upperInd = locateinLog(logFile , extract2, location )
    if upperInd == "Poison" or lowerInd == "Poison":
        raise ValueError(f"The Log file {logFile} did not terminate properly")
    atomOccupancies = {}
    with open(logFile , "r") as f:
        for idx, line in enumerate(f):
            if idx > lowerInd + 1  and idx < upperInd-1:
                occupancies = line.strip()
                occupancyLines = occupancies.split("    ")
                #print(occupancyLines)
                atomNums = occupancyLines[0].split("   ")
                atomInd = atomNums[0].split()
                print(atomInd)
                
                #print(occupancyLines[1:])
                nums = []
                for num in occupancyLines[1:]:
                    num = float(''.join(num.split()))
                    nums.append(num)
                
                #print(np.max(nums))
                atomOccupancies[int(atomInd[-1])] = [str(atomInd[0]) , np.max(nums)]
    #print(atomOccupancies)
    return atomOccupancies
def extractOccupanciesMull(logFile:str, extract1:str, extract2:str , location:str):

    lowerInd = locateinLog(logFile , extract1, location)
    upperInd = locateinLog(logFile , extract2, location )
    atomOccupancies = {}
    with open(logFile , "r") as f:
        for idx, line in enumerate(f):
            #print(lowerInd , upperInd)
            if idx > int(lowerInd) + 1  and idx < int(upperInd)-1:
                newLine = line.strip()
                #print(newLine.split())
                splits = newLine.split()
                try:
                    index = int(splits[1])
                    occupancy = float(splits[4])
                    atomName = str(splits[2])
                    atomOccupancies[index] = [atomName , occupancy]
                except ValueError:
                    atom =  list(atomOccupancies.keys())[-1]
                    atomOccupancies[atom][1] += float(splits[2])
                    
    return atomOccupancies
def getMeanDF(dfList):
    combined = pd.concat(dfList, axis=0)
    numericCols = combined.select_dtypes(include=['float64', 'float32']).columns
    meanDF = combined[numericCols].groupby(combined.index).mean()
    stringCols = combined.select_dtypes(exclude=['float64', 'float32']).columns
    stringDF = combined[stringCols].groupby(combined.index).first()
    finalDF = pd.concat([stringDF, meanDF], axis=1)
    return finalDF

def main(logDir , cationDir, anionDir , substrateCSV, outputDir , densityStr ):
    substrateScope = pd.read_csv(substrateCSV, encoding='utf-8')
    smileString = substrateScope["SMILES"]
    substrateName = substrateScope["ID"]
    neutralLogs = glob.glob(logDir + "/*.log")
    logsDicts = {}

    for log in neutralLogs:
        logList = []
        logName = log.split("/")[-1].split(".")[0]
        logList.append(log)
        if os.path.exists(cationDir + "/" + logName + "_cation.log"):
            logList.append(cationDir + "/" + logName + "_cation.log")
        else:
            print("Missing the cation .log File for " + logName)
            continue
        if os.path.exists(anionDir + "/" + logName + "_anion.log"):
            logList.append(anionDir + "/" + logName + "_anion.log")     
         
        else:
            print("Missing the anion .log File for " + logName)
            continue
        if len(logList) == 3:
            logsDicts[logName] = logList
    substrateDictMAST = { (name, smile): [] for name, smile in zip(substrateName, smileString) }
    for substrate, dirs in logsDicts.items():
        substrateDF = pd.DataFrame()
        atomCols = []
        for dir in dirs:
            name = dir.split("/")[-1]
            #print(name)
            if "ion" in name:
                name1 = name.split("_")[-1].split("ion")[0]
            else:
                name1 = "neut"
            if densityStr == "NBO7":
                occupanciesDict = extractOccupanciesNBO(dir,"Core      Valence    Rydberg" , "* Total *" , "earliest")
            elif densityStr == "mull":
                occupanciesDict = extractOccupanciesMull(dir,"Gross orbital populations:" , "Condensed to atoms (all electrons):" , "earliest")
            if not "atoms" in substrateDF.columns:
                substrateDF["atoms"] = list(map(int, occupanciesDict.keys()))
            substrateDF["pop_" + name1] = [values[-1] for values in occupanciesDict.values()]
            substrateDF["atom_" + name1] = [values[0] for values in occupanciesDict.values()]
            atomCols.append("atom_" + name1)
        if ((substrateDF[atomCols[0]] == substrateDF[atomCols[1]]) & (substrateDF[atomCols[1]] == substrateDF[atomCols[2]])).all():
            substrateDF = fukuiFunction(substrateDF)
            if "clust" in substrate:
                substrate_ = substrate.split("_clust")[0]
            else:
                substrate_ = substrate.split("_conf")[0]
            for key in substrateDictMAST:
                if substrate_ in key[0]:  # key[0] is the name
                    newKey= key + (substrate_,)
                    substrateDictMAST[newKey] = substrateDictMAST.pop(key)
                    substrateDictMAST[newKey].append(substrateDF)
                    break
        else:
            raise Exception(f"Atom Types on the 3 indexes for {substrate} do not line up")
    for identification , dfList in substrateDictMAST.items():
        logPaths = [log for log in neutralLogs if identification[0] in log]
        confList = list(identification[2:])
        if len(confList) == len(dfList):
            weightedList = []
            boltzmannDF = getBoltzmannWeightsGauss(logPaths, 298, 'electronic')
            boltzDict = boltzmannDF.set_index('logID')['boltzWeights'].to_dict()
            for index , conf in enumerate(confList):
                for key in boltzDict:
                    if conf in key:  # conf is a substring of the key
                        weight = boltzDict[key]
                        df = dfList[index]
                        numericDF = df.select_dtypes(include='number') * weight * len(confList)
                        nonNumericDF = df.select_dtypes(exclude='number')
                        finalDF = pd.concat([nonNumericDF, numericDF], axis=1)
                        weightedList.append(finalDF) 
            finalDF = getMeanDF(weightedList)
            if not os.path.exists(outputDir + "/" + str(identification[0])): 
                os.makedirs(outputDir + "/" + str(identification[0]))
            with open(outputDir + "/" + str(identification[0]) + "/identification.dat", "w") as f:
                f.write(f"{identification[0]},{identification[1]}")
            createCSV(finalDF ,outputDir + "/" + str(identification[0]), str(identification[0]))
        else:
            raise ValueError("Number of conformers does not match the number of Fukui dataframes.")


if __name__ == "__main__":
    neutralDir = str(sys.argv[1])
    cationDir = str(sys.argv[2])
    anionDir = str(sys.argv[3])
    substrateCSV = str(sys.argv[4])
    outputDir = str(sys.argv[5])
    densityStr = str(sys.argv[6])
    main(neutralDir , cationDir , anionDir , substrateCSV , outputDir , densityStr)

