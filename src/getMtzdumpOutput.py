import subprocess

def runMtzdump(commandString):
    process = subprocess.Popen(commandString, stdout=subprocess.PIPE, shell=True) #Run system command
    output = process.communicate() # interact with the process to get the data from the log.
    return output[0] #return the log file
