# script.py
import pandas as pd
def say_hello():
    print("Hello from Python!")

def add(a, b, c):
    return a + b + c

def getData():
    instanceData = pd.read_csv('Data/DEMSdata_FOCS_v1.csv')
    instanceData = instanceData[:10]
    #instanceData = instanceData.T
    return instanceData

def printdata(instanceData):
    print('How data is formatted:\n', instanceData)

def WriteToDataframe(data, names):
    #print('DATA :\n', data)
    instanceData = pd.DataFrame(data, names)
    instanceData = instanceData.T
    #instanceData = instanceData[:10]
    return instanceData