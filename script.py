# script.py
import pandas as pd
def say_hello():
    print("Hello from Python!")

def add(a, b, c):
    return a + b + c

def getData():
    instanceData = pd.read_excel('Data/ev_session_data_OR.xlsx')
    instanceData = instanceData[:10]
    return instanceData

def printdata(instanceData):
    print('How data is formatted:\n', instanceData)