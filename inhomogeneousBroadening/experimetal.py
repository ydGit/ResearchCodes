#!/usr/bin/python
from threading import Thread
from subprocess import call

def runJob(params):
    print "Starting a new thread: "
    call("./broadlinDipole2D "+params, shell=True)
    print "done."

parametersFile = open("parameters_.txt", 'r')
paramString = "0"
while (paramString != ''):
    paramString =  parametersFile.readline()
    paramString = paramString.rstrip('\n')
    if (len(paramString.split(' ')) == 11):
        newThread = Thread(target=runJob, args=(paramString, ))
        newThread.start()

parametersFile.close()
