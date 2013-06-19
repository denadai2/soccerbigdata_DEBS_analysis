import numpy
import string
import select
import sys
import matplotlib.pyplot as plt
import subprocess
import uuid
import datetime
import csv
import math
import copy
import thread
import pylab
import exceptions

fileDataset = None


keyTraj = None #players enqueued from the user
epsilonTraclus = None
bMatchEnd = False
gameInterruptions = []

# globals for speed  #
effectiveStart = None
varianceSpeedPlayers = {}
lLastSpeed = {}
lLastAvgTime = {}
######################

# globals for the speed plotting 
bPlotSpeed = False
TimeWindowSize = None
listOfLine = {}
ax = None
manager = None
timer = None
lAvgPlayers = {}
lLastTimeStampSpeed= 0
contSpeed = 0
######################

lastTimeStampStoreTraj = {}


def plotAreaAndMidfield(axSpec,fig):
    
    fieldX = [-52.447, 52.477, 52.477,52.477,52.477,-52.447,-52.447,-52.447,0.,0.]
    fieldY = [33.960,33.960,33.960,-33.960,-33.960,-33.960,-33.960,33.960,33.960,-33.960]
    
    area1x = [((52.477/2)-3.22),((52.477/2)-3.22),((52.477/2)+3.22),((52.477/2)+3.22)]
    area1y = [-33.960,(-33.960+4),(-33.960+4),-33.960]
    
    area2x = [((52.477/2)-3.22),((52.477/2)-3.22),((52.477/2)+3.22),((52.477/2)+3.22)]
    area2y = [33.960,(33.960-4),(33.960-4),33.960]
    
    lineMidFieldx = [0,(52.477)]
    lineMidFieldy = [0,0]
    
    
    circleMidField=plt.Circle(((52.477/2),0),4,color='b',fill=False)
    
    fig.gca().add_artist(circleMidField)
    axSpec.plot(area1x,area1y,'b')
    axSpec.plot(area2x,area2y,'b')
    axSpec.plot(lineMidFieldx,lineMidFieldy,'b')
    
    
def updateAvgVarStdDevPlayers(playerKey):
     
    buckets = []
    
    if varianceSpeedPlayers.has_key(playerKey):
        buckets = varianceSpeedPlayers[playerKey]
    else:
        return None

    lavg = []
    if lAvgPlayers.has_key(playerKey):
        lavg = lAvgPlayers[playerKey]
    lavg.append(buckets[len(buckets)-1][5])
    
    if (len(lavg) > 20):
        lavg[0] = "$"
        lavg.remove(lavg[0])   
    lAvgPlayers[playerKey] = lavg
    return buckets[len(buckets)-1][6], buckets[len(buckets)-1][7] 
  

def merge2VarianceBucketStar(buckets1, buckets2):
    
    n = buckets1[1] +  buckets2[1]
    av = (((float(buckets1[2]) * float(buckets1[1])) + (float(buckets2[2]) * float(buckets2[1]))) / float(n))
    Vij = float(float(buckets1[3]) + float(buckets2[3]) + (((float(buckets1[1]) *float(buckets2[1]))/(float(buckets1[1]) + float(buckets2[1])))  *  float(math.pow(((buckets1[2] - buckets2[2])),2))) )
    
    return [buckets1[0],n,av,Vij,"**"]
   
    
def merge2VarianceBucket(buckets1, buckets2):
    # timestamp, number of objects, average mu, variance, Bt-1*
    
    
    n = buckets1[1] +  buckets2[1]
    av = (((float(buckets1[2]) * float(buckets1[1])) + (float(buckets2[2]) * float(buckets2[1]))) / float(n))
    Vij = float(float(buckets1[3]) + float(buckets2[3]) + (((float(buckets1[1]) *float(buckets2[1]))/(float(buckets1[1]) + float(buckets2[1])))  *  float(math.pow(((buckets1[2] - buckets2[2])),2))) )
    
    newLastHistory = merge2VarianceBucketStar(buckets2,buckets2[4])

    return [buckets1[0],n,av,Vij,newLastHistory]

def loadSpeedPerformance(playerKey,speedVX,timestamp, T_window, ncomponent):

    bFullWindow = False
    global lLastTimeStampSpeed
    buckets = None
    if varianceSpeedPlayers.has_key(playerKey):
        buckets = copy.copy(varianceSpeedPlayers[playerKey])
    else:
        print "Start to store the speed performance of player: "+playerKey
        buckets = []
    
    # step 1 of the algorithm Motwani & Co.
    
    if (len(buckets) == 0):
        # timestamp, number of objects, average mu, variance, Bt-1*
        buckets.append([timestamp,1,speedVX,0,None])
    else:
        if (buckets[0][2] == speedVX):
            buckets[1][1] = buckets[1][1] +1
        else:
            buckets.insert(0,[timestamp,1,speedVX,0,None]) 
    
    # step 2 of the algorithm Motwani & Co.
    
    ind = len(buckets)-1
    oldestBucket = buckets[ind]
    
    diffSec = (timestamp - oldestBucket[0]) * math.pow(10,-12)
    
    if (diffSec > T_window and ind > 0):
        #it maintains Bm-1* and then m-1 becomes m
        ## bucket has expired
        bFullWindow = True
        newLastHistory = merge2VarianceBucketStar(oldestBucket[4],buckets[ind-1])
        buckets[ind-1][4] = newLastHistory
        
        lastBucketMerge = buckets[0]
        for i in xrange(1,(ind-1)):
            lastBucketMerge = merge2VarianceBucketStar(lastBucketMerge,buckets[i])
            buckets[i][4] = copy.copy(lastBucketMerge) 

        buckets.pop()
        ind = len(buckets)-1
        #print "window full"
    else:
        
        if(ind > 0):
            #it maintains Bm*
            lastBucketMerge = buckets[0]
            for i in xrange(1,(ind)):
                lastBucketMerge = merge2VarianceBucketStar(lastBucketMerge,buckets[i])
                buckets[i][4] = copy.copy(lastBucketMerge) 
                  
            buckets[ind][4] = copy.copy(lastBucketMerge)   
    
    
    #step3
    eps = 0.01 #error with the variance computation
    k = float(9.0/eps)
    if (len(buckets)>2):
        for i in xrange(3,(ind+1)):
            Vij = float(float(buckets[i][3]) + float(buckets[i-1][3]) + (((float(buckets[i][1]) *float(buckets[i-1][1]))/(float(buckets[i][1]) + float(buckets[i-1][1])))  *  float(math.pow(((buckets[i][2] - buckets[i-1][2])),2))) )
            ls = k*Vij
            if (ls<=buckets[i-1][4][3]):
                # merge 2 buckets
                bm = merge2VarianceBucket(buckets[i],buckets[i-1])
                buckets[i] = bm
                buckets[i-1] = "$"
                buckets.remove("$")
                
                break
            
            
            
        
    # output for the player
    avg = 0
    if bFullWindow and buckets[len(buckets)-1][1] > 1:
        # it takes in consideration the not yet expired point in the last bucket 

        nNotExp = (T_window * ncomponent)-len(buckets)
        avNotExp = buckets[len(buckets)-1][2]
        vNotExp = float(buckets[len(buckets)-1][3])/2.0
        Var = vNotExp + float(buckets[len(buckets)-1][4][3]) + (float((nNotExp*buckets[len(buckets)-1][4][1])/(nNotExp+buckets[len(buckets)-1][4][1])) * float(math.pow((avNotExp-float(buckets[len(buckets)-1][4][2])),2)))
        avg = buckets[len(buckets)-1][4][2]
        
        
        
    else:
        if len(buckets) > 1:
            Var = buckets[len(buckets)-1][4][3]
            avg = buckets[len(buckets)-1][4][2]
        else:    
            Var = buckets[len(buckets)-1][3]
            avg = buckets[len(buckets)-1][2]


    if (ncomponent > 1):
        actualPerformance = "-"
    else:            
    # keep track of the actual performance
        kmHSpeed = speedVX * 3.6     
        actualPerformance = ""
        
        if kmHSpeed < 1:
            actualPerformance = "stop"
        elif kmHSpeed <= 11:
            actualPerformance = "trot"
        elif kmHSpeed <= 14:
            actualPerformance = "low"
        elif kmHSpeed <= 17:
            actualPerformance = "medium"
        elif kmHSpeed <= 24:
            actualPerformance = "high"    
        elif kmHSpeed > 24:
            actualPerformance = "sprint"
        
        
    if len(buckets[len(buckets)-1]) == 5:
        buckets[len(buckets)-1].append(avg)
        buckets[len(buckets)-1].append(Var)
        buckets[len(buckets)-1].append(actualPerformance)

    else:
        buckets[len(buckets)-1][5] = avg
        buckets[len(buckets)-1][6] = Var
        buckets[len(buckets)-1][7] = actualPerformance
        
    varianceSpeedPlayers[playerKey] = buckets  
   

    lLastTimeStampSpeed = timestamp
    
    #print("Thanks to Motwani & Co.")


def loadInterruptions():
    

    firstLineRead = False
    tempTime = -1
    shift = 0
    with open('STOP_EVENTS/interruptions.csv', 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=';')
        for row in reader:
            row = string.split(row[0],'\t')
            if firstLineRead:
                
                if row[0].startswith('2'):
                    shift = MATCH_START_TIME
                else:
                    shift = SECOND_HALF_BEGIN
                    
                #case to handle zero
                if row[2] == '0':
                    int_struct_time = 0
                else:    
                    struct_time = datetime.datetime.strptime(row[2], "%H:%M:%S.%f")
                    struct_time = float(((struct_time.minute*60)+(struct_time.hour*(3600))+struct_time.second)*math.pow(10, 12))+((struct_time.microsecond)*math.pow(10, 6))
                    int_struct_time = int(struct_time)
                    
                if tempTime != -1:
                    gameInterruptions.append((tempTime, (shift+int_struct_time)))
                    tempTime = -1
                else:
                    tempTime = shift + int_struct_time
            else:
                firstLineRead = True


def format_seconds_to_mmss(seconds):
    minutes = seconds // 60
    seconds %= 60
    return "%02i:%02i" % (minutes, seconds)

def plotClusters(filename, wTraj, lPlayers,listTrajectories, bPlayerWithTheBall, filenameBall):
    
    figPl = plt.figure()
    axPl = plt.subplot(111)
    
    
    
    if (bPlayerWithTheBall):
    
        f = open(filenameBall)
        i=0
        numTrajRep = 0
        
        for line in iter(f):
            
            pointX = True
            listX = []
            listY = []
        
        
            lstdIn = line.replace('\n','')
            listPoint = string.split(lstdIn," ")
           
            if(i==0):
                numTrajRep = int(lstdIn)
                if int(numTrajRep) != 0:
                    
                    #field coordinates  
                    fieldX = [-52.447, 52.477, 52.477,52.477,52.477,-52.447,-52.447,-52.447,0.,0.]
                    fieldY = [33.960,33.960,33.960,-33.960,-33.960,-33.960,-33.960,33.960,33.960,-33.960]   
                    # plot the field   
                    axPl.plot(fieldX,fieldY)
                    plotAreaAndMidfield(axPl,figPl)
                    i+=1
                    
                
            else:
                numOfPoint = listPoint[1]
                numOfSubTrajInCluster = listPoint[2]
               
                if (int(numOfSubTrajInCluster) > 1):
                
                    for p in range(3,(int(numOfPoint)*2)+3):
                        if (pointX):
                            listX.append(listPoint[p])
                            pointX = False
                        else:
                            listY.append(listPoint[p])
                            pointX = True
                    
                    colRand = 'black'   
                    lines = axPl.plot(listX,listY,label = "Ball: "+str(numOfSubTrajInCluster)+': number of subtrajectories')
                    
                    p1 = float((listX[0]))
                    p2 = float((listY[0]))
                    
                    circle1=plt.Circle((p1,p2),1,color=colRand,label="start")
                    
                    p1 = float(listX[len(listX)-1])
                    p2 = float(listY[len(listY)-1])
                    
                    square = plt.Rectangle((p1,p2), 1 , 1 ,color=colRand,label="end")
                    
                    figPl.gca().add_artist(circle1)
                    figPl.gca().add_artist(square)
                    plt.setp(lines, color=colRand, linewidth=2.0) 
                    
    #### END PRINTING THE BALL IF NEEDED  
    
    f = open(filename)
    i=0
    numTrajRep = 0
    
    bShow = False
    
    for line in iter(f):
    
        pointX = True
        listX = []
        listY = []
        
        
        lstdIn = line.replace('\n','')
        listPoint = string.split(lstdIn," ")
        
        if(i==0):
            
            numTrajRep = int(lstdIn)
            
            if int(numTrajRep) != 0:
                
                if (not bPlayerWithTheBall):
                    #field coordinates  
                    
                    fieldX = [-52.447, 52.477, 52.477,52.477,52.477,-52.447,-52.447,-52.447,0.,0.]
                    fieldY = [33.960,33.960,33.960,-33.960,-33.960,-33.960,-33.960,33.960,33.960,-33.960]   
                    # plot the field   
                    axPl.plot(fieldX,fieldY)
                    plotAreaAndMidfield(axPl,figPl)
                    
                    printWholetrajectorySpec(listTrajectories, axPl,plotAreaAndMidfield)
                    
                i+=1
        
        else:
        
            bShow = True
            #idTraj = listPoint[0]
            numOfPoint = listPoint[1]
            numOfSubTrajInCluster = listPoint[2]
            
            if (int(numOfSubTrajInCluster) > 1):
                for p in range(3,(int(numOfPoint)*2)+3):
                    if (pointX):
                        listX.append(listPoint[p])
                        pointX = False
                    else:
                        listY.append(listPoint[p])
                        pointX = True
                
                colRand = numpy.random.rand(3,1)     
                lines = axPl.plot(listX,listY,label = str(numOfSubTrajInCluster)+': number of subtrajectories')
                
                p1 = float((listX[0]))
                p2 = float((listY[0]))
                
                circle1=plt.Circle((p1,p2),1,color=colRand,label="start")
                
                p1 = float(listX[len(listX)-1])
                p2 = float(listY[len(listY)-1])
                
                square = plt.Rectangle((p1,p2), 1 , 1 ,color=colRand,label="end")
                
                figPl.gca().add_artist(circle1)
                figPl.gca().add_artist(square)
                plt.setp(lines, color=colRand, linewidth=5.0)    

    if bShow:
        slPlayers = ""
        
        time = 0
        
        for key in lPlayers:
            time = max(time,int(lastTimeStampStoreTraj[key]))
            slPlayers = slPlayers +" "+ key
           
        time = (time-effectiveStart) * math.pow(10,-12)  
        strTime = format_seconds_to_mmss(time)  
        plt.suptitle("Trajectories for players: "+slPlayers+" at time (MM:SS): "+strTime)
        
    
        plt.legend( loc='upper left', numpoints = 1 )
        plt.show()
            
       
def findTeam(key):

#Dennis Dotterweich (Left Leg: 47, Right Leg:16)
#Niklas Waelzlein (Left Leg: 49, Right Leg: 88)
#Wili Sommer (Left Leg: 19, Right Leg: 52)
#Philipp Harlass (Left Leg: 53, Right Leg: 54)
#Roman Hartleb (Left Leg: 23, Right Leg: 24)
#Erik Engelhardt (Left Leg: 57, Right Leg: 58)
#Sandro Schneider (Left Leg: 59, Right Leg: 28
    if (key == "1647" or key == "4988"):
        return "Team A Defenders",2
    if (key == "1952" or key == "5354" or key == "2324"):
        return "Team A Midfielders",3
    if (key == "5758" or key == "2859"):
        return "Team A Strikers",2
    
    
#Kevin Baer (Left Leg: 63, Right Leg: 64)
#Luca Ziegler (Left Leg: 65, Right Leg: 66)
#Ben Mueller (Left Leg: 67, Right Leg: 68)
#Vale Reitstetter (Left Leg: 69, Right Leg: 38)
#Christopher Lee (Left Leg: 71, Right Leg: 40)
#Leon Heinze (Left Leg: 73, Right Leg: 74)
#Leo Langhans (Left Leg: 75, Right Leg: 44)
    if (key == "6364" or key == "6566"):
        return "Team B Defenders",2
    if (key == "6768" or key == "3869" or key == "4071"):
        return "Team B Midfielders",3
    if (key == "7374" or key == "4475"):
        return "Team B Strikers",2
    
 
    
    
    
def returnPairValue(sid):
    
    
#Dennis Dotterweich (Left Leg: 47, Right Leg:16)
#Niklas Waelzlein (Left Leg: 49, Right Leg: 88)
#Wili Sommer (Left Leg: 19, Right Leg: 52)
#Philipp Harlass (Left Leg: 53, Right Leg: 54)
#Roman Hartleb (Left Leg: 23, Right Leg: 24)
#Erik Engelhardt (Left Leg: 57, Right Leg: 58)
#Sandro Schneider (Left Leg: 59, Right Leg: 28

    
        
    if (sid == "47"):
        return "16"
    if (sid == "49"):
        return "88"
    if (sid == "19"):
        return "52"
    if (sid == "53"):
        return "54"
    if (sid == "23"):
        return "24"
    if (sid == "57"):
        return "58"
    if (sid == "59"):
        return "28"
    if (sid == "16"):
        return "47"
    if (sid == "88"):
        return "49"
    if (sid == "52"):
        return "19"
    if (sid == "54"):
        return "53"
    if (sid == "24"):
        return "23"
    if (sid == "58"):
        return "57"
    if (sid == "28"):
        return "59"
    

#Kevin Baer (Left Leg: 63, Right Leg: 64)
#Luca Ziegler (Left Leg: 65, Right Leg: 66)
#Ben Mueller (Left Leg: 67, Right Leg: 68)
#Vale Reitstetter (Left Leg: 69, Right Leg: 38)
#Christopher Lee (Left Leg: 71, Right Leg: 40)
#Leon Heinze (Left Leg: 73, Right Leg: 74)
#Leo Langhans (Left Leg: 75, Right Leg: 44)
    

    if (sid == "63"):
        return "64"
    if (sid == "65"):
        return "66"
    if (sid == "67"):
        return "68"
    if (sid == "69"):
        return "38"
    if (sid == "71"):
        return "40"
    if (sid == "73"):
        return "74"
    if (sid == "75"):
        return "44"
    if (sid == "64"):
        return "63"
    if (sid == "66"):
        return "65"
    if (sid == "68"):
        return "67"
    if (sid == "38"):
        return "69"
    if (sid == "40"):
        return "71"
    if (sid == "74"):
        return "73"
    if (sid == "44"):
        return "75"

def clusterBallPlayer(listTr):
    
    global keyTraj
    #CLUSTER OF THE TRAJECTORY FOR THE BALL
    lPlayersKeys =[]
    lT = []
    for kt in keyTraj:
        kt = kt.replace('\n','')
        if (kt == "ball"):
            if (listTr.has_key(kt)):
                lT.append(listTr[kt])
                lPlayersKeys.append(kt)
            else:
                print "the data about the ball don't exist"    
    
    if(len(lT)> 0):
        
        filename = str(uuid.uuid4())
        filename = "DATA_trajectories_ball/"+filename + '.tra'
        fnew = open(filename, 'w')
        fnew.write('2\n') # I write the number of dimension
        fnew.write(str(len(lT))+'\n') # I write the number of trajectory
        
        cont = 0
        for traj in lT:
            
            fnew.write(str(cont)+' ')
            fnew.write(str(len(traj))+' ')
            
            for point in traj:
                fnew.write(str(point[0])+' ')
                fnew.write(str(point[1])+' ')
                #fnew.write(str(point[2])+' ')
                
            fnew.write('\n')    
            cont+=1
        fnew.close()
        return filename
    
    return None
        

def printWholetrajectorySpec(listTrajectories, axSpec,figPlotSpec):
    
    
   
    for traj in listTrajectories:
        
        listX = [x[0] for x in traj]
        listY = [y[1] for y in traj]
        
        colRand = numpy.random.rand(3,1) 
        lines = axSpec.plot(listX,listY)
        plt.setp(lines, color=colRand , linewidth=0.5) 

        

    

def printWholetrajectory(listTrajectories):
    
    
    # ax = plt.subplot(111)
    figPlot = plt.figure()
    axPlot = figPlot.add_subplot(111,aspect='equal') 
    
    #field coordinatesCONST_FIELD_YLENGTH
        
    fieldX = [-52.447, 52.477, 52.477,52.477,52.477,-52.447,-52.447,-52.447,0.,0.]
    fieldY = [33.960,33.960,33.960,-33.960,-33.960,-33.960,-33.960,33.960,33.960,-33.960]
        
    # plot the field
        
    axPlot.plot(fieldX,fieldY)
    plotAreaAndMidfield(axPlot,figPlot)
    for traj in listTrajectories:
        
        listX = [x[0] for x in traj]
        listY = [y[1] for y in traj]
        
        colRand = numpy.random.rand(3,1) 
        circle1=plt.Circle((traj[0][0],traj[0][1]),3,color='b')
        square = plt.Rectangle((traj[len(traj)-1][0],traj[len(traj)-1][1]), 3 , 3 ,color='b')
        
        lines = axPlot.plot(listX,listY)
        
        
        
        plt.setp(lines, color=colRand , linewidth=2.0) 
        
        figPlot.gca().add_artist(circle1)
        figPlot.gca().add_artist(square)
        
    plt.show()


def startScanning(inputPath):
        
       
        global keyTraj
        global listOfLine
        global ax
        global manager
        global timer 
        global effectiveStart
        global TimeWindowSize
        global lastTimeStampStoreTraj
        global bMatchEnd
        global epsilonTraclus
        loadInterruptions()
        f = open(inputPath)
        
        lavgpoint = {}
        listTr = {}
        lFrequency = {}
        nTrStored = 0
        keyTraj = []
        bStartStoringData = False
        bMatchStarted = False
        bFirstIteration = True
        bAskClusterSpeed = False
        bClusterTrajectory = False
        sid = None
        cX = None
        cY = None
        cZ = None
        v = None
        lKeySpeedLogStarting = {}

        for line in iter(f):
              

            # CHOSE IF STORE SPEED PERFORMANCE OR TRAJECTORY  
            bAskTimeSize = False
            if (not bAskClusterSpeed):
                        print "Type 1 if you want cluster the trajectory\n2 if you want to check the speed performance"
                        while (not bAskClusterSpeed):
                            lfirstIn = sys.stdin.readline()
                            if (lfirstIn.startswith('1') ):
                                bClusterTrajectory = True
                                bAskClusterSpeed = True
                                
                                print "Insert epsilon (float) parameter for trajectory clustering: " 
                                while (not bAskTimeSize): 
                                    lfirstIn = sys.stdin.readline()
                                    
                                    try:
                                        epsilonTraclus =  float(lfirstIn)
                                        bAskTimeSize = True
                                    except exceptions.ValueError:
                                        print "Insert correct value: " 
                                
                                
                            else:
                                if (lfirstIn.startswith('2') ):
                                    bClusterTrajectory = False
                                    bAskClusterSpeed = True
                                    print "Insert time size window in minutes: " 
                                    while (not bAskTimeSize): 
                                        lfirstIn = sys.stdin.readline()
                                        
                                        try:
                                            TimeWindowSize =  float(lfirstIn)
                                            # it converts time window size in seconds
                                            TimeWindowSize = TimeWindowSize * 60
                                            bAskTimeSize = True
                                        except exceptions.ValueError:
                                            print "Insert correct value: " 

                                else:
                                    if (lfirstIn != ""):
                                        print "please insert a correct value 1: Clustering, 2: Speed performance!"
              
              
            # READ THE LINE FROM THE SENSORS
            tokens = string.split(line,',')
            sid = tokens[0]
            timestamp = int(tokens[1])
            ##########################
            #Check if the time-stamp is between the start and the end of the match 
            if (timestamp < MATCH_START_TIME ):
                
                if bFirstIteration:
                    print "The match will start in a few..."
                    bFirstIteration = False
               
                continue
            
            if (timestamp > MATCH_END_TIME):
                bMatchEnd = True
                  
              
            if (bStartStoringData):
                while sys.stdin in select.select([sys.stdin], [], [], 0)[0] or bMatchEnd:
                    
                    lstdIn = sys.stdin.readline()
                    
                    if (lstdIn.startswith('c') and bClusterTrajectory and bMatchStarted):
                        
                        #CLUSTER OF THE TRAJECTORY
                        lPlayersKeys =[]
                        lT = []
                        for kt in keyTraj:
                            kt = kt.replace('\n','')
                            if (kt != "ball"):
                                if (listTr.has_key(kt)):
                                    lT.append(listTr[kt])
                                    lPlayersKeys.append(kt)
                                else:
                                    print "the player with key "+kt+" doesn't exist"    

                        if(len(lT)> 0):
                            
                            filename = str(uuid.uuid4())
                            filename = "DATA_trajectories/"+filename + '.tra'
                            fnew = open(filename, 'w')
                            fnew.write('2\n') # I write the number of dimension
                            fnew.write(str(len(lT))+'\n') # I write the number of trajectory
                            
                            cont = 0
                            for traj in lT:
                                
                                fnew.write(str(cont)+' ')
                                fnew.write(str(len(traj))+' ')
                                
                                for point in traj:
                                    fnew.write(str(point[0])+' ')
                                    fnew.write(str(point[1])+' ')
                                    
                                fnew.write('\n')    
                                cont+=1
                            fnew.close()
                            
                           
                            filenameOut = filename + '.out'
                                
                            args = ("movebank/bin/traclus", filename, filenameOut, str(epsilonTraclus), str(len(lT)))                       
                            popen = subprocess.Popen(args, stdout=subprocess.PIPE)
                            popen.wait()
                            output = popen.stdout.read()
                            print output
                            print '\n'
                           
                            #os.remove(filename)
                            
                            if (lstdIn.startswith('cwithball')):
                                fn = clusterBallPlayer(listTr)
                                if (fn != None):
                                    filenameOutBall = fn + '.out'
                                    args = ("movebank/bin/traclus", fn, filenameOutBall, str(epsilonTraclus), str(1))                       
                                    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
                                    popen.wait()
                                    output = popen.stdout.read()
                                    print output
                                    print '\n'
                                    plotClusters(filenameOut,listTr,lPlayersKeys,lT,True,filenameOutBall)
                                    
                                else:
                                    plotClusters(filenameOut,listTr,lPlayersKeys,lT,False,None)
                                
                            else:
                                plotClusters(filenameOut,listTr,lPlayersKeys,lT,False,None)
                            
                            #os.remove(filenameOut)
                            #TRACLUS(lT)
                        else:
                            print "You have to insert at least one player:" 
                    
                    elif (lstdIn.startswith('bc') and bClusterTrajectory and bMatchStarted):
                        
                        #CLUSTER OF THE TRAJECTORY FOR THE BALL
                        lPlayersKeys =[]
                        lT = []
                        for kt in keyTraj:
                            kt = kt.replace('\n','')
                            if (kt == "ball"):
                                if (listTr.has_key(kt)):
                                    lT.append(listTr[kt])
                                    lPlayersKeys.append(kt)
                                else:
                                    print "the data about the ball don't exist"    

                        if(len(lT)> 0):
                            
                            filename = str(uuid.uuid4())
                            filename = "DATA_trajectories_ball/"+filename + '.tra'
                            fnew = open(filename, 'w')
                            fnew.write('2\n') # I write the number of dimension
                            fnew.write(str(len(lT))+'\n') # I write the number of trajectory
                            
                            cont = 0
                            for traj in lT:
                                
                                fnew.write(str(cont)+' ')
                                fnew.write(str(len(traj))+' ')
                                
                                for point in traj:
                                    fnew.write(str(point[0])+' ')
                                    fnew.write(str(point[1])+' ')
                                    #fnew.write(str(point[2])+' ')
                                    
                                fnew.write('\n')    
                                cont+=1
                            fnew.close()
                            
                           
                            filenameOut = filename + '.out'
                                
                            args = ("movebank/bin/traclus", filename, filenameOut, str(epsilonTraclus), str(len(lT)))                       
                            popen = subprocess.Popen(args, stdout=subprocess.PIPE)
                            popen.wait()
                            output = popen.stdout.read()
                            print output
                            print '\n'
                           
                            #os.remove(filename)
                            plotClusters(filenameOut,listTr,lPlayersKeys,lT,False,None)
                            
                            #os.remove(filenameOut)
                            #TRACLUS(lT)
                        else:
                            print "You have to insert at least one player:" 
                    
                    else: 
                        
                        #print the whole trajectory of a player
                        if(lstdIn.startswith('w') and bClusterTrajectory and bMatchStarted ):
                            
                            lT = []
                            for kt in keyTraj:
                                kt = kt.replace('\n','')
                                if (listTr.has_key(kt)):
                                    lT.append(listTr[kt])
                                else:
                                    print "the player with key "+kt+" doesn't exist"    
                            if(len(lT)> 0):
                                printWholetrajectory(lT)                        
                                #print the whole trajectory of a player
                            else:
                                print "You have to insert at least one player:"
                                
                        if(lstdIn.startswith('e') and bMatchStarted and not bPlotSpeed ):
                            
                            keyTraj = []   
                            print "All the players enqueued are deleted"
                            
                        if(lstdIn.startswith('p') and not bClusterTrajectory and bMatchStarted and not bPlotSpeed):
                            #show the plot with the speed
                            
                            thread.start_new_thread( startPlotting, () )
                                
                           
                            
                        if(lstdIn.startswith('bc') == False and lstdIn.startswith('c') == False and lstdIn.startswith('w') == False and lstdIn.startswith('e') == False and lstdIn.startswith('p') == False and bMatchStarted and not bPlotSpeed):

                            try:
                                keyTraj.index(lstdIn)
                                print "value already enqueued"
                            except ValueError:
                               
                                keyTraj.append(lstdIn)
                                if bClusterTrajectory:
                                    print "Type [bc] if you want to cluster the trajectories of the ball \nType [c] if you want to cluster the trajectories of the players \nType [w] if you want to plot the whole trajectories of the player\nType [e] if you want erase the player inserted otherwise put the player id to enqueue it\n"
                                else:
                                    print "Type [e] if you want erase the player inserted otherwise put the player id to enqueue it\n [p] if you want to plot the speed performance"
  
                                    ### END THE PART CONCERNING THE CHOICE FROM USER
            

            if (timestamp >= FIRST_HALF_END and timestamp <= SECOND_HALF_BEGIN ):
                
                if effectiveStart == MATCH_START_TIME:
                    print str(timestamp)
                    print "break beetwen first and second half"
                    effectiveStart = FIRST_HALF_END
                continue
            else:
                if (effectiveStart == FIRST_HALF_END):
                    effectiveStart = SECOND_HALF_BEGIN
                    
            ##########################
        

            # IT TAKES CARE OF INTERRUPTIONS
            for endevent in gameInterruptions:
                if (timestamp >=endevent[0] and timestamp<=endevent[1]):
                    #during an interruption
                    continue
                if (timestamp > endevent[1]):
                    gameInterruptions.pop()
                    break
            ###########################

                
            ### CHECK THE INTERSECTION OF THE PLAYER WITH THE FIELD
            
            ## check if the player's foot is in the field's rectangle or ball
            
            
            cX =  float(tokens[2])/1000.0
            cY =  float(tokens[3])/1000.0 

            cXg =  float(tokens[2])
            cYg =  float(tokens[3])
              
            if cXg < CONST_FIELD_MINX or cXg > CONST_FIELD_MAXX or cYg < CONST_FIELD_MINY or cYg > CONST_FIELD_MAXY:
                continue    
           # p = [cX,cY]
            
            #if (p[0]<-0.005 or p[0]> 52.489 or p[1]>33.965 or p[1]<-33.939):
                #print "point out of field"
            #    continue
            #else:
                #print "point into the field"
            
            #######################################################    


            if not bMatchStarted :
                bMatchStarted = True
                print("The referee put the match in play")
                effectiveStart = MATCH_START_TIME

            if not bClusterTrajectory:
                #uncomment to create a toydataset
                global fileDataset
                
                if fileDataset == None:
                    fileDataset = open("Dataset_firstHalf_2Minutes","w")
                
                fileDataset.write(line)
                
                

                storePoint = True
                #Every 2 seconds 200 Hz * 2 I store the point of the trajectory
                if (lFrequency.has_key(sid) == False):
                    lFrequency[sid] = 1
                else:
                    lFrequency[sid] += 1
                    if (lFrequency[sid] != 200):
                        storePoint = False
                    else:
                        lFrequency[sid] = 0
                        
                if (storePoint == False):
                        continue

                
                #COMPUTATION OF THE SPEED PERFORMANCE
                
                #vx =  int(tokens[7]) 
                # skip the ball and goalkeeper
                if (sid != "12" and 
                sid != "10" and 
                sid != "8" and 
                sid != "4" and
                sid != "14" and 
                sid != "13" and 
                sid != "97" and 
                sid != "98" and 
                sid != "61" and 
                sid != "62" and 
                sid != "99" and 
                sid != "100" and
                sid != "105" and
                sid != "106"):
                
                    v =  int(tokens[5])
                    pairValue = returnPairValue(sid)
                    if pairValue != None:
                        
                        lLastSpeed[sid] = v
                        key = str(min(int(pairValue),int(sid))) + str(max(int(pairValue),int(sid)))
                        if ( lLastSpeed.has_key(sid) and lLastSpeed.has_key(pairValue)):
    
                            
                                 
                            if (len(lKeySpeedLogStarting) == 14 and bStartStoringData == False):
                                bStartStoringData = True
                                print "Insert the number of player you want to check the speed performance: "
                            else:
                                if (not lKeySpeedLogStarting.has_key(key)):
                                    lKeySpeedLogStarting[key] = True
                                
       
                            v = (lLastSpeed[sid] + lLastSpeed[pairValue])/2
                            lLastSpeed.pop(sid)
                            lLastSpeed.pop(pairValue)
                           
                            speedVX =  v * (math.pow(10,-6)) #* vx #vx in m/s is derived by |v| * 1e-10 * vx
                            speedVX = math.fabs(speedVX)
                            loadSpeedPerformance(key,speedVX,timestamp,TimeWindowSize,1)
                            
                            TeamOfPlayer, ncomponent = findTeam(key)
                            loadSpeedPerformance(TeamOfPlayer,speedVX,timestamp,TimeWindowSize,ncomponent)
                ################################################################################
           
                ##END OF THE PART CONCERNING THE SPEED PERFORMANCE
           
            else: ## START OF THE PART CONCERNING THE TRAJECTORY CLUSTERING
                
                storePoint = True
                #Every 2 seconds 200 Hz * 2 I store the point of the trajectory
                if (lFrequency.has_key(sid) == False):
                    lFrequency[sid] = 1
                else:
                    lFrequency[sid] += 1
                    if (lFrequency[sid] != 200):
                        storePoint = False
                    else:
                        lFrequency[sid] = 0
                        
                if (storePoint == False):
                        continue
    
                
                # COMPUTATION OF THE TRAJECTORIES ###  
                # skip the goalkeepers and the referee
                if (sid != "14" and 
                    sid != "13" and 
                    sid != "97" and 
                    sid != "98" and 
                    sid != "61" and 
                    sid != "62" and 
                    sid != "99" and 
                    sid != "100" and
                    sid != "105" and
                    sid != "106" ):
                    
                  
                    
                    cX =  tokens[2]
                    cY =  tokens[3]
                    cZ =  tokens[4]
         
         
                    if (sid != "12" and sid != "10" and sid != "8" and sid != "4"):
                        
                        if (lavgpoint.has_key(sid) == False):
                            lavgpoint[sid] = [[float(cX)/1000.0,float(cY)/1000.0,float(cZ)/1000.0]]
                        else:
                            lavgpoint[sid].insert(0,[float(cX)/1000.0,float(cY)/1000.0,float(cZ)/1000.0])
                            
    
                        pairValue = returnPairValue(sid)
                        
                        if (lavgpoint.has_key(pairValue)):
                            if(len(lavgpoint[pairValue])>0):
                                
                                pointTraj = (((lavgpoint[pairValue][len(lavgpoint[pairValue])-1][0]+lavgpoint[sid][len(lavgpoint[sid])-1][0])/2.0),
                                             ((lavgpoint[pairValue][len(lavgpoint[pairValue])-1][1]+lavgpoint[sid][len(lavgpoint[sid])-1][1])/2.0),
                                             ((lavgpoint[pairValue][len(lavgpoint[pairValue])-1][2]+lavgpoint[sid][len(lavgpoint[sid])-1][2])/2.0))
                                
                                lavgpoint[sid].pop(len(lavgpoint[sid])-1)
                                lavgpoint[pairValue].pop(len(lavgpoint[pairValue])-1)
                                
                                key = str(min(int(pairValue),int(sid))) + str(max(int(pairValue),int(sid)))
                                
                                
                                if (listTr.has_key(key) == False):
                                    nTrStored = nTrStored +1
                                    print "Start to store trajectory for the player "+key
                                    listTr[key] = []
                                    listTr[key].append(pointTraj)
                                else:
                                    listTr[key].append(pointTraj)
                                    
                                 
                                lastTimeStampStoreTraj[key] = timestamp
                                    #print "trajectory for the player "+key+" has a length of "+str(len(listTr[key]))
                                
                                if (nTrStored == 15 and bStartStoringData == False):
                                    bStartStoringData = True
                                    print "Insert the number of player you want to cluster or insert [ball] if you want cluster the trajectory of the ball : "
                    else:
                        
                        ###HERE IT STORES THE BALL Trajectory
                        pointTraj = (float(cX)/1000.0,float(cY)/1000.0,float(cZ)/1000.0)
                        key = "ball"
                        if (listTr.has_key(key) == False):
                            nTrStored = nTrStored +1
                            print "Start to store trajectory for the ball "
                            listTr[key] = []
                            listTr[key].append(pointTraj)
                        else:
                            listTr[key].append(pointTraj)
                            
                        lastTimeStampStoreTraj[key] = timestamp          
                    ######################################################################
                        
    
                   

        # HERE IT FINISHED TO READ ALL THE DATA
        
        
        
        
        
        
        f.close()
        


def RealtimePloter(arg):
    
    global effectiveStart
    
    period = ""
    
    if (effectiveStart == FIRST_HALF_END):
        period = "interval"
    if (effectiveStart == MATCH_START_TIME ):
        period = "First half"
    if (effectiveStart == SECOND_HALF_BEGIN ):
        period = "Second half"            
    if bMatchEnd:
        period = "End match"
    
    actualTime = (lLastTimeStampSpeed - effectiveStart) * math.pow(10,-12)
    
    timeFormatted = format_seconds_to_mmss(actualTime) 
    
    
    global  ax, manager, keyTraj,contSpeed
    
    ax.set_title("Realtime Speed performance Plot - "+period+" Time: "+timeFormatted+" Time sliding window in seconds: "+str(TimeWindowSize))
   
    #print "real time plotting "+str(len(keyTraj))

    minV = sys.float_info.max
    maxV = sys.float_info.min
    
    
    minAxisX = None
    maxAxisX = None
    handles, labels = ax.get_legend_handles_labels()   
    cont = 0

    CurrentXAxis= pylab.arange(contSpeed,contSpeed+20,1)
    
    for kt in keyTraj:
        
            kt = kt.replace('\n','')
            var,actualPerformance = updateAvgVarStdDevPlayers(str(kt))
            
            minV = min(min(lAvgPlayers[str(kt)]),minV)
            maxV = max(max(lAvgPlayers[str(kt)]),maxV)
           
           
            
            #CurrentXAxis= pylab.arange(contSpeed,(contSpeed+len(lAvgPlayers[str(kt)])),1)

            
            #if (minAxisX == None):
            #    minAxisX = CurrentXAxis.min()
            #    maxAxisX = CurrentXAxis.max()
            #else:
            #    if minAxisX > CurrentXAxis.min():
            #        minAxisX = CurrentXAxis.min()
            #    if maxAxisX < CurrentXAxis.max():
            #        maxAxisX = CurrentXAxis.max()
             
            
            if (len(lAvgPlayers[str(kt)]) <20):
                nl = [0] * (20-len(lAvgPlayers[str(kt)]))
                nl = lAvgPlayers[str(kt)] + nl
            else:
                nl = lAvgPlayers[str(kt)]  
              
            
            
            if actualPerformance == "-":  
                labels[cont] = "player = "+kt+ " Variance = "+str(var)
                
            else:
                
                avgVel = 0   
                if (len(lAvgPlayers[str(kt)]) >0):   
                    avgVel = lAvgPlayers[str(kt)][(len(lAvgPlayers[str(kt)])-1)]
                
                
                if (lLastAvgTime.has_key(kt)):
                    avOld = lLastAvgTime[kt][0]
                    atOld = lLastAvgTime[kt][1]
                    
                    if (avOld != float(avgVel) and atOld != float(actualTime)):
                        
                        kmRan= (float(avgVel) * float(actualTime)) /1000.0
                        # fix the problem with measure error
                        if (kmRan < lLastAvgTime[kt][2]):
                            kmRan = lLastAvgTime[kt][2]
                        else:
                            lLastAvgTime[kt] = (float(avgVel),float(actualTime),kmRan)
                    else:
                        kmRan = lLastAvgTime[kt][2]
                else:
                    kmRan= (float(avgVel) * float(actualTime)) /1000.0
                    lLastAvgTime[kt] = (float(avgVel),float(actualTime),kmRan)
                
                labels[cont] = "player = "+kt+ " Variance = "+str(var)+" Actual performance = "+actualPerformance+ " Ran KM: "+str(kmRan)
           
            cont = cont +1   
            listOfLine[str(kt)][0].set_data(CurrentXAxis,pylab.array(nl))
           
            
    minV = minV - minV
    maxV = maxV + maxV
    
    ax.axis([CurrentXAxis.min(),CurrentXAxis.max(),minV,maxV])
    
    #ax.legend(handles, labels )
    pylab.legend(handles, labels,loc='upper left' )
   
            
    
    manager.canvas.draw()
    contSpeed = contSpeed + 1
    #manager.show()


def startPlotting():
    
    global ax ,listOfLine, manager, timer,lLastVarRegister,bPlotSpeed
    
    
    if len(keyTraj)>0:
        
        fig = pylab.figure(1)
        ax = fig.add_subplot(111)
        ax.grid(True)
        ax.set_title("Realtime Speed performance Plot - Waiting...")
        ax.set_xlabel("Time (real time, it could be differ from the time of match)")
        ax.set_ylabel("Total Speed Average m/s")
        ax.axis([0,1,-10,10])
        
        nShow = 0
        for kt in keyTraj:
            kt = kt.replace('\n','')
            if varianceSpeedPlayers.has_key(str(kt)):
                nShow = nShow + 1
                xAchse = pylab.arange(0,1,1)
                yAchse = pylab.arange(0,1,1)
                listOfLine[str(kt)] = ax.plot(xAchse,yAchse,'-',label = "player = "+kt) 
            else:
                print "The player "+ kt + " does not exist"
                
     
        if (nShow > 0 ):
            bPlotSpeed = True
            pylab.legend( loc='upper left' )
            manager = pylab.get_current_fig_manager()
            timer = fig.canvas.new_timer(interval=200)
            timer.add_callback(RealtimePloter, ())
            timer.start()
            pylab.show()
            
            lLastVarRegister = {}
            contSpeed = 0
            bPlotSpeed = False
        else:
            print "You have to insert at least one player:"
   
     
MATCH_START_TIME = 10753295594424116
FIRST_HALF_END =  12398000000000000              
SECOND_HALF_BEGIN =  13086639146403495
MATCH_END_TIME =  14879639146403495
MATCH_END_TIME = 14879639146403495

CONST_FIELD_MINX = 0
CONST_FIELD_MAXX = 52483
CONST_FIELD_MINY = -33960
CONST_FIELD_MAXY = 33960
CONST_FIELD_XLENGTH = 52483/2
CONST_FIELD_YLENGTH = 33960

#thread.start_new_thread( startScanning, ("DEBSDATA/full-game",) )                 
startScanning('DEBSDATA/full-game')



