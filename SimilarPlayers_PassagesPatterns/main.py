'''
Created on 07/mag/2013

@author: denadai2
'''
import csv
import sys
#Per vedere il possesso palla
import numpy
import scipy.spatial
import math
from pprint import pprint
from collections import namedtuple
from Game import Game

import dataset

db = dataset.connect('sqlite:///prova.sqlite')

CONST_LOGEVENTS = True
CONST_FIRSTMATCHSTART =     10753295594424116
#match end with the ball's error fix
CONST_FIRSTMATCHEND =       12398000000000000
CONST_SECONDMATCHSTART =    13086639146403495
CONST_SECONDMATCHEND =      14879639146403495
#distance from the ball to the player (mm)
CONST_DISTANCEGRANULARITY = 1000
#FIELD
CONST_FIELD_MINX = 0
CONST_FIELD_MAXX = 52483
CONST_FIELD_MINY = -33960
CONST_FIELD_MAXY = 33960
CONST_FIELD_XLENGTH = 52483/2
CONST_FIELD_YLENGTH = 33960
#GOAL
CONST_GOAL_MINX = 22578
CONST_GOAL_MAXX = 29898
CONST_GOAL_MAXZ = 2440

#METRICS WEIGHT
CONST_METRICWEIGHT_TACKLES_SF = 1.50
CONST_METRICWEIGHT_TOTALSHOTS = 0.50
CONST_METRICWEIGHT_SHOTSONGOAL = 0.75
CONST_METRICWEIGHT_CROSSES = 0.60
CONST_METRICWEIGHT_TACKLES = 0.40
CONST_METRICWEIGHT_PASSAGESDT = 0.35
CONST_METRICWEIGHT_PASSAGES = 0.15

#match sensor=>playerID
playersSensors = {}
playersSensors2 = {'Nick': (13, 14,97,98),
                  'Dennis': (47,16),
                  'Niklas': (49,88),
                  'Wili': (19,52),
                  'Philipp': (53,54),
                  'Roman': (23,24),
                  'Erik': (57,58),
                  'Sandro': (59,28),
                  
                  'Leon': (61,62,99,100),
                  'Kevin': (63,64),
                  'Luca': (65,66),
                  'Ben': (67,68),
                  'Vale': (69,38),
                  'Christopher': (71,40),
                  'Leon2': (73,74),
                  'Leo': (75,44),
                  
                  'Referee': (105,106)
                  }
teamPlayers = {   'Nick': 0,
                  'Dennis': 0,
                  'Niklas': 0,
                  'Wili': 0,
                  'Philipp': 0,
                  'Roman': 0,
                  'Erik': 0,
                  'Sandro': 0,
                  
                  'Leon': 1,
                  'Kevin': 1,
                  'Luca': 1,
                  'Ben': 1,
                  'Vale': 1,
                  'Christopher': 1,
                  'Leon2': 1,
                  'Leo': 1,
                  
                  'Referee': 2
                  }


#actions (duration, numberPassages, playersInvolved)
class _Action:
    duration = 0
    numberPassages = 0
    playersInvolved = []

    def __init__(self, duration, numberPassages, playersInvolved):
        self.duration = duration
        self.numberPassages = numberPassages
        self.playersInvolved = playersInvolved


actions = [_Action(0, 0, [])]
actionPointer = 0

    
class MatchTrackRow:
    ID = 0
    actualTime = 0
    position = ()
    #acceleration m/s^2
    acceleration = 0
    #velocity m/s
    velocity = 0

    def __init__(self, ID, actualTime, x, y, z, velocity, acceleration):
        self.ID = int(ID)
        self.actualTime = int(actualTime) - CONST_FIRSTMATCHSTART
        #position with corrected z (negative z represent errors)
        self.position = (float(x), float(y), max(float(z), 0))
        self.acceleration = float(acceleration) * math.pow(10, -6)
        self.velocity = float(velocity) * math.pow(10, -6)

    def isTheBall(self):
        if self.ID == 4 or self.ID == 8 or self.ID == 10 or self.ID == 12:
            return True
        else:
            return False


def initialize(game):

    #this will build the corrispondence sensorID => playerID
    for pKey in playersSensors2.keys():
        for iKey in playersSensors2[pKey]:
            playersSensors[iKey] = pKey

    #flush everything from the table
    db['dbscan'].delete()
    db['actions'].delete()
    db['actions_players'].delete()

    game.configure(teamPlayers)

if __name__ == '__main__':
    with open('full-game', 'rb') as csvfile:
        fileReader = csv.reader(csvfile)
        
        game = Game()
        #count the time elapsed from the last distance check
        lastDistanceCheckTimestamp = 0
        lastPlayer = ""
        lastHitPlayerName = ""
        #the two players with the minimum distance from the ball
        _PlayerNearTheBall = namedtuple("_PlayerNearTheBall", "ID distance timestamp")
        playersNearTheBall = [_PlayerNearTheBall("", sys.maxint, 0), _PlayerNearTheBall("", sys.maxint, 0)]

        '''SETUP'''
        initialize(game)
            
        #debug
        '''plt.rcParams['axes.unicode_minus'] = False
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.ylim([-34000, 34000])
        plt.xlim([0, 53000])'''
        
        '''READ THE DATASET'''
        #I read the whole dataset and every n seconds I search the player near the ball
        for row in fileReader:
            #read the row
            trackRow = MatchTrackRow(row[0], row[1], row[2], row[3], row[4], row[5], row[6])

            #if the ball is out from the field
            if trackRow.actualTime < 0 or (trackRow.actualTime + CONST_FIRSTMATCHSTART > CONST_FIRSTMATCHEND  and trackRow.actualTime < CONST_SECONDMATCHSTART - CONST_FIRSTMATCHSTART):
                continue
            #exclude the referee
            if trackRow.ID == 105 or trackRow.ID == 106:
                continue
            #exclude sensors out of the field
            if trackRow.position[0] < CONST_FIELD_MINX or trackRow.position[0] > CONST_FIELD_MAXX: #or trackRow.position[1] < CONST_FIELD_MINY or trackRow.position[1] > CONST_FIELD_MAXY
                #the ball went out, so:
                if trackRow.isTheBall() and game.lastBall_ID == trackRow.ID:
                    game.ball_is_out = True
                continue

            #check the interruptions
            if len(game.interruptions) > 0:
                if game.interruptions[0][0] < trackRow.actualTime < game.interruptions[0][1]:
                    game.gameWasInterrupted = True
                    continue
                else:
                    #If I have passed the time
                    if trackRow.actualTime > game.interruptions[0][1]:
                        game.interruptions.pop(0)

            #init some variables
            if lastDistanceCheckTimestamp == 0:
                lastDistanceCheckTimestamp = trackRow.actualTime
            timePassed = trackRow.actualTime - lastDistanceCheckTimestamp

  
            '''BallHit'''
            if trackRow.isTheBall():
                #debug
                '''if trackRow.ID == 4:
                    ax.plot(trackRow.position[0], trackRow.position[1], 'x')'''

                #Ball hit
                if trackRow.acceleration >= 55:
                    game.ball.setPosition(trackRow.position[0], trackRow.position[1], trackRow.position[2])

                if CONST_GOAL_MINX < trackRow.position[0] < CONST_GOAL_MAXX and (trackRow.position[1] > CONST_FIELD_MAXY or trackRow.position[1] < CONST_FIELD_MINY):
                    game.shotOnGoal = True
                    game.shotOnGoalTeam = teamPlayers[lastHitPlayerName]

                if trackRow.position[2] * math.pow(10, -3) > 2:
                    game.ballWasCrossed = True
                
            else:
                #if the ball is at almost 2 meters from the ground, anyone has the ball.
                if game.ball.position[2] * math.pow(10, -3) > 2:
                    continue

                #eucledian distance between ball and player
                distance = scipy.spatial.distance.euclidean(numpy.asarray(game.ball.position), numpy.asarray(trackRow.position))

                if distance < playersNearTheBall[0].distance:
                    #if the actual nearest player is different from the "new" nearest:
                    if playersNearTheBall[0].ID != "" and playersSensors[trackRow.ID] != playersSensors[playersNearTheBall[0].ID]:
                        playersNearTheBall[1] = playersNearTheBall[0]
                    playersNearTheBall[0] = _PlayerNearTheBall(trackRow.ID, distance, trackRow.actualTime)

                #debug
                '''if trackRow.ID == 71:
                    ax.plot(trackRow.position[0], trackRow.position[1], 'o')'''

                #check the distances every 0.3 seconds
                if timePassed * math.pow(10, -12) >= 0.3:

                    mostNearestPlayer = playersNearTheBall[0]
                    mostNearestPlayer2 = playersNearTheBall[1]

                    #debug
                    '''plt.show()
                    break'''

                    #if the minimum is near the ball
                    if mostNearestPlayer.distance < CONST_DISTANCEGRANULARITY:

                        '''Utilities'''
                        playerName = playersSensors[mostNearestPlayer.ID]
                        playerTeam = teamPlayers[playerName]

                        #if the player is changed and he possess the ball
                        if playersSensors[mostNearestPlayer.ID] != lastHitPlayerName and game.players[playersSensors[mostNearestPlayer.ID]].hasTheBall(trackRow.actualTime):
                            print playersSensors[mostNearestPlayer.ID] + " possessed the ball at second: " + str((3092000000000.0 + mostNearestPlayer.timestamp) * math.pow(10, -12))

                            if lastHitPlayerName == "":
                                lastHitPlayerName = playersSensors[mostNearestPlayer.ID]

                            if (game.ball_is_out or playerName != lastHitPlayerName or game.gameWasInterrupted) and not game.shotOnGoal:

                                if game.ball_is_out or playerTeam != teamPlayers[lastHitPlayerName]:

                                    '''Failed passages'''
                                    game.players[lastHitPlayerName].unsuccessfulPassages += 1 * CONST_METRICWEIGHT_PASSAGES

                                    '''tackles'''
                                    if not game.ball_is_out:
                                        distanceBetweenPlayers = abs(mostNearestPlayer.distance - mostNearestPlayer2.distance)
                                        #if the second nearest player is the last one who touched the ball, it means that he losed a tackle
                                        if distanceBetweenPlayers <= CONST_DISTANCEGRANULARITY * 2 and lastHitPlayerName == playersSensors[mostNearestPlayer2.ID]:
                                            game.players[lastHitPlayerName].lostTackles += 1 * CONST_METRICWEIGHT_TACKLES
                                            if trackRow.actualTime >= CONST_SECONDMATCHSTART:
                                                game.players[playerName].wonTackles += 1 * CONST_METRICWEIGHT_TACKLES_SF * CONST_METRICWEIGHT_TACKLES
                                            else:
                                                game.players[playerName].wonTackles += 1 * CONST_METRICWEIGHT_TACKLES

                                    '''Unsuccessful passages in the defensive third'''
                                    if trackRow.position[1] < CONST_FIELD_MINY + (CONST_FIELD_YLENGTH/2)/3 or trackRow.position[1] > CONST_FIELD_MAXY-(CONST_FIELD_YLENGTH/2)/3:
                                        game.players[lastHitPlayerName].unsuccessfulPassagesDefensiveThird += 1 * CONST_METRICWEIGHT_PASSAGESDT

                                    '''Unsuccessful crosses'''
                                    if game.ballWasCrossed:
                                        game.players[lastHitPlayerName].unsuccessfulCrosses += 1 * CONST_METRICWEIGHT_CROSSES

                                #If a player passed the ball to another player in the same team
                                elif playerTeam == teamPlayers[lastHitPlayerName] and not game.ball_is_out:
                                    game.insertPassage(lastHitPlayerName, playerName)

                                '''Successful passages/crosses'''
                                if playerTeam == teamPlayers[lastHitPlayerName] and not game.ball_is_out:
                                    game.players[lastHitPlayerName].successfulPassages += 1 * CONST_METRICWEIGHT_PASSAGES
                                    if trackRow.position[1] < CONST_FIELD_MINY+(CONST_FIELD_YLENGTH/2)/3 or trackRow.position[1] > CONST_FIELD_MAXY-(CONST_FIELD_YLENGTH/2)/3:
                                        game.players[lastHitPlayerName].successfulPassagesDefensiveThird += 1 * CONST_METRICWEIGHT_PASSAGESDT
                                    if game.ballWasCrossed:
                                        game.players[lastHitPlayerName].successfulCrosses += 1 * CONST_METRICWEIGHT_CROSSES

                            #if I changed the team with a shot, it's probable that I had a true shotOnGoal
                            if game.shotOnGoal and game.shotOnGoalTeam != playerTeam:
                                #if it is the goalkeeper
                                if len(playersSensors2[playerName]) == 4 or game.ball_is_out:
                                    print "## GOAL FAILED AT " + str((3092000000000.0 + mostNearestPlayer.timestamp) * math.pow(10, -12))
                                    game.players[lastHitPlayerName].shotsFailed += 1 * CONST_METRICWEIGHT_TOTALSHOTS
                                else:
                                    print "## GOAL DONE AT " + str((3092000000000.0 + mostNearestPlayer.timestamp) * math.pow(10, -12)) + " ora player "+str(playerName)
                                    game.players[lastHitPlayerName].shotsOnGoal += 1 * CONST_METRICWEIGHT_SHOTSONGOAL

                            game.shotOnGoal = False
                            game.shotOnGoalTeam = 0
                            lastHitPlayerName = playersSensors[mostNearestPlayer.ID]
                            game.ball_is_out = False
                            game.ballWasCrossed = False

                            #action insertion
                            l = len(actions[actionPointer].playersInvolved)

                            #If there is a new action
                            if l != 0 and teamPlayers[actions[actionPointer].playersInvolved[l-1]] != teamPlayers[playersSensors[mostNearestPlayer.ID]]:
                                actionPointer += 1
                                actions.insert(actionPointer, _Action(0, 0, []))
                                l = 0

                            actions[actionPointer].duration += (trackRow.actualTime - actions[actionPointer].duration)
                            if l == 0 or actions[actionPointer].playersInvolved[l-1] != playersSensors[mostNearestPlayer.ID]:
                                actions[actionPointer].numberPassages += 1
                                actions[actionPointer].playersInvolved.append(playersSensors[mostNearestPlayer.ID])

                        game.players[playersSensors[mostNearestPlayer.ID]].lastTouchTimestamp = trackRow.actualTime
                        game.gameWasInterrupted = False

                    playersNearTheBall[0] = _PlayerNearTheBall("", sys.maxint, "")
                    playersNearTheBall[1] = playersNearTheBall[0]
                    checkedPlayers = {}
                    lastDistanceCheckTimestamp = trackRow.actualTime

        #Save everything for the export
        DBSCANTable = db['dbscan']
        for p in teamPlayers.keys():
            player = game.players[p]
            DBSCANTable.insert(dict(name=player.ID,
                                    unsuccessfulPassages=player.unsuccessfulPassages,
                                    successfulPassages=player.successfulPassages,

                                    unsuccessfulPassagesDefensiveThird=player.unsuccessfulPassagesDefensiveThird,

                                    unsuccessfulCrosses=player.unsuccessfulCrosses,
                                    successfulCrosses=player.successfulCrosses,

                                    shotsOnGoal=player.shotsOnGoal,
                                    shotsFailed=player.shotsFailed,

                                    wonTackles=player.wonTackles,
                                    lostTackles=player.lostTackles))

        rowID = 1
        actionRowID = 1
        for action in actions:
            db['actions'].insert(dict(id=actionRowID,
                                          duration=action.duration,
                                          numberOfPassages=action.numberPassages))
            for p in action.playersInvolved:
                db['actions_players'].insert(dict(id=rowID,
                                                  actionID=actionRowID,
                                          playerInvolved=p))
                rowID += 1

            actionRowID += 1

        print "### PASSAGES MATRIX ###"

            
    pass