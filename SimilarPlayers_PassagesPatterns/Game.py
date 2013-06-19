__author__ = 'denadai2'

import csv
import datetime
import math
import numpy



class Ball:
    position = (0,0,0)
    positions = []

    def setPosition(self,x,y,z):
        self.position = (x,y,z)


class Player:
    ID = ""
    successfulPassages = 0
    unsuccessfulPassages = 0

    successfulPassagesDefensiveThird = 0
    unsuccessfulPassagesDefensiveThird = 0

    successfulCrosses = 0
    unsuccessfulCrosses = 0

    shotsOnGoal = 0
    shotsFailed = 0

    lostTackles = 0
    wonTackles = 0

    lastTouchTimestamp = 0

    def __init__(self, ID):
        self.ID = ID

    #if tha ball is possessed for more than 0.2 seconds, the Player is considered to have the ball
    def hasTheBall(self, actualTimestamp):
        diff = (actualTimestamp - self.lastTouchTimestamp) * math.pow(10, -12)
        if 0.2 <= diff < 1:
            return True

        return False


class Game:

    #(from,to)
    interruptions = []

    lastBall_ID = 4
    ball_is_out = False
    ball = Ball()
    ballWasCrossed = False
    #attemption of goal
    shotOnGoal = False
    shotOnGoalTeam = 0
    #gameWasInterrupted is true when there is an interruption between two actions
    gameWasInterrupted = False

    #players in the field
    players = {}

    #playerID => number (for the passages matrix)
    keyMap = {}
    #passages matrix
    passagesMatrix = []



    def configure(self, teamPlayers):
        self.passagesMatrix = numpy.zeros(shape=(len(teamPlayers), len(teamPlayers)))

        #inizialize players object
        count = 0
        for playerKey in teamPlayers.keys():
            self.players[playerKey] = Player(playerKey)
            self.keyMap[playerKey] = count
            count += 1

        #initialize interruptions
        self.loadInterruptions()

    def loadInterruptions(self):
        firstLineRead = False
        tempTime = -1
        with open('interruptions.csv', 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=';')
            for row in reader:
                if firstLineRead:
                    if row[2] == '0':
                        row[2] = "00:00:00.0"
                    struct_time = datetime.datetime.strptime(row[2], "%H:%M:%S.%f")
                    struct_time = float(struct_time.minute*60+struct_time.hour*60*60+struct_time.second)*math.pow(10, 12)+(struct_time.microsecond)*math.pow(10, 6)
                    if tempTime != -1:
                        self.interruptions.append((tempTime-3092000000000.0, struct_time-3092000000000.0))
                        tempTime = -1
                    else:
                        tempTime = struct_time
                else:
                    firstLineRead = True

    def insertPassage(self, playerIDFrom, playerIDTo):
        fromID = self.keyMap[playerIDFrom]
        toID = self.keyMap[playerIDTo]

        self.passagesMatrix[fromID][toID] += 1