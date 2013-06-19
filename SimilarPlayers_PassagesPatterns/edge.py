
import networkx as nx
import matplotlib.pyplot as plt
import dataset
from collections import namedtuple
from collections import defaultdict
from itertools import izip
from pprint import pprint
import numpy as np

def pairwise(iterable):
    "s -> (s0,s1), (s2,s3), (s4, s5), ..."
    a = iter(iterable)
    return izip(a, a)


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


db = dataset.connect('sqlite:///prova.sqlite')
_Action = namedtuple("_Action", "duration numberPassages playersInvolved")
actions=[_Action(0, 0, [])]
_Player = namedtuple("_Player", "efficiency")
players = {}

teamNumber = int(raw_input("Select the team you want to exploit (0-1):"))

for player in db['dbscan']:
    if player['name'] != 'Referee':
        players[player['name']] = _Player(float(player['successfulPassages'] / (player['successfulPassages'] + player['unsuccessfulPassages'])))

result = db.query('SELECT actionID, duration, numberOfPassages, playerInvolved FROM actions, actions_players WHERE actions.id = actions_players.actionID')
actionPointer = 0
for row in db['actions']:
   playersInvolved = db['actions_players'].find(actionID=row['id'])
   p = []
   for player in playersInvolved:
       p.append(player['playerInvolved'])

   actions.append(_Action(row['duration'], row['numberOfPassages'], p))

G = nx.Graph()
nodelist = []
#Add the players
for pKey in players.keys():
    print teamPlayers[pKey]
    if teamPlayers[pKey] == teamNumber:
        nodelist.append(pKey)
        G.add_node(pKey)

#Add the relations between players
adj_list = defaultdict(lambda: defaultdict(lambda: 0))
total = 0
num = 0

for action in actions:
    edges = []
    for i in range(0,(len(action.playersInvolved)-1)):
        if teamPlayers[action.playersInvolved[i]] == teamNumber:
            if adj_list[action.playersInvolved[i+1]][action.playersInvolved[i]] != 0:
                adj_list[action.playersInvolved[i+1]][action.playersInvolved[i]] += 1
            else:
                adj_list[action.playersInvolved[i]][action.playersInvolved[i+1]] += 1

pprint(adj_list.items())

pos=nx.circular_layout(G) # ..random_layout(G) # positions for all nodes
# nodes

total = 0
num = 0
for item in adj_list.keys():
    for item2 in adj_list[item].keys():
        total += adj_list[item][item2]
        num += 1

avg = total/num

for item in adj_list.keys():
    for item2 in adj_list[item].keys():
        if adj_list[item][item2] != 0:
            if teamPlayers[item] == teamNumber:
                G.add_edge(item, item2, weight=adj_list[item][item2])
                nx.draw_networkx_edges(G,pos,edgelist=[(item, item2)], width=(adj_list[item][item2]/avg)+0.1,edge_color='b')


total = 0
num = 0
for v in G.nodes():
    total += players[v].efficiency
    num += 1

avg = total / num

for v in G.nodes():
    nx.draw_networkx_nodes(G, pos, node_size=700, nodelist=[v], node_color=[(111*(players[v].efficiency/avg)/255,132*(players[v].efficiency/avg)/255,37*(players[v].efficiency/avg)/255)])
    print [(122*(players[v].efficiency/avg)/255,122*(players[v].efficiency/avg)/255,122*(players[v].efficiency/avg)/255)]

nx.draw_networkx_labels(G,pos,font_size=16,font_family='sans-serif')
G = nx.Graph()
nx.draw_networkx(G)
plt.show()