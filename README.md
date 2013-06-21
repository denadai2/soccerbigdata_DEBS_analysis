soccerbigdata_DEBS_analysis
===========================

The full dataset can be retrived here: http://www.orgs.ttu.edu/debs2013/index.php?goto=cfchallengedetails but we also split it in a 5m dataset, placed in the project. However for the passage patterns and similar players, our advice is to use the full dataset.


Similar players and passage patterns
===========================
The whole code is written in Python, so it is multi-platform. The main.py file is the analyzator of the dataset, that need to be placed in the same directory with name "full-game".
After the analyzation everything will be saved in a database and you can use kmeans.py for the k-means clustering, hierarchical.py for the agglomerative clustering and edge.py for the passage patterns.



Trajectories and speed variance
===========================
The dataset must be placed in the "DEBSDATA" directory with the name "full-game".
To start the project just run "TrajClustering_SpeedPerform.py"

You need thhese libraries:

numpy
string
select
sys
matplotlib.pyplot
subprocess
uuidthe 
datetime
csv
math
copy
thread
pylab
exceptions


through the stdin you can chose the several options to discover the functionality.

The part of trajectory clustering needs a c program, there already esist the compiled version of a program
in /movebank/bin/ (traclus). It is compiled on linux Ubuntu 64 bit. If it should not work recompile it.
You will find the make file in the source code directory /movebank/bin/traclus.

For every doubt or comments michele.linardi@alice.it/me@marcodena.it



