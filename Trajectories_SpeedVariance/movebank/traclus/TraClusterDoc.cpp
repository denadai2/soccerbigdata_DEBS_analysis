// TraClusterDoc.cpp : implementation of the CTraClusterDoc class
//
#include "TraClusterDoc.h"
#include "ClusterGen.h"
#include <climits>
#include <fstream>
#include <stdio.h>
#include <string.h>

using namespace std;


// CTraClusterDoc construction/destruction

CTraClusterDoc::CTraClusterDoc()
{
	// TODO: add one-time construction code here
	m_nTrajectories = 0;
	m_nClusters = 0;
	m_clusterRatio = 0.0;
}

CTraClusterDoc::~CTraClusterDoc()
{
	m_nTrajectories = 0;
	m_nClusters = 0;
	m_clusterRatio = 0.0;
}


// CTraClusterDoc commands

bool CTraClusterDoc::OnOpenDocument(char* inputFileName)
{
	int nDimensions = 2;		// default dimension = 2
	int nTrajectories = 0;
	int nTotalPoints = 0;
	int trajectoryId;
	int nPoints;
	float value;
	char buffer[1024000];
	char tmp_buf1[1024], tmp_buf2[1024];
	int offset = 0;

	ifstream istr(inputFileName);

	if(!istr)
	{
		fprintf(stderr, "Unable to open input file\n");
		return false;
	}

	istr.getline(buffer, sizeof(buffer));	// the number of dimensions
	sscanf(buffer, "%d", &nDimensions);
	m_nDimensions = nDimensions;
	istr.getline(buffer, sizeof(buffer));	// the number of trajectories
	sscanf(buffer, "%d", &nTrajectories);
	m_nTrajectories = nTrajectories;

	m_maxNPoints = INT_MIN;					// initialize for comparison

	// the trajectory Id, the number of points, the coordinate of a point ...
	for (int i = 0; i < nTrajectories; i++)
	{
		offset = 0;

		istr.getline(buffer, sizeof(buffer));	// each trajectory
		sscanf(buffer, "%d %d", &trajectoryId, &nPoints);
		sscanf(buffer, "%s %s", tmp_buf1, tmp_buf2);
		offset += (int)strlen(tmp_buf1) + (int)strlen(tmp_buf2) + 2;

		if (nPoints > m_maxNPoints) m_maxNPoints = nPoints;
		nTotalPoints += nPoints;

		CTrajectory* pTrajectoryItem = new CTrajectory(trajectoryId, nDimensions); 
		m_trajectoryList.push_back(pTrajectoryItem);

		for (int j = 0; j < nPoints; j++)
		{
			CMDPoint point(nDimensions);	// initialize the CMDPoint class for each point
			
			for (int k = 0; k < nDimensions; k++)
			{
				sscanf(buffer + offset, "%f", &value);
				sscanf(buffer + offset, "%s", tmp_buf1);
				offset += (int)strlen(tmp_buf1) + 1;

				point.SetCoordinate(k, value);
			}

			pTrajectoryItem->AddPointToArray(point);
		}
	}

	return true;
}

bool CTraClusterDoc::OnClusterGenerate(char* clusterFileName, float epsParam, int minLnsParam)
{
	CClusterGen generator(this);

	if (m_nTrajectories == 0)
	{
		fprintf(stderr, "Load a trajectory data set first\n");
		return false;
	}

	// FIRST STEP: Trajectory Partitioning
	if (!generator.PartitionTrajectory())
	{
		fprintf(stderr, "Unable to partition a trajectory\n");
		return false;
	}

	// SECOND STEP: Density-based Clustering
	if (!generator.PerformDBSCAN(epsParam, minLnsParam))
	{
		fprintf(stderr, "Unable to perform the DBSCAN algorithm\n");
		return false;
	}

	// THIRD STEP: Cluster Construction
	if (!generator.ConstructCluster())
	{
		fprintf(stderr, "Unable to construct a cluster\n");
		return false;
	}

	// Write the generated clusters
	// START ...
	ofstream outFile;
	outFile.open(clusterFileName);
	outFile << (int)m_clusterList.size() << endl;

	vector<CCluster*>::iterator iter;
	for (iter = m_clusterList.begin(); iter != m_clusterList.end(); iter++)
	{
		CCluster* pCluster = *iter;
		pCluster->WriteCluster(outFile);
	}

	outFile.close();
	// ... END

	return true;
}

bool CTraClusterDoc::OnEstimateParameter(float& epsParam, int& minLnsParam)
{
	CClusterGen generator(this);

	if (!generator.PartitionTrajectory())
	{
		fprintf(stderr, "Unable to partition a trajectory\n");
		return false;
	}

	if (!generator.EstimateParameterValue(epsParam, minLnsParam))
	{
		fprintf(stderr, "Unable to calculate the entropy\n");
		return false;
	}

	return true;
}
