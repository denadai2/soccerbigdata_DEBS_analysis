#ifndef __TRACLUSTERDOC_H__
#define __TRACLUSTERDOC_H__

// TraClusterDoc.h : interface of the CTraClusterDoc class
//

#include "Trajectory.h"
#include "Cluster.h"
#include "Common.h"

#include <vector>
using namespace std;

class CTraClusterDoc
{
public:
	CTraClusterDoc();
	virtual ~CTraClusterDoc();
// Attributes
public:
	int m_nDimensions;
	int m_nTrajectories;
	int m_nClusters;
	float m_clusterRatio;
	int m_maxNPoints;
	vector<CTrajectory*> m_trajectoryList;
	vector<CCluster*> m_clusterList;
public:
	bool OnOpenDocument(char* inputFileName);
	bool OnClusterGenerate(char* clusterFileName, float epsParam, int minLnsParam);
	bool OnEstimateParameter(float& epsParam, int& minLnsParam);
};

#endif
