#ifndef __CLUSTER_H__
#define __CLUSTER_H__

// CCluster command target

#include "MDPoint.h"

#include <fstream>
#include <vector>
using namespace std;

class CCluster
{
public:
	CCluster();
	CCluster(int id, int nDimensions);
	virtual ~CCluster();
private:
	int m_clusterId;		// the identifier of this cluster
	int m_nDimensions;		// the dimensionality of this cluster
	int m_nTrajectories;	// the number of trajectories belonging to this cluster
	int m_nDifferentSegment;	// the number of segment belonging to this cluster
	int m_nPoints;			// the number of points constituting a cluster 
	vector<CMDPoint> m_pointArray;	// the array of the cluster points
public:
	void SetId(int id) { m_clusterId = id; }
	const int GetId() const { return m_clusterId; }
	void SetDensity(int density) { m_nTrajectories = density; }
	void SetDifferentSegmentInCluster(int numberOfSegment) { m_nDifferentSegment = numberOfSegment; }
	const int GetDensity() const { return m_nTrajectories; }
	void AddPointToArray(CMDPoint point) { m_pointArray.push_back(point); m_nPoints++; }
	bool WriteCluster(ofstream& outFile);
};

#endif
