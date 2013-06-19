#ifndef __TRAJECTORY_H__
#define __TRAJECTORY_H__

// CTrajectory command target

#include "MDPoint.h"

#include <vector>
using namespace std;

class CTrajectory
{
friend class CClusterGen;
public:
	CTrajectory();
	CTrajectory(int id, int nDimensions);
	virtual ~CTrajectory();
private:
	int m_trajectoryId;		// the identifier of this trajectory
	int m_nDimensions;		// the dimensionality of this trajectory
	int m_nPoints;			// the number of points constituting a trajectory 
	vector<CMDPoint> m_pointArray;	// the array of the trajectory points
	int m_nPartitionPoints;		// the number of partition points in a trajectory
	vector<CMDPoint> m_partitionPointArray;	// the array of the partition points
public:
	void SetId(int id) { m_trajectoryId = id; }
	const int GetId() const { return m_trajectoryId; }
	void AddPointToArray(CMDPoint point) { m_pointArray.push_back(point); m_nPoints++; }
	void AddPartitionPointToArray(CMDPoint point) { m_partitionPointArray.push_back(point); m_nPartitionPoints++; }
};

#endif
