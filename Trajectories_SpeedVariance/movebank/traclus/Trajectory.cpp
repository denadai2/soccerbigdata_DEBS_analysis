// Trajectory.cpp : implementation file
//

#include "Trajectory.h"


// CTrajectory

CTrajectory::CTrajectory()
{
	m_trajectoryId = -1;
	m_nDimensions = 2;
	m_nPoints = 0;
	m_nPartitionPoints = 0;
}

CTrajectory::CTrajectory(int id, int nDimensions)
{
	m_trajectoryId = id;
	m_nDimensions = nDimensions;
	m_nPoints = 0;
	m_nPartitionPoints = 0;
}

CTrajectory::~CTrajectory()
{
	m_nPoints = 0;
}


// CTrajectory member functions
