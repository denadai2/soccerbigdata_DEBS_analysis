// Cluster.cpp : implementation file
//

#include "Cluster.h"
#include <iomanip>


// CCluster

CCluster::CCluster()
{
	m_clusterId = -1;
	m_nDimensions = 2;
	m_nTrajectories = 0;
	m_nPoints = 0;
	m_nDifferentSegment = 0;
}

CCluster::CCluster(int id, int nDimensions)
{
	m_clusterId = id;
	m_nDimensions = nDimensions;
	m_nTrajectories = 0;
	m_nPoints = 0;
	m_nDifferentSegment = 0;
}

CCluster::~CCluster()
{
	m_clusterId = -1;
	m_nTrajectories = 0;
	m_nPoints = 0;
	m_nDifferentSegment = 0;
}


// CCluster member functions

bool CCluster::WriteCluster(ofstream& outFile)
{
	outFile << (int) m_clusterId << ' ' << (int)m_pointArray.size() << ' ' << (int)m_nDifferentSegment << ' ';

	for (int i = 0; i < (int)m_pointArray.size(); i++)
	{
		outFile << fixed << setprecision(1);
		outFile << m_pointArray[i].GetCoordinate(0) << ' ' << m_pointArray[i].GetCoordinate(1) << ' ';
	}

	outFile << endl;

	return true;
}
