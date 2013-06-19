#ifndef __CLUSTERGEN_H__
#define __CLUSTERGEN_H__

#include "TraClusterDoc.h"
#include "Trajectory.h"
#include "Cluster.h"

#include <vector>
#include <list>
#include <set>
using namespace std;


// CClusterGen

typedef struct LineSegmentId {
	int trajectoryId;
	int order;
} LineSegmentId;

typedef struct CandidateClusterPoint {
	float orderingValue;
	int lineSegmentId;
	bool startPointFlag;
} CandidateClusterPoint;

typedef struct LineSegmentCluster {
	int lineSegmentClusterId;
	int nLineSegments;
	CMDPoint avgDirectionVector;
	float cosTheta, sinTheta;
	list<CandidateClusterPoint> candidatePointList;
	int nClusterPoints;
	vector<CMDPoint> clusterPointArray;
	int nTrajectories;
	vector<int> trajectoryIdList;
	bool enabled;
} LineSegmentCluster;

// used for performing the DBSCAN algorithm
const int UNCLASSIFIED = -2;
const int NOISE = -1;

// used for InsertClusterPoint() and ReplaceClusterPoint()
enum PointLocation {
	HEAD = 0,
	TAIL = 1
};

class CClusterGen
{
public:
	CClusterGen();
	CClusterGen(CTraClusterDoc* document);
	virtual ~CClusterGen();
public:
	CTraClusterDoc* m_document;
private:
	float m_epsParam;
	int m_minLnsParam;
	int m_nTotalLineSegments;
	int m_currComponentId;
	// the number of dense components discovered until now
	vector<int> m_componentIdArray;
	// the list of line segment clusters
	LineSegmentCluster* m_lineSegmentClusters;
	// programming trick: avoid frequent execution of the new and delete operations
	CMDPoint m_startPoint1, m_endPoint1, m_startPoint2, m_endPoint2;
	CMDPoint m_vector1, m_vector2;
	CMDPoint m_projectionPoint;
	float m_coefficient;

	vector<LineSegmentId> m_idArray;
	vector<CMDPoint> m_lineSegmentPointArray;
public:
	bool PartitionTrajectory();
	bool PerformDBSCAN(float eps, int minLns);
	bool ConstructCluster();
	bool EstimateParameterValue(float& eps, int& minLns);
private:
	void FindOptimalPartition(CTrajectory* pTrajectory);
	int ComputeModelCost(CTrajectory* pTrajectory, int startPIndex, int endPIndex);
	int ComputeEncodingCost(CTrajectory* pTrajectory, int startPIndex, int endPIndex);
	float MeasurePerpendicularDistance(CMDPoint* s1, CMDPoint* e1, CMDPoint* s2, CMDPoint* e2);
	float MeasureDistanceFromPointToLineSegment(CMDPoint* s, CMDPoint* e, CMDPoint* p);
	float MeasureDistanceFromPointToPoint(CMDPoint* point1, CMDPoint* point2);
	float ComputeVectorLength(CMDPoint* vector);
	float ComputeInnerProduct(CMDPoint* vector1, CMDPoint* vector2);
	float MeasureAngleDisntance(CMDPoint* s1, CMDPoint* e1, CMDPoint* s2, CMDPoint* e2);
	bool StoreClusterComponentIntoIndex();
	bool ExpandDenseComponent(int index, int componentId, float eps, int minDensity);
	void ComputeEPSNeighborhood(CMDPoint* startPoint, CMDPoint* endPoint, float eps, set<int>* result);
	float ComputeDistanceBetweenTwoLineSegments(CMDPoint* startPoint1, CMDPoint* endPoint1, CMDPoint* startPoint2, CMDPoint* endPoint2);
	void SubComputeDistanceBetweenTwoLineSegments(CMDPoint* startPoint1, CMDPoint* endPoint1, CMDPoint* startPoint2, CMDPoint* endPoint2, float& perpendicularDistance, float& parallelDistance, float& angleDistance);
	bool ConstructLineSegmentCluster();
	void RegisterAndUpdateLineSegmentCluster(int componentId, int lineSegmentId);
	void ComputeRepresentativeLines(LineSegmentCluster* clusterEntry);
	void ComputeAndRegisterClusterPoint(LineSegmentCluster* clusterEntry, float currValue, set<int>* lineSegments);
	void GetSweepPointOfLineSegment(LineSegmentCluster* clusterEntry, float currValue, int lineSegmentId, CMDPoint* sweepPoint);
	bool StoreLineSegmentCluster();
	inline void ExtractStartAndEndPoints(int index, CMDPoint* startPoint, CMDPoint* endPoint) {		// for speedup
		// compose the start and end points of the line segment
		for (int i = 0; i < m_document->m_nDimensions; i++)  {
			startPoint->SetCoordinate(i, m_lineSegmentPointArray[index].GetCoordinate(i));
			endPoint->SetCoordinate(i, m_lineSegmentPointArray[index].GetCoordinate(m_document->m_nDimensions + i));
		}
	}
};

#endif
