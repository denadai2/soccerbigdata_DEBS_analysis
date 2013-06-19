// ClusterGen.cpp : implementation file
//

#include "TraClusterDoc.h"
#include "ClusterGen.h"
#include <algorithm>
#include <cmath>
#include <climits>


// CClusterGen

CClusterGen::CClusterGen()
{
}

CClusterGen::CClusterGen(CTraClusterDoc* document)
{
	m_document = document;

	m_idArray.clear();
	m_lineSegmentPointArray.clear();
}

CClusterGen::~CClusterGen()
{
}



// CClusterGen message handlers

bool CClusterGen::PartitionTrajectory()
{
	vector<CTrajectory*>::iterator iter;
	for (iter = m_document->m_trajectoryList.begin(); iter != m_document->m_trajectoryList.end(); iter++)
		FindOptimalPartition(*iter);

	if (!StoreClusterComponentIntoIndex())
		return false;

	return true;
}

bool CClusterGen::StoreClusterComponentIntoIndex()
{
	int nDimensions = m_document->m_nDimensions;
	CMDPoint* startPoint;
	CMDPoint* endPoint;

	m_nTotalLineSegments = 0;

	vector<CTrajectory*>::iterator iter;
	for (iter = m_document->m_trajectoryList.begin(); iter != m_document->m_trajectoryList.end(); iter++)
	{
		CTrajectory* pTrajectory = *iter;
		for (int i = 0; i < pTrajectory->m_nPartitionPoints - 1; i++)
		{
			// convert an n-dimensional line segment into a 2n-dimensional point
			// i.e., the first n-dimension: the start point
			//       the last n-dimension: the end point
			startPoint = &(pTrajectory->m_partitionPointArray[i]);
			endPoint = &(pTrajectory->m_partitionPointArray[i + 1]);

			if (MeasureDistanceFromPointToPoint(startPoint, endPoint) < MIN_LINESEGMENT_LENGTH)
				continue;

			m_nTotalLineSegments++;

			CMDPoint lineSegmentPoint(nDimensions * 2);
			for (int j = 0; j < nDimensions; j++) 
			{
				lineSegmentPoint.SetCoordinate(j, startPoint->GetCoordinate(j));
				lineSegmentPoint.SetCoordinate(nDimensions + j, endPoint->GetCoordinate(j));
			}

			LineSegmentId id;
			id.trajectoryId = pTrajectory->GetId();
			id.order = i;

			m_idArray.push_back(id);
			m_lineSegmentPointArray.push_back(lineSegmentPoint);
		}
	}

	return true;
}

bool CClusterGen::PerformDBSCAN(float eps, int minLns)
{
	m_epsParam = eps;
	m_minLnsParam = minLns;

	m_currComponentId = 0;

	m_componentIdArray.resize(m_nTotalLineSegments, -1);
	for (int i = 0; i < m_nTotalLineSegments; i++)
		m_componentIdArray[i] = UNCLASSIFIED;

	for (int i = 0; i < m_nTotalLineSegments; i++)
	{
		if (m_componentIdArray[i] == UNCLASSIFIED)
			if (ExpandDenseComponent(i, m_currComponentId, eps, minLns))
				m_currComponentId++;
	}

	return true;
}

bool CClusterGen::ConstructCluster()
{
	// this step consists of two substeps
	// notice that the result of the previous substep is used in the following substeps

	if (!ConstructLineSegmentCluster())
		return false;

	if(!StoreLineSegmentCluster())
		return false;

	return true;
}


#include <math.h>

// NOTE: that this version does not guarantee the optimal partition

void CClusterGen::FindOptimalPartition(CTrajectory* pTrajectory)
{
	int nPoints = pTrajectory->m_nPoints;
	int startIndex = 0, length;
	int fullPartitionMDLCost, partialPartitionMDLCost;

	// add the start point of a trajectory
	pTrajectory->AddPartitionPointToArray(pTrajectory->m_pointArray[0]);

	for (;;)
	{
		fullPartitionMDLCost = partialPartitionMDLCost = 0;
		
		for (length = 1; startIndex + length < nPoints; length++)
		{
			// compute the total length of a trajectory
			fullPartitionMDLCost += ComputeModelCost(pTrajectory, startIndex + length - 1, startIndex + length);

			// compute the sum of (1) the length of a cluster component and
			//                    (2) the perpendicular and angle distances
			partialPartitionMDLCost = ComputeModelCost(pTrajectory, startIndex, startIndex + length) +
			                          ComputeEncodingCost(pTrajectory, startIndex, startIndex + length);

			if (fullPartitionMDLCost + MDL_COST_ADVANTAGE < partialPartitionMDLCost)
			{
				pTrajectory->AddPartitionPointToArray(pTrajectory->m_pointArray[startIndex + length - 1]);

				startIndex = startIndex + length - 1;
				length = 0;

				break;
			}
		}

		// if we reach at the end of a trajectory
		if (startIndex + length >= nPoints)
			break;
	}

	// add the end point of a trajectory
	pTrajectory->AddPartitionPointToArray(pTrajectory->m_pointArray[nPoints - 1]);

	return;
}

#define LOG2(_x)	(float)(log((float)(_x)) / log((float)2))

int CClusterGen::ComputeModelCost(CTrajectory* pTrajectory, int startPIndex, int endPIndex)
{
	CMDPoint* lineSegmentStart = &(pTrajectory->m_pointArray[startPIndex]);
	CMDPoint* lineSegmentEnd = &(pTrajectory->m_pointArray[endPIndex]);

	float distance = MeasureDistanceFromPointToPoint(lineSegmentStart, lineSegmentEnd);
	if (distance < 1.0) distance = 1.0;		// to take logarithm

	return (int)ceil(LOG2(distance));
}

int CClusterGen::ComputeEncodingCost(CTrajectory* pTrajectory, int startPIndex, int endPIndex)
{
	CMDPoint* clusterComponentStart;
	CMDPoint* clusterComponentEnd;
	CMDPoint* lineSegmentStart;
	CMDPoint* lineSegmentEnd;
	float perpendicularDistance;
	float angleDistance;
	int encodingCost = 0;

	clusterComponentStart = &(pTrajectory->m_pointArray[startPIndex]);
	clusterComponentEnd = &(pTrajectory->m_pointArray[endPIndex]);

	for (int i = startPIndex; i < endPIndex; i++)
	{
		lineSegmentStart = &(pTrajectory->m_pointArray[i]);
		lineSegmentEnd = &(pTrajectory->m_pointArray[i + 1]);

		perpendicularDistance = MeasurePerpendicularDistance(clusterComponentStart, clusterComponentEnd, lineSegmentStart, lineSegmentEnd);
		angleDistance = MeasureAngleDisntance(clusterComponentStart, clusterComponentEnd, lineSegmentStart, lineSegmentEnd);

		if (perpendicularDistance < 1.0) perpendicularDistance = 1.0;	// to take logarithm
		if (angleDistance < 1.0) angleDistance = 1.0;					// to take logarithm
		encodingCost += ((int)ceil(LOG2(perpendicularDistance)) + (int)ceil(LOG2(angleDistance)));
	}

	return encodingCost;
}

float CClusterGen::MeasurePerpendicularDistance(CMDPoint* s1, CMDPoint* e1, CMDPoint* s2, CMDPoint* e2)
{
	// we assume that the first line segment is longer than the second one
	float distance1;	// the distance from a start point to the cluster component
	float distance2;	// the distance from an end point to the cluster component

	distance1 = MeasureDistanceFromPointToLineSegment(s1, e1, s2);
	distance2 = MeasureDistanceFromPointToLineSegment(s1, e1, e2);

	// if the first line segment is exactly the same as the second one, the perpendicular distance should be zero
	if (distance1 == 0.0 && distance2 == 0.0) return 0.0;

	// return (d1^2 + d2^2) / (d1 + d2) as the perpendicular distance
	return ((pow(distance1, 2) + pow(distance2, 2)) / (distance1 + distance2));
}

float CClusterGen::MeasureDistanceFromPointToLineSegment(CMDPoint* s, CMDPoint* e, CMDPoint* p)
{
	int nDimensions = p->GetNDimensions();

	// NOTE: the variables m_vector1 and m_vector2 are declared as member variables

	// construct two vectors as follows
	// 1. the vector connecting the start point of the cluster component and a given point
	// 2. the vector representing the cluster component
	for (int i = 0; i < nDimensions; i++)
	{
		m_vector1.SetCoordinate(i, p->GetCoordinate(i) - s->GetCoordinate(i));
		m_vector2.SetCoordinate(i, e->GetCoordinate(i) - s->GetCoordinate(i));
	}

	// a coefficient (0 <= b <= 1)
	m_coefficient = ComputeInnerProduct(&m_vector1, &m_vector2) / ComputeInnerProduct(&m_vector2, &m_vector2);

	// the projection on the cluster component from a given point
	// NOTE: the variable m_projectionPoint is declared as a member variable

	for (int i = 0; i < nDimensions; i++)
		m_projectionPoint.SetCoordinate(i, s->GetCoordinate(i) + m_coefficient * m_vector2.GetCoordinate(i));

	// return the distance between the projection point and the given point
	return MeasureDistanceFromPointToPoint(p, &m_projectionPoint);
}

float CClusterGen::MeasureDistanceFromPointToPoint(CMDPoint* point1, CMDPoint* point2)
{
	int nDimensions = point1->GetNDimensions();
	float squareSum = 0.0;
	
	for (int i = 0; i < nDimensions; i++)
		squareSum += pow((point2->GetCoordinate(i) - point1->GetCoordinate(i)), 2);

	return sqrt(squareSum);
}

float CClusterGen::ComputeVectorLength(CMDPoint* vector)
{
	int nDimensions = vector->GetNDimensions();
	float squareSum = 0.0;

	for (int i = 0; i < nDimensions; i++)
		squareSum += pow(vector->GetCoordinate(i), 2);

	return sqrt(squareSum);
}

float CClusterGen::ComputeInnerProduct(CMDPoint* vector1, CMDPoint* vector2)
{
	int nDimensions = vector1->GetNDimensions();
	float innerProduct = 0.0;

	for (int i = 0; i < nDimensions; i++)
		innerProduct += (vector1->GetCoordinate(i) * vector2->GetCoordinate(i));

	return innerProduct;
}

float CClusterGen::MeasureAngleDisntance(CMDPoint* s1, CMDPoint* e1, CMDPoint* s2, CMDPoint* e2)
{
	int nDimensions = s1->GetNDimensions();
	
	// NOTE: the variables m_vector1 and m_vector2 are declared as member variables

	// construct two vectors representing the cluster component and a line segment, respectively
	for (int i = 0; i < nDimensions; i++)
	{
		m_vector1.SetCoordinate(i, e1->GetCoordinate(i) - s1->GetCoordinate(i));
		m_vector2.SetCoordinate(i, e2->GetCoordinate(i) - s2->GetCoordinate(i));
	}

	// we assume that the first line segment is longer than the second one
	// i.e., vectorLength1 >= vectorLength2
	float vectorLength1 = ComputeVectorLength(&m_vector1);
	float vectorLength2 = ComputeVectorLength(&m_vector2);

	// if one of two vectors is a point, the angle distance becomes zero
	if (vectorLength1 == 0.0 || vectorLength2 == 0.0) return 0.0;

	// compute the inner product of the two vectors
	float innerProduct = ComputeInnerProduct(&m_vector1, &m_vector2);

	// compute the angle between two vectors by using the inner product
	float cosTheta = innerProduct / (vectorLength1 * vectorLength2);
	// compensate the computation error (e.g., 1.00001)
	// cos(theta) should be in the range [-1.0, 1.0]
	// START ...
	if (cosTheta > 1.0) cosTheta = 1.0; 
	if (cosTheta < -1.0) cosTheta = -1.0;
	// ... END
	float sinTheta = sqrt(1 - pow(cosTheta, 2));
	// if 90 <= theta <= 270, the angle distance becomes the length of the line segment
	// if (cosTheta < -1.0) sinTheta = 1.0;

	return (vectorLength2 * sinTheta);
}

bool CClusterGen::ExpandDenseComponent(int index, int componentId, float eps, int minDensity)
{
	set <int> seeds, seedResult;
	set <int>::iterator iter;
	int currIndex;

	ExtractStartAndEndPoints(index, &m_startPoint1, &m_endPoint1);
	ComputeEPSNeighborhood(&m_startPoint1, &m_endPoint1, eps, &seeds);
	
	if ((int)seeds.size() < minDensity)			// not a core line segment
	{
		m_componentIdArray[index] = NOISE;
		return false;
	}
	else
	{
		for (iter = seeds.begin(); iter != seeds.end(); iter++)
			m_componentIdArray[*iter] = componentId;

		seeds.erase(seeds.find(index));

		while (!seeds.empty())
		{
			currIndex = *(seeds.begin());
			ExtractStartAndEndPoints(currIndex, &m_startPoint1, &m_endPoint1);
			ComputeEPSNeighborhood(&m_startPoint1, &m_endPoint1, eps, &seedResult);

			if ((int)seedResult.size() >= minDensity)
			{
				for (iter = seedResult.begin(); iter != seedResult.end(); iter++)
				{
					if (m_componentIdArray[*iter] == UNCLASSIFIED || m_componentIdArray[*iter] == NOISE)
					{
						if (m_componentIdArray[*iter] == UNCLASSIFIED)
							seeds.insert(*iter);

						m_componentIdArray[*iter] = componentId;
					}
				}
			}

			seeds.erase(seeds.find(currIndex));
		}

		return true;
	}
}

void CClusterGen::ComputeEPSNeighborhood(CMDPoint* startPoint, CMDPoint* endPoint, float eps, set<int>* result)
{
	result->clear();

	for (int j = 0; j < m_nTotalLineSegments; j++)
	{
		ExtractStartAndEndPoints(j, &m_startPoint2, &m_endPoint2);
		
		float distance = ComputeDistanceBetweenTwoLineSegments(startPoint, endPoint, &m_startPoint2, &m_endPoint2);
		// if the distance is below the threshold, this line segment belongs to the eps-neighborhood
		if (distance <= eps) result->insert(j);
	}

	return;
}

float CClusterGen::ComputeDistanceBetweenTwoLineSegments(CMDPoint* startPoint1, CMDPoint* endPoint1, CMDPoint* startPoint2, CMDPoint* endPoint2)
{
	float perpendicularDistance, parallelDistance, angleDistance;

	SubComputeDistanceBetweenTwoLineSegments(startPoint1, endPoint1, startPoint2, endPoint2, perpendicularDistance, parallelDistance, angleDistance);

	return (perpendicularDistance + parallelDistance + angleDistance);
}

void CClusterGen::SubComputeDistanceBetweenTwoLineSegments(CMDPoint* startPoint1, CMDPoint* endPoint1, CMDPoint* startPoint2, CMDPoint* endPoint2, float& perpendicularDistance, float& parallelDistance, float& angleDistance)
{
	float perDistance1, perDistance2;
	float parDistance1, parDistance2;
	float length1, length2;

	// the length of the first line segment
	length1 = MeasureDistanceFromPointToPoint(startPoint1, endPoint1);
	// the length of the second line segment
	length2 = MeasureDistanceFromPointToPoint(startPoint2, endPoint2);

	// compute the perpendicular distance and the parallel distance
	// START ...
	if (length1 > length2)
	{
		perDistance1 = MeasureDistanceFromPointToLineSegment(startPoint1, endPoint1, startPoint2);
		if (m_coefficient < 0.5) parDistance1 = MeasureDistanceFromPointToPoint(startPoint1, &m_projectionPoint);
		else parDistance1 = MeasureDistanceFromPointToPoint(endPoint1, &m_projectionPoint);

		perDistance2 = MeasureDistanceFromPointToLineSegment(startPoint1, endPoint1, endPoint2);
		if (m_coefficient < 0.5) parDistance2 = MeasureDistanceFromPointToPoint(startPoint1, &m_projectionPoint);
		else parDistance2 = MeasureDistanceFromPointToPoint(endPoint1, &m_projectionPoint);
	}
	else
	{
		perDistance1 = MeasureDistanceFromPointToLineSegment(startPoint2, endPoint2, startPoint1);
		if (m_coefficient < 0.5) parDistance1 = MeasureDistanceFromPointToPoint(startPoint2, &m_projectionPoint);
		else parDistance1 = MeasureDistanceFromPointToPoint(endPoint2, &m_projectionPoint);

		perDistance2 = MeasureDistanceFromPointToLineSegment(startPoint2, endPoint2, endPoint1);
		if (m_coefficient < 0.5) parDistance2 = MeasureDistanceFromPointToPoint(startPoint2, &m_projectionPoint);
		else parDistance2 = MeasureDistanceFromPointToPoint(endPoint2, &m_projectionPoint);
	}

	// compute the perpendicular distance; take (d1^2 + d2^2) / (d1 + d2)
	if (!(perDistance1 == 0.0 && perDistance2 == 0.0)) 
		perpendicularDistance = ((pow(perDistance1, 2) + pow(perDistance2, 2)) / (perDistance1 + perDistance2));
	else
		perpendicularDistance = 0.0;

	// compute the parallel distance; take the minimum
	parallelDistance = (parDistance1 < parDistance2) ? parDistance1 : parDistance2;
	// ... END

	// compute the angle distance
	// START ...
	// MeasureAngleDisntance() assumes that the first line segment is longer than the second one
	if (length1 > length2)
		angleDistance = MeasureAngleDisntance(startPoint1, endPoint1, startPoint2, endPoint2);
	else
		angleDistance = MeasureAngleDisntance(startPoint2, endPoint2, startPoint1, endPoint1);
	// ... END

	return;
}

bool CClusterGen::ConstructLineSegmentCluster()
{
	int nDimensions = m_document->m_nDimensions;

	m_lineSegmentClusters = new LineSegmentCluster [m_currComponentId];

	// initialize the list of line segment clusters
	// START ...
	for (int i = 0; i < m_currComponentId; i++)
	{
		m_lineSegmentClusters[i].lineSegmentClusterId = i;
		m_lineSegmentClusters[i].nLineSegments = 0;
		m_lineSegmentClusters[i].nClusterPoints = 0;
		m_lineSegmentClusters[i].nTrajectories = 0;
		m_lineSegmentClusters[i].enabled = false;
	}
	// ... END

	// accumulate the direction vector of a line segment
	for (int i = 0; i < m_nTotalLineSegments; i++)
	{
		int componentId = m_componentIdArray[i];
		if (componentId >= 0)
		{
			for (int j = 0; j < nDimensions; j++)
			{
				float difference = m_lineSegmentPointArray[i].GetCoordinate(nDimensions + j) - m_lineSegmentPointArray[i].GetCoordinate(j);
				float currSum = m_lineSegmentClusters[componentId].avgDirectionVector.GetCoordinate(j) + difference;
				m_lineSegmentClusters[componentId].avgDirectionVector.SetCoordinate(j, currSum);
			}
			m_lineSegmentClusters[componentId].nLineSegments++;
		}
	}

	// compute the average direction vector of a line segment cluster
	// START ...
	float vectorLength1, vectorLength2, innerProduct;
	float cosTheta, sinTheta;

	m_vector2.SetCoordinate(0, 1.0);
	m_vector2.SetCoordinate(1, 0.0);

	for (int i = 0; i < m_currComponentId; i++)
	{
		LineSegmentCluster* clusterEntry = &(m_lineSegmentClusters[i]);

		for (int j = 0; j < nDimensions; j++)
			clusterEntry->avgDirectionVector.SetCoordinate(j, clusterEntry->avgDirectionVector.GetCoordinate(j) / (float)clusterEntry->nLineSegments);

		vectorLength1 = ComputeVectorLength(&(clusterEntry->avgDirectionVector));
		vectorLength2 = 1.0;

		innerProduct = ComputeInnerProduct(&(clusterEntry->avgDirectionVector), &m_vector2);
		cosTheta = innerProduct / (vectorLength1 * vectorLength2);
		if (cosTheta > 1.0) cosTheta = 1.0; 
		if (cosTheta < -1.0) cosTheta = -1.0;
		sinTheta = sqrt(1 - pow(cosTheta, 2));

		if (clusterEntry->avgDirectionVector.GetCoordinate(1) < 0) sinTheta = -sinTheta;

		clusterEntry->cosTheta = cosTheta;
		clusterEntry->sinTheta = sinTheta;
	}
	// ... END

	// summarize the information about line segment clusters
	// the structure for summarization is as follows
	// [lineSegmentClusterId, nClusterPoints, clusterPointArray, nTrajectories, { trajectoryId, ... }]
	for (int i = 0; i < m_nTotalLineSegments; i++)
	{
		if (m_componentIdArray[i] >= 0)		// if the componentId < 0, it is a noise
			RegisterAndUpdateLineSegmentCluster(m_componentIdArray[i], i);
	}

	set <int, less<int> > trajectories;
	for (int i = 0; i < m_currComponentId; i++)
	{
		LineSegmentCluster* clusterEntry = &(m_lineSegmentClusters[i]);

		// a line segment cluster must have trajectories more than the minimum threshold
		if (clusterEntry->nTrajectories >= m_minLnsParam)
		{
			clusterEntry->enabled = true;

			// DEBUG: count the number of trajectories that belong to clusters
			for (int j = 0; j < (int)clusterEntry->trajectoryIdList.size(); j++)
				trajectories.insert(clusterEntry->trajectoryIdList[j]);
			// ... DEBUG

			ComputeRepresentativeLines(clusterEntry);
		}
		else
		{
			clusterEntry->candidatePointList.clear();
			clusterEntry->clusterPointArray.clear();
			clusterEntry->trajectoryIdList.clear();
		}
	}
	
	// DEBUG: compute the ratio of trajectories that belong to clusters
	m_document->m_clusterRatio = (float)trajectories.size() / (float)m_document->m_nTrajectories;

	return true;
}

#define GET_X_ROTATION(_x,_y,_cos,_sin) ((_x)*(_cos) + (_y)*(_sin))
#define GET_Y_ROTATION(_x,_y,_cos,_sin) (-(_x)*(_sin) + (_y)*(_cos))
#define GET_X_REV_ROTATION(_x,_y,_cos,_sin) ((_x)*(_cos) - (_y)*(_sin))
#define GET_Y_REV_ROTATION(_x,_y,_cos,_sin) ((_x)*(_sin) + (_y)*(_cos))

void CClusterGen::RegisterAndUpdateLineSegmentCluster(int componentId, int lineSegmentId)
{
	LineSegmentCluster* clusterEntry = &(m_lineSegmentClusters[componentId]);

	// the start and end values of the first dimension (e.g., the x value in the 2-dimension)
	// NOTE: this program code works only for the 2-dimensional data

	CMDPoint* aLineSegment = &(m_lineSegmentPointArray[lineSegmentId]);
	float orderingValue1 = GET_X_ROTATION(aLineSegment->GetCoordinate(0), aLineSegment->GetCoordinate(1), clusterEntry->cosTheta, clusterEntry->sinTheta); 
	float orderingValue2 = GET_X_ROTATION(aLineSegment->GetCoordinate(2), aLineSegment->GetCoordinate(3), clusterEntry->cosTheta, clusterEntry->sinTheta); 

	CandidateClusterPoint existingCandidatePoint, newCandidatePoint;
	int i, j;

	// sort the line segment points by the coordinate of the first dimension
	// simply use the insertion sort algorithm
	// START ...
	list<CandidateClusterPoint>::iterator iter1 = clusterEntry->candidatePointList.begin();
	for (i = 0; i < (int)clusterEntry->candidatePointList.size(); i++)
	{
		existingCandidatePoint = *iter1;
		if (existingCandidatePoint.orderingValue < orderingValue1)
			iter1++;
		else
			break;
	}

	newCandidatePoint.orderingValue = orderingValue1;
	newCandidatePoint.lineSegmentId = lineSegmentId;
	newCandidatePoint.startPointFlag = true;
	if (i == 0)
		clusterEntry->candidatePointList.push_front(newCandidatePoint);
	else if (i >= (int)clusterEntry->candidatePointList.size())
		clusterEntry->candidatePointList.push_back(newCandidatePoint);
	else
		clusterEntry->candidatePointList.insert(iter1, newCandidatePoint);

	list<CandidateClusterPoint>::iterator iter2 = clusterEntry->candidatePointList.begin();
	for (j = 0; j < (int)clusterEntry->candidatePointList.size(); j++)
	{
		existingCandidatePoint = *iter2;
		if (existingCandidatePoint.orderingValue < orderingValue2)
			iter2++;
		else
			break;
	}

	newCandidatePoint.orderingValue = orderingValue2;
	newCandidatePoint.lineSegmentId = lineSegmentId;
	newCandidatePoint.startPointFlag = false;
	if (j == 0)
		clusterEntry->candidatePointList.push_front(newCandidatePoint);
	else if (j >= (int)clusterEntry->candidatePointList.size())
		clusterEntry->candidatePointList.push_back(newCandidatePoint);
	else
		clusterEntry->candidatePointList.insert(iter2, newCandidatePoint);
	// ... END

	int trajectoryId = m_idArray[lineSegmentId].trajectoryId;

	// store the identifier of the trajectories that belong to this line segment cluster
	if (std::find(clusterEntry->trajectoryIdList.begin(), clusterEntry->trajectoryIdList.end(), trajectoryId) == clusterEntry->trajectoryIdList.end())
	{
		clusterEntry->trajectoryIdList.push_back(trajectoryId);
		clusterEntry->nTrajectories++;
	}

	return;
}

void CClusterGen::ComputeRepresentativeLines(LineSegmentCluster* clusterEntry)
{
	set <int> lineSegments, insertionList, deletionList;
	set <int>::iterator iter1, iter2, iter3;
	CandidateClusterPoint candidatePoint, nextCandidatePoint;
	float prevOrderingValue = 0.0;

	int nClusterPoints = 0;
	lineSegments.clear();

	// sweep the line segments in a line segment cluster
	list<CandidateClusterPoint>::iterator iter = clusterEntry->candidatePointList.begin();
	while (iter != clusterEntry->candidatePointList.end())
	{
		insertionList.clear();
		deletionList.clear();

		do {
			candidatePoint = *iter++;

			// check whether this line segment has begun or not
			iter1 = lineSegments.find(candidatePoint.lineSegmentId);
			if (iter1 == lineSegments.end())	{				// if there is no matched element,
				insertionList.insert(candidatePoint.lineSegmentId);		// this line segment begins at this point
				lineSegments.insert(candidatePoint.lineSegmentId);
			}
			else									// if there is a matched element,
				deletionList.insert(candidatePoint.lineSegmentId);		// this line segment ends at this point

			// check whether the next line segment begins or ends at the same point
			if (iter != clusterEntry->candidatePointList.end())
				nextCandidatePoint = *iter;
			else
				break;
		} while (candidatePoint.orderingValue == nextCandidatePoint.orderingValue);

		// check if a line segment is connected to another line segment in the same trajectory
		// if so, delete one of the line segments to remove duplicates
		for (iter2 = insertionList.begin(); iter2 != insertionList.end(); iter2++)
		{
			for (iter3 = deletionList.begin(); iter3 != deletionList.end(); iter3++)
			{
				if (*iter2 == *iter3)
				{
					lineSegments.erase(lineSegments.find(*iter3));
					deletionList.erase(iter3);		// now deleted
					break;
				}
			}
			for (iter3 = deletionList.begin(); iter3 != deletionList.end(); iter3++)
			{
				if (m_idArray[*iter2].trajectoryId == m_idArray[*iter3].trajectoryId)
				{
					lineSegments.erase(lineSegments.find(*iter3));
					deletionList.erase(iter3);		// now deleted
					break;
				}
			}
		}

		// if the current density exceeds a given threshold
		if ((int)(lineSegments.size()) >= m_minLnsParam)
		{
			if (abs(candidatePoint.orderingValue - prevOrderingValue) > ((double)MIN_LINESEGMENT_LENGTH / 1.414))
			{
				ComputeAndRegisterClusterPoint(clusterEntry, candidatePoint.orderingValue, &lineSegments);
				prevOrderingValue = candidatePoint.orderingValue;
				nClusterPoints++;
			}
		}

		// delete the line segment that is not connected to another line segment
		for (iter3 = deletionList.begin(); iter3 != deletionList.end(); iter3++)
			lineSegments.erase(lineSegments.find(*iter3));
	}

	if (nClusterPoints >= 2)
		clusterEntry->nClusterPoints = nClusterPoints;
	else
	{
		// there is no representative trend in this line segment cluster
		clusterEntry->enabled = false;
		clusterEntry->candidatePointList.clear();
		clusterEntry->clusterPointArray.clear();
		clusterEntry->trajectoryIdList.clear();
	}

	return;
}

void CClusterGen::ComputeAndRegisterClusterPoint(LineSegmentCluster* clusterEntry, float currValue, set<int>* lineSegments)
{
	int nDimensions = m_document->m_nDimensions;
	int nLineSegmentsInSet = (int)(lineSegments->size());
	CMDPoint clusterPoint(nDimensions);
	CMDPoint sweepPoint(nDimensions);
	
	set <int>::iterator iter;
	for (iter = lineSegments->begin(); iter != lineSegments->end(); iter++)
	{
		// get the sweep point of each line segment
		// this point is parallel to the current value of the sweeping direction
		GetSweepPointOfLineSegment(clusterEntry, currValue, *iter, &sweepPoint);
		// compute the average of all the sweep points
		for (int i = 0; i < nDimensions; i++)
			clusterPoint.SetCoordinate(i, clusterPoint.GetCoordinate(i) + (sweepPoint.GetCoordinate(i) / (float)nLineSegmentsInSet));
	}

	// NOTE: this program code works only for the 2-dimensional data
	float origX, origY;
	origX = GET_X_REV_ROTATION(clusterPoint.GetCoordinate(0), clusterPoint.GetCoordinate(1), clusterEntry->cosTheta, clusterEntry->sinTheta);
	origY = GET_Y_REV_ROTATION(clusterPoint.GetCoordinate(0), clusterPoint.GetCoordinate(1), clusterEntry->cosTheta, clusterEntry->sinTheta);
	clusterPoint.SetCoordinate(0, origX);
	clusterPoint.SetCoordinate(1, origY);

	// register the obtained cluster point (i.e., the average of all the sweep points)
	clusterEntry->clusterPointArray.push_back(clusterPoint);

	return;
}

void CClusterGen::GetSweepPointOfLineSegment(LineSegmentCluster* clusterEntry, float currValue, int lineSegmentId, CMDPoint* sweepPoint)
{
	CMDPoint* lineSegmentPoint = &(m_lineSegmentPointArray[lineSegmentId]);		// 2n-dimensional point
	float coefficient;

	// NOTE: this program code works only for the 2-dimensional data
	float newStartX, newEndX, newStartY, newEndY;
	newStartX = GET_X_ROTATION(lineSegmentPoint->GetCoordinate(0), lineSegmentPoint->GetCoordinate(1), clusterEntry->cosTheta, clusterEntry->sinTheta);
	newEndX   = GET_X_ROTATION(lineSegmentPoint->GetCoordinate(2), lineSegmentPoint->GetCoordinate(3), clusterEntry->cosTheta, clusterEntry->sinTheta);
	newStartY = GET_Y_ROTATION(lineSegmentPoint->GetCoordinate(0), lineSegmentPoint->GetCoordinate(1), clusterEntry->cosTheta, clusterEntry->sinTheta);
	newEndY   = GET_Y_ROTATION(lineSegmentPoint->GetCoordinate(2), lineSegmentPoint->GetCoordinate(3), clusterEntry->cosTheta, clusterEntry->sinTheta);

	coefficient = (currValue - newStartX) / (newEndX - newStartX);
	sweepPoint->SetCoordinate(0, currValue);
	sweepPoint->SetCoordinate(1, newStartY + coefficient * (newEndY - newStartY));

	return;
}

bool CClusterGen::StoreLineSegmentCluster()
{
	int currClusterId = 0;

	for (int i = 0; i < m_currComponentId; i++)
	{
		if (m_lineSegmentClusters[i].enabled)
		{
			// store the clusters finally identified
			// START ...
			CCluster* pClusterItem = new CCluster(currClusterId, m_document->m_nDimensions); 
			m_document->m_clusterList.push_back(pClusterItem);

			for (int j = 0; j < m_lineSegmentClusters[i].nClusterPoints; j++)
				pClusterItem->AddPointToArray(m_lineSegmentClusters[i].clusterPointArray[j]);

			pClusterItem->SetDensity(m_lineSegmentClusters[i].nTrajectories);
			pClusterItem->SetDifferentSegmentInCluster(m_lineSegmentClusters[i].nLineSegments);

			currClusterId++;	// increase the number of final clusters
			// ... END
		}
	}

	m_document->m_nClusters = currClusterId;

	return true;
}

bool CClusterGen::EstimateParameterValue(float& epsParam, int& minLnsParam)
{
	float entropy, minEntropy = (float)INT_MAX;
	float eps, minEps = (float)INT_MAX;
	int totalSize, minTotalSize = INT_MAX;
	set <int> seeds;

	int* EpsNeighborhoodSize = new int [m_nTotalLineSegments];

	for (eps = (float)20.0; eps <= (float)40.0; eps = eps + (float)1.0)
	{
		entropy = (float)0.0;
		totalSize = 0;
		seeds.clear();

		for (int i = 0; i < m_nTotalLineSegments; i++)
		{
			ExtractStartAndEndPoints(i, &m_startPoint1, &m_endPoint1);
			ComputeEPSNeighborhood(&m_startPoint1, &m_endPoint1, eps, &seeds);
			EpsNeighborhoodSize[i] = (int)seeds.size();
			totalSize += (int)seeds.size();
			seeds.clear();
		}

		for (int i = 0; i < m_nTotalLineSegments; i++)
			entropy += ((float)EpsNeighborhoodSize[i] / (float)totalSize) * LOG2((float)EpsNeighborhoodSize[i] / (float)totalSize);
		entropy = -entropy;

		if (entropy < minEntropy)
		{
			minEntropy = entropy;
			minTotalSize = totalSize;
			minEps = eps;
		}

		fprintf(stdout, ".");	fflush(stdout);
	} 

	delete EpsNeighborhoodSize;

	// setup output arguments
	epsParam    = minEps;
	minLnsParam = (int)ceil((float)minTotalSize / (float)m_nTotalLineSegments);

	fprintf(stdout, "\n");	fflush(stdout);

	return true;
}
