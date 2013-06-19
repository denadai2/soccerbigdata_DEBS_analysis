#include <stdio.h>
#include <iostream>

#include "TraClusterDoc.h"

int main(int argc, char** argv)
{
	if (argc != 2)
	{
		cout << "USAGE: traclus <input file>" << endl;
		cout << "       - <input file>: a trajectory file" << endl;
		return 0;
	}

	CTraClusterDoc* document = new CTraClusterDoc;
	
	if (!document->OnOpenDocument(argv[1]))
	{
		cout << "Cannot open a trajectory file" << endl;
		return 0;
	}

	float epsParam;
	int minLnsParam;

	if (!document->OnEstimateParameter(epsParam, minLnsParam))
	{
		cout << "Cannot estimate a parameter value" << endl;
		return 0;
	}

	cout << "The estimated values of parameters are as follows." << endl;
	cout << "- eps: around " << epsParam << endl;
	cout << "- MinLns: " << minLnsParam << "~" << minLnsParam + 2 << endl;

	return 1;
}
