#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "TraClusterDoc.h"

int main(int argc, char** argv)
{
	if (argc != 5)
	{
		cout << "USAGE: traclus <input file> <output file> [<eps> <minLns>]" << endl;
		cout << "       - <input file>: a trajectory file" << endl;
		cout << "       - <output file>: a cluster file" << endl;
		cout << "       - <eps>: the parameter epsilon (float)" << endl;
		cout << "       - <MinLns>: the parameter MinLnsi (integer)" << endl;
		return 0;
	}

	CTraClusterDoc* document = new CTraClusterDoc;
	
	if (!document->OnOpenDocument(argv[1]))
	{
		cout << "Cannot open a trajectory file" << endl;
		return 0;
	}

	if (!document->OnClusterGenerate(argv[2], atof(argv[3]), atoi(argv[4])))
	{
		cout << "Cannot generate a cluster file" << endl;
		return 0;
	}

	cout << "Clustering has been completed sucessfully" << endl;

	return 1;
}
