#ifndef FACE
#define FACE
#include <vector>
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "NodeClass.h"
#include "DesignVariableClass.h"

class Face {
	public:
	    int numNds;
	    int locNodes[6];
		int globNodes[6];
		bool onSurf;

	    Face();
		
		void setNode(int place, int locNd, int globNd);
		
		void sortedNodes(int srtNds[]);
		
		int getLowNd();

		//dup1

		void getAreaNormal(DiffDoub0& area, DiffDoub0 norm[], std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);

		//end dup
 
//skip 
 
//DiffDoub1 versions: 
		//dup1

		void getAreaNormal(DiffDoub1& area, DiffDoub1 norm[], std::vector<Node>& ndAr, std::vector<DesignVariable>& dvAr);

		//end dup
 
//end skip 
 
 
 
};

class FacePtList {
    public:
		std::list<int> fcList;

		FacePtList();

		void addFace(int newI);

		bool addIfAbsent(int newI, std::vector<Face>& globFaces);
};

#endif