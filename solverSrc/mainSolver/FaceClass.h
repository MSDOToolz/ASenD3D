#ifndef FACE
#define FACE
#include "DiffDoubClass.h"
#include "NodeClass.h"
#include "DesignVariableClass.h"

class Face {
	private:
	    int numNds;
	    int locNodes[6];
		int globNodes[6];
		bool onSurf;
		Face *next;
		
	public:
	    Face(int newNumNds);
		
		void setNode(int place, int locNd, int globNd);
		
		void setOnSurface(bool newOnSurf);
		
		void setNext(Face *newNext);
		
		void sortedNodes(int srtNds[]);
		
		int getNumNds();

		int* getLocNds();
		
		bool onSurface();
		
		int getLowNd();
		
		Face* getNext();

		//dup1

		void getAreaNormal(Doub& area, Doub norm[], Node* ndAr[], DesignVariable* dvAr[]);

		//end dup
 
//skip 
 
//DiffDoub versions: 
		//dup1

		void getAreaNormal(DiffDoub& area, DiffDoub norm[], Node* ndAr[], DesignVariable* dvAr[]);

		//end dup
 
//end skip 

};

class FacePt {
    public:
		Face* fc;
		FacePt* next;

		FacePt(Face* newFc);
};

class FaceList {
	private:
	    Face *first;
		Face *last;
		int len;
		
	public:
	    FaceList();
		
		void addFace(Face *newFc);
		
		void addIfAbsent(Face *newFc);
		
		Face* getFirst();
		
		~FaceList();
};

class FacePtList {
    private:
		FacePt* first;
		FacePt* last;
		int len;

    public:

		FacePtList();

		void addFace(FacePt* newFpt);

		bool addIfAbsent(Face* newFc);

		~FacePtList();
};

#endif