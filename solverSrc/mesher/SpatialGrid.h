#ifndef SPATIALGRID
#define SPATIALGRID
#include <list>
#include <vector>
#include "MeshNode.h"
#include "MeshElement.h"
#include "MeshFace.h"

//class MeshEnt {
//private:
//	MeshNode* nd;
//	MeshElement* el;
//	MeshFace* fc;
//	MeshEnt* next;
//
//public:
//	MeshEnt();
//
//	void setPt(MeshNode* newNd);
//
//	void setPt(MeshElement* newEl);
//
//	void setPt(MeshFace* newFc);
//
//	void setNext(MeshEnt* newNext);
//
//	MeshNode* getPt(MeshNode* pt);
//
//	MeshElement* getPt(MeshElement* pt);
//
//	MeshFace* getPt(MeshFace* pt);
//
//	MeshEnt* getNext();
//
//};
//
//class EntityList {
//private:
//	MeshEnt* first;
//	MeshEnt* last;
//	int length;
//
//public:
//	EntityList();
//
//	void addEntry(MeshNode* newNd);
//
//	void addEntry(MeshElement* newEl);
//
//	void addEntry(MeshFace* newFc);
//
//	int copyToArray(MeshEnt* outLst[], int maxLen);
//
//	~EntityList();
//};

class IntList {
public:
	std::list<int> iLst;

	IntList();

	int copy_to_vector(std::vector<int>& in_vec, int st_i, int max_len);
};

class SpatialGrid {
private:
	double xMin;
	double xSp;
	int xBins;
	double yMin;
	double ySp;
	int yBins;
	double zMin;
	double zSp;
	int zBins;
	std::vector<IntList> listAr;

public:
	SpatialGrid();

	void initialize(double xRange[], double xSpacing, double yRange[], double ySpacing, double zRange[], double zSpacing);

	void addEnt(int label, double crd[]);

	int getInXYZRange(std::vector<int>& outLst, int maxLen, double xRange[], double yRange[], double zRange[]);

	int getInRadius(std::vector<int>& outList, int maxLen, double pt[], double rad);
};

#endif
