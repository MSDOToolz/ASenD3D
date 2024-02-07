#ifndef SPATIALGRID
#define SPATIALGRID
#include "MeshNode.h"
#include "MeshElement.h"
#include "MeshFace.h"

class MeshEnt {
private:
	MeshNode* nd;
	MeshElement* el;
	MeshFace* fc;
	MeshEnt* next;

public:
	MeshEnt();

	void setPt(MeshNode* newNd);

	void setPt(MeshElement* newEl);

	void setPt(MeshFace* newFc);

	void setNext(MeshEnt* newNext);

	MeshNode* getPt(MeshNode* pt);

	MeshElement* getPt(MeshElement* pt);

	MeshFace* getPt(MeshFace* pt);

	MeshEnt* getNext();

};

class EntityList {
private:
	MeshEnt* first;
	MeshEnt* last;
	int length;

public:
	EntityList();

	void addEntry(MeshNode* newNd);

	void addEntry(MeshElement* newEl);

	void addEntry(MeshFace* newFc);

	int copyToArray(MeshEnt* outLst[], int maxLen);

	~EntityList();
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
	EntityList* listAr;

public:
	SpatialGrid();

	void initialize(double xRange[], double xSpacing, double yRange[], double ySpacing, double zRange[], double zSpacing);
	
	void addEnt(MeshNode* newNd, double crd[]);

	void addEnt(MeshElement* newEl, double crd[]);

	void addEnt(MeshFace* newFc, double crd[]);

	int getInXYZRange(MeshEnt* outLst[], int maxLen, double xRange[], double yRange[], double zRange[]);

	int getInRadius(MeshEnt* outList[], int maxLen, double pt[], double rad);

	~SpatialGrid();
};

#endif
