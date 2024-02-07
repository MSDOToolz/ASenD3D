#ifndef MESHNODE
#define MESHNODE

class MeshNode {
private:
	int label;
	double coord[3];
	MeshNode* next;

public:
	MeshNode(double newCrd[]);

	void setLabel(int newLab);

	void setNext(MeshNode* newNd);

	int getLabel();

	double* getCrd();

	MeshNode* getNext();
};

class MNList {
private:
	MeshNode* first;
	MeshNode* last;
	int length;

public:
	MNList();

	void addEnt(MeshNode* newNd);

	MeshNode* getFirst();

	int getLength();

	~MNList();
};

#endif