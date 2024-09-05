#ifndef FLUIDNODE
#define FLUIDNODE

class FluidNode {
private:
	int label;
	int dofIndex[5];
	int sortedRank;
	double coord[3];
	double displacement[3];
	double velocity[3];
	double flDen;
	double flVel[3];
	double flTemp;
	double denDot;
	double velDot[3];
	double tempDot;
	double prevDen;
	double prevVel[3];
	double prevTemp;
	double prevDenDot;
	double prevVelDot[3];
	double prevTempDot;
	double prevDenLF;
	double prevVelLF[3];
	double prevTempLF;
	FluidNode* nextNd;

public:
	FluidNode(int newLab);

	void setCrd(double newCrd[]);

	void setDofIndex(int dof, int index);

	void setSortedRank(int newRank);


};

#endif