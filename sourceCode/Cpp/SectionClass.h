#ifndef SECTION
#define SECTION
#include <string>

class Material {
	private:
	    std::string name;
		double density;
		double elastic[9];
		double stiffness[21];
		double conductivity[6];
		double expansion[6];
		double specHeat;
		double maxStress[9];
		double maxStrain[9];
		double maxStrEng;
		double maxMises;
		
	public:
	    Material(std::string newName);
		
		double getDensity();
		
		double* getElastic();

        double* getStiffMat();
		
};

class MaterialList {
	private:
	    Material *firstMat;
		Material *lastMat;
		int length;
		
	public:
	    MaterialList();
		
		void addMaterial(Material *newMat);
		
		int getLength();
		
		Material* getFirst();
};

class Layer {
	private:
	    std::string matName;
		Material *matPtr;
		double thickness;
		double angle;
		Layer *nextLay;
		
    public:
	    Layer(std::string newNm, double newThk, double newAng);
		
		std::string getMatName();
		
		double getThickness();
		
		double getAngle();
		
		Layer* getNext();
		
		void setNext(Layer *newNext);
		
		void setMatPtr(Material *newPtr);
};

class LayerList {
	private:
	    Layer *firstLay;
		Layer *lastLay;
		int length;
		
	public:
	    LayerList();
		
		void addLayer(std::string newNm, double newThk, double newAng);
		
		int getLength();
		
		Layer* getFirst();
};

class Section {
	private:
	    std::string type;
		std::string elSetName;
		std::string matName;
		Material *matPtr;
		double orientation[6];
		double zOffset;
		LayerList layers;
		double area;
		double areaMoment[5];
		double polarMoment;
		double stiffness[21];
		double mass[21];
		double expLoadCoef[6];
		double conductivity;
		double specHeat;
		Section *nextSection;
		
	public:
	    Section(std::string newType);
		
		Material* getMaterial();
		
		void setNext(Section *newNext);
};

class SectionList {
	private:
	    Node *firstSec;
		Node *lastSec;
		int length;
		
	public:
	    SectionList();
		
		void addSection(Section *newSec);
		
		int getLength();
		
		Section* getFirst();
};

#endif