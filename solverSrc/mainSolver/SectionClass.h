#ifndef SECTION
#define SECTION
#include <string>

class Material {
	private:
	    std::string name;
		double density;
		double modulus[3];
		double poissonRatio[3];
		double shearMod[3];
		double stiffness[36];
		double conductivity[6];
		double expansion[6];
		double specHeat;
		double damping[36];
		double maxTenStress[3];
		double maxCompStress[3];
		double maxShearStress[3];
		double maxTenStrain[3];
		double maxCompStrain[3];
		double maxShearStrain[3];
		double maxStrEng;
		double maxMises;
		Material *nextMat;
		
	public:
	    Material(std::string newName);
		
		void setDensity(double newDen);
		
		void setModulus(double newMod[]);
		
		void setPoissonRatio(double newPR[]);
		
		void setShearMod(double newMod[]);
		
		void setStiffness(int row, int col, double val);
		
		void setConductivity(double newCond[]);
		
		void setExpansion(double newExp[]);
		
		void setSpecHeat(double newSpecHeat);

		void setDamping(int row, int col, double val);
		
		void setMaxTenStress(double newMaxStr[]);
		
		void setMaxCompStress(double newMaxStr[]);
		
		void setMaxShearStress(double newMaxStr[]);
		
		void setMaxTenStrain(double newMaxStr[]);
		
		void setMaxCompStrain(double newMaxStr[]);
		
		void setMaxShearStrain(double newMaxStr[]);
		
		void setMaxStrEng(double newMax);
		
		void setMaxMises(double newMax);
		
		void setNext(Material *newNext);
		
		std::string getName();
		
		double getDensity();
		
		double* getModulus();
		
		double* getPoissonRatio();
		
		double* getShearMod();

        double* getStiffMat();

		double* getThermExp();

		double* getConductivity();

		double getSpecificHeat();

		double* getDamping();
		
		Material* getNext();	
		
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

		void destroy();
};

class Layer {
	private:
	    std::string matName;
		Material *matPtr;
		double thickness;
		double angle;
		Layer *nextLay;
		
    public:
	    Layer(std::string newNm);
		
		void setThickness(double newThk);
		
		void setAngle(double newAngle);
		
		std::string getMatName();
		
		Material* getMatPt();
		
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
		
		void addLayer(Layer *newLay);
		
		int getLength();
		
		Layer* getFirst();

		void destroy();
};

class Section {
	private:
	    std::string type;
		std::string elSetName;
		std::string matName;
		Material *matPtr;
		double orientation[9];
		double zOffset;
		LayerList layers;
		double area;
		double areaMoment[5];
		double polarMoment;
		double stiffness[36];
		double mass[36];
		double damping[36];
		double expLoadCoef[6];
		double conductivity;
		double specHeat;
		double massPerEl;
		double potCoef;
		double potExp;
		double dampCoef;
		double dampExp;
		Section *nextSection;
		
	public:
	    Section(std::string newType);
		
		void setElset(std::string newSet);
		
		void setMaterial(std::string newMat);
		
		void setMatPtr(Material *newMat);
		
		void setOrientation(double newOri[]);
		
		void setZOffset(double newZOff);
		
		void addLayer(Layer *newLay);
		
		void setArea(double newArea);
		
		void setAreaMoment(double newI[]);
		
		void setPolarMoment(double newJ);
		
		void setStiffness(int row, int col, double val);
		
		void setMass(int row, int col, double val);

		void setDamping(int row, int col, double val);
		
		void setExpLd(double newExpLd[]);
		
		void setConductivity(double newCond);
		
		void setSpecHeat(double specHeat);

		void setMassPerEl(double newMass);

		void setPotCoef(double newCoef);

		void setPotExp(double newExp);

		void setDampCoef(double newCoef);

		void setDampExp(double newExp);
		
		std::string getElset();
		
		std::string getMaterial();
		
		Material* getMatPtr();
		
		double* getOrientation();
		
		double getZOffset();

		int getNumLayers();
		
		Layer* getFirstLayer();
		
		double getArea();
		
		double* getAreaMoment();
		
		double getPolarMoment();
		
		double* getStiffMat();

		double* getMassMat();

		double* getDampMat();

		double* getExpLoad();

		double getConductivity();

		double getSpecificHeat();

		double getMassPerEl();

		double getPotCoef();

		double getPotExp();

		double getDampCoef();

		double getDampExp();
		
		Section* getNext();
		
		void setNext(Section *newNext);

		void destroy();
};

class SectionList {
	private:
		Section* firstSec;
		Section* lastSec;
		int length;
		
	public:
	    SectionList();
		
		void addSection(Section *newSec);
		
		int getLength();
		
		Section* getFirst();

		void destroy();
};

#endif
