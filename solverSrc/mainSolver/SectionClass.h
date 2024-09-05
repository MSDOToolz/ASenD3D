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

		~MaterialList();
};

class Fluid {
private:
	std::string name;
	double viscosity;
	double idealGas;
	double thermCond;
	double specHeat;
	Fluid* nextFl;

public:
	Fluid(std::string newName);

	void setViscosity(double newVis);

	void setIdealGas(double newIG);

	void setThermCond(double newTC);

	void setSpecHeat(double newSH);

	void setNext(Fluid* newNext);

	double getViscosity();

	double getIdealGas();

	double getThermCond();

	double getSpecHeat();

	Fluid* getNext();

};

class FluidList {
private:
	Fluid* firstFl;
	Fluid* lastFl;
	int length;

public:
	FluidList();

	void addFluid(Fluid* newFl);

	int getLength();

	Fluid* getFirst();

	~FluidList();
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

		~LayerList();
};

class Section {
	private:
	    std::string type;
		std::string elSetName;
		std::string matName;
		std::string flName;
		Material *matPtr;
		Fluid* flPtr;
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
		double condCoef;
		double radCoef;
		double denVisCoef;
		double tempVisCoef;
		double gradVisCoef;
		double enthCoef;
		double enthExp;
		double presCoef;
		double presExp;
		double refTemp;
		double refDen;
		double refGradV;
		double refEnth;
		Section *nextSection;
		
	public:
	    Section(std::string newType);
		
		void setElset(std::string newSet);
		
		void setMaterial(std::string newMat);

		void setFluid(std::string newFl);
		
		void setMatPtr(Material *newMat);

		void setFlPtr(Fluid* newFl);
		
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

		void setCondCoef(double newCoef);

		void setRadCoef(double newCoef);

		void setDenVisCoef(double newCoef);

		void setTempVisCoef(double newCoef);

		void setGradVisCoef(double newCoef);

		void setEnthCoef(double newCoef);

		void setEnthExp(double newExp);

		void setPresCoef(double newCoef);

		void setPresExp(double newExp);

		void setRefTemp(double newRT);

		void setRefDen(double newDen);

		void setRefGradV(double newGV);

		void setRefEnth(double newEnth);
		
		std::string getElset();
		
		std::string getMaterial();

		std::string getFluid();
		
		Material* getMatPtr();

		Fluid* getFluidPtr();
		
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

		double getCondCoef();

		double getRadCoef();

		double getRefTemp();

		double getDenVisCoef();

		double getTempVisCoef();

		double getGradVisCoef();

		double getEnthCoef();

		double getEnthExp();

		double getPresCoef();

		double getPresExp();

		double getRefTemp();

		double getRefDen();

		double getRefGradV();

		double getRefEnth();
		
		Section* getNext();
		
		void setNext(Section *newNext);
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

		int getMaxNumLayers();

		~SectionList();
};

#endif
