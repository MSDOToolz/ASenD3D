#ifndef SECTION
#define SECTION
#include <string>
#include <vector>
#include <list>

class Material {
	public:
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
		
	    Material();

		void setStiffness(int row, int col, double val);

		void setDamping(int row, int col, double val);
		
};

class Fluid {
public:
	std::string name;
	double viscosity;
	double idealGas;
	double thermCond;
	double specHeat;

	Fluid();

};

class Layer {
	public:
	    std::string matName;
		int matPtr;
		double thickness;
		double angle;
		
	    Layer();
};

class Section {
	public:
	    std::string type;
		std::string elSetName;
		std::string matName;
		std::string flName;
		int matPtr;
		int flPtr;
		double orientation[9];
		double zOffset;
		std::list<Layer> layers;
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
		double turbVisCoef;
		double gradVTurbCoef;
		double dissTurbCoef;
		double enthCoef;
		double enthExp;
		double presCoef;
		double presExp;
		double refTemp;
		double refDen;
		double refTurbE;
		double refEnth;
		
	    Section();
		
		void setOrientation(double newOri[]);
		
		void setAreaMoment(double newI[]);
		
		void setStiffness(int row, int col, double val);
		
		void setMass(int row, int col, double val);

		void setDamping(int row, int col, double val);
		
		void setExpLd(double newExpLd[]);

		int getLayerMatPtr(int layi);
		
};

#endif
