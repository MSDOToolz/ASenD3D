#include <string>

class Material {
	private:
	    string name;
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
	    Material(string newName="") {
			name = newName;
		}
		
		double getDensity() {
			return density;
		}
};