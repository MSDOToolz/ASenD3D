#include <string>
#include "MaterialClass.cpp"
#include "LayerClass.cpp"

class Section {
	private:
	    string type;
		string elSetName;
		string matName;
		Material *matPtr;
		double orientation[6];
		double zOffset;
		Layer *firstLayer;
		double area;
		double areaMoment[5];
		double polarMoment;
		double stiffness[21];
		double mass[21];
		double expLoadCoef[6];
		double conductivity;
		double specHeat;
		
	public:
	    Section(string newType="solid") {
			type = newType;
		}
};