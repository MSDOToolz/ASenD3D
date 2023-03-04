#include <string>

class Layer {
	private:
	    string matName;
		double thickness;
		double angle;
		Layer *nextLay;
		
    public:
	    Layer(string newNm, double newThk, double newAng) {
			matName = newNm;
			thickness = newThk;
			angle = newAng;
		}
		
		string getMatName() {
			return matName;
		}
		
		double getThickness() {
			return thickness;
		}
		
		double getAngle() {
			return angle;
		}
};