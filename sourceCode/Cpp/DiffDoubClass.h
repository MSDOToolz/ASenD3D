#ifndef DIFFDOUB
#define DIFFDOUB
#include <cmath>

class Doub {
    public:
	    double val;
		double dval;
	
	    Doub();
		
		bool diffType();
		
		void setVal(double inp);
		
		void setVal(Doub& inp);
		
		void add(Doub& inp);
		
		void sub(Doub& inp);
		
		void mult(Doub& inp);
		
		void dvd(Doub& inp);
		
		void sqt();
		
		void sqr();
		
		void sn();
		
		void cs();
		
		void asn();
		
		void atn();
		
		void neg();
};

class DiffDoub {
    public:
	    double val;
		double dval;
		double tmp;
		double tmp2;
	
	    DiffDoub();
		
		bool diffType();
		
		void setVal(double inp);
		
		void setVal(double inp, double dinp);

		void setVal(Doub& inp);
		
		void setVal(DiffDoub& inp);
		
		void add(DiffDoub& inp);
		
		void sub(DiffDoub& inp);
		
		void mult(DiffDoub& inp);
		
		void dvd(DiffDoub& inp);

		void sqt();
		
		void sqr();
		
		void sn();
		
		void cs();
		
		void asn();
		
		void atn();
		
		void neg();
};

class Diff2Doub {
    public:
	    double val;
		double dv1;
		double dv2;
		double dv12;
		double tmp;
		double tmp2;
		double tmp3;
		double tmp4;
	
	    Diff2Doub();
		
		bool diffType();
		
		void setVal(double inp);
		
		void setVal(double inp, double dinp);
		
		void setVal(double inp, double dinp1, double dinp2);
		
		void setVal(Diff2Doub& inp);
		
		void add(Diff2Doub& inp);
		
		void sub(Diff2Doub& inp);
		
		void mult(Diff2Doub& inp);
		
		void dvd(Diff2Doub& inp);

		void sqt();
		
		void sqr();
		
		void sn();
		
		void cs();
		
		void atn();
		
		void neg();
};

#endif