#ifndef diffdoub
#define diffdoub
#include <cmath>

class DiffDoub0 {
    public:
	    double val;
		double dval;
	
	    DiffDoub0();
		
		bool diff_type();
		
		void set_val(double inp);

		void set_val(double inp, double dinp);
		
		void set_val(DiffDoub0& inp);
		
		void add(DiffDoub0& inp);
		
		void sub(DiffDoub0& inp);
		
		void mult(DiffDoub0& inp);
		
		void dvd(DiffDoub0& inp);
		
		void sqt();
		
		void sqr();
		
		void sn();
		
		void cs();
		
		void asn();
		
		void atn();
		
		void neg();
};

class DiffDoub1 {
    public:
	    double val;
		double dval;
		double tmp;
		double tmp2;
	
	    DiffDoub1();
		
		bool diff_type();
		
		void set_val(double inp);
		
		void set_val(double inp, double dinp);

		void set_val(DiffDoub0& inp);
		
		void set_val(DiffDoub1& inp);
		
		void add(DiffDoub1& inp);
		
		void sub(DiffDoub1& inp);
		
		void mult(DiffDoub1& inp);
		
		void dvd(DiffDoub1& inp);

		void sqt();
		
		void sqr();
		
		void sn();
		
		void cs();
		
		void asn();
		
		void atn();
		
		void neg();
};

class DiffDoub2 {
    public:
	    double val;
		double dv1;
		double dv2;
		double dv12;
		double tmp;
		double tmp2;
		double tmp3;
		double tmp4;
	
	    DiffDoub2();
		
		bool diff_type();
		
		void set_val(double inp);
		
		void set_val(double inp, double dinp);
		
		void set_val(double inp, double dinp1, double dinp2);
		
		void set_val(DiffDoub2& inp);
		
		void add(DiffDoub2& inp);
		
		void sub(DiffDoub2& inp);
		
		void mult(DiffDoub2& inp);
		
		void dvd(DiffDoub2& inp);

		void sqt();
		
		void sqr();
		
		void sn();
		
		void cs();
		
		void atn();
		
		void neg();
};

#endif