#include "DiffDoubClass.h"
#include <cmath>

using namespace std;

 DiffDoub0::DiffDoub0() {
	val = 0.0;
	dval = 0.0;
}

bool DiffDoub0::diffType() {
	return false;
}

void DiffDoub0::setVal(double inp) {
	val = inp;
	return;
}

void DiffDoub0::setVal(double inp, double dinp) {
	val = inp;
	dval = dinp;
	return;
}

void DiffDoub0::setVal(DiffDoub0& inp) {
	val = inp.val;
	return;
}

void DiffDoub0::add(DiffDoub0& inp) {
	val+= inp.val;
	return;
}

void DiffDoub0::sub(DiffDoub0& inp) {
	val-= inp.val;
	return;
}

void DiffDoub0::mult(DiffDoub0& inp) {
	val*= inp.val;
	return;
}

void DiffDoub0::dvd(DiffDoub0& inp) {
	val/= inp.val;
	return;
}

void DiffDoub0::sqt() {
	val = sqrt(val);
	return;
}

void DiffDoub0::sqr() {
	val*= val;
	return;
}

void DiffDoub0::sn() {
	val = sin(val);
	return;
}

void DiffDoub0::cs() {
	val = cos(val);
	return;
}

void DiffDoub0::asn() {
	val = asin(val);
	return;
}

void DiffDoub0::atn() {
	val = atan(val);
	return;
}

void DiffDoub0::neg() {
	val = -val;
	return;
}



 DiffDoub1::DiffDoub1() {
	val = 0.0;
	dval = 0.0;
	tmp = 0.0;
	tmp2 = 0.0;
	return;
}

bool DiffDoub1::diffType() {
	return true;
}

void DiffDoub1::setVal(double inp) {
	val = inp;
	dval = 0.0;
	return;
}

void DiffDoub1::setVal(double inp, double dinp) {
	val = inp;
	dval = dinp;
	return;
}

void DiffDoub1::setVal(DiffDoub0& inp) {
	val = inp.val;
	dval = 0.0;
	return;
}

void DiffDoub1::setVal(DiffDoub1& inp) {
	val = inp.val;
	dval = inp.dval;
	return;
}

void DiffDoub1::add(DiffDoub1& inp) {
	val+= inp.val;
	dval+= inp.dval;
	return;
}

void DiffDoub1::sub(DiffDoub1& inp) {
	val-= inp.val;
	dval-= inp.dval;
	return;
}

void DiffDoub1::mult(DiffDoub1& inp) {
	tmp = val*inp.val;
	dval = dval*inp.val + val*inp.dval;
	val = tmp;
	return;
}

void DiffDoub1::dvd(DiffDoub1& inp) {
	tmp2 = 1.0/inp.val;
	tmp = val*tmp2;
	dval = tmp2*(dval - val*tmp2*inp.dval);
	val = tmp;
	return;
}

void DiffDoub1::sqt() {
	val = sqrt(val);
	dval = 0.5*dval/val;
	return;
}

void DiffDoub1::sqr() {
	tmp = val*val;
	dval = 2.0*val*dval;
	val = tmp;
	return;
}

void DiffDoub1::sn() {
	tmp = sin(val);
	dval = cos(val)*dval;
	val = tmp;
	return;
}

void DiffDoub1::cs() {
	tmp = cos(val);
	dval = -sin(val)*dval;
	val = tmp;
	return;
}

void DiffDoub1::asn() {
	tmp = asin(val);
	dval = dval/cos(tmp);
	val = tmp;
	return;
}

void DiffDoub1::atn() {
	tmp = atan(val);
	dval = dval/(1.0 + val*val);
	val = tmp;
	return;
}

void DiffDoub1::neg() {
	val = -val;
	dval = -dval;
	return;
}



 DiffDoub2::DiffDoub2() {
	val = 0.0;
	dv1 = 0.0;
	dv2 = 0.0;
	dv12 = 0.0;
	tmp = 0.0;
	tmp2 = 0.0;
	tmp3 = 0.0;
	tmp4 = 0.0;
	return;
}

bool DiffDoub2::diffType() {
	return true;
}

void DiffDoub2::setVal(double inp) {
	val = inp;
	dv1 = 0.0;
	dv2 = 0.0;
	dv12 = 0.0;
	return;
}

void DiffDoub2::setVal(double inp, double dinp) {
	val = inp;
	dv1 = dinp;
	dv2 = 0.0;
	dv12 = 0.0;
	return;
}

void DiffDoub2::setVal(double inp, double dinp1, double dinp2) {
	val = inp;
	dv1 = dinp1;
	dv2 = dinp2;
	dv12 = 0.0;
	return;
}

void DiffDoub2::setVal(DiffDoub2& inp) {
	val = inp.val;
	dv1 = inp.dv1;
	dv2 = inp.dv2;
	dv12 = inp.dv12;
	return;
}

void DiffDoub2::add(DiffDoub2& inp) {
	val+= inp.val;
	dv1+= inp.dv1;
	dv2+= inp.dv2;
	dv12+= inp.dv12;
	return;
}

void DiffDoub2::sub(DiffDoub2& inp) {
	val-= inp.val;
	dv1-= inp.dv1;
	dv2-= inp.dv2;
	dv12-= inp.dv12;
	return;
}

void DiffDoub2::mult(DiffDoub2& inp) {
	tmp = val*inp.val;
	tmp2 = dv1*inp.val + val*inp.dv1;
	tmp3 = dv2*inp.val + val*inp.dv2;
	dv12 = dv12*inp.val + dv1*inp.dv2 + dv2*inp.dv1 + val*inp.dv12;
	val = tmp;
	dv1 = tmp2;
	dv2 = tmp3;
	return;
}

void DiffDoub2::dvd(DiffDoub2& inp) {
	tmp = 1.0/inp.val; 
	tmp3 = tmp*(dv1 - val*tmp*inp.dv1); //dv1
	tmp4 = tmp*(dv2 - val*tmp*inp.dv2); //dv2
	dv12 = -tmp*inp.dv2*tmp3 + tmp*(dv12 - dv2*tmp*inp.dv1 + val*tmp*(tmp*inp.dv1*inp.dv2 - inp.dv12));
	val = val*tmp;
	dv1 = tmp3;
	dv2 = tmp4;
	return;
}

void DiffDoub2::sqt() {
	tmp = sqrt(val);
	tmp2 = 0.5*dv1/tmp;
	tmp3 = 0.5*dv2/tmp;
	dv12 = -0.25*dv1*dv2/(val*tmp) + 0.5*dv12/tmp;
	val = tmp;
	dv1 = tmp2;
	dv2 = tmp3;
	return;
}

void DiffDoub2::sqr() {
	tmp = val*val;
	tmp2 = 2.0*val*dv1;
	tmp3 = 2.0*val*dv2;
	dv12 = 2.0*(dv2*dv1 + val*dv12);
	val = tmp;
	dv1 = tmp2;
	dv2 = tmp3;
	return;
}

void DiffDoub2::sn() {
	tmp = sin(val);
	tmp2 = cos(val);
	tmp3 = tmp2*dv1;
	tmp4 = tmp2*dv2;
	dv12 = -tmp*dv1*dv2 + tmp2*dv12;
	val = tmp;
	dv1 = tmp3;
	dv2 = tmp4;
	return;
}

void DiffDoub2::cs() {
	tmp = cos(val);
	tmp2 = sin(val);
	tmp3 = -tmp2*dv1;
	tmp4 = -tmp2*dv2;
	dv12 = -tmp*dv1*dv2 - tmp2*dv12;
	val = tmp;
	dv1 = tmp3;
	dv2 = tmp4;
	return;
}

void DiffDoub2::atn() {
	tmp = 1.0/(1.0 + val*val);
	tmp2 = tmp*dv1;
	tmp3 = tmp*dv2;
	dv12 = tmp*(-tmp2*2.0*val*dv2 + dv12);
	val = atan(val);
	dv1 = tmp2;
	dv2 = tmp3;
	return;
}

void DiffDoub2::neg() {
	val = -val;
	dv1 = -dv1;
	dv2 = -dv2;
	dv12 = -dv12;
	return;
}

