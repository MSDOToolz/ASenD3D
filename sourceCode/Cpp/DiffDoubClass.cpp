#include "DiffDoubClass.h"
#include <cmath>

using namespace std;

Doub::Doub() {
}

bool Doub::diffType() {
	return false;
}

void Doub::setVal(double inp) {
	val = inp;
	return;
}

void Doub::setVal(Doub& inp) {
	val = inp.val;
	return;
}

void Doub::add(Doub& inp) {
	val+= inp.val;
	return;
}

void Doub::sub(Doub& inp) {
	val-= inp.val;
	return;
}

void Doub::mult(Doub& inp) {
	val*= inp.val;
	return;
}

void Doub::dvd(Doub& inp) {
	val/= inp.val;
	return;
}

void Doub::sqt() {
	val = sqrt(val);
	return;
}

void Doub::sqr() {
	val*= val;
	return;
}

void Doub::sn() {
	val = sin(val);
	return;
}

void Doub::cs() {
	val = cos(val);
	return;
}

void Doub::atn() {
	val = atan(val);
	return;
}

void Doub::neg() {
	val = -val;
	return;
}



DiffDoub::DiffDoub() {
	return;
}

bool DiffDoub::diffType() {
	return true;
}

void DiffDoub::setVal(double inp) {
	val = inp;
	dval = 0.0;
	return;
}

void DiffDoub::setVal(double inp, double dinp) {
	val = inp;
	dval = dinp;
	return;
}

void DiffDoub::setVal(DiffDoub& inp) {
	val = inp.val;
	dval = inp.dval;
	return;
}

void DiffDoub::add(DiffDoub& inp) {
	val+= inp.val;
	dval+= inp.dval;
	return;
}

void DiffDoub::sub(DiffDoub& inp) {
	val-= inp.val;
	dval-= inp.dval;
	return;
}

void DiffDoub::mult(DiffDoub& inp) {
	tmp = val*inp.val;
	dval = dval*inp.val + val*inp.dval;
	val = tmp;
	return;
}

void DiffDoub::dvd(DiffDoub& inp) {
	tmp2 = 1.0/inp.val;
	tmp = val*tmp2;
	dval = tmp2*(dval - val*tmp2*inp.dval);
	val = tmp;
	return;
}

void DiffDoub::sqt() {
	val = sqrt(val);
	dval = 0.5*dval/val;
	return;
}

void Doub::sqr() {
	tmp = val*val;
	dval = 2.0*val*dval;
	val = tmp;
	return;
}

void DiffDoub::sn() {
	tmp = sin(val);
	dval = cos(val)*dval;
	val = tmp;
	return;
}

void DiffDoub::cs() {
	tmp = cos(val);
	dval = -sin(val)*dval;
	val = tmp;
	return;
}

void DiffDoub::atn() {
	tmp = atan(val);
	dval = dval/(1.0 + val*val);
	val = tmp;
	return;
}

void Doub::neg() {
	val = -val;
	dval = -dval;
	return;
}



Diff2Doub::Diff2Doub() {
}

bool Diff2Doub::diffType() {
	return true;
}

void Diff2Doub::setVal(double inp) {
	val = inp;
	dv1 = 0.0;
	dv2 = 0.0;
	dv12 = 0.0;
	return;
}

void Diff2Doub::setVal(double inp, double dinp) {
	val = inp;
	dv1 = dinp;
	dv2 = 0.0;
	dv12 = 0.0;
	return;
}

void Diff2Doub::setVal(Diff2Doub& inp) {
	val = inp.val;
	dv1 = inp.dv1;
	dv2 = inp.dv2;
	dv12 = inp.dv12;
	return;
}

void Diff2Doub::add(Diff2Doub& inp) {
	val+= inp.val;
	dv1+= inp.dv1;
	dv2+= inp.dv2;
	dv12+= inp.dv12;
	return;
}

void Diff2Doub::sub(Diff2Doub& inp) {
	val-= inp.val;
	dv1-= inp.dv1;
	dv2-= inp.dv2;
	dv12-= inp.dv12;
	return;
}

void Diff2Doub::mult(Diff2Doub& inp) {
	tmp = val*inp.val;
	tmp2 = dv1*inp.val + val*inp.dv1;
	tmp3 = dv2*inp.val + val*inp.dv2;
	dv12 = dv12*inp.val + dv1*inp.dv2 + dv2*inp.dv1 + val*inp.dv12;
	val = tmp;
	dv1 = tmp2;
	dv2 = tmp3;
	return;
}

void Diff2Doub::dvd(Diff2Doub& inp) {
	tmp = 1.0/inp.val; 
	tmp3 = tmp*(dv1 - val*tmp*inp.dv1); //dv1
	tmp4 = tmp*(dv2 - val*tmp*inp.dv2); //dv2
	dv12 = -tmp*inp.dv2*tmp3 + tmp*(dv12 - dv2*tmp*inp.dv1 + val*tmp*(tmp*inp.dv1*inp.dv2 - inp.dv12));
	val = val*tmp;
	dv1 = tmp3;
	dv2 = tmp4;
	return;
}

void Diff2Doub::sqt() {
	tmp = sqrt(val);
	tmp2 = 0.5*dv1/tmp;
	tmp3 = 0.5*dv2/tmp;
	dv12 = -0.25*dv1*dv2/(val*tmp) + 0.5*dv12/tmp;
	val = tmp;
	dv1 = tmp2;
	dv2 = tmp3;
	return;
}

void Doub::sqr() {
	tmp = val*val;
	tmp2 = 2.0*val*dv1;
	tmp3 = 2.0*val*dv2;
	dv12 = 2.0*(dv2*dv1 + val*dv12);
	val = tmp;
	dv1 = tmp2;
	dv2 = tmp3;
	return;
}

void Diff2Doub::sn() {
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

void Diff2Doub::cs() {
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

void Diff2Doub::atn() {
	tmp = 1.0/(1.0 + val*val);
	tmp2 = tmp*dv1;
	tmp3 = tmp*dv2;
	dv12 = tmp*(-tmp2*2.0*val*dv2 + dv12);
	val = atan(val);
	dv1 = tmp2;
	dv2 = tmp3;
	return;
}

void Doub::neg() {
	val = -val;
	dv1 = -dv1;
	dv2 = -dv2;
	dv12 = -dv12;
	return;
}

