use crate::diff_doub::*;
use crate::fmath::*;

impl DiffDoub0 {
    pub fn diff_type(&mut self) -> bool {
        return  false;
    }

    pub fn set_val(&mut self, inp : f64) {
        self.val = inp;
        return;
    }

    pub fn set_val_2(&mut self, inp : f64, dinp : f64) {
        self.val = inp;
        self.dval = dinp;
        return;
    }

    pub fn set_val_dfd0(&mut self, inp : & DiffDoub0) {
        self.val = inp.val;
        return;
    }

    pub fn add(&mut self, inp : & DiffDoub0) {
        self.val +=  inp.val;
        return;
    }

    pub fn sub(&mut self, inp : & DiffDoub0) {
        self.val -=  inp.val;
        return;
    }

    pub fn mult(&mut self, inp : & DiffDoub0) {
        self.val *=  inp.val;
        return;
    }

    pub fn dvd(&mut self, inp : & DiffDoub0) {
        self.val /=  inp.val;
        return;
    }

    pub fn sqt(&mut self) {
        self.val = sqrt(self.val);
        return;
    }

    pub fn sqr(&mut self) {
        self.val *=  self.val;
        return;
    }

    pub fn sn(&mut self) {
        self.val = sin(self.val);
        return;
    }

    pub fn cs(&mut self) {
        self.val = cos(self.val);
        return;
    }

    pub fn asn(&mut self) {
        self.val = asin(self.val);
        return;
    }

    pub fn atn(&mut self) {
        self.val = atan(self.val);
        return;
    }

    pub fn neg(&mut self) {
        self.val = -self.val;
        return;
    }

}

impl DiffDoub1 {
    pub fn diff_type(&mut self) -> bool {
        return  true;
    }

    pub fn set_val(&mut self, inp : f64) {
        self.val = inp;
        self.dval = 0.0;
        return;
    }

    pub fn set_val_2(&mut self, inp : f64, dinp : f64) {
        self.val = inp;
        self.dval = dinp;
        return;
    }

    pub fn set_val_dfd0(&mut self, inp : & DiffDoub0) {
        self.val = inp.val;
        self.dval = 0.0;
        return;
    }

    pub fn set_val_dfd1(&mut self, inp : & DiffDoub1) {
        self.val = inp.val;
        self.dval = inp.dval;
        return;
    }

    pub fn add(&mut self, inp : & DiffDoub1) {
        self.val +=  inp.val;
        self.dval +=  inp.dval;
        return;
    }

    pub fn sub(&mut self, inp : & DiffDoub1) {
        self.val -=  inp.val;
        self.dval -=  inp.dval;
        return;
    }

    pub fn mult(&mut self, inp : & DiffDoub1) {
        self.tmp = self.val*inp.val;
        self.dval = self.dval*inp.val + self.val*inp.dval;
        self.val = self.tmp;
        return;
    }

    pub fn dvd(&mut self, inp : & DiffDoub1) {
        self.tmp2 = 1.0/inp.val;
        self.tmp = self.val*self.tmp2;
        self.dval = self.tmp2*(self.dval - self.val*self.tmp2*inp.dval);
        self.val = self.tmp;
        return;
    }

    pub fn sqt(&mut self) {
        self.val = sqrt(self.val);
        self.dval = 0.5*self.dval/self.val;
        return;
    }

    pub fn sqr(&mut self) {
        self.tmp = self.val*self.val;
        self.dval = 2.0*self.val*self.dval;
        self.val = self.tmp;
        return;
    }

    pub fn sn(&mut self) {
        self.tmp = sin(self.val);
        self.dval = cos(self.val)*self.dval;
        self.val = self.tmp;
        return;
    }

    pub fn cs(&mut self) {
        self.tmp = cos(self.val);
        self.dval = -sin(self.val)*self.dval;
        self.val = self.tmp;
        return;
    }

    pub fn asn(&mut self) {
        self.tmp = asin(self.val);
        self.dval = self.dval/cos(self.tmp);
        self.val = self.tmp;
        return;
    }

    pub fn atn(&mut self) {
        self.tmp = atan(self.val);
        self.dval = self.dval/(1.0 + self.val*self.val);
        self.val = self.tmp;
        return;
    }

    pub fn neg(&mut self) {
        self.val = -self.val;
        self.dval = -self.dval;
        return;
    }

}

impl DiffDoub2 {
    pub fn diff_type(&mut self) -> bool {
        return  true;
    }

    pub fn set_val(&mut self, inp : f64) {
        self.val = inp;
        self.dv1 = 0.0;
        self.dv2 = 0.0;
        self.dv12 = 0.0;
        return;
    }

    pub fn set_val_2(&mut self, inp : f64, dinp : f64) {
        self.val = inp;
        self.dv1 = dinp;
        self.dv2 = 0.0;
        self.dv12 = 0.0;
        return;
    }

    pub fn set_val_3(&mut self, inp : f64, dinp1 : f64, dinp2 : f64) {
        self.val = inp;
        self.dv1 = dinp1;
        self.dv2 = dinp2;
        self.dv12 = 0.0;
        return;
    }

    pub fn set_val_dfd2(&mut self, inp : & DiffDoub2) {
        self.val = inp.val;
        self.dv1 = inp.dv1;
        self.dv2 = inp.dv2;
        self.dv12 = inp.dv12;
        return;
    }

    pub fn add(&mut self, inp : & DiffDoub2) {
        self.val +=  inp.val;
        self.dv1 +=  inp.dv1;
        self.dv2 +=  inp.dv2;
        self.dv12 +=  inp.dv12;
        return;
    }

    pub fn sub(&mut self, inp : & DiffDoub2) {
        self.val -=  inp.val;
        self.dv1 -=  inp.dv1;
        self.dv2 -=  inp.dv2;
        self.dv12 -=  inp.dv12;
        return;
    }

    pub fn mult(&mut self, inp : & DiffDoub2) {
        self.tmp = self.val*inp.val;
        self.tmp2 = self.dv1*inp.val + self.val*inp.dv1;
        self.tmp3 = self.dv2*inp.val + self.val*inp.dv2;
        self.dv12 = self.dv12*inp.val + self.dv1*inp.dv2 + self.dv2*inp.dv1 + self.val*inp.dv12;
        self.val = self.tmp;
        self.dv1 = self.tmp2;
        self.dv2 = self.tmp3;
        return;
    }

    pub fn dvd(&mut self, inp : & DiffDoub2) {
        self.tmp = 1.0/inp.val;
        self.tmp3 = self.tmp*(self.dv1 - self.val*self.tmp*inp.dv1);//self.dv1
        self.tmp4 = self.tmp*(self.dv2 - self.val*self.tmp*inp.dv2);//self.dv2
        self.dv12 = -self.tmp*inp.dv2*self.tmp3 + self.tmp*(self.dv12 - self.dv2*self.tmp*inp.dv1 + self.val*self.tmp*(self.tmp*inp.dv1*inp.dv2 - inp.dv12));
        self.val = self.val*self.tmp;
        self.dv1 = self.tmp3;
        self.dv2 = self.tmp4;
        return;
    }

    pub fn sqt(&mut self) {
        self.tmp = sqrt(self.val);
        self.tmp2 = 0.5*self.dv1/self.tmp;
        self.tmp3 = 0.5*self.dv2/self.tmp;
        self.dv12 = -0.25*self.dv1*self.dv2/(self.val*self.tmp) + 0.5*self.dv12/self.tmp;
        self.val = self.tmp;
        self.dv1 = self.tmp2;
        self.dv2 = self.tmp3;
        return;
    }

    pub fn sqr(&mut self) {
        self.tmp = self.val*self.val;
        self.tmp2 = 2.0*self.val*self.dv1;
        self.tmp3 = 2.0*self.val*self.dv2;
        self.dv12 = 2.0*(self.dv2*self.dv1 + self.val*self.dv12);
        self.val = self.tmp;
        self.dv1 = self.tmp2;
        self.dv2 = self.tmp3;
        return;
    }

    pub fn sn(&mut self) {
        self.tmp = sin(self.val);
        self.tmp2 = cos(self.val);
        self.tmp3 = self.tmp2*self.dv1;
        self.tmp4 = self.tmp2*self.dv2;
        self.dv12 = -self.tmp*self.dv1*self.dv2 + self.tmp2*self.dv12;
        self.val = self.tmp;
        self.dv1 = self.tmp3;
        self.dv2 = self.tmp4;
        return;
    }

    pub fn cs(&mut self) {
        self.tmp = cos(self.val);
        self.tmp2 = sin(self.val);
        self.tmp3 = -self.tmp2*self.dv1;
        self.tmp4 = -self.tmp2*self.dv2;
        self.dv12 = -self.tmp*self.dv1*self.dv2 - self.tmp2*self.dv12;
        self.val = self.tmp;
        self.dv1 = self.tmp3;
        self.dv2 = self.tmp4;
        return;
    }

    pub fn atn(&mut self) {
        self.tmp = 1.0/(1.0 + self.val*self.val);
        self.tmp2 = self.tmp*self.dv1;
        self.tmp3 = self.tmp*self.dv2;
        self.dv12 = self.tmp*(-self.tmp2*2.0*self.val*self.dv2 + self.dv12);
        self.val = atan(self.val);
        self.dv1 = self.tmp2;
        self.dv2 = self.tmp3;
        return;
    }

    pub fn neg(&mut self) {
        self.val = -self.val;
        self.dv1 = -self.dv1;
        self.dv2 = -self.dv2;
        self.dv12 = -self.dv12;
        return;
    }

}


