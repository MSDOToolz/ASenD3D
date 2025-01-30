#include <string>
#include "ElementClass.h"
#include "constants.h"
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "DesignVariableClass.h"
#include "NodeClass.h"
#include "SectionClass.h"
#include "matrixFunctions.h"

using namespace std;

//dup1

void Element::get_gen_prop_dfd0(DiffDoub0& prop, string prop_key, vector<DesignVariable>& dv_ar) {
	int d_vind;
	DiffDoub0 dv_val;
	DiffDoub0 tmp;

	for (auto& dv : design_vars) {
		d_vind = dv.int_dat;
		DesignVariable& this_dv = dv_ar[d_vind];
		if (this_dv.category == prop_key) {
			this_dv.get_value_dfd0(dv_val);
			tmp.set_val(dv.doub_dat);
			dv_val.mult(tmp);
			prop.add(dv_val);
		}
	}

	return;
}

void Element::get_layer_thk_z_dfd0(vector<DiffDoub0>& lay_thk, vector<DiffDoub0>& lay_z, DiffDoub0& z_offset, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	//z_offset = 1: upper z surface is reference plane
	//z_offset = -1: lower z surface is reference plane
	int layi = 0;
	DiffDoub0 dv_val;
	DiffDoub0 coef;
	DiffDoub0 tot_thk;
	DiffDoub0 tmp;
	DiffDoub0 z_crd;
	DiffDoub0 z_next;
	DiffDoub0 z_mid;

	tot_thk.set_val(0.0);
	Section& this_sec = sec_ar[sect_ptr];
	for (auto& lay : this_sec.layers) {
		lay_thk[layi].set_val(lay.thickness);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "thickness" && this_dv.layer == layi) {
				this_dv.get_value_dfd0(dv_val);
				coef.set_val(dv.doub_dat);
				dv_val.mult(coef);
				lay_thk[layi].add(dv_val);
			}
		}
		tot_thk.add(lay_thk[layi]);
		layi++;
	}

	z_offset.set_val(this_sec.z_offset);
	for (auto& dv : design_vars) {
		DesignVariable& this_dv = dv_ar[dv.int_dat];
		if (this_dv.category == "zOffset") {
			this_dv.get_value_dfd0(dv_val);
			coef.set_val(dv.doub_dat);
			dv_val.mult(coef);
			z_offset.add(dv_val);
		}
	}

	tmp.set_val(1.0);
	tmp.add(z_offset);
	z_crd.set_val(-0.5);
	z_crd.mult(tot_thk);
	z_crd.mult(tmp);

	layi = 0;
	for (auto& lay : this_sec.layers) {
		z_next.set_val_dfd0(z_crd);
		z_next.add(lay_thk[layi]);
		tmp.set_val_dfd0(z_crd);
		tmp.add(z_next);
		z_mid.set_val(0.5);
		z_mid.mult(tmp);
		lay_z[layi].set_val_dfd0(z_mid);
		layi++;
		z_crd.set_val_dfd0(z_next);
	}

	return;
}

void Element::get_layer_q_dfd0(vector<DiffDoub0>& lay_q, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int layi;
	int dv_ind;
	int dv_comp;
	DiffDoub0 modulus_dv[3];
	DiffDoub0 poisson_dv[3];
	DiffDoub0 shear_mod_dv[3];
	DiffDoub0 smat[9];
	DiffDoub0 qmat[9];
	DiffDoub0 dv_val;
	DiffDoub0 coef;
	DiffDoub0 x_vec[3];
	DiffDoub0 b_vec[3];
	string dv_cat;

	i2 = 0;
	layi = 0;
	Section& this_sec = sec_ar[sect_ptr];
	for (auto& lay : this_sec.layers) {
		Material& this_mat = mat_ar[lay.mat_ptr];
		for (i1 = 0; i1 < 3; i1++) {
			modulus_dv[i1].set_val(this_mat.modulus[i1]);
			poisson_dv[i1].set_val(this_mat.poisson_ratio[i1]);
			shear_mod_dv[i1].set_val(this_mat.shear_mod[i1]);
		}
		for (auto& dv : design_vars) {
			dv_ind = dv.int_dat;
			DesignVariable& this_dv = dv_ar[dv_ind];
			if (this_dv.layer == layi) {
				dv_cat = this_dv.category;
				if (dv_cat == "modulus") {
					dv_comp = this_dv.component;
					coef.set_val(dv.doub_dat);
					this_dv.get_value_dfd0(dv_val);
					dv_val.mult(coef);
					modulus_dv[dv_comp - 1].add(dv_val);
				}
				else if (dv_cat == "poissonRatio") {
					dv_comp = this_dv.component;
					coef.set_val(dv.doub_dat);
					this_dv.get_value_dfd0(dv_val);
					dv_val.mult(coef);
					poisson_dv[dv_comp - 1].add(dv_val);
				}
				else if (dv_cat == "shearModulus") {
					dv_comp = this_dv.component;
					coef.set_val(dv.doub_dat);
					this_dv.get_value_dfd0(dv_val);
					dv_val.mult(coef);
					shear_mod_dv[dv_comp - 1].add(dv_val);
				}
			}
		}

		for (i1 = 0; i1 < 9; i1++) {
			smat[i1].set_val(0.0);
		}
		smat[0].set_val(1.0);
		smat[0].dvd(modulus_dv[0]);
		smat[1].set_val_dfd0(poisson_dv[0]);
		smat[1].neg();
		smat[1].dvd(modulus_dv[0]);
		smat[3].set_val_dfd0(smat[1]);
		smat[4].set_val(1.0);
		smat[4].dvd(modulus_dv[1]);
		smat[8].set_val(1.0);
		smat[8].dvd(shear_mod_dv[0]);
		get_det_inv_ar_dfd0(coef, qmat, smat, 3, 0, x_vec, b_vec);

		for (i1 = 0; i1 < 9; i1++) {
			lay_q[i2].set_val_dfd0(qmat[i1]);
			i2++;
		}

		layi++;
	}

	return;
}

void Element::get_layer_d_dfd0(vector<DiffDoub0>& lay_d, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int layi;
	int dv_ind;
	int dv_comp;
	DiffDoub0 damp_mat_dv[36];
	DiffDoub0 dv_val;
	DiffDoub0 coef;
	string dv_cat;

	i2 = 0;
	layi = 0;
	Section& this_sec = sec_ar[sect_ptr];
	for (auto& lay : this_sec.layers) {
		Material& this_mat = mat_ar[lay.mat_ptr];
		for (i1 = 0; i1 < 36; i1++) {
			damp_mat_dv[i1].set_val(this_mat.damping[i1]);
		}
		for (auto& dv : design_vars) {
			dv_ind = dv.int_dat;
			DesignVariable& this_dv = dv_ar[dv_ind];
			if (this_dv.layer == layi) {
				dv_cat = this_dv.category;
				if (dv_cat == "dampingMat") {
					dv_comp = this_dv.component;
					coef.set_val(dv.doub_dat);
					this_dv.get_value_dfd0(dv_val);
					dv_val.mult(coef);
					damp_mat_dv[dv_comp].add(dv_val);
				}
			}
		}

		for (i3 = 1; i3 < 6; i3++) {
			i5 = 6 * i3; // lower tri term
			i6 = i3; // upper tri term
			for (i4 = 0; i4 < i3; i4++) {
				damp_mat_dv[i5].set_val_dfd0(damp_mat_dv[i6]);
				i5++;
				i6 += 6;
			}
		}

		lay_d[i2].set_val_dfd0(damp_mat_dv[0]);
		lay_d[i2 + 1].set_val_dfd0(damp_mat_dv[1]);
		lay_d[i2 + 2].set_val_dfd0(damp_mat_dv[3]);
		lay_d[i2 + 3].set_val_dfd0(damp_mat_dv[6]);
		lay_d[i2 + 4].set_val_dfd0(damp_mat_dv[7]);
		lay_d[i2 + 5].set_val_dfd0(damp_mat_dv[9]);
		lay_d[i2 + 6].set_val_dfd0(damp_mat_dv[18]);
		lay_d[i2 + 7].set_val_dfd0(damp_mat_dv[19]);
		lay_d[i2 + 8].set_val_dfd0(damp_mat_dv[21]);

		layi++;
		i2 += 9;
	}
	return;
}

void Element::get_layer_angle_dfd0(vector<DiffDoub0>& lay_ang, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	int i2;
	int layi;
	int dv_ind;
	DiffDoub0 angle;
	DiffDoub0 dv_val;
	DiffDoub0 coef;
	string dv_cat;

	i2 = 0;
	layi = 0;
	Section& this_sec = sec_ar[sect_ptr];
	for (auto& lay : this_sec.layers) {
		angle.set_val(lay.angle);
		for (auto& dv : design_vars) {
			dv_ind = dv.int_dat;
			DesignVariable& this_dv = dv_ar[dv_ind];
			if (this_dv.layer == layi) {
				dv_cat = this_dv.category;
				if (dv_cat == "angle") {
					coef.set_val(dv.doub_dat);
					this_dv.get_value_dfd0(dv_val);
					dv_val.mult(coef);
					angle.add(dv_val);
				}
			}
		}
		lay_ang[layi].set_val_dfd0(angle);

		layi++;
	}

	return;
}

void Element::get_layer_th_exp_dfd0(vector<DiffDoub0>& lay_th_exp, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int layi;
	DiffDoub0 t_exp_dv[6];
	int dv_comp;
	DiffDoub0 tmp;
	DiffDoub0 dv_val;
	
	layi = 0;
	i2 = 0;
	Section& this_sec = sec_ar[sect_ptr];
	for (auto& lay : this_sec.layers) {
		Material& this_mat = mat_ar[lay.mat_ptr];
		for (i1 = 0; i1 < 6; i1++) {
			t_exp_dv[i1].set_val(this_mat.expansion[i1]);
		}
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "thermalExp" && this_dv.layer == layi) {
				dv_comp = this_dv.component - 1;
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(tmp);
				t_exp_dv[dv_comp].add(dv_val);
			}
		}
		lay_th_exp[i2].set_val_dfd0(t_exp_dv[0]);
		i2++;
		lay_th_exp[i2].set_val_dfd0(t_exp_dv[1]);
		i2++;
		lay_th_exp[i2].set_val_dfd0(t_exp_dv[3]);
		i2++;
		layi++;
	}
	return;
}

void Element::get_layer_einit_dfd0(vector<DiffDoub0>& lay_einit, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int layi;
	DiffDoub0 e0_dv[6];
	int dv_comp;
	DiffDoub0 tmp;
	DiffDoub0 dv_val;

	layi = 0;
	i2 = 0;
	Section& this_sec = sec_ar[sect_ptr];
	for (auto& lay : this_sec.layers) {
		for (i1 = 0; i1 < 6; i1++) {
			e0_dv[i1].set_val(0.0);
		}
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "initialStrain" && this_dv.layer == layi) {
				dv_comp = this_dv.component - 1;
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(tmp);
				e0_dv[dv_comp].add(dv_val);
			}
		}
		lay_einit[i2].set_val_dfd0(e0_dv[0]);
		i2++;
		lay_einit[i2].set_val_dfd0(e0_dv[1]);
		i2++;
		lay_einit[i2].set_val_dfd0(e0_dv[3]);
		i2++;
		layi++;
	}
	return;
}

void Element::get_layer_den_dfd0(vector<DiffDoub0>& layer_den, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int layi;
	double mat_den;
	DiffDoub0 den_dv;
	DiffDoub0 tmp;
	DiffDoub0 dv_val;

	layi = 0;
	Section& this_sec = sec_ar[sect_ptr];
	for (auto& lay : this_sec.layers) {
		mat_den = mat_ar[lay.mat_ptr].density;
		den_dv.set_val(mat_den);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "density" && this_dv.layer == layi) {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(tmp);
				den_dv.add(dv_val);
			}
		}
		layer_den[layi].set_val_dfd0(den_dv);
		layi++;
	}

	return;
}

void Element::get_layer_cond_dfd0(vector<DiffDoub0>& lay_cond, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int layi;
	DiffDoub0 cond_dv[6];
	DiffDoub0 tmp;
	DiffDoub0 dv_val;
	int dv_comp;

	i1 = 0;
	layi = 0;
	Section& this_sec = sec_ar[sect_ptr];
	for (auto& lay : this_sec.layers) {
		Material& this_mat = mat_ar[lay.mat_ptr];
		for (i2 = 0; i2 < 6; i2++) {
			cond_dv[i2].set_val(this_mat.conductivity[i2]);
		}
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "thermalCond" && this_dv.layer == layi) {
				dv_comp = this_dv.component - 1;
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(tmp);
				cond_dv[dv_comp].add(dv_val);
			}
		}
		lay_cond[i1].set_val_dfd0(cond_dv[0]);
		lay_cond[i1 + 1].set_val_dfd0(cond_dv[3]);
		lay_cond[i1 + 2].set_val_dfd0(cond_dv[4]);
		lay_cond[i1 + 3].set_val_dfd0(cond_dv[3]);
		lay_cond[i1 + 4].set_val_dfd0(cond_dv[1]);
		lay_cond[i1 + 5].set_val_dfd0(cond_dv[5]);
		lay_cond[i1 + 6].set_val_dfd0(cond_dv[4]);
		lay_cond[i1 + 7].set_val_dfd0(cond_dv[5]);
		lay_cond[i1 + 8].set_val_dfd0(cond_dv[2]);
		i1 += 9;
		layi++;
	}

	return;
}

void Element::get_layer_spec_heat_dfd0(vector<DiffDoub0>& lay_sh, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int layi;
	double mat_sh;
	DiffDoub0 sh_dv;
	DiffDoub0 tmp;
	DiffDoub0 dv_val;

	layi = 0;
	Section& this_sec = sec_ar[sect_ptr];
	for (auto& lay : this_sec.layers) {
		mat_sh = mat_ar[lay.mat_ptr].spec_heat;
		sh_dv.set_val(mat_sh);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "specHeat" && this_dv.layer == layi) {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(tmp);
				sh_dv.add(dv_val);
			}
		}
		lay_sh[layi].set_val_dfd0(sh_dv);
		layi++;
	}

	return;
}

void Element::transform_strain_dfd0(DiffDoub0 stn_new[], DiffDoub0 stn_orig[], DiffDoub0& angle) {
	DiffDoub0 angle_rad;
	DiffDoub0 a11;
	DiffDoub0 a12;
	DiffDoub0 a21;
	DiffDoub0 a22;
	DiffDoub0 tmp;
	DiffDoub0 te[9];


	angle_rad.set_val(r_pio180);
	angle_rad.mult(angle);
	a11.set_val_dfd0(angle_rad);
	a11.cs();
	a21.set_val_dfd0(angle_rad);
	a21.sn();
	a12.set_val_dfd0(a21);
	a12.neg();
	a22.set_val_dfd0(a11);

	te[0].set_val_dfd0(a11);
	te[0].sqr();
	te[1].set_val_dfd0(a12);
	te[1].sqr();
	te[2].set_val_dfd0(a11);
	te[2].mult(a12);
	te[3].set_val_dfd0(a21);
	te[3].sqr();
	te[4].set_val_dfd0(a22);
	te[4].sqr();
	te[5].set_val_dfd0(a22);
	te[5].mult(a21);
	te[6].set_val(2.0);
	te[6].mult(a11);
	te[6].mult(a21);
	te[7].set_val(2.0);
	te[7].mult(a12);
	te[7].mult(a22);
	te[8].set_val_dfd0(a11);
	te[8].mult(a22);
	tmp.set_val_dfd0(a12);
	tmp.mult(a21);
	te[8].add(tmp);

	mat_mul_ar_dfd0(stn_new, te, stn_orig, 3, 3, 1);

	return;
}

void Element::transform_q_dfd0(DiffDoub0 q_new[], DiffDoub0 q_orig[], DiffDoub0& angle) {
	int i1;
	DiffDoub0 angle_rad;
	DiffDoub0 a11;
	DiffDoub0 a12;
	DiffDoub0 a21;
	DiffDoub0 a22;
	DiffDoub0 tmp;
	DiffDoub0 coef;
	DiffDoub0 ts[9];
	DiffDoub0 te[9];
	DiffDoub0 te_inv[9];
	DiffDoub0 x_vec[3];
	DiffDoub0 b_vec[3];


	angle_rad.set_val(r_pio180);
	angle_rad.mult(angle);
	a11.set_val_dfd0(angle_rad);
	a11.cs();
	a21.set_val_dfd0(angle_rad);
	a21.sn();
	a12.set_val_dfd0(a21);
	a12.neg();
	a22.set_val_dfd0(a11);

	ts[0].set_val_dfd0(a11);
	ts[0].sqr();
	ts[1].set_val_dfd0(a12);
	ts[1].sqr();
	ts[2].set_val(2.0);
	ts[2].mult(a11);
	ts[2].mult(a12);
	ts[3].set_val_dfd0(a21);
	ts[3].sqr();
	ts[4].set_val_dfd0(a22);
	ts[4].sqr();
	ts[5].set_val(2.0);
	ts[5].mult(a22);
	ts[5].mult(a21);
	ts[6].set_val_dfd0(a11);
	ts[6].mult(a21);
	ts[7].set_val_dfd0(a12);
	ts[7].mult(a22);
	ts[8].set_val_dfd0(a11);
	ts[8].mult(a22);
	tmp.set_val_dfd0(a12);
	tmp.mult(a21);
	ts[8].add(tmp);

	for (i1 = 0; i1 < 9; i1++) {
		te[i1].set_val_dfd0(ts[i1]);
	}
	tmp.set_val(0.5);
	te[2].mult(tmp);
	te[5].mult(tmp);
	tmp.set_val(2.0);
	te[6].mult(tmp);
	te[7].mult(tmp);

	get_det_inv_ar_dfd0(coef, te_inv, te, 3, 0, x_vec, b_vec);

	mat_mul_ar_dfd0(te, q_orig, te_inv, 3, 3, 3);
	mat_mul_ar_dfd0(q_new, ts, te, 3, 3, 3);

	return;
}

void Element::get_solid_stiff_dfd0(vector<DiffDoub0>& cmat, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i4;
	int i5;
	int dv_ind;
	int dv_comp;
	DiffDoub0 modulus_dv[3];
	DiffDoub0 poisson_dv[3];
	DiffDoub0 shear_mod_dv[3];
	DiffDoub0 smat[36];
	DiffDoub0 ctmp[36];
	DiffDoub0 dv_val;
	DiffDoub0 coef;
	DiffDoub0 x_vec[6];
	DiffDoub0 b_vec[6];
	string dv_cat;


	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	if (this_mat.stiffness[0] > 0.0) {
		for (i1 = 0; i1 < 36; i1++) {
			cmat[i1].set_val(this_mat.stiffness[i1]);
		}
		for (auto& dv : design_vars) {
			dv_ind = dv.int_dat;
			DesignVariable& this_dv = dv_ar[dv_ind];
			if (this_dv.category == "stiffnessMat") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				cmat[dv_comp].add(dv_val);
			}
		}
		for (i1 = 1; i1 < 6; i1++) {
			i4 = 6 * i1;
			i5 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				cmat[i4].set_val_dfd0(cmat[i5]);
				i4++;
				i5 += 6;
			}
		}
	}
	else {
		for (i1 = 0; i1 < 3; i1++) {
			modulus_dv[i1].set_val(this_mat.modulus[i1]);
			poisson_dv[i1].set_val(this_mat.poisson_ratio[i1]);
			shear_mod_dv[i1].set_val(this_mat.shear_mod[i1]);
		}
		for (auto& dv : design_vars) {
			dv_ind = dv.int_dat;
			DesignVariable& this_dv = dv_ar[dv_ind];
			dv_cat = this_dv.category;
			if (dv_cat == "modulus") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				modulus_dv[dv_comp - 1].add(dv_val);
			}
			else if (dv_cat == "poissonRatio") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				poisson_dv[dv_comp - 1].add(dv_val);
			}
			else if (dv_cat == "shearModulus") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				shear_mod_dv[dv_comp - 1].add(dv_val);
			}
		}
		for (i1 = 0; i1 < 36; i1++) {
			smat[i1].set_val(0.0);
		}

		smat[0].set_val(1.0);
		smat[0].dvd(modulus_dv[0]);
		smat[1].set_val_dfd0(poisson_dv[0]);
		smat[1].neg();
		smat[1].dvd(modulus_dv[0]);
		smat[2].set_val_dfd0(poisson_dv[1]);
		smat[2].neg();
		smat[2].dvd(modulus_dv[0]);
		smat[6].set_val_dfd0(smat[1]);
		smat[7].set_val(1.0);
		smat[7].dvd(modulus_dv[1]);
		smat[8].set_val_dfd0(poisson_dv[2]);
		smat[8].neg();
		smat[8].dvd(modulus_dv[1]);
		smat[12].set_val_dfd0(smat[2]);
		smat[13].set_val_dfd0(smat[8]);
		smat[14].set_val(1.0);
		smat[14].dvd(modulus_dv[2]);
		smat[21].set_val(1.0);
		smat[21].dvd(shear_mod_dv[0]);
		smat[28].set_val(1.0);
		smat[28].dvd(shear_mod_dv[1]);
		smat[35].set_val(1.0);
		smat[35].dvd(shear_mod_dv[2]);

		get_det_inv_ar_dfd0(coef, ctmp, smat, 6, 0, x_vec, b_vec);
		ar_to_vec_dfd0(ctmp, cmat, 0, 36);
	}

	return;
}

void Element::get_abd_dfd0(vector<DiffDoub0>& cmat, vector<DiffDoub0>& lay_thk, vector<DiffDoub0>& lay_z, vector<DiffDoub0>& lay_q, vector<DiffDoub0>& lay_ang, vector<Section>& sec_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int num_lay;
	DiffDoub0 z_max;
	DiffDoub0 z_min;
	DiffDoub0 thk;
	DiffDoub0 qmat[9];
	DiffDoub0 tmp;
	DiffDoub0 tmp2;

	for (i1 = 0; i1 < 81; i1++) {
		cmat[i1].set_val(0.0);
	}

	num_lay = sec_ar[sect_ptr].layers.size();
	i2 = 0;
	for (i1 = 0; i1 < num_lay; i1++) {
		thk.set_val(0.5);
		thk.mult(lay_thk[i1]);
		z_min.set_val_dfd0(lay_z[i1]);
		z_min.sub(thk);
		z_max.set_val_dfd0(lay_z[i1]);
		z_max.add(thk);
		transform_q_dfd0(qmat, &lay_q[i2], lay_ang[i1]);
		
		// a matrix portion
		tmp.set_val_dfd0(z_max);
		tmp.sub(z_min);
		for (i3 = 0; i3 < 3; i3++) {
			i5 = 9 * i3;
			i6 = 3 * i3;
			for (i4 = 0; i4 < 3; i4++) {
				tmp2.set_val_dfd0(tmp);
				tmp2.mult(qmat[i6]);
				cmat[i5].add(tmp2);
				i5++;
				i6++;
			}
		}

		// b matrix portion
		tmp.set_val_dfd0(z_max);
		tmp.sqr();
		tmp2.set_val_dfd0(z_min);
		tmp2.sqr();
		tmp.sub(tmp2);
		tmp2.set_val(0.5);
		tmp.mult(tmp2);
		for (i3 = 0; i3 < 3; i3++) {
			i5 = 9 * i3 + 3;
			i6 = 3 * i3;
			for (i4 = 0; i4 < 3; i4++) {
				//i6 = i3*3 + i4;
				//i5 = i3*9 + (i4 + 3)
				tmp2.set_val_dfd0(tmp);
				tmp2.mult(qmat[i6]);
				cmat[i5].add(tmp2);
				i5++;
				i6++;
			}
		}

		// d matrix portion				
		tmp.set_val_dfd0(z_max);
		tmp.sqr();
		tmp.mult(z_max);
		tmp2.set_val_dfd0(z_min);
		tmp2.sqr();
		tmp2.mult(z_min);
		tmp.sub(tmp2);
		tmp2.set_val(r_1o3);
		tmp.mult(tmp2);
		for (i3 = 0; i3 < 3; i3++) {
			i5 = 9 * (i3 + 3) + 3;
			i6 = 3 * i3;
			for (i4 = 0; i4 < 3; i4++) {
				tmp2.set_val_dfd0(tmp);
				tmp2.mult(qmat[i6]);
				cmat[i5].add(tmp2);
				i5++;
				i6++;
			}
		}
		i2 += 9;
	}

	for (i1 = 1; i1 < 6; i1++) {
		i3 = 9 * i1;
		i4 = i1;
		for (i2 = 0; i2 < i1; i2++) {
			//i3 = i1*9 + i2
			//i4 = i2*9 + i1
			cmat[i3].set_val_dfd0(cmat[i4]);
			i3++;
			i4 += 9;
		}
	}

	tmp.set_val(1.0);
	tmp.mult(cmat[20]);
	cmat[60].set_val_dfd0(tmp);
	cmat[70].set_val_dfd0(tmp);
	cmat[80].set_val_dfd0(tmp);

	return;
}

void Element::get_beam_stiff_dfd0(vector<DiffDoub0>& cmat, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int dv_ind;
	int dv_comp;
	DiffDoub0 modulus_dv[3];
	DiffDoub0 shear_mod_dv[3];
	DiffDoub0 area_dv;
	DiffDoub0 idv[5];
	DiffDoub0 jdv;
	DiffDoub0 dv_val;
	DiffDoub0 coef;
	DiffDoub0 x_vec[6];
	DiffDoub0 b_vec[6];
	string dv_cat;

	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	if (this_sec.stiffness[0] > 0.0) {
		for (i1 = 0; i1 < 36; i1++) {
			cmat[i1].set_val(this_sec.stiffness[i1]);
		}
		for (auto& dv : design_vars) {
			dv_ind = dv.int_dat;
			DesignVariable& this_dv = dv_ar[dv_ind];
			if (this_dv.category == "stiffnessMat") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				cmat[dv_comp].add(dv_val);
			}
		}
		for (i1 = 1; i1 < 6; i1++) {
			i4 = 6 * i1;
			i5 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				cmat[i4].set_val_dfd0(cmat[i5]);
				i4++;
				i5 += 6;
			}
		}
	}
	else {
		for (i1 = 0; i1 < 3; i1++) {
			modulus_dv[i1].set_val(this_mat.modulus[i1]);
			shear_mod_dv[i1].set_val(this_mat.shear_mod[i1]);
		}
		area_dv.set_val(this_sec.area);
		for (i1 = 0; i1 < 5; i1++) {
			idv[i1].set_val(this_sec.area_moment[i1]);
		}
		jdv.set_val(this_sec.polar_moment);
		for (auto& dv : design_vars) {
			dv_ind = dv.int_dat;
			DesignVariable& this_dv = dv_ar[dv_ind];
			if (this_dv.category == "modulus") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				modulus_dv[dv_comp - 1].add(dv_val);
			}
			else if (this_dv.category == "shearModulus") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				shear_mod_dv[dv_comp - 1].add(dv_val);
			}
			else if (this_dv.category == "area") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				area_dv.add(dv_val);
			}
			else if (this_dv.category == "areaMoment") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				idv[dv_comp - 1].add(dv_val);
			}
			else if (this_dv.category == "polarMoment") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				jdv.add(dv_val);
			}
		}

		for (i1 = 0; i1 < 6; i1++) {
			i3 = 7 * i1;
			for (i2 = i1; i2 < 6; i2++) {
				cmat[i3].set_val(0.0);
				i3++;
			}
		}

		cmat[0].set_val_dfd0(modulus_dv[0]);
		cmat[0].mult(area_dv);
		cmat[4].set_val_dfd0(modulus_dv[0]);
		cmat[4].mult(idv[0]);
		cmat[5].set_val_dfd0(modulus_dv[0]);
		cmat[5].neg();
		cmat[5].mult(idv[1]);
		cmat[7].set_val_dfd0(shear_mod_dv[0]);
		cmat[7].mult(area_dv);
		cmat[9].set_val_dfd0(shear_mod_dv[0]);
		cmat[9].neg();
		cmat[9].mult(idv[0]);
		cmat[14].set_val_dfd0(shear_mod_dv[1]);
		cmat[14].mult(area_dv);
		cmat[15].set_val_dfd0(shear_mod_dv[2]);
		cmat[15].mult(idv[1]);
		cmat[21].set_val_dfd0(shear_mod_dv[0]);
		cmat[21].mult(jdv);
		cmat[28].set_val_dfd0(modulus_dv[0]);
		cmat[28].mult(idv[2]);
		cmat[29].set_val_dfd0(modulus_dv[0]);
		cmat[29].neg();
		cmat[29].mult(idv[4]);
		cmat[35].set_val_dfd0(modulus_dv[0]);
		cmat[35].mult(idv[3]);

		for (i1 = 1; i1 < 6; i1++) {
			i3 = 6 * i1;
			i4 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				cmat[i3].set_val_dfd0(cmat[i4]);
				i3++;
				i4 += 6;
			}
		}
	}

	return;
}

void Element::get_thermal_exp_dfd0(vector<DiffDoub0>& th_exp, vector<DiffDoub0>& einit, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	DesignVariable* this_dv;
	int dv_comp;
	DiffDoub0 dv_val;
	DiffDoub0 tmp;
	
	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	for (i1 = 0; i1 < 6; i1++) {
		th_exp[i1].set_val(this_mat.expansion[i1]);
		einit[i1].set_val(0.0);
	}

	for (auto& dv : design_vars) {
		DesignVariable& this_dv = dv_ar[dv.int_dat];
		if (this_dv.category == "thermalExp") {
			dv_comp = this_dv.component - 1;
			tmp.set_val(dv.doub_dat);
			this_dv.get_value_dfd0(dv_val);
			dv_val.mult(tmp);
			th_exp[dv_comp].add(dv_val);
		}
		else if (this_dv.category == "initialStrain") {
			dv_comp = this_dv.component - 1;
			tmp.set_val(dv.doub_dat);
			this_dv.get_value_dfd0(dv_val);
			dv_val.mult(tmp);
			einit[dv_comp].add(dv_val);
		}
	}

	return;
}

void Element::get_shell_exp_load_dfd0(vector<DiffDoub0>& exp_ld, vector<DiffDoub0>& e0_ld, vector<DiffDoub0>& lay_thk, vector<DiffDoub0>& lay_z, vector<DiffDoub0>& lay_q, vector<DiffDoub0>& lay_th_exp, vector<DiffDoub0>& lay_einit, vector<DiffDoub0>& lay_ang, vector<Section>& sec_ar) {
	int i1;
	int num_lay;
	int layi;
	int qi;
	int exi;
	DiffDoub0 sect_q[9];
	DiffDoub0 sect_te[3];
	DiffDoub0 sect_e0[3];
	DiffDoub0 qte_prod[3];
	DiffDoub0 qe0_prod[3];
	DiffDoub0 this_q[9];
	DiffDoub0 this_te[3];
	DiffDoub0 this_e0[3];
	DiffDoub0 z_min;
	DiffDoub0 z_max;
	DiffDoub0 tmp;
	DiffDoub0 tmp2;

	for (i1 = 0; i1 < 6; i1++) {
		exp_ld[i1].set_val(0.0);
		e0_ld[i1].set_val(0.0);
	}

	num_lay = sec_ar[sect_ptr].layers.size();
	qi = 0;
	exi = 0;
	for (layi = 0; layi < num_lay; layi++) {
		vec_to_ar_dfd0(this_q, lay_q, qi, qi + 9);
		transform_q_dfd0(sect_q, this_q, lay_ang[layi]);
		vec_to_ar_dfd0(this_te, lay_th_exp, exi, exi + 3);
		transform_strain_dfd0(sect_te, this_te, lay_ang[layi]);
		vec_to_ar_dfd0(this_e0, lay_einit, exi, exi + 3);
		transform_strain_dfd0(sect_e0, this_e0, lay_ang[layi]);
		mat_mul_ar_dfd0(qte_prod, sect_q, sect_te, 3, 3, 1);
		mat_mul_ar_dfd0(qe0_prod, sect_q, sect_e0, 3, 3, 1);
		
		for (i1 = 0; i1 < 3; i1++) {
			tmp.set_val_dfd0(qte_prod[i1]);
			tmp.mult(lay_thk[layi]);
			exp_ld[i1].add(tmp);
			tmp.set_val_dfd0(qe0_prod[i1]);
			tmp.mult(lay_thk[layi]);
			e0_ld[i1].add(tmp);
		}

		tmp.set_val(0.5);
		tmp.mult(lay_thk[layi]);
		z_min.set_val_dfd0(lay_z[layi]);
		z_min.sub(tmp);
		z_min.sqr();
		z_max.set_val_dfd0(lay_z[layi]);
		z_max.add(tmp);
		z_max.sqr();
		tmp.set_val(0.5);
		tmp2.set_val_dfd0(z_max);
		tmp2.sub(z_min);
		tmp.mult(tmp2); // tmp = 0.5*(z_max^2 - z_min^2)
		for (i1 = 0; i1 < 3; i1++) {
			tmp2.set_val_dfd0(qte_prod[i1]);
			tmp2.mult(tmp);
			exp_ld[i1 + 3].add(tmp2);
			tmp2.set_val_dfd0(qe0_prod[i1]);
			tmp2.mult(tmp);
			e0_ld[i1 + 3].add(tmp2);
		}
		qi += 9;
		exi += 3;
	}
	
	return;
}

void Element::get_beam_exp_load_dfd0(vector<DiffDoub0>& exp_ld, vector<DiffDoub0>& e0_ld, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	DiffDoub0 dv_val;
	string cat;
	string cat_list;
	DiffDoub0 tmp;
	int dv_comp;
	DiffDoub0 mod_dv[3];
	DiffDoub0 shr_mod_dv[3];
	DiffDoub0 te_coef_dv[6];
	DiffDoub0 e0_dv[6];
	DiffDoub0 area_dv;
	DiffDoub0 idv[5];
	DiffDoub0 qmat[9];
	DiffDoub0 qte[3];
	DiffDoub0 qe0[3];
	DiffDoub0 dedgu[18];
	DiffDoub0 tmp_exp[6];
	
	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	if (this_sec.exp_load_coef[0] > 0.0) {
		for (i1 = 0; i1 < 6; i1++) {
			exp_ld[i1].set_val(this_sec.exp_load_coef[i1]);
			e0_ld[i1].set_val(0.0);
		}
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			cat = this_dv.category;
			if (cat == "thermalExp") {
				dv_comp = this_dv.component - 1;
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(tmp);
				exp_ld[dv_comp].add(dv_val);
			}
			else if (cat == "initialStrain") {
				dv_comp = this_dv.component - 1;
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(tmp);
				e0_ld[dv_comp].add(dv_val);
			}
		}
	}
	else {
		for (i1 = 0; i1 < 3; i1++) {
			mod_dv[i1].set_val(this_mat.modulus[i1]);
			shr_mod_dv[i1].set_val(this_mat.shear_mod[i1]);
			te_coef_dv[i1].set_val(this_mat.expansion[i1]);
			te_coef_dv[i1 + 3].set_val(this_mat.expansion[i1]);
			e0_dv[i1].set_val(0.0);
			e0_dv[i1 + 3].set_val(0.0);
		}
		area_dv.set_val(this_sec.area);
		for (i1 = 0; i1 < 5; i1++) {
			idv[i1].set_val(this_sec.area_moment[i1]);
		}
		cat_list = "modulus shearModulus thermalExp initialStrain area areaMoment";
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			cat = this_dv.category;
			i2 = cat_list.find(cat);
			if (i2 > -1) {
				dv_comp = this_dv.component - 1;
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(tmp);
				if (cat == "modulus") {
					mod_dv[dv_comp].add(dv_val);
				}
				else if (cat == "shearModulus") {
					shr_mod_dv[dv_comp].add(dv_val);
				}
				else if (cat == "thermalExp") {
					te_coef_dv[dv_comp].add(dv_val);
				}
				else if (cat == "initialStrain") {
					e0_dv[dv_comp].add(dv_val);
				}
				else if (cat == "area") {
					area_dv.add(dv_val);
				}
				else if (cat == "areaMoment") {
					idv[dv_comp].add(dv_val);
				}
			}
		}
		for (i1 = 1; i1 < 8; i1++) {
			qmat[i1].set_val(0.0);
		}
		qmat[0].set_val_dfd0(mod_dv[0]);
		qmat[4].set_val_dfd0(shr_mod_dv[0]);
		qmat[8].set_val_dfd0(shr_mod_dv[1]);
		te_coef_dv[1].set_val_dfd0(te_coef_dv[3]);
		te_coef_dv[2].set_val_dfd0(te_coef_dv[4]);
		mat_mul_ar_dfd0(qte, qmat, te_coef_dv, 3, 3, 1);
		e0_dv[1].set_val_dfd0(e0_dv[3]);
		e0_dv[2].set_val_dfd0(e0_dv[4]);
		mat_mul_ar_dfd0(qe0, qmat, e0_dv, 3, 3, 1);
		for (i1 = 0; i1 < 18; i1++) {
			dedgu[i1].set_val(0.0);
		}
		dedgu[0].set_val_dfd0(area_dv);
		dedgu[4].set_val_dfd0(area_dv);
		dedgu[8].set_val_dfd0(area_dv);
		dedgu[10].set_val_dfd0(idv[0]);
		dedgu[10].neg();
		dedgu[11].set_val_dfd0(idv[1]);
		dedgu[12].set_val_dfd0(idv[0]);
		dedgu[15].set_val_dfd0(idv[1]);
		dedgu[15].neg();

		mat_mul_ar_dfd0(tmp_exp, dedgu, qte, 6, 3, 1);
		ar_to_vec_dfd0(tmp_exp, exp_ld, 0, 6);
		mat_mul_ar_dfd0(tmp_exp, dedgu, qe0, 6, 3, 1);
		ar_to_vec_dfd0(tmp_exp, e0_ld, 0, 6);
	}

	return;
}

void Element::get_density_dfd0(DiffDoub0& den, int layer, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int layi;
	string cat;
	int dv_lay;
	DiffDoub0 coef;
	DiffDoub0 dv_val;

	Section& this_sec = sec_ar[sect_ptr];
	if (type == 3 || type == 41) {
		Material& this_mat = mat_ar[this_sec.get_layer_mat_ptr(layer)];
		den.set_val(this_mat.density);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			cat = this_dv.category;
			dv_lay = this_dv.layer;
			if (cat == "density" && dv_lay == layer) {
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				den.add(dv_val);
			}
		}
	} else {
		Material& this_mat = mat_ar[this_sec.mat_ptr];
		den.set_val(this_mat.density);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			cat = this_dv.category;
			if (cat == "density") {
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				den.add(dv_val);
			}
		}
	}

	return;
}

void Element::get_shell_mass_dfd0(vector<DiffDoub0>& mmat, vector<DiffDoub0>& lay_thk, vector<DiffDoub0>& lay_z, vector<DiffDoub0>& lay_den, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int layi;
	DiffDoub0 tmp;
	DiffDoub0 tmp2;
	DiffDoub0 z_min;
	DiffDoub0 z_min2;
	DiffDoub0 z_max;
	DiffDoub0 z_max2;

	for (i1 = 0; i1 < 36; i1++) {
		mmat[i1].set_val(0.0);
	}

	i1 = sec_ar[sect_ptr].layers.size();
	for (layi = 0; layi < i1; layi++) {
		tmp.set_val_dfd0(lay_den[layi]);
		tmp.mult(lay_thk[layi]);
		mmat[0].add(tmp);
		mmat[7].add(tmp);
		mmat[14].add(tmp);

		tmp.set_val(0.5);
		tmp.mult(lay_thk[layi]);
		z_min.set_val_dfd0(lay_z[layi]);
		z_min.sub(tmp);
		z_min2.set_val_dfd0(z_min);
		z_min2.sqr();
		z_max.set_val_dfd0(lay_z[layi]);
		z_max.add(tmp);
		z_max2.set_val_dfd0(z_max);
		z_max2.sqr();
		tmp.set_val_dfd0(z_max2);
		tmp.sub(z_min2);  // tmp == z_max^2 - z_min^2
		tmp2.set_val(0.5);
		tmp2.mult(lay_den[layi]);
		tmp2.mult(tmp); // tmp2 = 0.5*rho*(z_max^2 - z_min^2)
		mmat[4].add(tmp2);
		mmat[24].add(tmp2);
		mmat[9].sub(tmp2);
		mmat[19].sub(tmp2);

		z_max2.mult(z_max);
		z_min2.mult(z_min);
		tmp.set_val_dfd0(z_max2);
		tmp.sub(z_min2); // tmp == z_max^3 - z_min^3
		tmp2.set_val(r_1o3);
		tmp2.mult(lay_den[layi]);
		tmp2.mult(tmp);
		mmat[21].add(tmp2);
		mmat[28].add(tmp2);
		mmat[35].add(tmp2);
	}

	return;
}

void Element::get_beam_mass_dfd0(vector<DiffDoub0>& mmat, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int dv_comp;

	DesignVariable* this_dv;
	DiffDoub0 dv_val;
	DiffDoub0 coef;

	string d_cat;
	Material* mat_pt;
	double* area_mom;
	DiffDoub0 den_dv;
	DiffDoub0 area_dv;
	DiffDoub0 idv[5];
	DiffDoub0 jdv;
	DiffDoub0 tmp;

	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	if (this_sec.mass[0] > 0.0) {
		for (i1 = 0; i1 < 36; i1++) {
			mmat[i1].set_val(this_sec.mass[i1]);
		}
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "massMat") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				mmat[dv_comp].add(dv_val);
			}
		}
		for (i1 = 1; i1 < 6; i1++) {
			i3 = 6 * i1; // lower tri term
			i4 = i1; // upper tri term
			for (i2 = 0; i2 < i1; i2++) {
				mmat[i3].set_val_dfd0(mmat[i4]);
				i3++;
				i4 += 6;
			}
		}
	}
	else {
		for (i1 = 0; i1 < 36; i1++) {
			mmat[i1].set_val(0.0);
		}
		for (i1 = 0; i1 < 5; i1++) {
			idv[i1].set_val(this_sec.area_moment[i1]);
		}
		area_dv.set_val(this_sec.area);
		jdv.set_val(this_sec.polar_moment);
		den_dv.set_val(this_mat.density);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			d_cat = this_dv.category;
			if (d_cat == "density") {
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				den_dv.add(dv_val);
			}
			else if (d_cat == "area") {
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				area_dv.add(dv_val);
			}
			else if (d_cat == "areaMoment") {
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				dv_comp = this_dv.component;
				idv[dv_comp - 1].add(dv_val);
			}
			else if (d_cat == "polarMoment") {
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				jdv.add(dv_val);
			}
		}
		tmp.set_val_dfd0(den_dv);
		tmp.mult(area_dv);
		mmat[0].set_val_dfd0(tmp);
		mmat[7].set_val_dfd0(tmp);
		mmat[14].set_val_dfd0(tmp);
		tmp.set_val_dfd0(den_dv);
		tmp.mult(idv[0]);
		mmat[4].set_val_dfd0(tmp);
		mmat[9].set_val_dfd0(tmp);
		mmat[9].neg();
		tmp.set_val_dfd0(den_dv);
		tmp.mult(idv[1]);
		mmat[5].set_val_dfd0(tmp);
		mmat[5].neg();
		mmat[15].set_val_dfd0(tmp);
		tmp.set_val_dfd0(den_dv);
		tmp.mult(jdv);
		mmat[21].set_val_dfd0(tmp);
		tmp.set_val_dfd0(den_dv);
		tmp.mult(idv[2]);
		mmat[28].set_val_dfd0(tmp);
		tmp.set_val_dfd0(den_dv);
		tmp.mult(idv[4]);
		mmat[29].set_val_dfd0(tmp);
		tmp.set_val_dfd0(den_dv);
		tmp.mult(idv[3]);
		mmat[35].set_val_dfd0(tmp);

		for (i1 = 3; i1 < 6; i1++) {
			i3 = 6 * i1; // lower tri term
			i4 = i1; // upper tri term
			for (i2 = 0; i2 < i1; i2++) {
				mmat[i3].set_val_dfd0(mmat[i4]);
				i3++;
				i4 += 6;
			}
		}
	}

	return;
}

void Element::get_solid_damp_dfd0(vector<DiffDoub0>& dmat, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	DesignVariable* this_dv;
	DiffDoub0 temp;
	DiffDoub0 dv_val;
	int dv_comp;

	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	for (i1 = 0; i1 < 36; i1++) {
		dmat[i1].set_val(this_mat.damping[i1]);
	}

	for (auto& dv : design_vars) {
		DesignVariable& this_dv = dv_ar[dv.int_dat];
		if (this_dv.category == "dampingMat") {
			this_dv.get_value_dfd0(dv_val);
			dv_comp = this_dv.component;
			temp.set_val(dv.doub_dat);
			dv_val.mult(temp);
			dmat[dv_comp].add(dv_val);
		}
	}

	for (i1 = 1; i1 < 6; i1++) {
		i3 = 6 * i1; // lower tri term
		i4 = i1; // upper tri term
		for (i2 = 0; i2 < i1; i2++) {
			dmat[i3].set_val_dfd0(dmat[i4]);
			i3++;
			i4 += 6;
		}
	}

	return;
}

void Element::get_shell_damp_dfd0(vector<DiffDoub0>& dmat, vector<DiffDoub0>& lay_thk, vector<DiffDoub0>& lay_z, vector<DiffDoub0>& lay_d, vector<DiffDoub0>& lay_ang, vector<Section>& sec_ar) {
	get_abd_dfd0(dmat, lay_thk, lay_z, lay_d, lay_ang, sec_ar);
	dmat[60].set_val(0.0);
	dmat[70].set_val(0.0);
	dmat[80].set_val(0.0);
	return;
}

void Element::get_beam_damp_dfd0(vector<DiffDoub0>& dmat, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int dv_ind;
	int dv_comp;
	DiffDoub0 dmat_dv[36];
	DiffDoub0 area_dv;
	DiffDoub0 idv[5];
	DiffDoub0 jdv;
	DiffDoub0 dv_val;
	DiffDoub0 coef;
	DiffDoub0 x_vec[6];
	DiffDoub0 b_vec[6];
	string dv_cat;

	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	if (this_sec.damping[0] > 0.0) {
		for (i1 = 0; i1 < 36; i1++) {
			dmat[i1].set_val(this_sec.damping[i1]);
		}
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "dampingMat") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				dmat[dv_comp].add(dv_val);
			}
		}
		for (i1 = 1; i1 < 6; i1++) {
			i4 = 6 * i1;
			i5 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				dmat[i4].set_val_dfd0(dmat[i5]);
				i4++;
				i5 += 6;
			}
		}
	}
	else {
		for (i1 = 0; i1 < 36; i1++) {
			dmat_dv[i1].set_val(this_mat.damping[i1]);
		}
		area_dv.set_val(this_sec.area);
		for (i1 = 0; i1 < 5; i1++) {
			idv[i1].set_val(this_sec.area_moment[i1]);
		}
		jdv.set_val(this_sec.polar_moment);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "dampingMat") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				dmat_dv[dv_comp].add(dv_val);
			}
			else if (this_dv.category == "area") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				area_dv.add(dv_val);
			}
			else if (this_dv.category == "areaMoment") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				idv[dv_comp - 1].add(dv_val);
			}
			else if (this_dv.category == "polarMoment") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(coef);
				jdv.add(dv_val);
			}
		}

		for (i1 = 0; i1 < 6; i1++) {
			i3 = 7 * i1;
			for (i2 = i1; i2 < 6; i2++) {
				dmat[i3].set_val(0.0);
				i3++;
			}
		}

		dmat[0].set_val_dfd0(dmat_dv[0]);
		dmat[0].mult(area_dv);
		dmat[4].set_val_dfd0(dmat_dv[0]);
		dmat[4].mult(idv[0]);
		dmat[5].set_val_dfd0(dmat_dv[0]);
		dmat[5].neg();
		dmat[5].mult(idv[1]);
		dmat[7].set_val_dfd0(dmat_dv[21]);
		dmat[7].mult(area_dv);
		dmat[9].set_val_dfd0(dmat_dv[21]);
		dmat[9].neg();
		dmat[9].mult(idv[0]);
		dmat[14].set_val_dfd0(dmat_dv[28]); // 28
		dmat[14].mult(area_dv);
		dmat[15].set_val_dfd0(dmat_dv[35]); // 35
		dmat[15].mult(idv[1]);
		dmat[21].set_val_dfd0(dmat_dv[21]); // 21
		dmat[21].mult(jdv);
		dmat[28].set_val_dfd0(dmat_dv[0]);
		dmat[28].mult(idv[2]);
		dmat[29].set_val_dfd0(dmat_dv[0]);
		dmat[29].neg();
		dmat[29].mult(idv[4]);
		dmat[35].set_val_dfd0(dmat_dv[0]);
		dmat[35].mult(idv[3]);

		for (i1 = 1; i1 < 6; i1++) {
			i3 = 6 * i1;
			i4 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				dmat[i3].set_val_dfd0(dmat[i4]);
				i3++;
				i4 += 6;
			}
		}
	}


	return;
}

void Element::get_conductivity_dfd0(vector<DiffDoub0>& t_cond, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	DiffDoub0 cond_dv[6];
	DesignVariable* this_dv;
	DiffDoub0 temp;
	DiffDoub0 dv_val;
	int dv_comp;

	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	for (i1 = 0; i1 < 6; i1++) {
		cond_dv[i1].set_val(this_mat.conductivity[i1]);
	}

	for (auto& dv : design_vars) {
		DesignVariable& this_dv = dv_ar[dv.int_dat];
		if (this_dv.category == "thermalCond") {
			dv_comp = this_dv.component - 1;
			temp.set_val(dv.doub_dat);
			this_dv.get_value_dfd0(dv_val);
			dv_val.mult(temp);
			cond_dv[dv_comp].add(dv_val);
		}
	}

	t_cond[0] = cond_dv[0];
	t_cond[1] = cond_dv[3];
	t_cond[2] = cond_dv[4];
	t_cond[3] = cond_dv[3];
	t_cond[4] = cond_dv[1];
	t_cond[5] = cond_dv[5];
	t_cond[6] = cond_dv[4];
	t_cond[7] = cond_dv[5];
	t_cond[8] = cond_dv[2];

	return;
}

void Element::get_shell_cond_dfd0(vector<DiffDoub0>& t_cond, vector<DiffDoub0>& lay_thk, vector<DiffDoub0>& lay_ang, vector<DiffDoub0>& lay_cond, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int num_lay = sec_ar[sect_ptr].layers.size();
	DiffDoub0 cond_dv[6];
	DiffDoub0 layer_mat[9];
	DiffDoub0 al_mat[9];
	DiffDoub0 al_t[9];
	DiffDoub0 tmp[9];
	DiffDoub0 this_cond[9];

	for (i1 = 0; i1 < 9; i1++) {
		t_cond[i1].set_val(0.0);
		al_mat[i1].set_val(0.0);
	}
	al_mat[8].set_val(1.0);

	for (i1 = 0; i1 < num_lay; i1++) {
		al_mat[0].set_val(r_pio180);
		al_mat[0].mult(lay_ang[i1]);
		al_mat[0].cs();
		al_mat[4].set_val_dfd0(al_mat[0]);
		al_mat[3].set_val(r_pio180);
		al_mat[3].mult(lay_ang[i1]);
		al_mat[3].sn();
		al_mat[1].set_val_dfd0(al_mat[3]);
		al_mat[1].neg();
		transpose_ar_dfd0(al_t, al_mat, 3, 3);
		i2 = 9 * i1;
		vec_to_ar_dfd0(this_cond, lay_cond, i2, i2 + 9);
		mat_mul_ar_dfd0(tmp, this_cond, al_t, 3, 3, 3);
		mat_mul_ar_dfd0(layer_mat, al_mat, tmp, 3, 3, 3);
		for (i2 = 0; i2 < 9; i2++) {
			layer_mat[i2].mult(lay_thk[i1]);
			t_cond[i2].add(layer_mat[i2]);
		}
	}

	return;
}

void Element::get_beam_cond_dfd0(vector<DiffDoub0>& t_cond, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	DesignVariable* this_dv;
	string cat;
	DiffDoub0 cond_dv;
	DiffDoub0 area_dv;
	DiffDoub0 tmp;
	DiffDoub0 dv_val;
	
	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	if (this_sec.conductivity > 0.0) {
		cond_dv.set_val(this_sec.conductivity);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "thermCond") {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(tmp);
				cond_dv.add(dv_val);
			}
		}
	}
	else {
		cond_dv.set_val(this_mat.conductivity[0]);
		area_dv.set_val(this_sec.area);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			cat = this_dv.category;
			if (cat == "thermCond") {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(tmp);
				cond_dv.add(dv_val);
			}
			else if (cat == "area") {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(tmp);
				area_dv.add(dv_val);
			}
		}
		cond_dv.mult(area_dv);
	}

	for (i1 = 1; i1 < 8; i1++) {
		t_cond[i1].set_val(0.0);
	}

	t_cond[0].set_val_dfd0(cond_dv);
	t_cond[4].set_val_dfd0(cond_dv);
	t_cond[8].set_val_dfd0(cond_dv);

	return;
}

void Element::get_specific_heat_dfd0(DiffDoub0& spec_heat, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	DiffDoub0 spec_heat_dv;
	DiffDoub0 temp;
	DiffDoub0 dv_val;
	int dv_comp;

	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	spec_heat.set_val(this_mat.spec_heat);

	for (auto& dv : design_vars) {
		DesignVariable& this_dv = dv_ar[dv.int_dat];
		if (this_dv.category == "specHeat") {
			temp.set_val(dv.doub_dat);
			this_dv.get_value_dfd0(dv_val);
			dv_val.mult(temp);
			spec_heat.add(dv_val);
		}
	}
	
	return;
}

void Element::get_shell_spec_heat_dfd0(DiffDoub0& spec_heat, vector<DiffDoub0>& lay_thk, vector<DiffDoub0>& lay_sh, vector<DiffDoub0>& lay_den, vector<Section>& sec_ar) {
	int i1;
	int num_lay;
	DiffDoub0 tmp;

	spec_heat.set_val(0.0);
	num_lay = sec_ar[sect_ptr].layers.size();
	for (i1 = 0; i1 < num_lay; i1++) {
		tmp.set_val_dfd0(lay_den[i1]);
		tmp.mult(lay_sh[i1]);
		tmp.mult(lay_thk[i1]);
		spec_heat.add(tmp);
	}

	return;
}

void Element::get_beam_spec_heat_dfd0(DiffDoub0& spec_heat, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	string cat;
	DiffDoub0 density_dv;
	DiffDoub0 area_dv;
	DiffDoub0 tmp;
	DiffDoub0 dv_val;

	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	double sec_sh = this_sec.spec_heat;
	double mat_sh = this_mat.spec_heat;
	double mat_den = this_mat.density;

	if (sec_sh > 0.0) {
		spec_heat.set_val(sec_sh);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "specHeat") {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(tmp);
				spec_heat.add(dv_val);
			}
		}
	}
	else {
		spec_heat.set_val(mat_sh);
		density_dv.set_val(mat_den);
		area_dv.set_val(this_sec.area);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			cat = this_dv.category;
			if (cat == "specHeatDV") {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(tmp);
				spec_heat.add(dv_val);
			}
			else if (cat == "density") {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(tmp);
				density_dv.add(dv_val);
			}
			else if (cat == "area") {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd0(dv_val);
				dv_val.mult(tmp);
				area_dv.add(dv_val);
			}
		}
		spec_heat.mult(density_dv);
		spec_heat.mult(area_dv);
	}

	return;
}

void Element::get_nd_crds_dfd0(vector<DiffDoub0>& x_glob, vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	DiffDoub0 nd_crd[3];
	
	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		this_nd.get_crd_dfd0(nd_crd,dv_ar);
		x_glob[i1].set_val_dfd0(nd_crd[0]);
		x_glob[i1+num_nds].set_val_dfd0(nd_crd[1]);
		x_glob[i1+2*num_nds].set_val_dfd0(nd_crd[2]);
	}
	
	return;
}

void Element::get_loc_ori_dfd0(vector<DiffDoub0>& loc_ori, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int dv_ind;
	DesignVariable *this_dvpt;
	string dv_cat;
	DiffDoub0 rot[3];
	DiffDoub0 dv_val;
	DiffDoub0 coef;
	DiffDoub0 ori_copy[9];
	DiffDoub0 tmp_ori[9];
	
	Section& this_sec = sec_ar[sect_ptr];
	for (i1 = 0; i1 < 9; i1++) {
		ori_copy[i1].set_val(this_sec.orientation[i1]);
	}
	
	rot[0].set_val(0.0);
	rot[1].set_val(0.0);
    rot[2].set_val(0.0);
	for (auto& dv : design_vars) {
		DesignVariable& this_dv = dv_ar[dv.int_dat];
		dv_cat = this_dv.category;
		if(dv_cat == "orientation") {
			i1 = this_dv.component - 1;
			coef.set_val(dv.doub_dat);
			this_dv.get_value_dfd0(dv_val);
			dv_val.mult(coef);
			rot[i1].add(dv_val);
		}
	}
	
	rotate_orient_dfd0(tmp_ori, ori_copy, rot);
	ar_to_vec_dfd0(tmp_ori, loc_ori, 0, 9);
	
	return;
}

void Element::correct_orient_dfd0(vector<DiffDoub0>& loc_ori, vector<DiffDoub0>& x_glob) {
	int i1;
	DiffDoub0 v1[3];
	DiffDoub0 v2[3];
	DiffDoub0 v3[3];
	DiffDoub0 rot[3];
	DiffDoub0 ori_copy[9];
	DiffDoub0 tmp_ori[9];
	
	DiffDoub0 dp;
	DiffDoub0 magv3;
	DiffDoub0 mag_cp;
	DiffDoub0 theta;
	DiffDoub0 tmp;
	
	for (i1 = 0; i1 < 9; i1++) {
		ori_copy[i1].set_val_dfd0(loc_ori[i1]);
	}
	
	if(type == 3 || type == 41) {
		if(type == 3) {
			v1[0].set_val_dfd0(x_glob[1]);
			v1[0].sub(x_glob[0]);
			v1[1].set_val_dfd0(x_glob[4]);
			v1[1].sub(x_glob[3]);
			v1[2].set_val_dfd0(x_glob[7]);
			v1[2].sub(x_glob[6]);
			
			v2[0].set_val_dfd0(x_glob[2]);
			v2[0].sub(x_glob[0]);
			v2[1].set_val_dfd0(x_glob[5]);
			v2[1].sub(x_glob[3]);
			v2[2].set_val_dfd0(x_glob[8]);
			v2[2].sub(x_glob[6]);		
		} else {
			v1[0].set_val_dfd0(x_glob[2]);
			v1[0].sub(x_glob[0]);
			v1[1].set_val_dfd0(x_glob[6]);
			v1[1].sub(x_glob[4]);
			v1[2].set_val_dfd0(x_glob[10]);
			v1[2].sub(x_glob[8]);
			
			v2[0].set_val_dfd0(x_glob[3]);
			v2[0].sub(x_glob[1]);
			v2[1].set_val_dfd0(x_glob[7]);
			v2[1].sub(x_glob[5]);
			v2[2].set_val_dfd0(x_glob[11]);
			v2[2].sub(x_glob[9]);	
		}
		cross_prod_dfd0(v3, v1, v2);
		
		dp.set_val_dfd0(v3[0]);
		dp.mult(loc_ori[6]);
		tmp.set_val_dfd0(v3[1]);
		tmp.mult(loc_ori[7]);
		dp.add(tmp);
		tmp.set_val_dfd0(v3[2]);
		tmp.mult(loc_ori[8]);
		dp.add(tmp);
		
		if(dp.val < 0.0) {
			v3[0].neg();
			v3[1].neg();
			v3[2].neg();
		}
		
		magv3.set_val_dfd0(v3[0]);
		magv3.sqr();
		tmp.set_val_dfd0(v3[1]);
		tmp.sqr();
		magv3.add(tmp);
		tmp.set_val_dfd0(v3[2]);
		tmp.sqr();
		magv3.add(tmp);
		magv3.sqt();
		
		cross_prod_dfd0(rot,&loc_ori[6],v3);
		mag_cp.set_val_dfd0(rot[0]);
		mag_cp.sqr();
		tmp.set_val_dfd0(rot[1]);
		tmp.sqr();
		mag_cp.add(tmp);
		tmp.set_val_dfd0(rot[2]);
		tmp.sqr();
		mag_cp.add(tmp);
		mag_cp.sqt();
		if (mag_cp.val < 1e-12) {
			return;
		}
		
		theta.set_val_dfd0(mag_cp);
		theta.dvd(magv3);
		theta.asn();
		
		tmp.set_val_dfd0(theta);
		tmp.dvd(mag_cp);
		
		rot[0].mult(tmp);
		rot[1].mult(tmp);
		rot[2].mult(tmp);
		
		rotate_orient_dfd0(tmp_ori, ori_copy, rot);
		ar_to_vec_dfd0(tmp_ori, loc_ori, 0, 9);
	}
	
	return;
}

void Element::get_frc_fld_const_dfd0(vector<DiffDoub0>& coef, vector<DiffDoub0>& exp, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	DiffDoub0 dv_val;
	DiffDoub0 tmp;
	string cat;

	Section& this_sec = sec_ar[sect_ptr];
	coef[0].set_val(this_sec.pot_coef);
	coef[1].set_val(this_sec.damp_coef);
	exp[0].set_val(this_sec.pot_exp);
	exp[1].set_val(this_sec.damp_exp);

	for (auto& dv : design_vars) {
		DesignVariable& this_dv = dv_ar[dv.int_dat];
		cat = this_dv.category;
		if (cat == "potFldCoef") {
			this_dv.get_value_dfd0(dv_val);
			tmp.set_val(dv.doub_dat);
			dv_val.mult(tmp);
			coef[0].add(dv_val);
		}
		else if (cat == "dampFldCoef") {
			this_dv.get_value_dfd0(dv_val);
			tmp.set_val(dv.doub_dat);
			dv_val.mult(tmp);
			coef[1].add(dv_val);
		}
	}

	return;
}

void Element::get_thrm_fld_const_dfd0(vector<DiffDoub0>& coef, DiffDoub0& ref_t, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	DiffDoub0 dv_val;
	DiffDoub0 tmp;
	string cat;

	Section& this_sec = sec_ar[sect_ptr];
	coef[0].set_val(this_sec.cond_coef);
	coef[1].set_val(this_sec.rad_coef);
	ref_t.set_val(this_sec.ref_temp);

	for (auto& dv : design_vars) {
		DesignVariable& this_dv = dv_ar[dv.int_dat];
		cat = this_dv.category;
		if (cat == "condCoef") {
			this_dv.get_value_dfd0(dv_val);
			tmp.set_val(dv.doub_dat);
			dv_val.mult(tmp);
			coef[0].add(dv_val);
		}
		else if (cat == "radCoef") {
			this_dv.get_value_dfd0(dv_val);
			tmp.set_val(dv.doub_dat);
			dv_val.mult(tmp);
			coef[1].add(dv_val);
		}
	}

	return;
}

void Element::get_mass_per_el_dfd0(DiffDoub0& mass_per_el, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	DiffDoub0 dv_val;
	DiffDoub0 tmp;
	string cat;

	Section& this_sec = sec_ar[sect_ptr];
	mass_per_el.set_val(this_sec.mass_per_el);

	for(auto& dv : design_vars) {
		DesignVariable& this_dv = dv_ar[dv.int_dat];
		cat = this_dv.category;
		if (cat == "massPerEl") {
			this_dv.get_value_dfd0(dv_val);
			tmp.set_val(dv.doub_dat);
			dv_val.mult(tmp);
			mass_per_el.add(dv_val);
		}
	}

	return;
}

//end dup
 
//skip 
 
//DiffDoub1 versions: 
//dup1

void Element::get_gen_prop_dfd1(DiffDoub1& prop, string prop_key, vector<DesignVariable>& dv_ar) {
	int d_vind;
	DiffDoub1 dv_val;
	DiffDoub1 tmp;

	for (auto& dv : design_vars) {
		d_vind = dv.int_dat;
		DesignVariable& this_dv = dv_ar[d_vind];
		if (this_dv.category == prop_key) {
			this_dv.get_value_dfd1(dv_val);
			tmp.set_val(dv.doub_dat);
			dv_val.mult(tmp);
			prop.add(dv_val);
		}
	}

	return;
}

void Element::get_layer_thk_z_dfd1(vector<DiffDoub1>& lay_thk, vector<DiffDoub1>& lay_z, DiffDoub1& z_offset, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	//z_offset = 1: upper z surface is reference plane
	//z_offset = -1: lower z surface is reference plane
	int layi = 0;
	DiffDoub1 dv_val;
	DiffDoub1 coef;
	DiffDoub1 tot_thk;
	DiffDoub1 tmp;
	DiffDoub1 z_crd;
	DiffDoub1 z_next;
	DiffDoub1 z_mid;

	tot_thk.set_val(0.0);
	Section& this_sec = sec_ar[sect_ptr];
	for (auto& lay : this_sec.layers) {
		lay_thk[layi].set_val(lay.thickness);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "thickness" && this_dv.layer == layi) {
				this_dv.get_value_dfd1(dv_val);
				coef.set_val(dv.doub_dat);
				dv_val.mult(coef);
				lay_thk[layi].add(dv_val);
			}
		}
		tot_thk.add(lay_thk[layi]);
		layi++;
	}

	z_offset.set_val(this_sec.z_offset);
	for (auto& dv : design_vars) {
		DesignVariable& this_dv = dv_ar[dv.int_dat];
		if (this_dv.category == "zOffset") {
			this_dv.get_value_dfd1(dv_val);
			coef.set_val(dv.doub_dat);
			dv_val.mult(coef);
			z_offset.add(dv_val);
		}
	}

	tmp.set_val(1.0);
	tmp.add(z_offset);
	z_crd.set_val(-0.5);
	z_crd.mult(tot_thk);
	z_crd.mult(tmp);

	layi = 0;
	for (auto& lay : this_sec.layers) {
		z_next.set_val_dfd1(z_crd);
		z_next.add(lay_thk[layi]);
		tmp.set_val_dfd1(z_crd);
		tmp.add(z_next);
		z_mid.set_val(0.5);
		z_mid.mult(tmp);
		lay_z[layi].set_val_dfd1(z_mid);
		layi++;
		z_crd.set_val_dfd1(z_next);
	}

	return;
}

void Element::get_layer_q_dfd1(vector<DiffDoub1>& lay_q, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int layi;
	int dv_ind;
	int dv_comp;
	DiffDoub1 modulus_dv[3];
	DiffDoub1 poisson_dv[3];
	DiffDoub1 shear_mod_dv[3];
	DiffDoub1 smat[9];
	DiffDoub1 qmat[9];
	DiffDoub1 dv_val;
	DiffDoub1 coef;
	DiffDoub1 x_vec[3];
	DiffDoub1 b_vec[3];
	string dv_cat;

	i2 = 0;
	layi = 0;
	Section& this_sec = sec_ar[sect_ptr];
	for (auto& lay : this_sec.layers) {
		Material& this_mat = mat_ar[lay.mat_ptr];
		for (i1 = 0; i1 < 3; i1++) {
			modulus_dv[i1].set_val(this_mat.modulus[i1]);
			poisson_dv[i1].set_val(this_mat.poisson_ratio[i1]);
			shear_mod_dv[i1].set_val(this_mat.shear_mod[i1]);
		}
		for (auto& dv : design_vars) {
			dv_ind = dv.int_dat;
			DesignVariable& this_dv = dv_ar[dv_ind];
			if (this_dv.layer == layi) {
				dv_cat = this_dv.category;
				if (dv_cat == "modulus") {
					dv_comp = this_dv.component;
					coef.set_val(dv.doub_dat);
					this_dv.get_value_dfd1(dv_val);
					dv_val.mult(coef);
					modulus_dv[dv_comp - 1].add(dv_val);
				}
				else if (dv_cat == "poissonRatio") {
					dv_comp = this_dv.component;
					coef.set_val(dv.doub_dat);
					this_dv.get_value_dfd1(dv_val);
					dv_val.mult(coef);
					poisson_dv[dv_comp - 1].add(dv_val);
				}
				else if (dv_cat == "shearModulus") {
					dv_comp = this_dv.component;
					coef.set_val(dv.doub_dat);
					this_dv.get_value_dfd1(dv_val);
					dv_val.mult(coef);
					shear_mod_dv[dv_comp - 1].add(dv_val);
				}
			}
		}

		for (i1 = 0; i1 < 9; i1++) {
			smat[i1].set_val(0.0);
		}
		smat[0].set_val(1.0);
		smat[0].dvd(modulus_dv[0]);
		smat[1].set_val_dfd1(poisson_dv[0]);
		smat[1].neg();
		smat[1].dvd(modulus_dv[0]);
		smat[3].set_val_dfd1(smat[1]);
		smat[4].set_val(1.0);
		smat[4].dvd(modulus_dv[1]);
		smat[8].set_val(1.0);
		smat[8].dvd(shear_mod_dv[0]);
		get_det_inv_ar_dfd1(coef, qmat, smat, 3, 0, x_vec, b_vec);

		for (i1 = 0; i1 < 9; i1++) {
			lay_q[i2].set_val_dfd1(qmat[i1]);
			i2++;
		}

		layi++;
	}

	return;
}

void Element::get_layer_d_dfd1(vector<DiffDoub1>& lay_d, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int layi;
	int dv_ind;
	int dv_comp;
	DiffDoub1 damp_mat_dv[36];
	DiffDoub1 dv_val;
	DiffDoub1 coef;
	string dv_cat;

	i2 = 0;
	layi = 0;
	Section& this_sec = sec_ar[sect_ptr];
	for (auto& lay : this_sec.layers) {
		Material& this_mat = mat_ar[lay.mat_ptr];
		for (i1 = 0; i1 < 36; i1++) {
			damp_mat_dv[i1].set_val(this_mat.damping[i1]);
		}
		for (auto& dv : design_vars) {
			dv_ind = dv.int_dat;
			DesignVariable& this_dv = dv_ar[dv_ind];
			if (this_dv.layer == layi) {
				dv_cat = this_dv.category;
				if (dv_cat == "dampingMat") {
					dv_comp = this_dv.component;
					coef.set_val(dv.doub_dat);
					this_dv.get_value_dfd1(dv_val);
					dv_val.mult(coef);
					damp_mat_dv[dv_comp].add(dv_val);
				}
			}
		}

		for (i3 = 1; i3 < 6; i3++) {
			i5 = 6 * i3; // lower tri term
			i6 = i3; // upper tri term
			for (i4 = 0; i4 < i3; i4++) {
				damp_mat_dv[i5].set_val_dfd1(damp_mat_dv[i6]);
				i5++;
				i6 += 6;
			}
		}

		lay_d[i2].set_val_dfd1(damp_mat_dv[0]);
		lay_d[i2 + 1].set_val_dfd1(damp_mat_dv[1]);
		lay_d[i2 + 2].set_val_dfd1(damp_mat_dv[3]);
		lay_d[i2 + 3].set_val_dfd1(damp_mat_dv[6]);
		lay_d[i2 + 4].set_val_dfd1(damp_mat_dv[7]);
		lay_d[i2 + 5].set_val_dfd1(damp_mat_dv[9]);
		lay_d[i2 + 6].set_val_dfd1(damp_mat_dv[18]);
		lay_d[i2 + 7].set_val_dfd1(damp_mat_dv[19]);
		lay_d[i2 + 8].set_val_dfd1(damp_mat_dv[21]);

		layi++;
		i2 += 9;
	}
	return;
}

void Element::get_layer_angle_dfd1(vector<DiffDoub1>& lay_ang, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	int i2;
	int layi;
	int dv_ind;
	DiffDoub1 angle;
	DiffDoub1 dv_val;
	DiffDoub1 coef;
	string dv_cat;

	i2 = 0;
	layi = 0;
	Section& this_sec = sec_ar[sect_ptr];
	for (auto& lay : this_sec.layers) {
		angle.set_val(lay.angle);
		for (auto& dv : design_vars) {
			dv_ind = dv.int_dat;
			DesignVariable& this_dv = dv_ar[dv_ind];
			if (this_dv.layer == layi) {
				dv_cat = this_dv.category;
				if (dv_cat == "angle") {
					coef.set_val(dv.doub_dat);
					this_dv.get_value_dfd1(dv_val);
					dv_val.mult(coef);
					angle.add(dv_val);
				}
			}
		}
		lay_ang[layi].set_val_dfd1(angle);

		layi++;
	}

	return;
}

void Element::get_layer_th_exp_dfd1(vector<DiffDoub1>& lay_th_exp, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int layi;
	DiffDoub1 t_exp_dv[6];
	int dv_comp;
	DiffDoub1 tmp;
	DiffDoub1 dv_val;
	
	layi = 0;
	i2 = 0;
	Section& this_sec = sec_ar[sect_ptr];
	for (auto& lay : this_sec.layers) {
		Material& this_mat = mat_ar[lay.mat_ptr];
		for (i1 = 0; i1 < 6; i1++) {
			t_exp_dv[i1].set_val(this_mat.expansion[i1]);
		}
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "thermalExp" && this_dv.layer == layi) {
				dv_comp = this_dv.component - 1;
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(tmp);
				t_exp_dv[dv_comp].add(dv_val);
			}
		}
		lay_th_exp[i2].set_val_dfd1(t_exp_dv[0]);
		i2++;
		lay_th_exp[i2].set_val_dfd1(t_exp_dv[1]);
		i2++;
		lay_th_exp[i2].set_val_dfd1(t_exp_dv[3]);
		i2++;
		layi++;
	}
	return;
}

void Element::get_layer_einit_dfd1(vector<DiffDoub1>& lay_einit, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int layi;
	DiffDoub1 e0_dv[6];
	int dv_comp;
	DiffDoub1 tmp;
	DiffDoub1 dv_val;

	layi = 0;
	i2 = 0;
	Section& this_sec = sec_ar[sect_ptr];
	for (auto& lay : this_sec.layers) {
		for (i1 = 0; i1 < 6; i1++) {
			e0_dv[i1].set_val(0.0);
		}
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "initialStrain" && this_dv.layer == layi) {
				dv_comp = this_dv.component - 1;
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(tmp);
				e0_dv[dv_comp].add(dv_val);
			}
		}
		lay_einit[i2].set_val_dfd1(e0_dv[0]);
		i2++;
		lay_einit[i2].set_val_dfd1(e0_dv[1]);
		i2++;
		lay_einit[i2].set_val_dfd1(e0_dv[3]);
		i2++;
		layi++;
	}
	return;
}

void Element::get_layer_den_dfd1(vector<DiffDoub1>& layer_den, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int layi;
	double mat_den;
	DiffDoub1 den_dv;
	DiffDoub1 tmp;
	DiffDoub1 dv_val;

	layi = 0;
	Section& this_sec = sec_ar[sect_ptr];
	for (auto& lay : this_sec.layers) {
		mat_den = mat_ar[lay.mat_ptr].density;
		den_dv.set_val(mat_den);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "density" && this_dv.layer == layi) {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(tmp);
				den_dv.add(dv_val);
			}
		}
		layer_den[layi].set_val_dfd1(den_dv);
		layi++;
	}

	return;
}

void Element::get_layer_cond_dfd1(vector<DiffDoub1>& lay_cond, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int layi;
	DiffDoub1 cond_dv[6];
	DiffDoub1 tmp;
	DiffDoub1 dv_val;
	int dv_comp;

	i1 = 0;
	layi = 0;
	Section& this_sec = sec_ar[sect_ptr];
	for (auto& lay : this_sec.layers) {
		Material& this_mat = mat_ar[lay.mat_ptr];
		for (i2 = 0; i2 < 6; i2++) {
			cond_dv[i2].set_val(this_mat.conductivity[i2]);
		}
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "thermalCond" && this_dv.layer == layi) {
				dv_comp = this_dv.component - 1;
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(tmp);
				cond_dv[dv_comp].add(dv_val);
			}
		}
		lay_cond[i1].set_val_dfd1(cond_dv[0]);
		lay_cond[i1 + 1].set_val_dfd1(cond_dv[3]);
		lay_cond[i1 + 2].set_val_dfd1(cond_dv[4]);
		lay_cond[i1 + 3].set_val_dfd1(cond_dv[3]);
		lay_cond[i1 + 4].set_val_dfd1(cond_dv[1]);
		lay_cond[i1 + 5].set_val_dfd1(cond_dv[5]);
		lay_cond[i1 + 6].set_val_dfd1(cond_dv[4]);
		lay_cond[i1 + 7].set_val_dfd1(cond_dv[5]);
		lay_cond[i1 + 8].set_val_dfd1(cond_dv[2]);
		i1 += 9;
		layi++;
	}

	return;
}

void Element::get_layer_spec_heat_dfd1(vector<DiffDoub1>& lay_sh, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int layi;
	double mat_sh;
	DiffDoub1 sh_dv;
	DiffDoub1 tmp;
	DiffDoub1 dv_val;

	layi = 0;
	Section& this_sec = sec_ar[sect_ptr];
	for (auto& lay : this_sec.layers) {
		mat_sh = mat_ar[lay.mat_ptr].spec_heat;
		sh_dv.set_val(mat_sh);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "specHeat" && this_dv.layer == layi) {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(tmp);
				sh_dv.add(dv_val);
			}
		}
		lay_sh[layi].set_val_dfd1(sh_dv);
		layi++;
	}

	return;
}

void Element::transform_strain_dfd1(DiffDoub1 stn_new[], DiffDoub1 stn_orig[], DiffDoub1& angle) {
	DiffDoub1 angle_rad;
	DiffDoub1 a11;
	DiffDoub1 a12;
	DiffDoub1 a21;
	DiffDoub1 a22;
	DiffDoub1 tmp;
	DiffDoub1 te[9];


	angle_rad.set_val(r_pio180);
	angle_rad.mult(angle);
	a11.set_val_dfd1(angle_rad);
	a11.cs();
	a21.set_val_dfd1(angle_rad);
	a21.sn();
	a12.set_val_dfd1(a21);
	a12.neg();
	a22.set_val_dfd1(a11);

	te[0].set_val_dfd1(a11);
	te[0].sqr();
	te[1].set_val_dfd1(a12);
	te[1].sqr();
	te[2].set_val_dfd1(a11);
	te[2].mult(a12);
	te[3].set_val_dfd1(a21);
	te[3].sqr();
	te[4].set_val_dfd1(a22);
	te[4].sqr();
	te[5].set_val_dfd1(a22);
	te[5].mult(a21);
	te[6].set_val(2.0);
	te[6].mult(a11);
	te[6].mult(a21);
	te[7].set_val(2.0);
	te[7].mult(a12);
	te[7].mult(a22);
	te[8].set_val_dfd1(a11);
	te[8].mult(a22);
	tmp.set_val_dfd1(a12);
	tmp.mult(a21);
	te[8].add(tmp);

	mat_mul_ar_dfd1(stn_new, te, stn_orig, 3, 3, 1);

	return;
}

void Element::transform_q_dfd1(DiffDoub1 q_new[], DiffDoub1 q_orig[], DiffDoub1& angle) {
	int i1;
	DiffDoub1 angle_rad;
	DiffDoub1 a11;
	DiffDoub1 a12;
	DiffDoub1 a21;
	DiffDoub1 a22;
	DiffDoub1 tmp;
	DiffDoub1 coef;
	DiffDoub1 ts[9];
	DiffDoub1 te[9];
	DiffDoub1 te_inv[9];
	DiffDoub1 x_vec[3];
	DiffDoub1 b_vec[3];


	angle_rad.set_val(r_pio180);
	angle_rad.mult(angle);
	a11.set_val_dfd1(angle_rad);
	a11.cs();
	a21.set_val_dfd1(angle_rad);
	a21.sn();
	a12.set_val_dfd1(a21);
	a12.neg();
	a22.set_val_dfd1(a11);

	ts[0].set_val_dfd1(a11);
	ts[0].sqr();
	ts[1].set_val_dfd1(a12);
	ts[1].sqr();
	ts[2].set_val(2.0);
	ts[2].mult(a11);
	ts[2].mult(a12);
	ts[3].set_val_dfd1(a21);
	ts[3].sqr();
	ts[4].set_val_dfd1(a22);
	ts[4].sqr();
	ts[5].set_val(2.0);
	ts[5].mult(a22);
	ts[5].mult(a21);
	ts[6].set_val_dfd1(a11);
	ts[6].mult(a21);
	ts[7].set_val_dfd1(a12);
	ts[7].mult(a22);
	ts[8].set_val_dfd1(a11);
	ts[8].mult(a22);
	tmp.set_val_dfd1(a12);
	tmp.mult(a21);
	ts[8].add(tmp);

	for (i1 = 0; i1 < 9; i1++) {
		te[i1].set_val_dfd1(ts[i1]);
	}
	tmp.set_val(0.5);
	te[2].mult(tmp);
	te[5].mult(tmp);
	tmp.set_val(2.0);
	te[6].mult(tmp);
	te[7].mult(tmp);

	get_det_inv_ar_dfd1(coef, te_inv, te, 3, 0, x_vec, b_vec);

	mat_mul_ar_dfd1(te, q_orig, te_inv, 3, 3, 3);
	mat_mul_ar_dfd1(q_new, ts, te, 3, 3, 3);

	return;
}

void Element::get_solid_stiff_dfd1(vector<DiffDoub1>& cmat, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i4;
	int i5;
	int dv_ind;
	int dv_comp;
	DiffDoub1 modulus_dv[3];
	DiffDoub1 poisson_dv[3];
	DiffDoub1 shear_mod_dv[3];
	DiffDoub1 smat[36];
	DiffDoub1 ctmp[36];
	DiffDoub1 dv_val;
	DiffDoub1 coef;
	DiffDoub1 x_vec[6];
	DiffDoub1 b_vec[6];
	string dv_cat;


	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	if (this_mat.stiffness[0] > 0.0) {
		for (i1 = 0; i1 < 36; i1++) {
			cmat[i1].set_val(this_mat.stiffness[i1]);
		}
		for (auto& dv : design_vars) {
			dv_ind = dv.int_dat;
			DesignVariable& this_dv = dv_ar[dv_ind];
			if (this_dv.category == "stiffnessMat") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				cmat[dv_comp].add(dv_val);
			}
		}
		for (i1 = 1; i1 < 6; i1++) {
			i4 = 6 * i1;
			i5 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				cmat[i4].set_val_dfd1(cmat[i5]);
				i4++;
				i5 += 6;
			}
		}
	}
	else {
		for (i1 = 0; i1 < 3; i1++) {
			modulus_dv[i1].set_val(this_mat.modulus[i1]);
			poisson_dv[i1].set_val(this_mat.poisson_ratio[i1]);
			shear_mod_dv[i1].set_val(this_mat.shear_mod[i1]);
		}
		for (auto& dv : design_vars) {
			dv_ind = dv.int_dat;
			DesignVariable& this_dv = dv_ar[dv_ind];
			dv_cat = this_dv.category;
			if (dv_cat == "modulus") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				modulus_dv[dv_comp - 1].add(dv_val);
			}
			else if (dv_cat == "poissonRatio") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				poisson_dv[dv_comp - 1].add(dv_val);
			}
			else if (dv_cat == "shearModulus") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				shear_mod_dv[dv_comp - 1].add(dv_val);
			}
		}
		for (i1 = 0; i1 < 36; i1++) {
			smat[i1].set_val(0.0);
		}

		smat[0].set_val(1.0);
		smat[0].dvd(modulus_dv[0]);
		smat[1].set_val_dfd1(poisson_dv[0]);
		smat[1].neg();
		smat[1].dvd(modulus_dv[0]);
		smat[2].set_val_dfd1(poisson_dv[1]);
		smat[2].neg();
		smat[2].dvd(modulus_dv[0]);
		smat[6].set_val_dfd1(smat[1]);
		smat[7].set_val(1.0);
		smat[7].dvd(modulus_dv[1]);
		smat[8].set_val_dfd1(poisson_dv[2]);
		smat[8].neg();
		smat[8].dvd(modulus_dv[1]);
		smat[12].set_val_dfd1(smat[2]);
		smat[13].set_val_dfd1(smat[8]);
		smat[14].set_val(1.0);
		smat[14].dvd(modulus_dv[2]);
		smat[21].set_val(1.0);
		smat[21].dvd(shear_mod_dv[0]);
		smat[28].set_val(1.0);
		smat[28].dvd(shear_mod_dv[1]);
		smat[35].set_val(1.0);
		smat[35].dvd(shear_mod_dv[2]);

		get_det_inv_ar_dfd1(coef, ctmp, smat, 6, 0, x_vec, b_vec);
		ar_to_vec_dfd1(ctmp, cmat, 0, 36);
	}

	return;
}

void Element::get_abd_dfd1(vector<DiffDoub1>& cmat, vector<DiffDoub1>& lay_thk, vector<DiffDoub1>& lay_z, vector<DiffDoub1>& lay_q, vector<DiffDoub1>& lay_ang, vector<Section>& sec_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int i6;
	int num_lay;
	DiffDoub1 z_max;
	DiffDoub1 z_min;
	DiffDoub1 thk;
	DiffDoub1 qmat[9];
	DiffDoub1 tmp;
	DiffDoub1 tmp2;

	for (i1 = 0; i1 < 81; i1++) {
		cmat[i1].set_val(0.0);
	}

	num_lay = sec_ar[sect_ptr].layers.size();
	i2 = 0;
	for (i1 = 0; i1 < num_lay; i1++) {
		thk.set_val(0.5);
		thk.mult(lay_thk[i1]);
		z_min.set_val_dfd1(lay_z[i1]);
		z_min.sub(thk);
		z_max.set_val_dfd1(lay_z[i1]);
		z_max.add(thk);
		transform_q_dfd1(qmat, &lay_q[i2], lay_ang[i1]);
		
		// a matrix portion
		tmp.set_val_dfd1(z_max);
		tmp.sub(z_min);
		for (i3 = 0; i3 < 3; i3++) {
			i5 = 9 * i3;
			i6 = 3 * i3;
			for (i4 = 0; i4 < 3; i4++) {
				tmp2.set_val_dfd1(tmp);
				tmp2.mult(qmat[i6]);
				cmat[i5].add(tmp2);
				i5++;
				i6++;
			}
		}

		// b matrix portion
		tmp.set_val_dfd1(z_max);
		tmp.sqr();
		tmp2.set_val_dfd1(z_min);
		tmp2.sqr();
		tmp.sub(tmp2);
		tmp2.set_val(0.5);
		tmp.mult(tmp2);
		for (i3 = 0; i3 < 3; i3++) {
			i5 = 9 * i3 + 3;
			i6 = 3 * i3;
			for (i4 = 0; i4 < 3; i4++) {
				//i6 = i3*3 + i4;
				//i5 = i3*9 + (i4 + 3)
				tmp2.set_val_dfd1(tmp);
				tmp2.mult(qmat[i6]);
				cmat[i5].add(tmp2);
				i5++;
				i6++;
			}
		}

		// d matrix portion				
		tmp.set_val_dfd1(z_max);
		tmp.sqr();
		tmp.mult(z_max);
		tmp2.set_val_dfd1(z_min);
		tmp2.sqr();
		tmp2.mult(z_min);
		tmp.sub(tmp2);
		tmp2.set_val(r_1o3);
		tmp.mult(tmp2);
		for (i3 = 0; i3 < 3; i3++) {
			i5 = 9 * (i3 + 3) + 3;
			i6 = 3 * i3;
			for (i4 = 0; i4 < 3; i4++) {
				tmp2.set_val_dfd1(tmp);
				tmp2.mult(qmat[i6]);
				cmat[i5].add(tmp2);
				i5++;
				i6++;
			}
		}
		i2 += 9;
	}

	for (i1 = 1; i1 < 6; i1++) {
		i3 = 9 * i1;
		i4 = i1;
		for (i2 = 0; i2 < i1; i2++) {
			//i3 = i1*9 + i2
			//i4 = i2*9 + i1
			cmat[i3].set_val_dfd1(cmat[i4]);
			i3++;
			i4 += 9;
		}
	}

	tmp.set_val(1.0);
	tmp.mult(cmat[20]);
	cmat[60].set_val_dfd1(tmp);
	cmat[70].set_val_dfd1(tmp);
	cmat[80].set_val_dfd1(tmp);

	return;
}

void Element::get_beam_stiff_dfd1(vector<DiffDoub1>& cmat, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int dv_ind;
	int dv_comp;
	DiffDoub1 modulus_dv[3];
	DiffDoub1 shear_mod_dv[3];
	DiffDoub1 area_dv;
	DiffDoub1 idv[5];
	DiffDoub1 jdv;
	DiffDoub1 dv_val;
	DiffDoub1 coef;
	DiffDoub1 x_vec[6];
	DiffDoub1 b_vec[6];
	string dv_cat;

	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	if (this_sec.stiffness[0] > 0.0) {
		for (i1 = 0; i1 < 36; i1++) {
			cmat[i1].set_val(this_sec.stiffness[i1]);
		}
		for (auto& dv : design_vars) {
			dv_ind = dv.int_dat;
			DesignVariable& this_dv = dv_ar[dv_ind];
			if (this_dv.category == "stiffnessMat") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				cmat[dv_comp].add(dv_val);
			}
		}
		for (i1 = 1; i1 < 6; i1++) {
			i4 = 6 * i1;
			i5 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				cmat[i4].set_val_dfd1(cmat[i5]);
				i4++;
				i5 += 6;
			}
		}
	}
	else {
		for (i1 = 0; i1 < 3; i1++) {
			modulus_dv[i1].set_val(this_mat.modulus[i1]);
			shear_mod_dv[i1].set_val(this_mat.shear_mod[i1]);
		}
		area_dv.set_val(this_sec.area);
		for (i1 = 0; i1 < 5; i1++) {
			idv[i1].set_val(this_sec.area_moment[i1]);
		}
		jdv.set_val(this_sec.polar_moment);
		for (auto& dv : design_vars) {
			dv_ind = dv.int_dat;
			DesignVariable& this_dv = dv_ar[dv_ind];
			if (this_dv.category == "modulus") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				modulus_dv[dv_comp - 1].add(dv_val);
			}
			else if (this_dv.category == "shearModulus") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				shear_mod_dv[dv_comp - 1].add(dv_val);
			}
			else if (this_dv.category == "area") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				area_dv.add(dv_val);
			}
			else if (this_dv.category == "areaMoment") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				idv[dv_comp - 1].add(dv_val);
			}
			else if (this_dv.category == "polarMoment") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				jdv.add(dv_val);
			}
		}

		for (i1 = 0; i1 < 6; i1++) {
			i3 = 7 * i1;
			for (i2 = i1; i2 < 6; i2++) {
				cmat[i3].set_val(0.0);
				i3++;
			}
		}

		cmat[0].set_val_dfd1(modulus_dv[0]);
		cmat[0].mult(area_dv);
		cmat[4].set_val_dfd1(modulus_dv[0]);
		cmat[4].mult(idv[0]);
		cmat[5].set_val_dfd1(modulus_dv[0]);
		cmat[5].neg();
		cmat[5].mult(idv[1]);
		cmat[7].set_val_dfd1(shear_mod_dv[0]);
		cmat[7].mult(area_dv);
		cmat[9].set_val_dfd1(shear_mod_dv[0]);
		cmat[9].neg();
		cmat[9].mult(idv[0]);
		cmat[14].set_val_dfd1(shear_mod_dv[1]);
		cmat[14].mult(area_dv);
		cmat[15].set_val_dfd1(shear_mod_dv[2]);
		cmat[15].mult(idv[1]);
		cmat[21].set_val_dfd1(shear_mod_dv[0]);
		cmat[21].mult(jdv);
		cmat[28].set_val_dfd1(modulus_dv[0]);
		cmat[28].mult(idv[2]);
		cmat[29].set_val_dfd1(modulus_dv[0]);
		cmat[29].neg();
		cmat[29].mult(idv[4]);
		cmat[35].set_val_dfd1(modulus_dv[0]);
		cmat[35].mult(idv[3]);

		for (i1 = 1; i1 < 6; i1++) {
			i3 = 6 * i1;
			i4 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				cmat[i3].set_val_dfd1(cmat[i4]);
				i3++;
				i4 += 6;
			}
		}
	}

	return;
}

void Element::get_thermal_exp_dfd1(vector<DiffDoub1>& th_exp, vector<DiffDoub1>& einit, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	DesignVariable* this_dv;
	int dv_comp;
	DiffDoub1 dv_val;
	DiffDoub1 tmp;
	
	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	for (i1 = 0; i1 < 6; i1++) {
		th_exp[i1].set_val(this_mat.expansion[i1]);
		einit[i1].set_val(0.0);
	}

	for (auto& dv : design_vars) {
		DesignVariable& this_dv = dv_ar[dv.int_dat];
		if (this_dv.category == "thermalExp") {
			dv_comp = this_dv.component - 1;
			tmp.set_val(dv.doub_dat);
			this_dv.get_value_dfd1(dv_val);
			dv_val.mult(tmp);
			th_exp[dv_comp].add(dv_val);
		}
		else if (this_dv.category == "initialStrain") {
			dv_comp = this_dv.component - 1;
			tmp.set_val(dv.doub_dat);
			this_dv.get_value_dfd1(dv_val);
			dv_val.mult(tmp);
			einit[dv_comp].add(dv_val);
		}
	}

	return;
}

void Element::get_shell_exp_load_dfd1(vector<DiffDoub1>& exp_ld, vector<DiffDoub1>& e0_ld, vector<DiffDoub1>& lay_thk, vector<DiffDoub1>& lay_z, vector<DiffDoub1>& lay_q, vector<DiffDoub1>& lay_th_exp, vector<DiffDoub1>& lay_einit, vector<DiffDoub1>& lay_ang, vector<Section>& sec_ar) {
	int i1;
	int num_lay;
	int layi;
	int qi;
	int exi;
	DiffDoub1 sect_q[9];
	DiffDoub1 sect_te[3];
	DiffDoub1 sect_e0[3];
	DiffDoub1 qte_prod[3];
	DiffDoub1 qe0_prod[3];
	DiffDoub1 this_q[9];
	DiffDoub1 this_te[3];
	DiffDoub1 this_e0[3];
	DiffDoub1 z_min;
	DiffDoub1 z_max;
	DiffDoub1 tmp;
	DiffDoub1 tmp2;

	for (i1 = 0; i1 < 6; i1++) {
		exp_ld[i1].set_val(0.0);
		e0_ld[i1].set_val(0.0);
	}

	num_lay = sec_ar[sect_ptr].layers.size();
	qi = 0;
	exi = 0;
	for (layi = 0; layi < num_lay; layi++) {
		vec_to_ar_dfd1(this_q, lay_q, qi, qi + 9);
		transform_q_dfd1(sect_q, this_q, lay_ang[layi]);
		vec_to_ar_dfd1(this_te, lay_th_exp, exi, exi + 3);
		transform_strain_dfd1(sect_te, this_te, lay_ang[layi]);
		vec_to_ar_dfd1(this_e0, lay_einit, exi, exi + 3);
		transform_strain_dfd1(sect_e0, this_e0, lay_ang[layi]);
		mat_mul_ar_dfd1(qte_prod, sect_q, sect_te, 3, 3, 1);
		mat_mul_ar_dfd1(qe0_prod, sect_q, sect_e0, 3, 3, 1);
		
		for (i1 = 0; i1 < 3; i1++) {
			tmp.set_val_dfd1(qte_prod[i1]);
			tmp.mult(lay_thk[layi]);
			exp_ld[i1].add(tmp);
			tmp.set_val_dfd1(qe0_prod[i1]);
			tmp.mult(lay_thk[layi]);
			e0_ld[i1].add(tmp);
		}

		tmp.set_val(0.5);
		tmp.mult(lay_thk[layi]);
		z_min.set_val_dfd1(lay_z[layi]);
		z_min.sub(tmp);
		z_min.sqr();
		z_max.set_val_dfd1(lay_z[layi]);
		z_max.add(tmp);
		z_max.sqr();
		tmp.set_val(0.5);
		tmp2.set_val_dfd1(z_max);
		tmp2.sub(z_min);
		tmp.mult(tmp2); // tmp = 0.5*(z_max^2 - z_min^2)
		for (i1 = 0; i1 < 3; i1++) {
			tmp2.set_val_dfd1(qte_prod[i1]);
			tmp2.mult(tmp);
			exp_ld[i1 + 3].add(tmp2);
			tmp2.set_val_dfd1(qe0_prod[i1]);
			tmp2.mult(tmp);
			e0_ld[i1 + 3].add(tmp2);
		}
		qi += 9;
		exi += 3;
	}
	
	return;
}

void Element::get_beam_exp_load_dfd1(vector<DiffDoub1>& exp_ld, vector<DiffDoub1>& e0_ld, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	DiffDoub1 dv_val;
	string cat;
	string cat_list;
	DiffDoub1 tmp;
	int dv_comp;
	DiffDoub1 mod_dv[3];
	DiffDoub1 shr_mod_dv[3];
	DiffDoub1 te_coef_dv[6];
	DiffDoub1 e0_dv[6];
	DiffDoub1 area_dv;
	DiffDoub1 idv[5];
	DiffDoub1 qmat[9];
	DiffDoub1 qte[3];
	DiffDoub1 qe0[3];
	DiffDoub1 dedgu[18];
	DiffDoub1 tmp_exp[6];
	
	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	if (this_sec.exp_load_coef[0] > 0.0) {
		for (i1 = 0; i1 < 6; i1++) {
			exp_ld[i1].set_val(this_sec.exp_load_coef[i1]);
			e0_ld[i1].set_val(0.0);
		}
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			cat = this_dv.category;
			if (cat == "thermalExp") {
				dv_comp = this_dv.component - 1;
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(tmp);
				exp_ld[dv_comp].add(dv_val);
			}
			else if (cat == "initialStrain") {
				dv_comp = this_dv.component - 1;
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(tmp);
				e0_ld[dv_comp].add(dv_val);
			}
		}
	}
	else {
		for (i1 = 0; i1 < 3; i1++) {
			mod_dv[i1].set_val(this_mat.modulus[i1]);
			shr_mod_dv[i1].set_val(this_mat.shear_mod[i1]);
			te_coef_dv[i1].set_val(this_mat.expansion[i1]);
			te_coef_dv[i1 + 3].set_val(this_mat.expansion[i1]);
			e0_dv[i1].set_val(0.0);
			e0_dv[i1 + 3].set_val(0.0);
		}
		area_dv.set_val(this_sec.area);
		for (i1 = 0; i1 < 5; i1++) {
			idv[i1].set_val(this_sec.area_moment[i1]);
		}
		cat_list = "modulus shearModulus thermalExp initialStrain area areaMoment";
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			cat = this_dv.category;
			i2 = cat_list.find(cat);
			if (i2 > -1) {
				dv_comp = this_dv.component - 1;
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(tmp);
				if (cat == "modulus") {
					mod_dv[dv_comp].add(dv_val);
				}
				else if (cat == "shearModulus") {
					shr_mod_dv[dv_comp].add(dv_val);
				}
				else if (cat == "thermalExp") {
					te_coef_dv[dv_comp].add(dv_val);
				}
				else if (cat == "initialStrain") {
					e0_dv[dv_comp].add(dv_val);
				}
				else if (cat == "area") {
					area_dv.add(dv_val);
				}
				else if (cat == "areaMoment") {
					idv[dv_comp].add(dv_val);
				}
			}
		}
		for (i1 = 1; i1 < 8; i1++) {
			qmat[i1].set_val(0.0);
		}
		qmat[0].set_val_dfd1(mod_dv[0]);
		qmat[4].set_val_dfd1(shr_mod_dv[0]);
		qmat[8].set_val_dfd1(shr_mod_dv[1]);
		te_coef_dv[1].set_val_dfd1(te_coef_dv[3]);
		te_coef_dv[2].set_val_dfd1(te_coef_dv[4]);
		mat_mul_ar_dfd1(qte, qmat, te_coef_dv, 3, 3, 1);
		e0_dv[1].set_val_dfd1(e0_dv[3]);
		e0_dv[2].set_val_dfd1(e0_dv[4]);
		mat_mul_ar_dfd1(qe0, qmat, e0_dv, 3, 3, 1);
		for (i1 = 0; i1 < 18; i1++) {
			dedgu[i1].set_val(0.0);
		}
		dedgu[0].set_val_dfd1(area_dv);
		dedgu[4].set_val_dfd1(area_dv);
		dedgu[8].set_val_dfd1(area_dv);
		dedgu[10].set_val_dfd1(idv[0]);
		dedgu[10].neg();
		dedgu[11].set_val_dfd1(idv[1]);
		dedgu[12].set_val_dfd1(idv[0]);
		dedgu[15].set_val_dfd1(idv[1]);
		dedgu[15].neg();

		mat_mul_ar_dfd1(tmp_exp, dedgu, qte, 6, 3, 1);
		ar_to_vec_dfd1(tmp_exp, exp_ld, 0, 6);
		mat_mul_ar_dfd1(tmp_exp, dedgu, qe0, 6, 3, 1);
		ar_to_vec_dfd1(tmp_exp, e0_ld, 0, 6);
	}

	return;
}

void Element::get_density_dfd1(DiffDoub1& den, int layer, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int layi;
	string cat;
	int dv_lay;
	DiffDoub1 coef;
	DiffDoub1 dv_val;

	Section& this_sec = sec_ar[sect_ptr];
	if (type == 3 || type == 41) {
		Material& this_mat = mat_ar[this_sec.get_layer_mat_ptr(layer)];
		den.set_val(this_mat.density);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			cat = this_dv.category;
			dv_lay = this_dv.layer;
			if (cat == "density" && dv_lay == layer) {
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				den.add(dv_val);
			}
		}
	} else {
		Material& this_mat = mat_ar[this_sec.mat_ptr];
		den.set_val(this_mat.density);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			cat = this_dv.category;
			if (cat == "density") {
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				den.add(dv_val);
			}
		}
	}

	return;
}

void Element::get_shell_mass_dfd1(vector<DiffDoub1>& mmat, vector<DiffDoub1>& lay_thk, vector<DiffDoub1>& lay_z, vector<DiffDoub1>& lay_den, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int layi;
	DiffDoub1 tmp;
	DiffDoub1 tmp2;
	DiffDoub1 z_min;
	DiffDoub1 z_min2;
	DiffDoub1 z_max;
	DiffDoub1 z_max2;

	for (i1 = 0; i1 < 36; i1++) {
		mmat[i1].set_val(0.0);
	}

	i1 = sec_ar[sect_ptr].layers.size();
	for (layi = 0; layi < i1; layi++) {
		tmp.set_val_dfd1(lay_den[layi]);
		tmp.mult(lay_thk[layi]);
		mmat[0].add(tmp);
		mmat[7].add(tmp);
		mmat[14].add(tmp);

		tmp.set_val(0.5);
		tmp.mult(lay_thk[layi]);
		z_min.set_val_dfd1(lay_z[layi]);
		z_min.sub(tmp);
		z_min2.set_val_dfd1(z_min);
		z_min2.sqr();
		z_max.set_val_dfd1(lay_z[layi]);
		z_max.add(tmp);
		z_max2.set_val_dfd1(z_max);
		z_max2.sqr();
		tmp.set_val_dfd1(z_max2);
		tmp.sub(z_min2);  // tmp == z_max^2 - z_min^2
		tmp2.set_val(0.5);
		tmp2.mult(lay_den[layi]);
		tmp2.mult(tmp); // tmp2 = 0.5*rho*(z_max^2 - z_min^2)
		mmat[4].add(tmp2);
		mmat[24].add(tmp2);
		mmat[9].sub(tmp2);
		mmat[19].sub(tmp2);

		z_max2.mult(z_max);
		z_min2.mult(z_min);
		tmp.set_val_dfd1(z_max2);
		tmp.sub(z_min2); // tmp == z_max^3 - z_min^3
		tmp2.set_val(r_1o3);
		tmp2.mult(lay_den[layi]);
		tmp2.mult(tmp);
		mmat[21].add(tmp2);
		mmat[28].add(tmp2);
		mmat[35].add(tmp2);
	}

	return;
}

void Element::get_beam_mass_dfd1(vector<DiffDoub1>& mmat, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int dv_comp;

	DesignVariable* this_dv;
	DiffDoub1 dv_val;
	DiffDoub1 coef;

	string d_cat;
	Material* mat_pt;
	double* area_mom;
	DiffDoub1 den_dv;
	DiffDoub1 area_dv;
	DiffDoub1 idv[5];
	DiffDoub1 jdv;
	DiffDoub1 tmp;

	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	if (this_sec.mass[0] > 0.0) {
		for (i1 = 0; i1 < 36; i1++) {
			mmat[i1].set_val(this_sec.mass[i1]);
		}
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "massMat") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				mmat[dv_comp].add(dv_val);
			}
		}
		for (i1 = 1; i1 < 6; i1++) {
			i3 = 6 * i1; // lower tri term
			i4 = i1; // upper tri term
			for (i2 = 0; i2 < i1; i2++) {
				mmat[i3].set_val_dfd1(mmat[i4]);
				i3++;
				i4 += 6;
			}
		}
	}
	else {
		for (i1 = 0; i1 < 36; i1++) {
			mmat[i1].set_val(0.0);
		}
		for (i1 = 0; i1 < 5; i1++) {
			idv[i1].set_val(this_sec.area_moment[i1]);
		}
		area_dv.set_val(this_sec.area);
		jdv.set_val(this_sec.polar_moment);
		den_dv.set_val(this_mat.density);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			d_cat = this_dv.category;
			if (d_cat == "density") {
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				den_dv.add(dv_val);
			}
			else if (d_cat == "area") {
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				area_dv.add(dv_val);
			}
			else if (d_cat == "areaMoment") {
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				dv_comp = this_dv.component;
				idv[dv_comp - 1].add(dv_val);
			}
			else if (d_cat == "polarMoment") {
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				jdv.add(dv_val);
			}
		}
		tmp.set_val_dfd1(den_dv);
		tmp.mult(area_dv);
		mmat[0].set_val_dfd1(tmp);
		mmat[7].set_val_dfd1(tmp);
		mmat[14].set_val_dfd1(tmp);
		tmp.set_val_dfd1(den_dv);
		tmp.mult(idv[0]);
		mmat[4].set_val_dfd1(tmp);
		mmat[9].set_val_dfd1(tmp);
		mmat[9].neg();
		tmp.set_val_dfd1(den_dv);
		tmp.mult(idv[1]);
		mmat[5].set_val_dfd1(tmp);
		mmat[5].neg();
		mmat[15].set_val_dfd1(tmp);
		tmp.set_val_dfd1(den_dv);
		tmp.mult(jdv);
		mmat[21].set_val_dfd1(tmp);
		tmp.set_val_dfd1(den_dv);
		tmp.mult(idv[2]);
		mmat[28].set_val_dfd1(tmp);
		tmp.set_val_dfd1(den_dv);
		tmp.mult(idv[4]);
		mmat[29].set_val_dfd1(tmp);
		tmp.set_val_dfd1(den_dv);
		tmp.mult(idv[3]);
		mmat[35].set_val_dfd1(tmp);

		for (i1 = 3; i1 < 6; i1++) {
			i3 = 6 * i1; // lower tri term
			i4 = i1; // upper tri term
			for (i2 = 0; i2 < i1; i2++) {
				mmat[i3].set_val_dfd1(mmat[i4]);
				i3++;
				i4 += 6;
			}
		}
	}

	return;
}

void Element::get_solid_damp_dfd1(vector<DiffDoub1>& dmat, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	DesignVariable* this_dv;
	DiffDoub1 temp;
	DiffDoub1 dv_val;
	int dv_comp;

	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	for (i1 = 0; i1 < 36; i1++) {
		dmat[i1].set_val(this_mat.damping[i1]);
	}

	for (auto& dv : design_vars) {
		DesignVariable& this_dv = dv_ar[dv.int_dat];
		if (this_dv.category == "dampingMat") {
			this_dv.get_value_dfd1(dv_val);
			dv_comp = this_dv.component;
			temp.set_val(dv.doub_dat);
			dv_val.mult(temp);
			dmat[dv_comp].add(dv_val);
		}
	}

	for (i1 = 1; i1 < 6; i1++) {
		i3 = 6 * i1; // lower tri term
		i4 = i1; // upper tri term
		for (i2 = 0; i2 < i1; i2++) {
			dmat[i3].set_val_dfd1(dmat[i4]);
			i3++;
			i4 += 6;
		}
	}

	return;
}

void Element::get_shell_damp_dfd1(vector<DiffDoub1>& dmat, vector<DiffDoub1>& lay_thk, vector<DiffDoub1>& lay_z, vector<DiffDoub1>& lay_d, vector<DiffDoub1>& lay_ang, vector<Section>& sec_ar) {
	get_abd_dfd1(dmat, lay_thk, lay_z, lay_d, lay_ang, sec_ar);
	dmat[60].set_val(0.0);
	dmat[70].set_val(0.0);
	dmat[80].set_val(0.0);
	return;
}

void Element::get_beam_damp_dfd1(vector<DiffDoub1>& dmat, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int i3;
	int i4;
	int i5;
	int dv_ind;
	int dv_comp;
	DiffDoub1 dmat_dv[36];
	DiffDoub1 area_dv;
	DiffDoub1 idv[5];
	DiffDoub1 jdv;
	DiffDoub1 dv_val;
	DiffDoub1 coef;
	DiffDoub1 x_vec[6];
	DiffDoub1 b_vec[6];
	string dv_cat;

	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	if (this_sec.damping[0] > 0.0) {
		for (i1 = 0; i1 < 36; i1++) {
			dmat[i1].set_val(this_sec.damping[i1]);
		}
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "dampingMat") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				dmat[dv_comp].add(dv_val);
			}
		}
		for (i1 = 1; i1 < 6; i1++) {
			i4 = 6 * i1;
			i5 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				dmat[i4].set_val_dfd1(dmat[i5]);
				i4++;
				i5 += 6;
			}
		}
	}
	else {
		for (i1 = 0; i1 < 36; i1++) {
			dmat_dv[i1].set_val(this_mat.damping[i1]);
		}
		area_dv.set_val(this_sec.area);
		for (i1 = 0; i1 < 5; i1++) {
			idv[i1].set_val(this_sec.area_moment[i1]);
		}
		jdv.set_val(this_sec.polar_moment);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "dampingMat") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				dmat_dv[dv_comp].add(dv_val);
			}
			else if (this_dv.category == "area") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				area_dv.add(dv_val);
			}
			else if (this_dv.category == "areaMoment") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				idv[dv_comp - 1].add(dv_val);
			}
			else if (this_dv.category == "polarMoment") {
				dv_comp = this_dv.component;
				coef.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(coef);
				jdv.add(dv_val);
			}
		}

		for (i1 = 0; i1 < 6; i1++) {
			i3 = 7 * i1;
			for (i2 = i1; i2 < 6; i2++) {
				dmat[i3].set_val(0.0);
				i3++;
			}
		}

		dmat[0].set_val_dfd1(dmat_dv[0]);
		dmat[0].mult(area_dv);
		dmat[4].set_val_dfd1(dmat_dv[0]);
		dmat[4].mult(idv[0]);
		dmat[5].set_val_dfd1(dmat_dv[0]);
		dmat[5].neg();
		dmat[5].mult(idv[1]);
		dmat[7].set_val_dfd1(dmat_dv[21]);
		dmat[7].mult(area_dv);
		dmat[9].set_val_dfd1(dmat_dv[21]);
		dmat[9].neg();
		dmat[9].mult(idv[0]);
		dmat[14].set_val_dfd1(dmat_dv[28]); // 28
		dmat[14].mult(area_dv);
		dmat[15].set_val_dfd1(dmat_dv[35]); // 35
		dmat[15].mult(idv[1]);
		dmat[21].set_val_dfd1(dmat_dv[21]); // 21
		dmat[21].mult(jdv);
		dmat[28].set_val_dfd1(dmat_dv[0]);
		dmat[28].mult(idv[2]);
		dmat[29].set_val_dfd1(dmat_dv[0]);
		dmat[29].neg();
		dmat[29].mult(idv[4]);
		dmat[35].set_val_dfd1(dmat_dv[0]);
		dmat[35].mult(idv[3]);

		for (i1 = 1; i1 < 6; i1++) {
			i3 = 6 * i1;
			i4 = i1;
			for (i2 = 0; i2 < i1; i2++) {
				dmat[i3].set_val_dfd1(dmat[i4]);
				i3++;
				i4 += 6;
			}
		}
	}


	return;
}

void Element::get_conductivity_dfd1(vector<DiffDoub1>& t_cond, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	DiffDoub1 cond_dv[6];
	DesignVariable* this_dv;
	DiffDoub1 temp;
	DiffDoub1 dv_val;
	int dv_comp;

	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	for (i1 = 0; i1 < 6; i1++) {
		cond_dv[i1].set_val(this_mat.conductivity[i1]);
	}

	for (auto& dv : design_vars) {
		DesignVariable& this_dv = dv_ar[dv.int_dat];
		if (this_dv.category == "thermalCond") {
			dv_comp = this_dv.component - 1;
			temp.set_val(dv.doub_dat);
			this_dv.get_value_dfd1(dv_val);
			dv_val.mult(temp);
			cond_dv[dv_comp].add(dv_val);
		}
	}

	t_cond[0] = cond_dv[0];
	t_cond[1] = cond_dv[3];
	t_cond[2] = cond_dv[4];
	t_cond[3] = cond_dv[3];
	t_cond[4] = cond_dv[1];
	t_cond[5] = cond_dv[5];
	t_cond[6] = cond_dv[4];
	t_cond[7] = cond_dv[5];
	t_cond[8] = cond_dv[2];

	return;
}

void Element::get_shell_cond_dfd1(vector<DiffDoub1>& t_cond, vector<DiffDoub1>& lay_thk, vector<DiffDoub1>& lay_ang, vector<DiffDoub1>& lay_cond, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int i2;
	int num_lay = sec_ar[sect_ptr].layers.size();
	DiffDoub1 cond_dv[6];
	DiffDoub1 layer_mat[9];
	DiffDoub1 al_mat[9];
	DiffDoub1 al_t[9];
	DiffDoub1 tmp[9];
	DiffDoub1 this_cond[9];

	for (i1 = 0; i1 < 9; i1++) {
		t_cond[i1].set_val(0.0);
		al_mat[i1].set_val(0.0);
	}
	al_mat[8].set_val(1.0);

	for (i1 = 0; i1 < num_lay; i1++) {
		al_mat[0].set_val(r_pio180);
		al_mat[0].mult(lay_ang[i1]);
		al_mat[0].cs();
		al_mat[4].set_val_dfd1(al_mat[0]);
		al_mat[3].set_val(r_pio180);
		al_mat[3].mult(lay_ang[i1]);
		al_mat[3].sn();
		al_mat[1].set_val_dfd1(al_mat[3]);
		al_mat[1].neg();
		transpose_ar_dfd1(al_t, al_mat, 3, 3);
		i2 = 9 * i1;
		vec_to_ar_dfd1(this_cond, lay_cond, i2, i2 + 9);
		mat_mul_ar_dfd1(tmp, this_cond, al_t, 3, 3, 3);
		mat_mul_ar_dfd1(layer_mat, al_mat, tmp, 3, 3, 3);
		for (i2 = 0; i2 < 9; i2++) {
			layer_mat[i2].mult(lay_thk[i1]);
			t_cond[i2].add(layer_mat[i2]);
		}
	}

	return;
}

void Element::get_beam_cond_dfd1(vector<DiffDoub1>& t_cond, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	DesignVariable* this_dv;
	string cat;
	DiffDoub1 cond_dv;
	DiffDoub1 area_dv;
	DiffDoub1 tmp;
	DiffDoub1 dv_val;
	
	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	if (this_sec.conductivity > 0.0) {
		cond_dv.set_val(this_sec.conductivity);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "thermCond") {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(tmp);
				cond_dv.add(dv_val);
			}
		}
	}
	else {
		cond_dv.set_val(this_mat.conductivity[0]);
		area_dv.set_val(this_sec.area);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			cat = this_dv.category;
			if (cat == "thermCond") {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(tmp);
				cond_dv.add(dv_val);
			}
			else if (cat == "area") {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(tmp);
				area_dv.add(dv_val);
			}
		}
		cond_dv.mult(area_dv);
	}

	for (i1 = 1; i1 < 8; i1++) {
		t_cond[i1].set_val(0.0);
	}

	t_cond[0].set_val_dfd1(cond_dv);
	t_cond[4].set_val_dfd1(cond_dv);
	t_cond[8].set_val_dfd1(cond_dv);

	return;
}

void Element::get_specific_heat_dfd1(DiffDoub1& spec_heat, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	DiffDoub1 spec_heat_dv;
	DiffDoub1 temp;
	DiffDoub1 dv_val;
	int dv_comp;

	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	spec_heat.set_val(this_mat.spec_heat);

	for (auto& dv : design_vars) {
		DesignVariable& this_dv = dv_ar[dv.int_dat];
		if (this_dv.category == "specHeat") {
			temp.set_val(dv.doub_dat);
			this_dv.get_value_dfd1(dv_val);
			dv_val.mult(temp);
			spec_heat.add(dv_val);
		}
	}
	
	return;
}

void Element::get_shell_spec_heat_dfd1(DiffDoub1& spec_heat, vector<DiffDoub1>& lay_thk, vector<DiffDoub1>& lay_sh, vector<DiffDoub1>& lay_den, vector<Section>& sec_ar) {
	int i1;
	int num_lay;
	DiffDoub1 tmp;

	spec_heat.set_val(0.0);
	num_lay = sec_ar[sect_ptr].layers.size();
	for (i1 = 0; i1 < num_lay; i1++) {
		tmp.set_val_dfd1(lay_den[i1]);
		tmp.mult(lay_sh[i1]);
		tmp.mult(lay_thk[i1]);
		spec_heat.add(tmp);
	}

	return;
}

void Element::get_beam_spec_heat_dfd1(DiffDoub1& spec_heat, vector<Section>& sec_ar, vector<Material>& mat_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	string cat;
	DiffDoub1 density_dv;
	DiffDoub1 area_dv;
	DiffDoub1 tmp;
	DiffDoub1 dv_val;

	Section& this_sec = sec_ar[sect_ptr];
	Material& this_mat = mat_ar[this_sec.mat_ptr];
	double sec_sh = this_sec.spec_heat;
	double mat_sh = this_mat.spec_heat;
	double mat_den = this_mat.density;

	if (sec_sh > 0.0) {
		spec_heat.set_val(sec_sh);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			if (this_dv.category == "specHeat") {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(tmp);
				spec_heat.add(dv_val);
			}
		}
	}
	else {
		spec_heat.set_val(mat_sh);
		density_dv.set_val(mat_den);
		area_dv.set_val(this_sec.area);
		for (auto& dv : design_vars) {
			DesignVariable& this_dv = dv_ar[dv.int_dat];
			cat = this_dv.category;
			if (cat == "specHeatDV") {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(tmp);
				spec_heat.add(dv_val);
			}
			else if (cat == "density") {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(tmp);
				density_dv.add(dv_val);
			}
			else if (cat == "area") {
				tmp.set_val(dv.doub_dat);
				this_dv.get_value_dfd1(dv_val);
				dv_val.mult(tmp);
				area_dv.add(dv_val);
			}
		}
		spec_heat.mult(density_dv);
		spec_heat.mult(area_dv);
	}

	return;
}

void Element::get_nd_crds_dfd1(vector<DiffDoub1>& x_glob, vector<Node>& nd_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	DiffDoub1 nd_crd[3];
	
	for (i1 = 0; i1 < num_nds; i1++) {
		Node& this_nd = nd_ar[nodes[i1]];
		this_nd.get_crd_dfd1(nd_crd,dv_ar);
		x_glob[i1].set_val_dfd1(nd_crd[0]);
		x_glob[i1+num_nds].set_val_dfd1(nd_crd[1]);
		x_glob[i1+2*num_nds].set_val_dfd1(nd_crd[2]);
	}
	
	return;
}

void Element::get_loc_ori_dfd1(vector<DiffDoub1>& loc_ori, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	int i1;
	int dv_ind;
	DesignVariable *this_dvpt;
	string dv_cat;
	DiffDoub1 rot[3];
	DiffDoub1 dv_val;
	DiffDoub1 coef;
	DiffDoub1 ori_copy[9];
	DiffDoub1 tmp_ori[9];
	
	Section& this_sec = sec_ar[sect_ptr];
	for (i1 = 0; i1 < 9; i1++) {
		ori_copy[i1].set_val(this_sec.orientation[i1]);
	}
	
	rot[0].set_val(0.0);
	rot[1].set_val(0.0);
    rot[2].set_val(0.0);
	for (auto& dv : design_vars) {
		DesignVariable& this_dv = dv_ar[dv.int_dat];
		dv_cat = this_dv.category;
		if(dv_cat == "orientation") {
			i1 = this_dv.component - 1;
			coef.set_val(dv.doub_dat);
			this_dv.get_value_dfd1(dv_val);
			dv_val.mult(coef);
			rot[i1].add(dv_val);
		}
	}
	
	rotate_orient_dfd1(tmp_ori, ori_copy, rot);
	ar_to_vec_dfd1(tmp_ori, loc_ori, 0, 9);
	
	return;
}

void Element::correct_orient_dfd1(vector<DiffDoub1>& loc_ori, vector<DiffDoub1>& x_glob) {
	int i1;
	DiffDoub1 v1[3];
	DiffDoub1 v2[3];
	DiffDoub1 v3[3];
	DiffDoub1 rot[3];
	DiffDoub1 ori_copy[9];
	DiffDoub1 tmp_ori[9];
	
	DiffDoub1 dp;
	DiffDoub1 magv3;
	DiffDoub1 mag_cp;
	DiffDoub1 theta;
	DiffDoub1 tmp;
	
	for (i1 = 0; i1 < 9; i1++) {
		ori_copy[i1].set_val_dfd1(loc_ori[i1]);
	}
	
	if(type == 3 || type == 41) {
		if(type == 3) {
			v1[0].set_val_dfd1(x_glob[1]);
			v1[0].sub(x_glob[0]);
			v1[1].set_val_dfd1(x_glob[4]);
			v1[1].sub(x_glob[3]);
			v1[2].set_val_dfd1(x_glob[7]);
			v1[2].sub(x_glob[6]);
			
			v2[0].set_val_dfd1(x_glob[2]);
			v2[0].sub(x_glob[0]);
			v2[1].set_val_dfd1(x_glob[5]);
			v2[1].sub(x_glob[3]);
			v2[2].set_val_dfd1(x_glob[8]);
			v2[2].sub(x_glob[6]);		
		} else {
			v1[0].set_val_dfd1(x_glob[2]);
			v1[0].sub(x_glob[0]);
			v1[1].set_val_dfd1(x_glob[6]);
			v1[1].sub(x_glob[4]);
			v1[2].set_val_dfd1(x_glob[10]);
			v1[2].sub(x_glob[8]);
			
			v2[0].set_val_dfd1(x_glob[3]);
			v2[0].sub(x_glob[1]);
			v2[1].set_val_dfd1(x_glob[7]);
			v2[1].sub(x_glob[5]);
			v2[2].set_val_dfd1(x_glob[11]);
			v2[2].sub(x_glob[9]);	
		}
		cross_prod_dfd1(v3, v1, v2);
		
		dp.set_val_dfd1(v3[0]);
		dp.mult(loc_ori[6]);
		tmp.set_val_dfd1(v3[1]);
		tmp.mult(loc_ori[7]);
		dp.add(tmp);
		tmp.set_val_dfd1(v3[2]);
		tmp.mult(loc_ori[8]);
		dp.add(tmp);
		
		if(dp.val < 0.0) {
			v3[0].neg();
			v3[1].neg();
			v3[2].neg();
		}
		
		magv3.set_val_dfd1(v3[0]);
		magv3.sqr();
		tmp.set_val_dfd1(v3[1]);
		tmp.sqr();
		magv3.add(tmp);
		tmp.set_val_dfd1(v3[2]);
		tmp.sqr();
		magv3.add(tmp);
		magv3.sqt();
		
		cross_prod_dfd1(rot,&loc_ori[6],v3);
		mag_cp.set_val_dfd1(rot[0]);
		mag_cp.sqr();
		tmp.set_val_dfd1(rot[1]);
		tmp.sqr();
		mag_cp.add(tmp);
		tmp.set_val_dfd1(rot[2]);
		tmp.sqr();
		mag_cp.add(tmp);
		mag_cp.sqt();
		if (mag_cp.val < 1e-12) {
			return;
		}
		
		theta.set_val_dfd1(mag_cp);
		theta.dvd(magv3);
		theta.asn();
		
		tmp.set_val_dfd1(theta);
		tmp.dvd(mag_cp);
		
		rot[0].mult(tmp);
		rot[1].mult(tmp);
		rot[2].mult(tmp);
		
		rotate_orient_dfd1(tmp_ori, ori_copy, rot);
		ar_to_vec_dfd1(tmp_ori, loc_ori, 0, 9);
	}
	
	return;
}

void Element::get_frc_fld_const_dfd1(vector<DiffDoub1>& coef, vector<DiffDoub1>& exp, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	DiffDoub1 dv_val;
	DiffDoub1 tmp;
	string cat;

	Section& this_sec = sec_ar[sect_ptr];
	coef[0].set_val(this_sec.pot_coef);
	coef[1].set_val(this_sec.damp_coef);
	exp[0].set_val(this_sec.pot_exp);
	exp[1].set_val(this_sec.damp_exp);

	for (auto& dv : design_vars) {
		DesignVariable& this_dv = dv_ar[dv.int_dat];
		cat = this_dv.category;
		if (cat == "potFldCoef") {
			this_dv.get_value_dfd1(dv_val);
			tmp.set_val(dv.doub_dat);
			dv_val.mult(tmp);
			coef[0].add(dv_val);
		}
		else if (cat == "dampFldCoef") {
			this_dv.get_value_dfd1(dv_val);
			tmp.set_val(dv.doub_dat);
			dv_val.mult(tmp);
			coef[1].add(dv_val);
		}
	}

	return;
}

void Element::get_thrm_fld_const_dfd1(vector<DiffDoub1>& coef, DiffDoub1& ref_t, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	DiffDoub1 dv_val;
	DiffDoub1 tmp;
	string cat;

	Section& this_sec = sec_ar[sect_ptr];
	coef[0].set_val(this_sec.cond_coef);
	coef[1].set_val(this_sec.rad_coef);
	ref_t.set_val(this_sec.ref_temp);

	for (auto& dv : design_vars) {
		DesignVariable& this_dv = dv_ar[dv.int_dat];
		cat = this_dv.category;
		if (cat == "condCoef") {
			this_dv.get_value_dfd1(dv_val);
			tmp.set_val(dv.doub_dat);
			dv_val.mult(tmp);
			coef[0].add(dv_val);
		}
		else if (cat == "radCoef") {
			this_dv.get_value_dfd1(dv_val);
			tmp.set_val(dv.doub_dat);
			dv_val.mult(tmp);
			coef[1].add(dv_val);
		}
	}

	return;
}

void Element::get_mass_per_el_dfd1(DiffDoub1& mass_per_el, vector<Section>& sec_ar, vector<DesignVariable>& dv_ar) {
	DiffDoub1 dv_val;
	DiffDoub1 tmp;
	string cat;

	Section& this_sec = sec_ar[sect_ptr];
	mass_per_el.set_val(this_sec.mass_per_el);

	for(auto& dv : design_vars) {
		DesignVariable& this_dv = dv_ar[dv.int_dat];
		cat = this_dv.category;
		if (cat == "massPerEl") {
			this_dv.get_value_dfd1(dv_val);
			tmp.set_val(dv.doub_dat);
			dv_val.mult(tmp);
			mass_per_el.add(dv_val);
		}
	}

	return;
}

//end dup
 
//end skip 
 
 
 
