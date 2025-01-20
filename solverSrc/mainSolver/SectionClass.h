#ifndef section
#define section
#include <string>
#include <vector>
#include <list>

class Material {
	public:
	    std::string name;
		double density;
		double modulus[3];
		double poisson_ratio[3];
		double shear_mod[3];
		double stiffness[36];
		double conductivity[6];
		double expansion[6];
		double spec_heat;
		double damping[36];
		double max_ten_stress[3];
		double max_comp_stress[3];
		double max_shear_stress[3];
		double max_ten_strain[3];
		double max_comp_strain[3];
		double max_shear_strain[3];
		double max_str_eng;
		double max_mises;
		
	    Material();

		void set_stiffness(int row, int col, double val);

		void set_damping(int row, int col, double val);
		
};

class Fluid {
public:
	std::string name;
	double viscosity;
	double ideal_gas;
	double therm_cond;
	double spec_heat;

	Fluid();

};

class Layer {
	public:
	    std::string mat_name;
		int mat_ptr;
		double thickness;
		double angle;
		
	    Layer();
};

class Section {
	public:
	    std::string type;
		std::string el_set_name;
		std::string mat_name;
		std::string fl_name;
		int mat_ptr;
		int fl_ptr;
		double orientation[9];
		double z_offset;
		std::list<Layer> layers;
		double area;
		double area_moment[5];
		double polar_moment;
		double stiffness[36];
		double mass[36];
		double damping[36];
		double exp_load_coef[6];
		double conductivity;
		double spec_heat;
		double mass_per_el;
		double pot_coef;
		double pot_exp;
		double damp_coef;
		double damp_exp;
		double cond_coef;
		double rad_coef;
		double den_vis_coef;
		double temp_vis_coef;
		double turb_vis_coef;
		double grad_vturb_coef;
		double diss_turb_coef;
		double enth_coef;
		double enth_exp;
		double pres_coef;
		double pres_exp;
		double ref_temp;
		double ref_den;
		double ref_turb_e;
		double ref_enth;
		
	    Section();
		
		void set_orientation(double new_ori[]);
		
		void set_area_moment(double new_i[]);
		
		void set_stiffness(int row, int col, double val);
		
		void set_mass(int row, int col, double val);

		void set_damping(int row, int col, double val);
		
		void set_exp_ld(double new_exp_ld[]);

		int get_layer_mat_ptr(int layi);
		
};

#endif
