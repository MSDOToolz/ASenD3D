#ifndef objective
#define objective
#include "ListEntClass.h"
#include "SetClass.h"
#include "NodeClass.h"
#include "ElementClass.h"
#include "DesignVariableClass.h"

class ObjectiveTerm {
	public:
	    std::string category;
		std::string optr;
		double active_time[2];
		int component;
		int layer;
		double coef;
		double expnt;
		std::string el_set_name;
		int el_set_ptr;
		std::string nd_set_name;
		int nd_set_ptr;
		std::string tgt_tag;
		std::list<double> tgt_vals;
		double value;

		std::vector<double> q_vec;
		std::vector<double> el_vol_vec;
		std::vector<double> tgt_vec;
		std::vector<double> err_norm_vec;

		int q_len;

		SparseMat d_qd_u;
		SparseMat d_qd_v;
		SparseMat d_qd_a;
		SparseMat d_qd_t;
		SparseMat d_qd_tdot;
		SparseMat d_qd_d;
		SparseMat d_vd_d;
		
	    ObjectiveTerm();
		
		void set_active_time(double new_at[]);

		void allocate_obj(std::vector<Set>& nd_sets, std::vector<Set>& el_sets);

		void allocate_obj_grad();

		double get_power_norm();

		void d_power_normd_u(std::vector<double>& d_ld_u, std::vector<double>& d_ld_v, std::vector<double>& d_ld_a, std::vector<double>& d_ld_t, std::vector<double>& d_ld_tdot);

		void d_power_normd_d(std::vector<double>& d_ld_d);

		double get_vol_integral();

		void d_vol_integrald_u(std::vector<double>& d_ld_u, std::vector<double>& d_ld_v, std::vector<double>& d_ld_a, std::vector<double>& d_ld_t, std::vector<double>& d_ld_tdot);

		void d_vol_integrald_d(std::vector<double>& d_ld_d);

		double get_vol_average();

		void d_vol_averaged_u(std::vector<double>& d_ld_u, std::vector<double>& d_ld_v, std::vector<double>& d_ld_a, std::vector<double>& d_ld_t, std::vector<double>& d_ld_tdot);

		void d_vol_averaged_d(std::vector<double>& d_ld_d);

		void get_obj_val(double time, bool n_lgeom, std::vector<Node>& nd_ar, std::vector<Element>& el_ar, std::vector<Set>& nd_sets, std::vector<Set>& el_sets, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar, DiffDoub0StressPrereq& st_pre);

		void getd_ld_u(std::vector<double>& d_ld_u, std::vector<double>& d_ld_v, std::vector<double>& d_ld_a, std::vector<double>& d_ld_t, std::vector<double>& d_ld_tdot, double time, bool n_lgeom, std::vector<Node>& nd_ar,  std::vector<Element>& el_ar, std::vector<Set>& nd_sets, std::vector<Set>& el_sets, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar, DiffDoub0StressPrereq& st_pre);

		void getd_ld_d(std::vector<double>& d_ld_d, double time, bool n_lgeom, std::vector<Node>& nd_ar, std::vector<Element>& el_ar, std::vector<Set>& nd_sets, std::vector<Set>& el_sets, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar, DiffDoub1StressPrereq& st_pre);
};

class Objective {
	public:
		std::vector<ObjectiveTerm> terms;
		
	    Objective();

		void clear_values();

		void calculate_terms(double time, bool n_lgeom, std::vector<Node>& nd_ar,  std::vector<Element>& el_ar, std::vector<Set>& nd_sets, std::vector<Set>& el_sets, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar, DiffDoub0StressPrereq& st_pre);

		void calculated_ld_u(std::vector<double>& d_ld_u, std::vector<double>& d_ld_v, std::vector<double>& d_ld_a, std::vector<double>& d_ld_t, std::vector<double>& d_ld_tdot, double time, bool n_lgeom, std::vector<Node>& nd_ar, std::vector<Element>& el_ar, std::vector<Set>& nd_sets, std::vector<Set>& el_sets, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar, DiffDoub0StressPrereq& st_pre);

		void calculated_ld_d(std::vector<double>& d_ld_d, double time, bool n_lgeom, std::vector<Node>& nd_ar,  std::vector<Element>& el_ar, std::vector<Set>& nd_sets, std::vector<Set>& el_sets, std::vector<Section>& sec_ar, std::vector<Material>& mat_ar, std::vector<DesignVariable>& dv_ar, DiffDoub1StressPrereq& st_pre);

};

#endif
