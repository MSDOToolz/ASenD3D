#ifndef face
#define face
#include <vector>
#include "DiffDoubClass.h"
#include "ListEntClass.h"
#include "NodeClass.h"
#include "DesignVariableClass.h"

class Face {
	public:
	    int num_nds;
	    int loc_nodes[6];
		int glob_nodes[6];
		bool on_surf;

	    Face();
		
		void set_node(int place, int loc_nd, int glob_nd);
		
		void sorted_nodes(int srt_nds[]);
		
		int get_low_nd();

		//dup1

		void get_area_normal_dfd0(DiffDoub0& area, DiffDoub0 norm[], std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);

		//end dup
 
//skip 
 
//DiffDoub1 versions: 
		//dup1

		void get_area_normal_dfd1(DiffDoub1& area, DiffDoub1 norm[], std::vector<Node>& nd_ar, std::vector<DesignVariable>& dv_ar);

		//end dup
 
//end skip 
 
 
 
};

class FacePtList {
    public:
		std::list<int> fc_list;

		FacePtList();

		void add_face(int new_i);

		bool add_if_absent(int new_i, std::vector<Face>& glob_faces);
};

#endif