use crate::mesher::*;
use crate::constants::*;
use crate::fmath::*;
use crate::utilities::*;
use crate::cpp_str::CppStr;
use crate::mesh_node::*;
use crate::mesh_element::*;
use crate::mesh_face::*;

use std::fs::File;
use std::io::Write;

impl Mesher {

    pub fn read_input(&mut self, file_name : &mut CppStr) {
        let mut file_line = CppStr::new();
        let mut headings  = vec![CppStr::new(); 4];
        let mut hd_ld_space : [usize; 4] = [0usize; 4];
        let mut data  = vec![CppStr::new(); 3];
        let mut data_len : usize = 0usize;
        //let mut self.new_nd : usize;
        
        self.nd_ct = 0;
        self.fc_ct = 0;
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
                if headings[0].s == "nodes" && data_len == 3 {
                    self.nd_ct += 1usize;
                }
                else if headings[0].s == "faces" && data_len == 3 {
                    self.fc_ct += 1usize;
                }
                else if headings[0].s == "globProjWt" && data_len == 1 {
                    self.glob_proj_wt = CppStr::stod(&mut data[0]);
                }
                else if headings[0].s == "maxNumEls" && data_len == 1 {
                    self.max_num_els = CppStr::stoi(&mut data[0]);
                }
            }
        }
        else {
            panic!("Error: Could not open input file {} for unstructured 3D mesh generation.", file_name.s);
        }
        
        let mut rad : usize =  1;
        while 12*rad*rad < self.nd_ct {
            rad += 1usize;
        }
        
        self.nd_cap = 15 * rad * rad * rad;
        self.el_cap = 6 * self.nd_cap;
        if self.max_num_els == MAX_INT {
            self.max_num_els = self.el_cap;
        }
        else if self.el_cap > self.max_num_els {
            self.el_cap = self.max_num_els;
        }
        self.fc_cap = 2 * self.el_cap;
        
        self.nodes = vec![MeshNode::new(); self.nd_cap];
        self.elements = vec![MeshElement::new(); self.el_cap];
        self.faces = vec![MeshFace::new(); self.fc_cap];
        
        self.nd_ct = 0;
        self.fc_ct = 0;
        
        if let Ok(lines) = read_lines(file_name.s.clone()) {
            for line in lines.map_while(Result::ok) {
                file_line.s = line;
                read_input_line(&mut file_line, &mut  headings, &mut  hd_ld_space, &mut  data, &mut  data_len);
                if headings[0].s == "nodes" && data_len == 3 {
                    self.nodes[self.nd_ct].coord[0] = CppStr::stod(&mut data[0]);
                    self.nodes[self.nd_ct].coord[1] = CppStr::stod(&mut data[1]);
                    self.nodes[self.nd_ct].coord[2] = CppStr::stod(&mut data[2]);
                    self.nd_ct += 1usize;
                }
                else if headings[0].s == "faces" && data_len == 3 {
                    self.faces[self.fc_ct].nodes[0] = CppStr::stoi(&mut data[0]);
                    self.faces[self.fc_ct].nodes[1] = CppStr::stoi(&mut data[1]);
                    self.faces[self.fc_ct].nodes[2] = CppStr::stoi(&mut data[2]);
                    self.fc_ct += 1usize;
                }
            }
        }
        
        return;
    }

    pub fn init_boundary_normals(&mut self) {
        let mut lst_len : usize;
        let mut cent : [f64; 3] = [0f64; 3];
        let mut srch_dir : usize;
        let mut x_range : [f64; 2] = [0f64; 2];
        let mut y_range : [f64; 2] = [0f64; 2];
        let mut z_range : [f64; 2] = [0f64; 2];
        let mut vec : [f64; 3] = [0f64; 3];
        let mut int_out : [f64; 3] = [0f64; 3];
        let mut dp : f64;
        let mut intersects : bool;
        let mut int_ct : i32;
        let mut fci2 : usize;
        for fci in 0..self.fc_ct {
            //let mut this_fc = &mut self.faces[fci];
            self.faces[fci].get_centroid(&mut cent,  &mut self.nodes);
            //let mut norm_dir = &mut this_fc.norm_dir;
            srch_dir = 0;
            if fabs(self.faces[fci].norm_dir[1]) > fabs(self.faces[fci].norm_dir[0]) {
                srch_dir = 1;
            }
            if fabs(self.faces[fci].norm_dir[2]) > fabs(self.faces[fci].norm_dir[1]) {
                srch_dir = 2;
            }
            if srch_dir == 0 {
                vec[0] = 1.0;
                vec[1] = 7.52893558402e-4;
                vec[2] = 5.39890329009e-4;
                x_range[0] = 1.0;
                x_range[1] = 0.0;
                y_range[0] = cent[1] - self.avg_proj;
                y_range[1] = cent[1] + self.avg_proj;
                z_range[0] = cent[2] - self.avg_proj;
                z_range[1] = cent[2] + self.avg_proj;
            }
            else if srch_dir == 1 {
                vec[0] = 7.52893558402e-4;
                vec[1] = 1.0;
                vec[2] = 5.39890329009e-4;
                x_range[0] = cent[0] - self.avg_proj;
                x_range[1] = cent[0] + self.avg_proj;
                y_range[0] = 1.0;
                y_range[1] = 0.0;
                z_range[0] = cent[2] - self.avg_proj;
                z_range[1] = cent[2] + self.avg_proj;
            }
            else {
                vec[0] = 7.52893558402e-4;
                vec[1] = 5.39890329009e-4;
                vec[2] = 1.0;
                x_range[0] = cent[0] - self.avg_proj;
                x_range[1] = cent[0] + self.avg_proj;
                y_range[0] = cent[1] - self.avg_proj;
                y_range[1] = cent[1] + self.avg_proj;
                z_range[0] = 1.0;
                z_range[1] = 0.0;
            }
            lst_len = self.face_grid.get_in_xyzrange(&mut self.grid_out1,  self.g_olen, &mut  x_range, &mut  y_range, &mut  z_range);
            int_ct = 1;
            for i1 in 0..lst_len {
                fci2 = self.grid_out1[i1];
                //let mut this_fc2 = &mut self.faces[fci2];
                if fci != fci2 {
                    intersects = self.faces[fci2].get_intersection(&mut int_out, &mut  cent, &mut  vec, &mut  self.nodes);
                    if intersects && int_out[0] > 0.0 {
                        int_ct  *=  -1;
                    }
                }
            }
            // if in_ct > 0, number of intersections is even, vec points out of the surface
            dp = vec[0] * self.faces[fci].norm_dir[0] + vec[1] * self.faces[fci].norm_dir[1] + vec[2] * self.faces[fci].norm_dir[2];
            if (int_ct > 0 && dp > 0.0) || (int_ct < 0 && dp < 0.0) {
                self.faces[fci].norm_dir[0]  *=  -1.0;
                self.faces[fci].norm_dir[1]  *=  -1.0;
                self.faces[fci].norm_dir[2]  *=  -1.0;
            }
            //
            // if cent[0] > 9.99 {
            //     dp = dp;
            // }
            //
        }
        return;
    }

    pub fn prep(&mut self) {
        self.num_bound_nds = self.nd_ct;
        
        //allocate the grid output list
        let num_faces : usize =  self.fc_ct;
        if num_faces > self.num_bound_nds {
            self.g_olen = 2 * num_faces;
        }
        else {
            self.g_olen = 2 * self.num_bound_nds;
        }
        self.grid_out1 = vec![0usize; self.g_olen];
        self.grid_out2 = vec![0usize; self.g_olen];
        
        //initialize grids
        let mut x_range : [f64; 2] = [ 1.0e+100,-1.0e+100 ];
        let mut y_range : [f64; 2] = [ 1.0e+100,-1.0e+100 ];
        let mut z_range : [f64; 2] = [ 1.0e+100,-1.0e+100 ];
        
        for i1 in 0..self.nd_ct {
            let this_nd = &self.nodes[i1];
            let crd = &this_nd.coord;
            if crd[0] < x_range[0] {
                x_range[0] = crd[0];
            }
            if crd[0] > x_range[1] {
                x_range[1] = crd[0];
            }
            if crd[1] < y_range[0] {
                y_range[0] = crd[1];
            }
            if crd[1] > y_range[1] {
                y_range[1] = crd[1];
            }
            if crd[2] < z_range[0] {
                z_range[0] = crd[2];
            }
            if crd[2] > z_range[1] {
                z_range[1] = crd[2];
            }
        }
        
        let mut spacing : f64;
        self.avg_proj = 0.0;
        self.max_proj = 0.0;
        self.max_edge_len = 0.0;
        for i1 in 0..self.fc_ct {
            //let mut this_fc = &self.faces[i1];
            spacing = self.faces[i1].proj_dist;
            self.avg_proj  +=  spacing;
            if spacing > self.max_proj {
                self.max_proj = spacing;
            }
            spacing = self.faces[i1].get_longest_edge_len(&mut self.nodes);
            if spacing > self.max_edge_len {
                self.max_edge_len = spacing;
            }
        }
        self.avg_proj  /=  num_faces as f64;
        spacing = self.avg_proj;
        self.node_grid.initialize(&mut x_range,  spacing, &mut  y_range,  spacing, &mut  z_range,  spacing);
        self.element_grid.initialize(&mut x_range,  spacing, &mut  y_range,  spacing, &mut  z_range,  spacing);
        self.face_grid.initialize(&mut x_range,  spacing, &mut  y_range,  spacing, &mut  z_range,  spacing);
        
        //add boundary self.nodes and self.faces to their grids
        let mut cent : [f64; 3] = [0f64; 3];
        for i1 in 0..self.nd_ct {
            let this_nd = &self.nodes[i1];
            self.node_grid.add_ent(i1, &this_nd.coord);
        }
        
        for i1 in 0..self.fc_ct {
            //let mut this_fc = &self.faces[i1];
            self.faces[i1].get_centroid(&mut cent,  &mut self.nodes);
            self.face_grid.add_ent(i1, &mut  cent);
        }
        
        // check boundary completeness
        
        let mut lst_len : usize;
        let mut lst_nd : [usize; 3] = [0usize; 3];
        let mut shared : [bool; 3] = [false; 3];
        let mut num_shared : usize;
        let mut neighb_ct : usize;
        let mut fi1 : usize;
        for i2 in 0..self.fc_ct {
            //let mut this_fc = &self.faces[i2];
            self.faces[i2].get_centroid(&mut cent,  &mut self.nodes);
            lst_len = self.face_grid.get_in_radius(&mut self.grid_out1,  self.g_olen, &mut  cent,  1.01*self.max_edge_len);
            neighb_ct = 0;
            for i1 in 0..lst_len {
                fi1 = self.grid_out1[i1];
                //let mut lst_fc = &self.faces[self.grid_out1[i1]];
                //fc_nds = lst_fc->get_nd_pt();
                //cout << fc_nds[0]->get_label() << ", " << fc_nds[1]->get_label() << ", " << fc_nds[2]->get_label() << endl;
                num_shared = self.faces[fi1].get_shared_nodes(&mut lst_nd, &mut  shared,  i2, &self.faces);
                if num_shared == 2 {
                    neighb_ct += 1usize;
                }
            }
            if neighb_ct != 3 {
                //cout << "neighbCt " << neighb_ct;
                //fc_nds = this_fc->get_nd_pt();
                //cout << fc_nds[0]->get_label() << ", " << fc_nds[1]->get_label() << ", " << fc_nds[2]->get_label() << endl;
                panic!("Error: Invalid boundary surface input to 3D unstructured mesh generator.\nBoundary must be a completely closed surface of triangular faces.");
            }
        }
        
        // check face normal directions
        
        self.init_boundary_normals();
        
        return;
    }

    pub fn check_new_el(&mut self) -> bool {
        let mut i3 : usize;
        let nei : usize = self.el_ct;
        let mut n1 : usize;
        let mut n2 : usize;
        let mut lst_len : usize;
        let mut cent : [f64; 3] = [0f64; 3];
        let el_edges : [usize; 12] = [ 0,1,0,2,0,3,1,2,1,3,2,3 ];
        let face_edges : [usize; 6] = [ 0,1,0,2,1,2 ];
        
        
        let mut vec : [f64; 3] = [0f64; 3];
        let mut out_p : [f64; 3] = [0f64; 3];
        let mut intersects : bool;
        let el_nds = &self.elements[nei].nodes;
        
        
        let mut fc_nd_out : [usize; 3] = [0usize; 3];
        let mut shared : [bool; 3] = [false; 3];
        
        self.elements[nei].get_centroid(&mut cent,  &mut self.nodes);
        lst_len = self.face_grid.get_in_radius(&mut self.grid_out2,  self.g_olen, &mut  cent,  1.01 * self.max_edge_len);
        for i1 in 0..lst_len {
            let this_fc = &self.faces[self.grid_out2[i1]];
            
            for i2 in 0..4 {
                // any face of new element already exists and is closed
                i3 = this_fc.get_shared_nodes(&mut fc_nd_out, &mut  shared,  i2, &self.new_el_fcs);
                let fc_els = &this_fc.elements;
                if i3 == 3 && fc_els[1] < MAX_INT {
                    return  false;
                }
                // any edge of new element intersects existing edge
                intersects = this_fc.edges_intersect(i2,  1.0e-6 * self.avg_proj, &mut  self.nodes, &self.new_el_fcs);
                if intersects {
                    return  false;
                }
            }
            
            // any edge of new element intersects existing face
            i3 = 0;
            for _i2 in 0..6 {
                n1 = el_edges[i3];
                n2 = el_edges[i3 + 1];
                let pt = &self.nodes[el_nds[n1]].coord;
                let pt2 = &self.nodes[el_nds[n2]].coord;
                vec[0] = pt2[0] - pt[0];
                vec[1] = pt2[1] - pt[1];
                vec[2] = pt2[2] - pt[2];
                intersects = this_fc.get_intersection(&mut out_p, pt, &mut  vec, &self.nodes);
                if intersects && out_p[0] < 0.9999999999 && out_p[0] > 0.0000000001 {
                    return  false;
                }
                i3  +=  2;
            }
            
            // any edge of existing face intersects face of new element
            
            let fc_nds = &this_fc.nodes;
            i3 = 0;
            for _i2 in 0..3 {
                n1 = face_edges[i3];
                n2 = face_edges[i3 + 1];
                let pt = &self.nodes[fc_nds[n1]].coord;
                let pt2 = &self.nodes[fc_nds[n2]].coord;
                vec[0] = pt2[0] - pt[0];
                vec[1] = pt2[1] - pt[1];
                vec[2] = pt2[2] - pt[2];
                for i4 in 0..4 {
                    intersects = self.new_el_fcs[i4].get_intersection(&mut out_p, pt, &mut  vec, &self.nodes);
                    if intersects && out_p[0] < 0.9999999999 && out_p[0] > 0.0000000001 {
                        return  false;
                    }
                }
                i3  +=  2;
            }
        }
        
        // any node of new element in an existing element
        
        lst_len = self.element_grid.get_in_radius(&mut self.grid_out2,  self.g_olen, &mut  cent,  1.01 * self.max_edge_len);
        for i1 in 0..lst_len {
            let this_el = &self.elements[self.grid_out2[i1]];
            for i2 in 0..4 {
                let pt = & self.nodes[el_nds[i2]].coord;
                intersects = this_el.point_in(pt, &self.nodes);
                if intersects {
                    return  false;
                }
            }
        }
        
        // any existing node in the new element
        
        lst_len = self.node_grid.get_in_radius(&mut self.grid_out2,  self.g_olen, &mut  cent,  1.01 * self.max_edge_len);
        for i1 in 0..lst_len {
            let this_nd = &self.nodes[self.grid_out2[i1]];
            let pt = &this_nd.coord;
            intersects = self.elements[nei].point_in(pt, &self.nodes);
            if intersects {
                return  false;
            }
        }
        
        return  true;
    }

    pub fn add_face_if_absent(&mut self, new_el : usize) -> bool {
        let mut cent : [f64; 3] = [0f64; 3];
        let new_face = &self.faces[self.fc_ct];
        new_face.get_centroid(&mut cent,  &mut  self.nodes);
        
        let mut this_nds : [usize; 3] = [0usize; 3];
        
        let mut shared : [bool; 3] = [false; 3];
        let mut num_shared : usize;
        
        let lst_len : usize =  self.face_grid.get_in_radius(&mut self.grid_out2,  self.g_olen, &mut  cent,  0.1*self.avg_proj);
        let mut fi1 : usize;
        for i1 in 0..lst_len {
            fi1 = self.grid_out2[i1];
            let this_fc = &self.faces[fi1];
            num_shared = this_fc.get_shared_nodes(&mut this_nds, &mut  shared,  self.fc_ct,  &self.faces);
            if num_shared == 3 {
                self.faces[fi1].elements[1] = new_el;
                return  false;
            }
        }
        
        self.face_grid.add_ent(self.fc_ct, &mut  cent);
        self.fc_ct += 1usize;
        
        return  true;
    }

    pub fn adopt_connected_nd(&mut self, fc_i : usize, tgt_pt : &mut [f64], srch_rad : f64) -> bool {
        //let mut this_fc = &self.faces[fc_i];
        let mut face_nds : [usize; 3] = [0usize; 3];
        let mut shared : [bool; 3] = [false; 3];
        let mut num_shared : usize;
        let mut un_shared : usize;
        
        let mut d_vec : [f64; 3] = [0f64; 3];
        let mut dist : f64;
        let mut dp : f64;
        let mut cent : [f64; 3] = [0f64; 3];
        //let mut new_el_nds = &mut self.elements[self.el_ct].nodes;
        let ec : usize = self.el_ct;
        let mut el_check : bool;
        let mut _fc_added : bool;
        
        
        self.new_el_fcs[0].copy_data(&mut self.faces[fc_i]);
        
        //let mut norm_dir = &self.faces[fc_i].norm_dir;
        //let mut this_fc_nds = &self.faces[fc_i].nodes;
        //let mut this_fc_els = &self.faces[fc_i].elements;
        self.faces[fc_i].get_centroid(&mut cent,  &mut  self.nodes);
        let lst_len : usize =  self.face_grid.get_in_radius(&mut self.grid_out1,  self.g_olen, tgt_pt,  1.01 * self.max_edge_len);
        for i1 in 0..lst_len {
            let list_fc = &self.faces[self.grid_out1[i1]];
            num_shared = list_fc.get_shared_nodes(&mut face_nds, &mut  shared,  fc_i, &self.faces);
            if num_shared == 2 {
                un_shared = MAX_INT;
                for i2 in 0..3 {
                    if !shared[i2] {
                        un_shared = face_nds[i2];
                    }
                }
                let crd = &self.nodes[un_shared].coord;
                d_vec[0] = crd[0] - tgt_pt[0];
                d_vec[1] = crd[1] - tgt_pt[1];
                d_vec[2] = crd[2] - tgt_pt[2];
                dist = sqrt(d_vec[0] * d_vec[0] + d_vec[1] * d_vec[1] + d_vec[2] * d_vec[2]);
                if dist < srch_rad {
                    d_vec[0] = crd[0] - cent[0];
                    d_vec[1] = crd[1] - cent[1];
                    d_vec[2] = crd[2] - cent[2];
                    dp = d_vec[0] * self.faces[fc_i].norm_dir[0] + d_vec[1] * self.faces[fc_i].norm_dir[1] + d_vec[2] * self.faces[fc_i].norm_dir[2];
                    if dp > 1.0e-3*self.avg_proj {
                        self.elements[ec].nodes[0] = self.faces[fc_i].nodes[0];
                        self.elements[ec].nodes[1] = self.faces[fc_i].nodes[1];
                        self.elements[ec].nodes[2] = self.faces[fc_i].nodes[2];
                        self.elements[ec].nodes[3] = un_shared;
                        
                        let mut new_fc = &mut self.new_el_fcs[1];
                        new_fc.nodes[0] = self.faces[fc_i].nodes[0];
                        new_fc.nodes[1] = self.faces[fc_i].nodes[1];
                        new_fc.nodes[2] = un_shared;
                        new_fc.elements[0] = self.el_ct;
                        new_fc.elements[1] = MAX_INT;
                        new_fc.init_norm_dir(&mut self.nodes);
                        
                        new_fc = &mut self.new_el_fcs[2];
                        new_fc.nodes[0] = self.faces[fc_i].nodes[1];
                        new_fc.nodes[1] = self.faces[fc_i].nodes[2];
                        new_fc.nodes[2] = un_shared;
                        new_fc.elements[0] = self.el_ct;
                        new_fc.elements[1] = MAX_INT;
                        new_fc.init_norm_dir(&mut self.nodes);
                        
                        new_fc = &mut self.new_el_fcs[3];
                        new_fc.nodes[0] = self.faces[fc_i].nodes[2];
                        new_fc.nodes[1] = self.faces[fc_i].nodes[0];
                        new_fc.nodes[2] = un_shared;
                        new_fc.elements[0] = self.el_ct;
                        new_fc.elements[1] = MAX_INT;
                        new_fc.init_norm_dir(&mut self.nodes);
                        
                        el_check = self.check_new_el();
                        if el_check {
                            self.elements[self.el_ct].get_centroid(&mut cent,  &mut  self.nodes);
                            self.element_grid.add_ent(self.el_ct, &mut  cent);
                            
                            self.faces[fc_i].elements[1] = self.el_ct;
                            
                            self.new_el_fcs[1].norm_dir_from_el_cent(&mut cent, &mut  self.nodes);
                            self.faces[self.fc_ct].copy_data(&mut self.new_el_fcs[1]);
                            _fc_added = self.add_face_if_absent(self.el_ct);
                            
                            self.new_el_fcs[2].norm_dir_from_el_cent(&mut cent, &mut  self.nodes);
                            self.faces[self.fc_ct].copy_data(&mut self.new_el_fcs[2]);
                            _fc_added = self.add_face_if_absent(self.el_ct);
                            
                            self.new_el_fcs[3].norm_dir_from_el_cent(&mut cent, &mut  self.nodes);
                            self.faces[self.fc_ct].copy_data(&mut self.new_el_fcs[3]);
                            _fc_added = self.add_face_if_absent(self.el_ct);
                            
                            self.el_ct += 1usize;
                            return  true;
                        }
                    }
                }
            }
        }
        
        return  false;
    }

    pub fn adopt_any_nd(&mut self, fc_i : usize, tgt_pt : &mut [f64], srch_rad : f64) -> bool {
        //let mut this_fc = &self.faces[fc_i];
                
        let mut ndi : usize;
                
        let mut d_vec : [f64; 3] = [0f64; 3];
        let mut dist : f64;
        let mut dp : f64;
        let mut cent : [f64; 3] = [0f64; 3];
        //let mut new_el_nds = &mut self.elements[self.el_ct].nodes;
        let ec : usize = self.el_ct;
        let mut el_check : bool;
        let mut _fc_added : bool;
        
        self.new_el_fcs[0].copy_data(&self.faces[fc_i]);
        
        //let mut norm_dir = &this_fc.norm_dir;
        //let mut this_fc_nds = &this_fc.nodes;
        //let mut this_fc_els = &this_fc.elements;
        self.faces[fc_i].get_centroid(&mut cent,  &mut  self.nodes);
        let lst_len : usize =  self.node_grid.get_in_radius(&mut self.grid_out1, self.g_olen, tgt_pt, srch_rad);
        for i1 in 0..lst_len {
            //list_nd = self.grid_out1[i1]->get_pt(list_nd);
            ndi = self.grid_out1[i1];
            let list_nd = &self.nodes[ndi];
            if ndi != self.faces[fc_i].nodes[0] && ndi != self.faces[fc_i].nodes[1] && ndi != self.faces[fc_i].nodes[2] {
                let crd = &list_nd.coord;
                d_vec[0] = crd[0] - tgt_pt[0];
                d_vec[1] = crd[1] - tgt_pt[1];
                d_vec[2] = crd[2] - tgt_pt[2];
                dist = sqrt(d_vec[0] * d_vec[0] + d_vec[1] * d_vec[1] + d_vec[2] * d_vec[2]);
                if dist < srch_rad {
                    d_vec[0] = crd[0] - cent[0];
                    d_vec[1] = crd[1] - cent[1];
                    d_vec[2] = crd[2] - cent[2];
                    dp = d_vec[0] * self.faces[fc_i].norm_dir[0] + d_vec[1] * self.faces[fc_i].norm_dir[1] + d_vec[2] * self.faces[fc_i].norm_dir[2];
                    if dp > 1.0e-3*self.avg_proj {
                        self.elements[ec].nodes[0] = self.faces[fc_i].nodes[0];
                        self.elements[ec].nodes[1] = self.faces[fc_i].nodes[1];
                        self.elements[ec].nodes[2] = self.faces[fc_i].nodes[2];
                        self.elements[ec].nodes[3] = ndi;
                        
                        let new_fc = &mut self.new_el_fcs[1];
                        new_fc.nodes[0] = self.faces[fc_i].nodes[0];
                        new_fc.nodes[1] = self.faces[fc_i].nodes[1];
                        new_fc.nodes[2] = ndi;
                        new_fc.elements[0] = self.el_ct;
                        new_fc.elements[1] = MAX_INT;
                        new_fc.init_norm_dir(&mut self.nodes);
                        
                        let new_fc = &mut self.new_el_fcs[2];
                        new_fc.nodes[0] = self.faces[fc_i].nodes[1];
                        new_fc.nodes[1] = self.faces[fc_i].nodes[2];
                        new_fc.nodes[2] = ndi;
                        new_fc.elements[0] = self.el_ct;
                        new_fc.elements[1] = MAX_INT;
                        new_fc.init_norm_dir(&mut self.nodes);
                        
                        let new_fc = &mut self.new_el_fcs[3];
                        new_fc.nodes[0] = self.faces[fc_i].nodes[2];
                        new_fc.nodes[1] = self.faces[fc_i].nodes[0];
                        new_fc.nodes[2] = ndi;
                        new_fc.elements[0] = self.el_ct;
                        new_fc.elements[1] = MAX_INT;
                        new_fc.init_norm_dir(&mut self.nodes);
                        
                        el_check = self.check_new_el();
                        if el_check {
                            //self.elements.add_ent(self.new_el);
                            self.elements[self.el_ct].get_centroid(&mut cent,  &mut  self.nodes);
                            self.element_grid.add_ent(self.el_ct, &mut  cent);
                            
                            self.faces[fc_i].elements[1] = self.el_ct;
                            
                            self.new_el_fcs[1].norm_dir_from_el_cent(&mut cent, &mut  self.nodes);
                            self.faces[self.fc_ct].copy_data(&mut self.new_el_fcs[1]);
                            _fc_added = self.add_face_if_absent(self.el_ct);
                            
                            self.new_el_fcs[2].norm_dir_from_el_cent(&mut cent, &mut  self.nodes);
                            self.faces[self.fc_ct].copy_data(&mut self.new_el_fcs[2]);
                            _fc_added = self.add_face_if_absent(self.el_ct);
                            
                            self.new_el_fcs[3].norm_dir_from_el_cent(&mut cent, &mut  self.nodes);
                            self.faces[self.fc_ct].copy_data(&mut self.new_el_fcs[3]);
                            _fc_added = self.add_face_if_absent(self.el_ct);
                            
                            self.el_ct += 1usize;
                            
                            return  true;
                        }
                    }
                }
            }
        }
        
        return  false;
    }

    pub fn create_new_nd(&mut self, fc_i : usize, tgt_pt : &mut [f64]) -> bool {
        self.nodes[self.nd_ct].coord[0] = tgt_pt[0];
        self.nodes[self.nd_ct].coord[1] = tgt_pt[1];
        self.nodes[self.nd_ct].coord[2] = tgt_pt[2];
        //let mut this_fc = &self.faces[fc_i];
        let el_check : bool;
        let mut _fc_added : bool;
        let mut cent : [f64; 3] = [0f64; 3];
        
        //let mut this_fc_nds = &this_fc.nodes;
        //let mut this_fc_els = &this_fc.elements;
        //let mut new_el_nds = &self.elements[self.el_ct].nodes;
        let ec : usize = self.el_ct;
        self.elements[ec].nodes[0] = self.faces[fc_i].nodes[0];
        self.elements[ec].nodes[1] = self.faces[fc_i].nodes[1];
        self.elements[ec].nodes[2] = self.faces[fc_i].nodes[2];
        self.elements[ec].nodes[3] = self.nd_ct;
        self.elements[self.el_ct].get_centroid(&mut cent,  &mut  self.nodes);
        
        self.new_el_fcs[0].copy_data(&self.faces[fc_i]);
        
        let mut new_fc = &mut self.new_el_fcs[1];
        new_fc.nodes[0] = self.faces[fc_i].nodes[0];
        new_fc.nodes[1] = self.faces[fc_i].nodes[1];
        new_fc.nodes[2] = self.nd_ct;
        new_fc.elements[0] = self.el_ct;
        new_fc.elements[1] = MAX_INT;
        new_fc.init_norm_dir(&mut self.nodes);
        new_fc.norm_dir_from_el_cent(&mut cent, &mut  self.nodes);
        
        new_fc = &mut self.new_el_fcs[2];
        new_fc.nodes[0] = self.faces[fc_i].nodes[1];
        new_fc.nodes[1] = self.faces[fc_i].nodes[2];
        new_fc.nodes[2] = self.nd_ct;
        new_fc.elements[0] = self.el_ct;
        new_fc.elements[1] = MAX_INT;
        new_fc.init_norm_dir(&mut self.nodes);
        new_fc.norm_dir_from_el_cent(&mut cent, &mut  self.nodes);
        
        new_fc = &mut self.new_el_fcs[3];
        new_fc.nodes[0] = self.faces[fc_i].nodes[2];
        new_fc.nodes[1] = self.faces[fc_i].nodes[0];
        new_fc.nodes[2] = self.nd_ct;
        new_fc.elements[0] = self.el_ct;
        new_fc.elements[1] = MAX_INT;
        new_fc.init_norm_dir(&mut self.nodes);
        new_fc.norm_dir_from_el_cent(&mut cent, &mut  self.nodes);
        
        el_check = self.check_new_el();
        if el_check {
            self.node_grid.add_ent(self.nd_ct, tgt_pt);
            
            self.element_grid.add_ent(self.el_ct, &mut  cent);
            
            self.faces[fc_i].elements[1] = self.el_ct;
            
            self.faces[self.fc_ct].copy_data(&mut self.new_el_fcs[1]);
            _fc_added = self.add_face_if_absent(self.el_ct);
            
            self.faces[self.fc_ct].copy_data(&mut self.new_el_fcs[2]);
            _fc_added = self.add_face_if_absent(self.el_ct);
            
            self.faces[self.fc_ct].copy_data(&mut self.new_el_fcs[3]);
            _fc_added = self.add_face_if_absent(self.el_ct);
            
            self.nd_ct += 1usize;
            self.el_ct += 1usize;
            
            return  true;
        }
        
        return  false;
    }

    pub fn generate_mesh(&mut self) -> bool {
        let mut fc_i : usize;
        
        let mut fc_cent : [f64; 3] = [0f64; 3];
        let mut fc_proj : f64;
        
        let mut proj : f64;
        let mut tgt_pt : [f64; 3] = [0f64; 3];
        let mut srch_rad : f64;
        let mut edge_closed : bool;
        
        let mut el_added : bool =  true;
        while el_added {
            el_added = false;
            fc_i = 0;
            while fc_i < self.fc_ct && self.el_ct < self.max_num_els {
                //let mut this_fc = &self.faces[fc_i];
                //let mut fc_els = &this_fc.elements;
                if self.faces[fc_i].elements[1] == MAX_INT {
                    self.faces[fc_i].get_centroid(&mut fc_cent,  &mut  self.nodes);
                    fc_proj = self.faces[fc_i].proj_dist;
                    proj = self.glob_proj_wt * self.avg_proj + (1.0 - self.glob_proj_wt) * fc_proj;
                    //let mut fc_norm = &self.faces[fc_i].norm_dir;
                    tgt_pt[0] = fc_cent[0] + 0.5 * proj * self.faces[fc_i].norm_dir[0];
                    tgt_pt[1] = fc_cent[1] + 0.5 * proj * self.faces[fc_i].norm_dir[1];
                    tgt_pt[2] = fc_cent[2] + 0.5 * proj * self.faces[fc_i].norm_dir[2];
                    srch_rad = 0.75 * proj;
                    edge_closed = self.adopt_connected_nd(fc_i, &mut  tgt_pt,  srch_rad);
                    if !edge_closed {
                        edge_closed = self.adopt_any_nd(fc_i, &mut  tgt_pt,  srch_rad);
                    }
                    if !edge_closed {
                        tgt_pt[0] = fc_cent[0] + proj * self.faces[fc_i].norm_dir[0];
                        tgt_pt[1] = fc_cent[1] + proj * self.faces[fc_i].norm_dir[1];
                        tgt_pt[2] = fc_cent[2] + proj * self.faces[fc_i].norm_dir[2];
                        edge_closed = self.create_new_nd(fc_i, &mut  tgt_pt);
                    }
                    if !edge_closed {
                        tgt_pt[0] = fc_cent[0] + 0.5 * proj * self.faces[fc_i].norm_dir[0];
                        tgt_pt[1] = fc_cent[1] + 0.5 * proj * self.faces[fc_i].norm_dir[1];
                        tgt_pt[2] = fc_cent[2] + 0.5 * proj * self.faces[fc_i].norm_dir[2];
                        edge_closed = self.create_new_nd(fc_i, &mut  tgt_pt);
                    }
                    if edge_closed {
                        el_added = true;
                    }
                }
                fc_i += 1usize;
            }
        }
        
        if self.el_ct >= self.max_num_els {
            return  false;
        }
        
        return  true;
    }

    pub fn distribute_nodes(&mut self) {
        let mut i1 : usize;
        let mut i3 : usize;
        let mut i4 : usize;
        let mut i5 : usize;
        let mut i6 : usize;
        let mut i7 : usize;
        let mut i8 : usize;
        let mut i9 : usize;
        let dim : usize =  self.nd_ct * 3;
        
        let mut e_vol : f64;
        let mut avg_wt : f64;
        
        let mut res : f64;
        let mut alpha : f64;
        let mut beta : f64;
        let mut r_next : f64;
        let mut dp : f64;
        
        let mut el_mat = vec![0f64; 144];
        let mut x_vec = vec![0f64; dim];
        let mut h_vec = vec![0f64; dim];
        let mut z_vec = vec![0f64; dim];
        let mut g_vec = vec![0f64; dim];
        let mut w_vec = vec![0f64; dim];
        let mut d_mat = vec![0f64; dim];
        let mut p_mat = vec![0f64; dim];
        let mut p_inv = vec![0f64; dim];
        
        i1 = self.el_ct;
        let mut el_wt = vec![0f64; i1];
        
        // make static element matrix
        for i1 in 0..144 {
            el_mat[i1] = 0.0;
        }
        
        i3 = 0;
        for i1 in 0..12 {
            if i3 < 9 {
                i4 = 0;
            }
            else {
                i4 = i3 - 9;
            }
            i5 = i1 * 12;
            i6 = (i1 + 1) * 12;
            while i4 <= (i3 + 9) {
                if i4 >= i5 && i4 < i6 {
                    el_mat[i4] = -1.0;
                }
                i4  +=  3;
            }
            el_mat[i3] = 3.0;
            i3  +=  13;
        }
        
        // determine element weights
        i1 = 0;
        avg_wt = 0.0;
        for i2 in 0..self.el_ct {
            let this_el = &self.elements[i2];
            e_vol = this_el.get_volume(&mut self.nodes);
            el_wt[i1] = e_vol;
            avg_wt  +=  e_vol;
            i1 += 1usize;
        }
        avg_wt  /=  i1 as f64;
        for i2 in 0..i1 {
            el_wt[i2]  /=  avg_wt;
        }
        
        // initialize matrices
        for i1 in 0..dim {
            d_mat[i1] = 0.0;
            p_mat[i1] = 30.0;
            g_vec[i1] = 0.0;
        }
        
        //this_nd = self.nodes.get_first();
        for i1 in 0..self.num_bound_nds {
            let this_nd = &self.nodes[i1];
            let crd = &this_nd.coord;
            for i3 in 0..3 {
                i4 = i1 * 3 + i3;
                d_mat[i4] = 100000.0;
                p_mat[i4]  +=  100000.0;
                g_vec[i4] = -100000.0 * crd[i3];
            }
        }
        
        for i1 in 0..dim {
            p_inv[i1] = 1.0 / p_mat[i1];
        }
        
        // intialize vectors
        
        res = 0.0;
        for i1 in 0..dim {
            x_vec[i1] = 0.0;
            w_vec[i1] = p_inv[i1] * g_vec[i1];
            h_vec[i1] = -w_vec[i1];
            res  +=  w_vec[i1] * g_vec[i1];
        }
        
        i1 = 0;
        while i1 < dim && res > 1.0e-12 {
            for i2 in 0..dim {
                z_vec[i2] = 0.0;
            }
            for i2 in 0..self.el_ct {
                let this_el = &self.elements[i2];
                let el_nds = &this_el.nodes;
                i7 = 0;//index in el_mat
                for i3 in 0..4 {
                    for i4 in 0..3 {
                        i8 = 3*el_nds[i3] + i4;// global row
                        for i5 in 0..4 {
                            for i6 in 0..3 {
                                i9 = 3 * el_nds[i5] + i6;//global col
                                z_vec[i8]  +=  el_wt[i2] * el_mat[i7] * h_vec[i9];
                                i7 += 1usize;
                            }
                        }
                    }
                }
            }
            for i2 in 0..dim {
                z_vec[i2]  +=  d_mat[i2] * h_vec[i2];
            }
            dp = 0.0;
            for i2 in 0..dim {
                dp  +=  z_vec[i2] * h_vec[i2];
            }
            alpha = res / dp;
            r_next = 0.0;
            for i2 in 0..dim {
                x_vec[i2]  +=  alpha * h_vec[i2];
                g_vec[i2]  +=  alpha * z_vec[i2];
                w_vec[i2] = p_inv[i2] * g_vec[i2];
                r_next  +=  g_vec[i2] * w_vec[i2];
            }
            beta = r_next / res;
            for i2 in 0..dim {
                h_vec[i2] = -w_vec[i2] + beta * h_vec[i2];
            }
            res = r_next;
            i1 += 1usize;
        }
        
        for i1 in self.num_bound_nds..self.nd_ct {
            let this_nd = &mut self.nodes[i1];
            let crd = &mut this_nd.coord;
            for i3 in 0..3 {
                i4 = 3 * i1 + i3;
                crd[i3] = x_vec[i4];
            }
        }
        
        return;
    }

    pub fn write_output(&mut self, file_name : &mut CppStr) {
        
        let mut out_file = match File::create(&file_name.s) {
            Err(_why) => panic!("could not open file {}", file_name.s),
            Ok(file) => file,
        };
        
        let _ = out_file.write(format!("{}", "nodes:\n").as_bytes());
        for i1 in 0..self.nd_ct {
            let this_nd = &self.nodes[i1];
            let crd = &this_nd.coord;
            let _ = out_file.write(format!("{}{}{}{}{}{}{}", "  - [" , crd[0] , ", " , crd[1] , ", " , crd[2] , "]\n").as_bytes());
        }
        
        let _ = out_file.write(format!("{}", "elements:\n").as_bytes());
        for i1 in 0..self.el_ct {
            let this_el = &self.elements[i1];
            let el_nds = &this_el.nodes;
            let _ = out_file.write(format!("{}{}", "  - [" , el_nds[0]).as_bytes());
            for i1 in 1..4 {
                let _ = out_file.write(format!("{}{}", ", " , el_nds[i1]).as_bytes());
            }
            let _ = out_file.write(format!("{}", "]\n").as_bytes());
        }
        
        return;
    }

    pub fn print_current_mesh(&mut self) {
        
        println!("{}", "Current Mesh:" );
        println!("{}", "Nodes:" );
        
        //this_nd = self.nodes.get_first();
        for i1 in 0..self.nd_ct {
            let this_nd = &self.nodes[i1];
            let crd = &this_nd.coord;
            println!("{}{}{}{}{}", crd[0] , ", " , crd[1] , ", " , crd[2] );
        }
        println!("{}", "Elements:" );
        //this_el = self.elements.get_first();
        for i1 in 0..self.el_ct {
            let this_el = &self.elements[i1];
            let el_nds = &this_el.nodes;
            print!("{}{}", el_nds[0] , ", ");
            print!("{}{}", el_nds[1] , ", ");
            print!("{}{}", el_nds[2] , ", ");
            println!("{}", el_nds[3] );
        }
        
        println!("{}", "Faces:" );
        //this_fc = self.faces.get_first();
        for i1 in 0..self.fc_ct {
            let this_fc = &self.faces[i1];
            print!("{}", "  nodes: ");
            let el_nds = &this_fc.nodes;
            print!("{}{}", el_nds[0] , ", ");
            print!("{}{}", el_nds[1] , ", ");
            print!("{}{}", el_nds[2] , ", ");
            let fc_els = &this_fc.elements;
            print!("{}", "element pt: ");
            println!("{}{}{}", fc_els[0] , ", " , fc_els[1] );
        }
        
        return;
    }

}


