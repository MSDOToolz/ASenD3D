use crate::fluid_domain::*;

impl FluidDomain {
    pub fn update_reference(&mut self) {
        let mut i1 : usize;
        let mut i2 : usize;

        i1 = 0;
        for sd in self.sub_domains.iter_mut() {
            i2 = match self.cs_map.get(&sd.cell_set_name.s) {
                None => panic!("Error: cell set {} named in subdomain does not exist.", sd.cell_set_name.s),
                Some(x) => *x,
            };
            for ci in self.cell_sets[i2].labels.iter() {
                self.cells[*ci].sub_dom_pt = i1;
            }
            i2 = 0;
            for fl in self.fluids.iter() {
                if fl.name.s == sd.fluid_name.s {
                    sd.fluid_ptr = i2;
                }
                i2 += 1;
            }
            i1 += 1;
        }

        for cnst in self.constraints.const_vec.iter_mut() {
            for tm in cnst.terms.iter_mut() {
                tm.ns_ptr = match self.ns_map.get(&tm.node_set.s) {
                    None => panic!("Error: node set {} named for fluid constraint does not exist.", tm.node_set.s),
                    Some(x) => *x,
                };
            }
        }

        for ld in self.loads.iter_mut() {
            ld.set_pt = match self.cs_map.get(&ld.cell_set.s) {
                None => panic!("Error: cell set {} named in load does not exist.", ld.cell_set.s),
                Some(x) => *x,
            };
        }

    }

    pub fn reorder_nodes(&mut self, block_dim : usize) {
        let mut nd1 : usize;
        let mut nd2 : usize;

        let num_nodes = self.nodes.len();
        // build nodal connectivity

        let mut nodal_conn = vec![Set::new(); num_nodes];
        
        for c in self.cells.iter_mut() {
            for i in 0..4 {
                nd1 = c.nodes[i];
                for j in i+1..4 {
                    nd2 = c.nodes[j];
                    nodal_conn[nd1].add_if_absent(nd2);
                    nodal_conn[nd2].add_if_absent(nd1);
                    self.nodes[nd1].add_conn_nd(nd2);
                    self.nodes[nd2].add_conn_nd(nd1);
                }
                self.nodes[nd1].add_cell(c.label);
            }
        }


    }
}