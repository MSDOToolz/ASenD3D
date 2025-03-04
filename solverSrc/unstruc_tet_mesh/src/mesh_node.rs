

#[derive(Clone)]
pub struct MeshNode {
    pub label : usize,
    pub coord : [f64; 3],
}

impl MeshNode {
    pub fn new() -> MeshNode {
        MeshNode {
            label : 0,
            coord : [0.0; 3],
        }
    }
}


