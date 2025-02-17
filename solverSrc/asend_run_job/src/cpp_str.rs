use crate::constants::*;

#[derive(Clone)]
pub struct CppStr {
	pub s : String,
}

impl CppStr {
	pub fn new() -> CppStr {
		CppStr {
			s : String::new(),
		}
	}

	pub fn from(in_st : & str) -> CppStr {
		CppStr {
			s : String::from(in_st),
		}
	}

	pub fn to_string(&mut self) -> String {
		self.s.clone()
	}

	pub fn len(&mut self) -> usize {
		self.s.len()
	}

	pub fn find(&mut self, substr :& str) -> usize {
		match self.s.find(substr) {
			None => max_int,
			Some(x) => x as usize,
		}
	}
	
	pub fn find_first_of(&mut self, substr :& str) -> usize {
		let mut ind : usize = 0usize;
		let mut res_i : usize;
		for c_str in self.s.chars() {
			res_i = match substr.find(c_str) {
				None => max_int,
				Some(x) => x as usize,
			};
			if (res_i < max_int) {
				return ind;
			}
			ind += 1usize;
		}
		max_int
	}
	
	pub fn find_first_not_of(&mut self, substr :& str) -> usize {
		let mut ind : usize = 0usize;
		let mut res_i : usize;
		for c_str in self.s.chars() {
			res_i = match substr.find(c_str) {
				None => max_int,
				Some(x) => x as usize,
			};
			if (res_i == max_int) {
				return ind;
			}
			ind += 1usize;
		}
		max_int
	}

	pub fn substr(&mut self, st_c : usize, st_len : usize) -> CppStr {
		if ((st_c + st_len) > self.s.len()) {
			let sub = match self.s.get(st_c..) {
				None => "".to_string(),
				Some(x) => x.to_string(),
			};
			return CppStr {
				s : sub.clone(),
			};
		}
		else {
			let sub = match self.s.get(st_c..(st_c+st_len)) {
				None => "".to_string(),
				Some(x) => x.to_string(),
			};
			return CppStr {
				s : sub.clone(),
			};
		}
	}
	
	pub fn is_int(&mut self) -> bool {
		match self.s.parse::<usize>() {
			Err(_why) => false,
			Ok(_x) => true,
		}
	}
	
	pub fn stoi(&mut self) -> usize {
		match self.s.parse::<usize>() {
			Err(_why) => 0usize,
			Ok(x) => x,
		}
	}
	
	pub fn is_doub(&mut self) -> bool {
		match self.s.parse::<f64>() {
			Err(_why) => false,
			Ok(_x) => true,
		}
	}
	
	pub fn stod(&mut self) -> f64 {
		match self.s.parse::<f64>() {
			Err(_why) => 0.0f64,
			Ok(x) => x,
		}
	}

}