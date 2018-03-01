use std::vec::Vec;
use std::collections::HashMap;
use std::f64::consts::PI;

struct Occupations {
    n_electrons: u32,
    k_pf: f64,
    ee_pf: f64,
    occ: Vec<(i32, i32, i32)>
}

fn build_occupations(r_s: f64, k2_max: f64) -> Occupations { 
    let k_limits = 2*k2_max.sqrt().ceil() as i32;
    let mut occ: Vec<(i32, i32, i32)> = Vec::new();
    let mut unocc: Vec<(i32, i32, i32)> = Vec::new();
    for x in -k_limits..k_limits {
        for y in -k_limits..k_limits {
            for z in -k_limits..k_limits {
                let kpoint:(i32, i32, i32) = (x,y,z);
                if ((x*x + y*y + z*z) as f64) < k2_max {
                    occ.push(kpoint);
                } else {
                    unocc.push(kpoint);
                }
            };
        };
    };

    let mut k_map :HashMap<(i32, i32, i32),bool> = HashMap::new();
    for &k_pnt in occ.iter() {
        k_map.insert(k_pnt,true);
    }
    for &k_pnt in unocc.iter() {
        k_map.insert(k_pnt,false);
    }

    let n_electrons = occ.len() as u32;
    let volume = 4.0/3.0*PI*(n_electrons as f64)*r_s.powi(3);
    let length = volume.powf(1.0/3.0);

    let k_pf = (8.0*PI)*length.powi(-2) as f64;
    let ee_pf = -2.0/(length*PI);

    Occupations{
        n_electrons: n_electrons,
        k_pf: k_pf/n_electrons as f64,
        ee_pf: ee_pf/n_electrons as f64,
        occ: occ
    }
}

fn get_energy(occ: &Occupations) -> f64 {
    let mut ke: f64 = 0.0;
    let mut ee: f64 = 0.0;
    for index_0 in 0..occ.n_electrons {
        let k_point = occ.occ[index_0 as usize];
        ke += (k_point.0.pow(2) + k_point.1.pow(2) + k_point.2.pow(2)) as f64;
        for index_1 in index_0+1..occ.n_electrons {
            let j_point = occ.occ[index_1 as usize];
            let dx2 = (k_point.0 - j_point.0).pow(2);
            let dy2 = (k_point.1 - j_point.1).pow(2);
            let dz2 = (k_point.2 - j_point.2).pow(2);
            ee += ((dx2+dy2+dz2) as f64).powi(-1);
        }
    }
    (ke*occ.k_pf+ee*occ.ee_pf) as f64
}

#[no_mangle]
pub extern fn get_heg_info(r_s: f64, k2_max: f64) -> f64 {
    let electrons = build_occupations(r_s,k2_max);
    get_energy(&electrons)
}
