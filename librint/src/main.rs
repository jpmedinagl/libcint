use std::io;

use librint::cint_bas::CINTcgto_cart;
use librint::cint1e::cint1e_ovlp_cart;
use librint::cint1e::CINTOpt;

pub const ATM_SLOTS: usize = 6;
pub const BAS_SLOTS: usize = 8;

fn main() -> io::Result<()> {
    const natm: usize = 2;
    const nbas: usize = 2;
    
    // let mut atm: Box<[i32]> = vec![0; natm * ATM_SLOTS].into_boxed_slice();
    // let mut bas: Box<[i32]> = vec![0; nbas * BAS_SLOTS].into_boxed_slice();
    // let mut env: Box<[f32]> = vec![0.0; 10000].into_boxed_slice();

    // let mut file = File::open("/home/juanpmg/enzyme/h2/basis.txt")?;

    // read_arrays(&mut file, &mut atm, &mut bas, &mut env)?;

    let mut atm_arr: [i32; natm * ATM_SLOTS] = [1, 20, 1, 23, 0, 0, 1, 24, 1, 27, 0, 0];
    let mut bas_arr: [i32; nbas * BAS_SLOTS] = [0, 0, 3, 1, 0, 28, 31, 0, 1, 0, 3, 1, 0, 28, 31, 0];
    let mut env_arr: [f64; 34] = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
        -1.5117809, 0., 0., 0., 1.5117809, 0., 3.42525091, 0.62391373, 0.1688554, 0.98170675, 0.94946401, 0.29590645];

    let mut shls_arr: [i32; 4] = [0, 0, 0, 0];

	println!("buf");
    for i in 0..nbas {
        for j in 0..nbas {
            shls_arr[0] = i as i32;
            shls_arr[1] = j as i32;
            
            let di = CINTcgto_cart(i, &bas_arr);
            let dj = CINTcgto_cart(j, &bas_arr);

            let mut buf = vec![0.0; (di * dj) as usize];

            cint1e_ovlp_cart(&mut buf, &mut shls_arr, &mut atm_arr, natm as i32, &mut bas_arr, nbas as i32, &mut env_arr);

            for i in 0..(34) {
                println!("{} ", buf[i]);
            }
        }
    }
    Ok(())
}