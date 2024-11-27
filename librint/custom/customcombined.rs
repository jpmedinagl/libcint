#![feature(autodiff)]

use std::autodiff::autodiff;

#[inline(never)]
unsafe extern "C" fn square_(src: &f64, dest: &mut f64) {
    *dest = *src * *src;
}

#[no_mangle]
#[autodiff(dsquare, Reverse, Active, Duplicated)]
unsafe extern "C" fn square(x: f64) -> f64 {
    let mut y = 0.0;
    square_(&x, &mut y);
    y
}

static mut AUGMENT: i32 = 0;

unsafe extern "C" fn augment_square_(src: &f64, d_src: &f64, dest: &mut f64, d_dest: &mut f64) {
    unsafe {
        AUGMENT += 1;
    }
    // intentionally incorrect for debugging
    *dest = 7.0;
    *d_dest = 11.0;
}

static mut GRADIENT: i32 = 0;

unsafe extern "C" fn gradient_square_(src: &f64, d_src: &mut f64, dest: &f64, d_dest: &f64, tape: *mut std::ffi::c_void) {
    unsafe {
        GRADIENT += 1;
    }
    // intentionally incorrect for debugging
    *d_src = 13.0;
}


extern "C" {
    fn __enzyme_autodiff(arg1: *const (), arg2: f64) -> f64;
}


// fn dsquare(x: f64) -> f64 {
//     unsafe {
//         return __enzyme_autodiff(square as *const (), x);
//     }
// }


unsafe extern "C" fn f1_impl(p1: &f64, p2: &mut f64) { /* logic */ }
#[repr(C)]
struct register {
    f1: unsafe extern "C" fn(&f64, &mut f64),
    f2: unsafe extern "C" fn(&f64, &f64, &mut f64, &mut f64),
    f3: unsafe extern "C" fn(&f64, &mut f64, &f64, &f64, *mut std::ffi::c_void)
}

#[no_mangle] static __enzyme_register_gradient_square: register = register { 
    f1: square_,
    f2: augment_square_,
    f3: gradient_square_
};


fn main() {
    unsafe {
        let res: f64 = dsquare(3.0);
        println!("res={} augment={} gradient={}", res, AUGMENT, GRADIENT);
        assert!((res - 13.0).abs() < 1e-10);
        assert_eq!(AUGMENT, 1);
        assert_eq!(GRADIENT, 1);
    }
}