#![feature(autodiff)]

use std::autodiff::autodiff;

// Function to compute the square of a number
fn square(out: &mut f64, input: f64) {
    *out = input * input;
}

// Custom derivative implementation using autodiff
#[autodiff(dsquare, Reverse, Duplicated, Const)]
fn square_ad(
    out: &mut f64,  // Output
    input: f64,     // Input
) {
    square(out, input);
}

fn main() {
    let mut out = 0.0;

    // Input value
    let input = 3.0;
    let mut dout = 0.0;

    // Call autodiff derivative function
    dsquare(&mut out, &mut dout, input, 1.0);

    // Debugging output
    println!("Derivative of square at {}: {} {}", input, out, dout);

    // Assertions to check correctness
    // assert!((out - 6.0).abs() < 1e-10); // Derivative of x^2 at 3 is 2 * 3 = 6
}
