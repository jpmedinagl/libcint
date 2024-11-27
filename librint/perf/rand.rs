use std::hint::black_box;

#[no_mangle]
pub fn sum() -> i32 {
    let mut sum = 0;
    for i in 0..1000 {
        sum += std::hint::black_box(1);
    }
    return sum;
}

fn main() {
    println!("Sum: {}", black_box(sum()));
}