#![feature(autodiff)]

#[no_mangle]
#[autodiff(diff, Reverse)]
fn eprintfunc() {
    eprintln!("eprintln error {}", 0);
}

fn main() {
    diff();
}