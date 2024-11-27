use std::io;

fn another() {
    println!("within greet");
    return;
}

#[no_mangle]
fn greet() {
    println!("Hello, world!");
    another();
    return;
}

fn main() -> io::Result<()> {
    greet();
    Ok(())
}