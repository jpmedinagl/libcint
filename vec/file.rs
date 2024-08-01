#![feature(autodiff)]

fn func(
    env: &mut Vec<f64>
) -> f64 {
    return env[1];
}

#[no_mangle]
#[autodiff(diff, Reverse, Const, Duplicated, Active)]
fn wrapper(
    env1: Vec<f64>,
    env2: Vec<f64>,
) -> f64 {
    let mut env = Vec::new();
    env.extend(env1);
    env.extend(env2);
    return func(env);
}

fn main() {
    let env1 = vec![1.0; 1];
    let env2 = vec![2.0; 1];
    let denv2 = vec![0.0; 1];

    let _ = diff(env1, env2, denv2, 1.0);
}