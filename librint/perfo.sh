perf record -F 40000 -e cycles:u -g ./target/release/pscf
perf script -i perf.data | inferno-collapse-perf | flamelens