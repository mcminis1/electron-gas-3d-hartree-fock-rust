#!/bin/bash

cargo +nightly build --target wasm32-unknown-unknown --release
cp target/wasm32-unknown-unknown/release/electron_gas_3d_hartree_fock_rust.wasm .