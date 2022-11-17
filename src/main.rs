#![feature(generic_const_exprs)]

use crate::galois::{GFElement, EGF};

mod galois;

fn main() {
    //вариант 5
    const P: usize = 5;
    const M: usize = 3;

    type GF5 = GFElement<P>;

    println!("addition table");

    print!("+ ");
    (0..P).for_each(|v| print!("{v} "));
    println!();

    for y in 0..P {
        let y: GF5 = y.into();
        print!("{y} ");

        for x in 0..P {
            let x: GF5 = x.into();
            print!("{} ", x + y);
        }
        println!()
    }

    println!("\nmultiplication table");

    print!("* ");
    (0..P).for_each(|v| print!("{v} "));
    println!();

    for y in 0..P {
        let y: GF5 = y.into();
        print!("{y} ");

        for x in 0..P {
            let x: GF5 = x.into();
            print!("{} ", x * y);
        }
        println!()
    }

    let egf: EGF<P, M> = EGF::new([1, 0, 3, 2].map(GF5::from));

    let el = egf.construct_from_digits([1, 2, 0]);

    println!("{}", el.as_polynomial());

    println!("primitive: {:?}", egf.primitive().as_polynomial());

    println!(
        "{:?} = {:?} ^ {}",
        el.as_polynomial(),
        egf.primitive().as_polynomial(),
        el.primitive_power()
    );

    let table = egf.build_log_table();

    println!("log table: ");

    println!("Infinity -> 0");

    for (i, entry) in table.iter().enumerate() {
        println!(
            "{i: >3} -> {}",
            entry
                .as_ref()
                .map(ToString::to_string)
                .unwrap_or_else(|| "Infinity".to_string())
        )
    }
}
