
#[macro_use]
extern crate nom;
#[macro_use]
extern crate clap;
extern crate num;

pub mod field;
#[macro_use]
pub mod mono;
#[macro_use]
pub mod poly;
pub mod ideal;

pub use field::*;
pub use mono::*;
pub use poly::*;
pub use ideal::*;

use clap::AppSettings;

fn main() {

    let matches = clap_app! ( grobner =>
        //(about: "compute Gröbner bases over the rationals")
        (version: crate_version!())
        (author: crate_authors!())
        (long_about: "compute Gröbner bases for quotients of tensor algebras over simple fields")
        (after_help:
"NOTE:
    Polynomials should be entered between single quotes, for example:
        grobner basis 'y^2 - x' 'yx - y'

    This program was developed as part of 18.821.")
        (setting: AppSettings::SubcommandRequired)
        (@arg ITERATIONS: -i --iterations +takes_value default_value("5")
            {|s| s.parse::<usize>().map(|_| ()).map_err(|_| "expected nonnegative integer".into())}
            "iterations to run when computing a basis")
        (@arg CHAR: -f --field +takes_value default_value("0")
            {|s| s.parse::<usize>().map(|_| ()).map_err(|_| "expected nonnegative integer".into())}
            "coefficient field; 0 for rationals, p positive for mod p")
        (@subcommand basis =>
            (about: "compute a basis for Q<...>/(p1,...,pn)")
            (@arg IDEAL: +required ... {|s| Poly::from_str(&s).map(|_| ())}
                "polynomials which generate the ideal of relations")
        )
        (@subcommand reduce =>
            (about: "reduce a polynomial in Q<...>/(p1,...,pn)")
            (@arg POLY: -p --poly +required +takes_value {|s| Poly::from_str(&s).map(|_| ())}
                "the polynomial to reduce")
            (@arg IDEAL: +required ... {|s| Poly::from_str(&s).map(|_| ())}
                "polynomials which generate the ideal of relations")
        )
    ).get_matches();

    let iterations = matches.value_of("ITERATIONS").unwrap().parse().unwrap();

    let prime = matches.value_of("CHAR").unwrap().parse().unwrap();

    unsafe {
        field::set_prime(prime);
    }

    macro_rules! subcommands {
        ($new_prime:expr) => {
            if let Some(matches) = matches.subcommand_matches("basis") {
                let mut ideal = Ideal::new(matches.values_of("IDEAL")
                        .unwrap()
                        .map($new_prime));
                ideal.simplify();
                ideal.expand_overlaps(iterations);
                for r in &ideal.0 {
                    println!("{}", r.into_poly());
                }
            } else if let Some(matches) = matches.subcommand_matches("reduce") {
                let mut ideal = Ideal::new(matches.values_of("IDEAL")
                        .unwrap()
                        .map($new_prime));
                ideal.simplify();
                ideal.expand_overlaps(iterations);
                let poly = $new_prime(&matches.value_of("POLY").unwrap());
                println!("a = {}", ideal);
                println!("{} = {} (mod a)", poly, poly.reduce(&ideal));
            }
        };
    }

    if prime == 0 {
        subcommands!(|s: &str| Poly::from_str(&s).unwrap());
    } else {
        subcommands!(|s: &str| Poly::from_str(&s).unwrap().change_field::<Finite>());
    }


}
