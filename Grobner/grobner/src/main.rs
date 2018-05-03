
use std::collections::BTreeMap;

#[macro_use]
extern crate nom;
#[macro_use]
extern crate clap;
extern crate num;

use clap::AppSettings;

use num::traits::*;
use num::rational::Rational;

use std::mem;
use std::ops::{Add, AddAssign, SubAssign, MulAssign};
use std::fmt::{Result as FResult, Display, Formatter};

macro_rules! __monomial_inner {
    ($var:ident) => {stringify!($var).as_bytes()[0]}
}

macro_rules! monomial {
    [] => {$crate::Monomial::one()};
    [$($var:ident)*] => {{
        $crate::Monomial::new(&[$(__monomial_inner!($var),)*])
    }}
}

macro_rules! __poly_coeff {
    () => {
        ::num::rational::Rational::new(1, 1);
    };
    ($a:expr) => {
        ::num::rational::Rational::new($a, 1);
    };
    ($a:expr, $b:expr) => {
        ::num::rational::Rational::new($a, $b);
    }
}

macro_rules! poly {
    [$(($($c:expr)*) $($var:ident)* $(+)*)*] => {{
        $crate::Poly::new(
            vec![$((monomial![$($var)*], __poly_coeff!($($c)*)),)*]
        )
    }}
}

/// A single monomial, represented as
/// a vector of `u8`-indexed variables.
#[derive(Clone, Hash, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Monomial(Vec<u8>);

impl Monomial {

    /// Create a new empty monomial.
    pub fn one() -> Monomial {
        Monomial(vec![])
    }

    /// Create a new monomial from an iterator.
    pub fn new<I>(vars: I) -> Monomial
    where I: IntoIterator<Item = u8>, {
        Monomial(vars.into_iter().collect())
    }

    /// The degree of this monomial, i.e. the
    /// length of the backing vector.
    pub fn degree(&self) -> usize {
        self.0.len()
    }


    pub fn into_poly(self) -> Poly {
        Poly(vec![(self, 1.into())])
    }


    /// Reduce this monomial with respect to a relation,
    /// doing a single pass.
    pub fn reduce(&self, r: &Relation) -> Poly {
        let mut res = poly![(1)];

        if self.degree() < r.m.degree() {
            return res
        }

        let mut c = 0;
        loop {
            if &self.0[c..c + r.m.degree()] == &r.m.0[..] {
                res *= r.phi.clone();
                c += r.m.degree();
            } else {
                res *= Monomial(vec![self.0[c]]);
                c += 1;
            }

            if c + r.m.degree() > self.degree() {
                res *= Monomial(self.0[c..].into());
                break res;
            }
        }
    }

    pub fn divides(&self, dividend: &Monomial) -> bool {
        if self.degree() > dividend.degree() {
            return false
        }

        for sub in dividend.0.windows(self.degree()) {
            if self.0 == sub {
                return true
            }
        }

        return false
    }
}

impl Display for Monomial {
    fn fmt(&self, f: &mut Formatter) -> FResult {
        if self.0.is_empty() {
            return write!(f, "1");
        }

        let mut buf = Vec::new();
        let mut last = self.0[0];
        let mut count = 0;
        for &v in &self.0 {
            if v != last {
                buf.push(::std::str::from_utf8(&[last]).unwrap().to_string());
                if count > 1 {
                    buf.push(format!("^{}", count))
                }
                last = v;
                count = 1;
            } else {
                count += 1
            }
        }
        buf.push(::std::str::from_utf8(&[last]).unwrap().to_string());
        if count > 1 {
            buf.push(format!("^{}", count))
        }
        write!(f, "{}", buf.join(""))
    }
}

/// A polynomial with rational coefficients.
// this vec is in sorted order; monomials are inserted.
// this also makes multiplication not waste space.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Poly(Vec<(Monomial, Rational)>);

impl Poly {

    /// Create an empty polynomial
    pub fn zero() -> Poly {
        Poly(Vec::new())
    }

    fn insert_monomials(&mut self, mut monomials: Vec<(Monomial, Rational)>) {
        self.0.clear();
        monomials.sort_by(|a, b| a.0.cmp(&b.0));
        let mut total = 0.into();
        let mut last = Monomial::one();
        for (m, c) in monomials {
            if m != last {
                if total != 0.into() {
                    self.0.push((last, total));
                    total = 0.into();
                }
                last = m;
            }
            total += c;
        }
        if total != 0.into() {
            self.0.push((last, total));
        }
    }

    pub fn new<I>(monomials: I) -> Poly
    where I: IntoIterator<Item=(Monomial, Rational)> {
        let mut p = Poly::zero();
        p.insert_monomials(monomials.into_iter().collect());
        p
    }

    /// Check whether this is the zero polynomial
    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    /// Reduce a polynomial with respect to an ideal
    pub fn reduce(&self, ideal: &Ideal) -> Poly {
        fn reduce_rel(poly: &Poly, rel: &Relation) -> Poly {
            let mut res = Poly::zero();
            for &(ref k, ref v) in &poly.0 {
                let mut reduced = k.reduce(rel);
                reduced *= *v;
                res += reduced;
            }
            res
        }

        let mut last = self.clone();
        let mut res = last.clone();

        loop {
            for rel in &ideal.0 {
                res = reduce_rel(&res, rel);
            }
            if res == last {
                break res;
            } else {
                last = res.clone();
            }
        }
    }

    pub fn from_str(s: &str) -> Result<Poly, String> {

        use std::str;
        use std::str::FromStr;
        use nom::*;
        use std::ops::*;

        pub fn one_alpha<T>(input:T) -> IResult<T, T> where
            T: Slice<Range<usize>> + Slice<RangeFrom<usize>> + Slice<RangeTo<usize>>,
            T: InputIter+InputLength,
            <T as InputIter>::Item: AsChar {
            let input_length = input.input_len();
            if input_length == 0 {
                return IResult::Incomplete(Needed::Unknown);
            }

            let (_, next) = input.iter_indices().next().unwrap();

            if !next.is_alpha() {
                return IResult::Error(error_position!(ErrorKind::Alpha, input))
            } else {
                return IResult::Done(input.slice(1..), input.slice(0..1))
            }
        }

        named!(parse_lit<isize>,
            map_res!(
                map_res!(
                    ws!(digit),
                    str::from_utf8
                ), FromStr::from_str
            )
        );

        named!(parse_integral<Rational>,
            map!(
                parse_lit, Rational::from_integer
            )
        );

        named!(parse_fractional<Rational>, map!(
            do_parse!(
                numer: parse_lit >>
                ws!(tag!("/")) >>
                denom: parse_lit >>
                (numer, denom)
            ), |(n, d)| Rational::new(n, d)
        ));

        named!(parse_rational<Rational>, alt!(
           parse_fractional | parse_integral
        ));

        named!(parse_power<isize>, do_parse!(
            ws!(tag!("^")) >>
            power: ws!(parse_lit) >>
            (power)
        ));

        named!(parse_var<(u8, usize)>, do_parse!(
            var: ws!(one_alpha) >>
            power: opt!(parse_power) >>
            (var[0], power.unwrap_or(1) as usize)
        ));

        named!(parse_monomial<Monomial>, map!(
            fold_many0!(
                ws!(parse_var), Vec::new(), |mut acc: Vec<_>, (var, pow)| {
                    for _ in 0..pow {
                        acc.push(var)
                    }
                    acc
                }
            ), Monomial::new
        ));

        named!(parse_scaled_monomial<(Monomial, Rational)>, do_parse!(
           c: ws!(opt!(parse_rational)) >>
           m: ws!(parse_monomial) >>
           (m, c.unwrap_or(1.into()))
        ));

        named!(parse_signed_monomial<(Monomial, Rational, bool)>, do_parse!(
            sign: ws!(one_of!("+-")) >>
            rest: parse_scaled_monomial >>
            (rest.0, rest.1, sign == '+')
        ));

        named!(parse_polynomial_rest<Poly>, fold_many0!(
            parse_signed_monomial, Poly::zero(), |mut poly: Poly, m: (Monomial, Rational, bool)| {
                let (m, mut c, s) = m;
                if !s { c = -c }
                poly += (m, c);
                poly
            }
        ));

        named!(parse_polynomial<Poly>, do_parse!(
            leading: parse_scaled_monomial >>
            rest: parse_polynomial_rest >>
            ({let mut rest = rest; rest += leading; rest})
        ));

        parse_polynomial(s.as_bytes()).to_result().map_err(|s| format!("{}", s))
    }
}

// Addition is O(n log n). could be made better
// if we assign ids to all monomials of a given degree.

impl AddAssign<(Monomial, Rational)> for Poly {
    fn add_assign(&mut self, rhs: (Monomial, Rational)) {
        if rhs.1 == 0.into() {
            return
        }
        match self.0.binary_search_by_key(&&rhs.0, |&(ref k, _)| k) {
            Ok(idx) => {
                self.0[idx].1 += rhs.1;
                if self.0[idx].1 == 0.into() {
                    self.0.remove(idx);
                }
            },
            Err(idx) => self.0.insert(idx, rhs),
        }
    }
}

impl AddAssign<Monomial> for Poly {
    fn add_assign(&mut self, rhs: Monomial) {
        *self += (rhs, 1.into());
    }
}

impl AddAssign<Poly> for Poly {
    fn add_assign(&mut self, rhs: Poly) {
        for (m, c) in rhs.0 {
            *self += (m, c)
        }
    }
}

impl SubAssign<(Monomial, Rational)> for Poly {
    fn sub_assign(&mut self, rhs: (Monomial, Rational)) {
        *self += (rhs.0, -rhs.1);
    }
}

impl SubAssign<Monomial> for Poly {
    fn sub_assign(&mut self, rhs: Monomial) {
        *self -= (rhs, 1.into());
    }
}

impl SubAssign<Poly> for Poly {
    fn sub_assign(&mut self, rhs: Poly) {
        for (m, c) in rhs.0 {
            *self -= (m, c)
        }
    }
}

impl MulAssign<Rational> for Poly {
    fn mul_assign(&mut self, rhs: Rational) {
        if rhs == 0.into() {
            self.0.clear();
        } else {
            for &mut (_, ref mut v) in &mut self.0 {
                *v *= rhs;
            }
        }
    }
}

impl MulAssign<(Monomial, Rational)> for Poly {
    fn mul_assign(&mut self, rhs: (Monomial, Rational)) {
        if rhs.0.degree() == 0 {
            *self *= rhs.1;
            return;
        }

        // we can just append to each monomial, since
        // this is guaranteed to preserve ordering.
        for &mut (ref mut k, ref mut v) in &mut self.0 {
            k.0.extend_from_slice(&(rhs.0).0);
            *v *= rhs.1;
        }
    }
}

impl MulAssign<Monomial> for Poly {
    fn mul_assign(&mut self, rhs: Monomial) {
        *self *= (rhs, 1.into());
    }
}

impl MulAssign<Poly> for Poly {
    fn mul_assign(&mut self, mut rhs: Poly) {
        use std::mem;
        let mut monomials = Vec::new();
        for &mut (ref mut m1, ref c1) in &mut self.0 {
            for &mut (ref mut m2, ref c2) in &mut rhs.0 {
                if *c1 != 0.into() && *c2 != 0.into() {
                    let mut m = m1.clone();
                    m.0.extend_from_slice(&m2.0);
                    monomials.push((m, *c1 * c2))
                }
            }
        }
        self.insert_monomials(monomials);
    }
}

impl Display for Poly {
    fn fmt(&self, f: &mut Formatter) -> FResult {
        if self.is_zero() {
            return write!(f, "0");
        }
        let mut parts = Vec::new();
        let mut first = true;
        for &(ref m, mut c) in self.0.iter().rev() {
            if !first {
                parts.push(if c < 0.into() {
                    "-".into()
                } else {
                    "+".into()
                })
            }
            parts.push({

                if first {
                    first = false;
                } else {
                    c = c.abs();
                }

                if c == 1.into() {
                    format!("{}", m)
                } else if c == (-1).into() {
                    format!("-{}", m)
                } else {
                    if m.degree() == 0 {
                        format!("{}", c)
                    } else {
                        format!("{}{}", c, m)
                    }
                }
            })
        }

        write!(f, "{}", parts.join(" "))
    }
}

#[derive(Clone, Debug)]
struct Relation {
    m: Monomial,
    phi: Poly,
}

impl Relation {

    pub fn new(mut poly: Poly) -> Relation {
        let (k, v) = poly.0.pop().unwrap();
        poly *= -v.recip();
        Relation{ m: k, phi: poly }
    }

    pub fn into_poly(&self) -> Poly {
        let mut poly = self.phi.clone();
        poly *= Rational::new(-1, 1);
        poly += self.m.clone();
        poly
    }

}

impl Display for Relation {
    fn fmt(&self, f: &mut Formatter) -> FResult {
        write!(f, "{} = {}", self.m, self.phi)
    }
}

#[derive(Clone, Debug)]
struct Ideal(Vec<Relation>);

impl Ideal {

    pub fn new<I>(polys: I) -> Ideal
    where I: IntoIterator<Item = Poly> {
        let mut ideal = Ideal(Vec::new());
        for mut p in polys {
            if p.is_zero() {
                continue;
            }

            ideal.0.push(Relation::new(p));
        }

        ideal
    }


    /// Simplify an ideal by reducing relations
    /// with respect to other relations.
    pub fn simplify(&mut self) {
        loop {
            let mut changed = false;
            for i in 0..self.0.len() {
                for j in 0..self.0.len() {
                    if i == j {
                        continue;
                    }
                    let (a, b): (&Relation, &mut Relation) = unsafe {
                        let a = &*(&self.0[j] as *const _);
                        let b = &mut *(&mut self.0[j] as *mut _);
                        (a, b)
                    };

                    if b.phi.is_zero() || a.m.degree() >= b.m.degree() ||
                        !a.m.divides(&b.m) {
                        continue;
                    }

                    changed = true;

                    let red = b.m.reduce(&a);
                    b.phi -= red;

                    if b.phi.is_zero() {
                        continue;
                    }

                    let (k, v) = b.phi.0.pop().unwrap();
                    b.phi *= v.recip();
                    b.m = k;
                }
            }
            if !changed {
                break;
            }
        }
    }

    /// Expand all overlaps, checking for
    /// consistency.
    pub fn expand_overlaps(&mut self, iterations: usize) {
        // run through all pairs of relations
        let mut new_rels = Vec::new();
        use std::collections::HashSet;
        let mut checked_pairs = HashSet::new();
        for _ in 0..iterations {
            for i in 0..self.0.len() {
                for j in 0..self.0.len() {
                    let r1 = &self.0[i];
                    let r2 = &self.0[j];

                    if !checked_pairs.contains(&(r1.m.clone(), r2.m.clone())) {
                        checked_pairs.insert((r1.m.clone(), r2.m.clone()));
                    } else {
                        continue
                    }

                    // check if a tail of the first is a head of the second
                    for k in 0..r1.m.degree() {
                        let tail = &r1.m.0[k..];
                        if r2.m.0.starts_with(tail) {
                            // we've found an overlap!

                            // convert the overlapped monomials into polys:

                            let mut p1 = Monomial::new(r1.m.0[..k].iter().cloned()).into_poly();
                            p1 *= r2.phi.clone();

                            let mut p2 = r1.phi.clone();
                            p2 *= Monomial::new(r2.m.0[r1.m.degree() - k..].iter().cloned());


                            p1 = p1.reduce(&self);
                            p2 = p2.reduce(&self);

                            // the reduced forms are inconsistent, so
                            // we need to add another relation.
                            if p1 != p2 {

                                p1 -= p2;
                                new_rels.push(Relation::new(p1));
                            }
                        }
                    }
                }
            }
            if new_rels.is_empty() {
                // everything is consistent!
                break;
            } else {
                for x in new_rels.drain(..) {
                    self.0.push(x);
                }
            }
        }
    }
}

impl Display for Ideal {
    fn fmt(&self, f: &mut Formatter) -> FResult {
        write!(f, "({})", self.0.iter().map(|x| format!("{}", x)).collect::<Vec<_>>().join(", "))
    }
}

fn main() {

    let matches = clap_app! ( grobner =>
        (about: "compute Gröbner bases of quotients of tensor algebras over Q")
        (setting: AppSettings::SubcommandRequired)
        (@arg ITERATIONS: -i --iterations +takes_value default_value("20")
            {|s| s.parse::<usize>().map(|_| ()).map_err(|_| "expected nonnegative integer".into())}
            "iterations to run when computing a basis")
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

    if let Some(matches) = matches.subcommand_matches("basis") {
        let mut ideal = Ideal::new(matches.values_of("IDEAL")
                        .unwrap()
                        .map(|s| Poly::from_str(&s).unwrap()));
        ideal.simplify();
        ideal.expand_overlaps(iterations);
        for r in &ideal.0 {
            println!("{}", r.into_poly());
        }
    } else if let Some(matches) = matches.subcommand_matches("reduce") {
        let mut ideal = Ideal::new(matches.values_of("IDEAL")
            .unwrap()
            .map(|s| Poly::from_str(&s).unwrap()));
        ideal.simplify();
        ideal.expand_overlaps(iterations);
        let poly = Poly::from_str(&matches.value_of("POLY").unwrap()).unwrap();
        println!("{} = {} (mod a)", poly, poly.reduce(&ideal));
    }

}
