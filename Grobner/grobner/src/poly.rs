use *;

use num::traits::*;
use num::Rational;

use std::ops::*;
use std::fmt::{ Result as FResult, Display, Formatter };

#[allow(unused_macros)]
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

#[macro_export]
macro_rules! poly {
    [$(($($c:expr)*) $($var:ident)* $(+)*)*] => {{
        $crate::Poly::new(
            vec![$((monomial![$($var)*], __poly_coeff!($($c)*)),)*]
        )
    }}
}

/// A polynomial with rational coefficients.
// this vec is in sorted order; monomials are inserted.
// this also makes multiplication not waste space.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Poly<F>(pub Vec<(Monomial, F)>);

impl<F: Field> Poly<F> {
    /// Create an empty polynomial
    pub fn zero() -> Poly<F> {
        Poly(Vec::new())
    }

    fn insert_monomials(&mut self, mut monomials: Vec<(Monomial, F)>) {
        self.0.clear();
        monomials.sort_by(|a, b| a.0.cmp(&b.0));
        let mut total = F::zero();
        let mut last = Monomial::one();
        for (m, c) in monomials {
            if m != last {
                if total != F::zero() {
                    self.0.push((last, total));
                    total = F::zero();
                }
                last = m;
            }
            total = total + c;
        }
        if total != F::zero() {
            self.0.push((last, total));
        }
    }

    pub fn new<I>(monomials: I) -> Poly<F>
        where I: IntoIterator<Item = (Monomial, F)> {
        let mut p = Poly::zero();
        p.insert_monomials(monomials.into_iter().collect());
        p
    }

    /// Check whether this is the zero polynomial
    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    /// Reduce a polynomial with respect to an ideal
    pub fn reduce(&self, ideal: &Ideal<F>) -> Poly<F> {
        fn reduce_rel<F: Field>(poly: &Poly<F>, rel: &Relation<F>) -> Poly<F> {
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
}

impl Poly<Rational> {

    pub fn change_field<F: Field>(self) -> Poly<F> {
        Poly(self.0.into_iter()
            .map(|(m, c)| (m, F::from_rational(c)))
            .filter(|&(_, c)| !c.is_zero())
            .collect::<Vec<_>>())
    }

    pub fn from_str(s: &str) -> Result<Poly<Rational>, String> {

        use nom::*;
        use std::ops::*;
        use std::str::{ self, FromStr };

        pub fn one_alpha<T>(input:T) -> IResult<T, T> where
            T: Slice<Range<usize>> + Slice<RangeFrom<usize>> + Slice<RangeTo<usize>>,
            T: InputIter+InputLength,
            <T as InputIter>::Item: AsChar {
            let input_length = input.input_len();
            if input_length == 0 {
                return IResult::Incomplete(Needed::Unknown)
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
           complete!(parse_fractional) | complete!(parse_integral)
        ));

        named!(parse_power<isize>, do_parse!(
            ws!(tag!("^")) >>
            power: ws!(parse_lit) >>
            (power)
        ));

        named!(parse_var<(u8, usize)>, do_parse!(
            var: ws!(one_alpha) >>
            power: opt!(complete!(parse_power)) >>
            (var[0], power.unwrap_or(1) as usize)
        ));

        named!(parse_monomial<Monomial>, map!(
            fold_many0!(
                ws!(complete!(parse_var)), Vec::new(), |mut acc: Vec<_>, (var, pow)| {
                    for _ in 0..pow {
                        acc.push(var)
                    }
                    acc
                }
            ), Monomial::new
        ));

        named!(parse_scaled_monomial<(Monomial, Rational)>, alt!(
            complete!(do_parse!(
                c: ws!(opt!(parse_rational)) >>
                m: ws!(parse_monomial) >>
                (m, c.unwrap_or(1.into()))
            )) |
            complete!(do_parse!(
                c: ws!(parse_rational) >>
                (Monomial::one(), c)
            ))
        ));

        named!(parse_signed_monomial<(Monomial, Rational, bool)>, do_parse!(
            sign: ws!(one_of!("+-")) >>
            rest: parse_scaled_monomial >>
            (rest.0, rest.1, sign == '+')
        ));

        named!(parse_polynomial_rest<Poly<Rational>>, fold_many0!(
            parse_signed_monomial, Poly::zero(), |mut poly: Poly<Rational>, m: (Monomial, Rational, bool)| {
                let (m, mut c, s) = m;
                if !s { c = -c }
                poly += (m, c);
                poly
            }
        ));

        named!(parse_polynomial<Poly<Rational>>, do_parse!(
            leading: parse_scaled_monomial >>
            rest: parse_polynomial_rest >>
            ({let mut rest = rest; rest += leading; rest})
        ));

        let res = parse_polynomial(s.as_bytes());

        res.to_result().map_err(|s| format!("{}", s))
    }
}

// Addition is O(n log n). could be made better
// if we assign ids to all monomials of a given degree.

impl<F: Field> AddAssign<(Monomial, F)> for Poly<F> {
    fn add_assign(&mut self, rhs: (Monomial, F)) {
        if rhs.1 == F::zero() {
            return
        }
        match self.0.binary_search_by_key(&&rhs.0, |&(ref k, _)| k) {
            Ok(idx) => {
                self.0[idx].1 = self.0[idx].1 + rhs.1;
                if self.0[idx].1 == F::zero() {
                    self.0.remove(idx);
                }
            },
            Err(idx) => self.0.insert(idx, rhs),
        }
    }
}

impl<F: Field> AddAssign<Monomial> for Poly<F> {
    fn add_assign(&mut self, rhs: Monomial) {
        *self += (rhs, F::one());
    }
}

impl<F: Field> AddAssign<Poly<F>> for Poly<F> {
    fn add_assign(&mut self, rhs: Poly<F>) {
        for (m, c) in rhs.0 {
            *self += (m, c)
        }
    }
}

impl<F: Field> SubAssign<(Monomial, F)> for Poly<F> {
    fn sub_assign(&mut self, rhs: (Monomial, F)) {
        *self += (rhs.0, -rhs.1);
    }
}

impl<F: Field> SubAssign<Monomial> for Poly<F> {
    fn sub_assign(&mut self, rhs: Monomial) {
        *self -= (rhs, F::one());
    }
}

impl<F: Field> SubAssign<Poly<F>> for Poly<F> {
    fn sub_assign(&mut self, rhs: Poly<F>) {
        for (m, c) in rhs.0 {
            *self -= (m, c)
        }
    }
}

impl<F: Field> MulAssign<F> for Poly<F> {
    fn mul_assign(&mut self, rhs: F) {
        if rhs == F::zero() {
            self.0.clear();
        } else {
            for &mut (_, ref mut v) in &mut self.0 {
                *v = *v * rhs;
            }
        }
    }
}

impl<F: Field> DivAssign<F> for Poly<F> {
    fn div_assign(&mut self, rhs: F) {
        if rhs == F::zero() {
            panic!("Division by zero")
        } else {
            for &mut (_, ref mut v) in &mut self.0 {
                *v = *v / rhs;
            }
        }
    }
}

impl<F: Field> MulAssign<(Monomial, F)> for Poly<F> {
    fn mul_assign(&mut self, rhs: (Monomial, F)) {
        if rhs.0.degree() == 0 {
            *self *= rhs.1;
            return;
        }

        // we can just append to each monomial, since
        // this is guaranteed to preserve ordering.
        for &mut (ref mut k, ref mut v) in &mut self.0 {
            k.0.extend_from_slice(&(rhs.0).0);
            *v = *v * rhs.1;
        }
    }
}

impl<F: Field> MulAssign<Monomial> for Poly<F> {
    fn mul_assign(&mut self, rhs: Monomial) {
        *self *= (rhs, F::one());
    }
}

impl<F: Field> MulAssign<Poly<F>> for Poly<F> {
    fn mul_assign(&mut self, mut rhs: Poly<F>) {
        let mut monomials = Vec::new();
        for &mut (ref mut m1, ref c1) in &mut self.0 {
            for &mut (ref mut m2, ref c2) in &mut rhs.0 {
                if *c1 != F::zero() && *c2 != F::zero() {
                    let mut m = m1.clone();
                    m.0.extend_from_slice(&m2.0);
                    monomials.push((m, *c1 * *c2))
                }
            }
        }
        self.insert_monomials(monomials);
    }
}

impl<F: Field + Signed + Display> Display for Poly<F> {
    fn fmt(&self, f: &mut Formatter) -> FResult {
        if self.is_zero() {
            return write!(f, "0");
        }
        let mut parts = Vec::new();
        let mut first = true;
        for &(ref m, mut c) in self.0.iter().rev() {
            if !first {
                parts.push(if c.is_negative() {
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

                if c == F::one() {
                    format!("{}", m)
                } else if c == -F::one() {
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