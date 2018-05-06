use *;

use std::cmp::{ PartialOrd, Ord, Ordering };
use std::fmt::{ Result as FResult, Display, Formatter };

#[allow(unused_macros)]
macro_rules! __monomial_inner {
    ($var:ident) => {stringify!($var).as_bytes()[0]}
}

#[macro_export]
macro_rules! monomial {
    [] => {$crate::Monomial::one()};
    [$($var:ident)*] => {{
        $crate::Monomial::new(&[$(__monomial_inner!($var),)*])
    }}
}

/// A single monomial, represented as
/// a vector of `u8`-indexed variables.
#[derive(Clone, Hash, Debug, PartialEq, Eq)]
pub struct Monomial(pub Vec<u8>);

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


    pub fn into_poly<F: Field>(self) -> Poly<F> {
        Poly(vec![(self, F::one())])
    }


    /// Reduce this monomial with respect to a relation,
    /// doing a single pass.
    pub fn reduce<F: Field>(&self, r: &Relation<F>) -> Poly<F> {
        let mut res = Monomial::one().into_poly();

        if self.degree() < r.m.degree() {
            return self.clone().into_poly()
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

impl PartialOrd for Monomial {

    fn partial_cmp(&self, other: &Monomial) -> Option<Ordering> {
        Some(<Self as Ord>::cmp(self, other))
    }
}

impl Ord for Monomial {

    fn cmp(&self, other: &Monomial) -> Ordering {
        match (self.degree(), other.degree()) {
            (a, b) if a < b => Ordering::Less,
            (a, b) if a > b => Ordering::Greater,
            _ => self.0.cmp(&other.0)
        }
    }
}