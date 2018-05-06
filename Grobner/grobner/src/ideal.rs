use *;

use num::traits::*;
use std::fmt::{ Result as FResult, Display, Formatter };

#[derive(Clone, Debug)]
pub struct Relation<F> {
    pub m: Monomial,
    pub phi: Poly<F>,
}

impl<F: Field> Relation<F> {

    pub fn new(mut poly: Poly<F>) -> Relation<F> {
        let (k, v) = poly.0.pop().unwrap();
        poly /= -v;
        Relation { m: k, phi: poly }
    }

    pub fn into_poly(&self) -> Poly<F> {
        let mut poly = self.phi.clone();
        poly *= -F::one();
        poly += self.m.clone();
        poly
    }

}

impl<F: Field + Signed + Display> Display for Relation<F> {
    fn fmt(&self, f: &mut Formatter) -> FResult {
        write!(f, "{} = {}", self.m, self.phi)
    }
}

#[derive(Clone, Debug)]
pub struct Ideal<F>(pub Vec<Relation<F>>);

impl<F: Field> Ideal<F> {

    pub fn new<I>(polys: I) -> Ideal<F>
        where I: IntoIterator<Item = Poly<F>> {
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
                    let (a, b): (&Relation<F>, &mut Relation<F>) = unsafe {
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

                    //if i == j {
                    //    continue;
                    //}

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

impl<F: Field + Signed + Display> Display for Ideal<F> {
    fn fmt(&self, f: &mut Formatter) -> FResult {
        write!(f, "({})", self.0.iter().map(|x| format!("{}", x)).collect::<Vec<_>>().join(", "))
    }
}