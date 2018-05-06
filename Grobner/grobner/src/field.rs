
use num::traits::*;
use num::Rational;
use std::ops::*;
use std::fmt::{ Result as FResult, Display, Formatter };

pub trait Field: Num + Neg<Output = Self> + Copy {

    fn recip(self) -> Self {
        Self::one() / self
    }

    fn from_rational(q: Rational) -> Self;
}

impl Field for Rational {
    fn from_rational(q: Rational) -> Self {
        q
    }
}

// this part is horrible

static mut PRIME: usize = 0;

#[inline]
pub fn prime() -> usize {
    unsafe { PRIME }
}

pub unsafe fn set_prime(p: usize) {
    if PRIME == 0 {
        PRIME = p
    } else {
        panic!("double-set of PRIME")
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Finite(usize);

impl Finite {

    #[inline]
    pub fn new(x: usize) -> Finite {
        Finite(x % prime())
    }

    #[inline]
    pub fn into_raw(self) -> usize {
        self.0
    }
}

impl Add for Finite {
    type Output = Finite;
    #[inline]
    fn add(self, rhs: Finite) -> Self::Output {
        Finite::new(self.0 + rhs.0)
    }
}

impl Neg for Finite {
    type Output = Finite;
    #[inline]
    fn neg(self) -> Self::Output {
        Finite(prime() - self.0)
    }
}

impl Sub for Finite {
    type Output = Finite;
    #[inline]
    fn sub(self, rhs: Finite) -> Self::Output {
        self + (-rhs)
    }
}

impl Mul for Finite {
    type Output = Finite;
    #[inline]
    fn mul(self, rhs: Finite) -> Self::Output {
        Finite::new(self.0 * rhs.0)
    }
}

impl Div for Finite {
    type Output = Finite;
    fn div(self, rhs: Finite) -> Self::Output {
        if rhs.0 == 0 {
            panic!("attempted division by zero")
        }
        // extended Euclidean algorithm,
        // for solving at = m mod n
        let mut t = 0;
        let mut r = prime() as isize;
        let mut next_t = self.0 as isize;
        let mut next_r = rhs.0 as isize;

        while next_r != 0 {
            let quot = r / next_r;

            let temp_t = t - quot * next_t;
            t = next_t;
            next_t = temp_t;

            let temp_r = r - quot * next_r;
            r = next_r;
            next_r = temp_r;
        }

        if r > 1 {
            panic!("attempted division by nonzero singular element")
        }

        if t < 0 {
            t += prime() as isize;
        }

        Finite(t as usize)
    }
}

impl Rem for Finite {
    type Output = Finite;
    fn rem(self, _: Finite) -> Self::Output {
        unimplemented!()
    }
}

impl Zero for Finite {
    #[inline]
    fn zero() -> Self {
        Finite(0)
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.0 == 0
    }
}

impl One for Finite {
    #[inline]
    fn one() -> Self {
        Finite(1)
    }
}

impl Num for Finite {
    type FromStrRadixErr = <usize as Num>::FromStrRadixErr;
    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        <usize as Num>::from_str_radix(str, radix).map(Finite::new)
    }
}

impl Field for Finite {
    #[inline]
    fn from_rational(q: Rational) -> Self {
        Finite::new(*q.numer() as usize) / Finite::new(*q.denom() as usize)
    }
}

impl Signed for Finite {
    fn abs(&self) -> Self {
        *self
    }

    fn abs_sub(&self, other: &Self) -> Self {
        *self - *other
    }

    fn signum(&self) -> Self {
        if self.0 == 0 {
            Finite::zero()
        } else {
            Finite::one()
        }
    }

    fn is_positive(&self) -> bool {
        self.signum().0 != 0
    }

    fn is_negative(&self) -> bool {
        false
    }
}

impl Display for Finite {
    fn fmt(&self, f: &mut Formatter) -> FResult {
        write!(f, "{}", self.0)
    }
}