use std::{
    fmt::Display,
    ops::{Add, Div, Mul, Sub},
};

#[derive(Copy, Clone, Debug)]
pub struct GFElement<const P: usize> {
    _value: usize,
}

impl<const P: usize> GFElement<P> {
    pub const fn modulo() -> usize {
        P
    }

    pub fn value(&self) -> usize {
        self._value
    }

    pub fn inv(&self) -> Option<Self> {
        if self.value() == 0 {
            None
        } else {
            // p is expected to be prime
            Some(Self::from(self.value().pow((P - 2) as u32)))
        }
    }
}

impl<const P: usize> From<usize> for GFElement<P> {
    fn from(value: usize) -> Self {
        Self { _value: value % P }
    }
}

impl<const P: usize> Mul for GFElement<P> {
    type Output = Self;

    fn mul(self, rhs: GFElement<P>) -> Self::Output {
        Self::from(self.value() * rhs.value())
    }
}

impl<const P: usize> Add for GFElement<P> {
    type Output = Self;

    fn add(self, rhs: GFElement<P>) -> Self::Output {
        Self::from(self.value() + rhs.value())
    }
}

impl<const P: usize> Sub for GFElement<P> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::from(self.value() + Self::modulo() - rhs.value())
    }
}

impl<const P: usize> Div for GFElement<P> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self * rhs
            .inv()
            .expect("computing inverse for zero element in division")
    }
}

impl<const P: usize> Display for GFElement<P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.value())
    }
}

#[derive(Clone, Debug)]

pub struct EGF<const P: usize, const M: usize> {
    /// polynomial is big-endian - lowest powers come last
    polynomial: Vec<GFElement<P>>,
}

pub struct EGFElement<'f, const P: usize, const M: usize> {
    /// polynomial is big-endian - lowest powers come last
    _value: [GFElement<P>; M],
    field: &'f EGF<P, M>,
}

impl<const P: usize, const M: usize> EGF<P, M> {
    /// Construct Extended Galois Field over elements mod P with base polynomial.
    ///
    /// Polynomial is expected to be in big-endian form (lowest powers come last)
    pub fn new(polynomial: Vec<GFElement<P>>) -> Self {
        assert_eq!(polynomial.len(), M + 1);
        Self { polynomial }
    }

    /// Constrict EGF element from digits
    ///
    /// Digits are expected to be in big-endian form
    pub fn construct_element(&self, value: [GFElement<P>; M]) -> EGFElement<P, M> {
        EGFElement {
            _value: value,
            field: self,
        }
    }

    pub fn construct_from_digits(&self, value: [usize; M]) -> EGFElement<P, M> {
        self.construct_element(value.map(<GFElement<P>>::from))
    }

    pub const fn prime_base(&self) -> usize {
        P
    }

    pub const fn power(&self) -> usize {
        M
    }
}

impl<'f, const P: usize, const M: usize> EGFElement<'f, P, M> {
    pub fn into_digits(self) -> [GFElement<P>; M] {
        self._value
    }

    pub fn as_polynomial(&self) -> String {
        let mut items = self
            ._value
            .iter()
            .rev()
            .enumerate()
            .filter_map(|(power, digit)| {
                let digit = digit.value();

                Some(match (power, digit) {
                    (_, 0) => return None,
                    (0, d) => d.to_string(),
                    (1, 1) => "x".to_string(),
                    (1, d) => format!("{d}*x"),
                    (power, 1) => format!("x^{power}"),
                    (power, digit) => format!("{digit}*x^{power}"),
                })
            })
            .collect::<Vec<_>>();

        items.reverse();

        match items.len() {
            0 => "0".to_string(),
            1 => items.pop().unwrap(),
            _ => items
                .into_iter()
                .reduce(|mut acc, value| {
                    acc.push_str(" + ");
                    acc.push_str(&value);
                    acc
                })
                .unwrap(),
        }
    }
}

impl<'f, const P: usize, const M: usize> Add for EGFElement<'f, P, M> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut buf = self._value;
        for digit in 0..buf.len() {
            buf[digit] = buf[digit] + rhs._value[digit];
        }

        self.field.construct_element(buf)
    }
}

impl<'f, const P: usize, const M: usize> Sub for EGFElement<'f, P, M> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut buf = self._value;
        for digit in 0..buf.len() {
            buf[digit] = buf[digit] - rhs._value[digit];
        }

        self.field.construct_element(buf)
    }
}

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

    let egf: EGF<P, M> = EGF::new([1, 0, 3, 2].map(GF5::from).to_vec());

    let el = egf.construct_from_digits([1, 2, 0]);

    println!("{}", el.as_polynomial())
}
