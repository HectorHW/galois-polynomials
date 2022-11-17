#![feature(generic_const_exprs)]
use std::{
    collections::{HashMap, HashSet},
    fmt::Display,
    hash::Hash,
    ops::{Add, AddAssign, Div, Mul, Sub},
};

#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
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

impl<const P: usize> AddAssign for GFElement<P> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
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

    #[allow(clippy::suspicious_arithmetic_impl)]
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

#[derive(Clone, Debug, PartialEq, Eq)]

pub struct EGF<const P: usize, const M: usize>
where
    [(); M + 1]: Sized,
{
    /// polynomial is big-endian - lowest powers come last
    polynomial: [GFElement<P>; M + 1],

    _primitive: [GFElement<P>; M],
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EGFElement<'f, const P: usize, const M: usize>
where
    [(); M + 1]: Sized,
{
    /// polynomial is big-endian - lowest powers come last
    _value: [GFElement<P>; M],
    field: &'f EGF<P, M>,
}

impl<'f, const P: usize, const M: usize> Hash for EGFElement<'f, P, M>
where
    [(); M + 1]: Sized,
{
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self._value.hash(state);
    }
}

impl<const P: usize, const M: usize> EGF<P, M>
where
    [(); 2 * M + 1]: Sized,
    [(); 2 * M - 1]: Sized,
    [(); M + 1]: Sized,
{
    /// Construct Extended Galois Field over elements mod P with base polynomial.
    ///
    /// Polynomial is expected to be in big-endian form (lowest powers come last)
    pub fn new(polynomial: [GFElement<P>; M + 1]) -> Self {
        let mut value = Self {
            polynomial,
            _primitive: [Default::default(); M],
        };
        value.find_primitive();
        value
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

    pub fn primitive(&self) -> EGFElement<P, M> {
        self.construct_element(self._primitive)
    }

    fn find_primitive(&mut self) {
        let mut element: [GFElement<P>; M] = [Default::default(); M];
        element[M - 1] = 2.into();

        let element_power = P.pow(M as u32) - 1;

        let mut terms = HashSet::with_capacity(element_power);

        loop {
            terms.insert(element);

            let mut exponent = element;

            let mut power = 1;

            while power < element_power {
                exponent = (self.construct_element(exponent) * self.construct_element(element))
                    .into_digits();

                if !terms.insert(exponent) {
                    terms.clear();
                    break;
                }
                power += 1;
            }

            if power == element_power {
                self._primitive = element;
                return;
            }

            element = positional_inc(element);
            if is_zero(&element) {
                break;
            }
        }

        panic!("failed to find primitive element after exhaustive search")
    }
}

impl<'f, const P: usize, const M: usize> EGFElement<'f, P, M>
where
    [(); 2 * M + 1]: Sized,
    [(); 2 * M - 1]: Sized,
    [(); M + 1]: Sized,
{
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

    pub fn primitive_power(&self) -> usize {
        assert!(
            !is_zero(&self._value),
            "trying to find primitive power of zero"
        );

        let mut exp = 1;
        let mut element = self.field.primitive();

        loop {
            if &element == self {
                break;
            }
            exp += 1;
            element = element * self.field.primitive();
        }
        exp
    }
}

impl<'f, const P: usize, const M: usize> EGFElement<'f, P, M>
where
    [(); 2 * M + 1]: Sized,
    [(); 2 * M - 1]: Sized,
    [(); M + 1]: Sized,
{
    fn pow(self, power: usize) -> Self {
        if power == 0 {
            let mut buf = [Default::default(); M];
            buf[M - 1] = 1usize;
            return self.field.construct_from_digits(buf);
        }
        if power == 1 {
            return self;
        }

        if power % 2 == 0 {
            self.pow(power / 2) * self.pow(power / 2)
        } else {
            self * self.pow(power - 1)
        }
    }
}

impl<'f, const P: usize, const M: usize> Add for EGFElement<'f, P, M>
where
    [(); 2 * M + 1]: Sized,
    [(); 2 * M - 1]: Sized,
    [(); M + 1]: Sized,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut buf = self._value;
        for digit in 0..buf.len() {
            buf[digit] += rhs._value[digit];
        }

        self.field.construct_element(buf)
    }
}

impl<'f, const P: usize, const M: usize> Sub for EGFElement<'f, P, M>
where
    [(); 2 * M + 1]: Sized,
    [(); 2 * M - 1]: Sized,
    [(); M + 1]: Sized,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut buf = self._value;
        for digit in 0..buf.len() {
            buf[digit] = buf[digit] - rhs._value[digit];
        }

        self.field.construct_element(buf)
    }
}

fn multiply_by_digit<const P: usize, const M: usize>(
    mut value: [GFElement<P>; M],
    digit: GFElement<P>,
) -> [GFElement<P>; M] {
    for i in 0..M {
        value[i] = value[i] * digit;
    }
    value
}

fn sub_array<const P: usize, const M: usize>(
    mut a: [GFElement<P>; M],
    b: [GFElement<P>; M],
) -> [GFElement<P>; M] {
    for i in 0..M {
        a[i] = a[i] - b[i];
    }
    a
}

fn shift_by_power<const P: usize, const M: usize>(
    mut value: [GFElement<P>; M],
    power: usize,
) -> [GFElement<P>; M] {
    for i in 0..(M - power) {
        value[i] = value[i + power];
    }
    for i in (M - power)..M {
        value[i] = Default::default();
    }
    value
}

fn find_power<const P: usize, const M: usize>(value: [GFElement<P>; M]) -> usize {
    for i in 0..M {
        if value[i].value() != 0 {
            return M - 1 - i;
        }
    }

    0
}

fn positional_inc<const P: usize, const M: usize>(
    mut value: [GFElement<P>; M],
) -> [GFElement<P>; M] {
    let mut carry = true;

    for i in (0..M).rev() {
        if carry {
            let new_digit = value[i] + 1.into();
            carry = (value[i].value() + 1) >= P;
            value[i] = new_digit;
        }
    }

    value
}

fn is_zero<const P: usize, const M: usize>(value: &[GFElement<P>; M]) -> bool {
    for digit in value {
        if digit.value() != 0 {
            return false;
        }
    }
    true
}

fn euclidian_divrem<const P: usize, const M: usize>(
    mut a: [GFElement<P>; M],
    b: [GFElement<P>; M],
) -> ([GFElement<P>; M], [GFElement<P>; M]) {
    let mut quotient = [<GFElement<P>>::default(); M];

    let divisor_power = find_power(b);

    for power in (0..(M - divisor_power)).rev() {
        let divisor = shift_by_power(b, power);
        let idx = M - 1 - (power + divisor_power);
        while a[idx].value() != 0 {
            a = sub_array(a, divisor);
            quotient[idx] += <GFElement<P>>::from(1);
        }
    }
    (quotient, a)
}

impl<'f, const P: usize, const M: usize> Mul for EGFElement<'f, P, M>
where
    [(); 2 * M + 1]: Sized,
    [(); 2 * M - 1]: Sized,
    [(); M + 1]: Sized,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut powers: HashMap<usize, GFElement<P>> = Default::default();
        for digit1 in 0..self._value.len() {
            let power1 = M - digit1 - 1;

            for digit2 in 0..rhs._value.len() {
                let power2 = M - digit2 - 1;

                *powers.entry(power1 + power2).or_default() +=
                    self._value[digit1] * rhs._value[digit2];
            }
        }

        let mut multiplication_result = [Default::default(); { 2 * M - 1 }];

        for (power, power_value) in powers {
            multiplication_result[power] = power_value;
        }

        multiplication_result.reverse();

        let mut upcast_modulo = [Default::default(); { 2 * M - 1 }];

        upcast_modulo[(M - 1 - 1)..].copy_from_slice(&self.field.polynomial[..]);

        let (_d, r) = euclidian_divrem(multiplication_result, upcast_modulo);

        let mut downcast_q = [Default::default(); M];
        downcast_q.copy_from_slice(&r[(M - 1)..]);
        Self {
            _value: downcast_q,
            field: self.field,
        }
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

    let egf: EGF<P, M> = EGF::new([1, 0, 3, 2].map(GF5::from));

    let el = egf.construct_from_digits([1, 2, 0]);

    println!("{}", el.as_polynomial());

    println!("primitive: {:?}", egf.primitive().as_polynomial());

    println!(
        "{:?} = {:?} ^ {}",
        el.as_polynomial(),
        egf.primitive().as_polynomial(),
        el.primitive_power()
    )
}

#[cfg(test)]
mod tests {
    use crate::{positional_inc, GFElement, EGF};

    #[test]
    fn multiplication_works_without_overflow() {
        let field: EGF<2, 3> = EGF::new([1, 0, 1, 1].map(<GFElement<2>>::from));

        let q1 = field.construct_from_digits([0, 1, 0]);
        let q2 = field.construct_from_digits([0, 1, 1]);

        assert_eq!((q1 * q2).into_digits(), [1, 1, 0].map(<GFElement<2>>::from));
    }

    #[test]
    fn multiplication_by_zero_gives_zero() {
        let field: EGF<2, 3> = EGF::new([1, 0, 1, 1].map(<GFElement<2>>::from));

        let q1 = field.construct_from_digits([0, 1, 0]);
        let q2 = field.construct_from_digits([0, 0, 0]);

        assert_eq!((q1 * q2).into_digits(), [0, 0, 0].map(<GFElement<2>>::from));
    }

    #[test]
    fn multiplication_by_one_gives_same_element() {
        let field: EGF<2, 3> = EGF::new([1, 0, 1, 1].map(<GFElement<2>>::from));

        let q1 = field.construct_from_digits([0, 1, 0]);
        let q2 = field.construct_from_digits([0, 0, 1]);

        assert_eq!((q1 * q2).into_digits(), q1.into_digits());
    }

    #[test]
    fn multiplication_wraps_by_polynomial() {
        let field: EGF<2, 3> = EGF::new([1, 0, 1, 1].map(<GFElement<2>>::from));

        let q1 = field.construct_from_digits([1, 1, 0]);
        let q2 = field.construct_from_digits([1, 0, 1]);

        assert_eq!(q1 * q2, field.construct_from_digits([0, 1, 1]));
    }

    #[test]
    fn power_computation() {
        let field: EGF<2, 3> = EGF::new([1, 0, 1, 1].map(<GFElement<2>>::from));

        let q1 = field.construct_from_digits([1, 1, 0]);

        assert_eq!(q1.pow(3), field.construct_from_digits([1, 1, 1]));
    }

    #[test]
    fn add_one_increases() {
        let vector = [0, 1, 0].map(<GFElement<2>>::from);
        assert_eq!(positional_inc(vector), [0, 1, 1].map(<GFElement<2>>::from));
    }

    #[test]
    fn add_one_carrying() {
        let vector = [0, 1, 1].map(<GFElement<2>>::from);
        assert_eq!(positional_inc(vector), [1, 0, 0].map(<GFElement<2>>::from));
    }

    #[test]
    fn add_one_wraps() {
        let vector = [1, 1, 1].map(<GFElement<2>>::from);
        assert_eq!(positional_inc(vector), [0, 0, 0].map(<GFElement<2>>::from));
    }

    #[test]
    fn add_one_works_with_other_bases() {
        let vector = [1, 1, 4].map(<GFElement<5>>::from);
        assert_eq!(positional_inc(vector), [1, 2, 0].map(<GFElement<5>>::from));
    }
}
