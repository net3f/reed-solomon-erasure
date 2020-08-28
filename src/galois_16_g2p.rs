//! GF(2^16) interface to g2p crate implementation

use crate::galois_8;
use g2p::{g2p, GaloisField};
use std::ops::{Add, Div, Mul, Sub};

// let one: GF16 = 1.into();
// let a: GF16 = 5.into();
// let b: GF16 = 4.into();
// let c: GF16 = 7.into();
// assert_eq!(a + c, 2.into());
// assert_eq!(a - c, 2.into());
// assert_eq!(a * b, c);
// assert_eq!(a / c, one / b);
// assert_eq!(b / b, one);
g2p::g2p!(GF2p16, 16);
/// An element of `GF(2^16)`.
type Element = GF2p16;
//struct Element(pub GF2p16);

impl Default for GF2p16 {
    fn default() -> Self {
        GF2p16(0)
    }
}

/// The field GF(2^16).
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq)]
pub struct Field;

impl crate::Field for Field {
    const ORDER: usize = 65536;

    type Elem = u16;

    fn add(a: u16, b: u16) -> u16 {
        (GF2p16(a) + GF2p16(b)).0
    }

    fn mul(a: u16, b: u16) -> u16 {
        (GF2p16(a) * GF2p16(b)).0
    }

    fn div(a: u16, b: u16) -> u16 {
        (GF2p16(a) / GF2p16(b)).0
    }

    fn zero() -> u16 {
        0
    }

    fn one() -> u16 {
        1
    }

    fn exp(elem: u16, n: usize) -> u16 {
        GF2p16(elem).pow(n).0
    }

    fn nth_internal(n: usize) -> u16 {
        GF2p16::GENERATOR.pow(n).0
    }
}

/// Type alias of ReedSolomon over GF(2^8).
pub type ReedSolomon = crate::ReedSolomon<Field>;

/// Type alias of ShardByShard over GF(2^8).
pub type ShardByShard<'a> = crate::ShardByShard<'a, Field>;

//#[derive(Debug, Copy, Clone, PartialEq, Eq)]
impl GF2p16 {
    // Create the zero element.
    fn zero() -> Self {
        GF2p16::ZERO
    }

    // A constant element evaluating to `n`.
    fn constant(n: u16) -> GF2p16 {
        GF2p16(n)
    }

    // Whether this is the zero element.
    fn is_zero(&self) -> bool {
        self == &GF2p16::ZERO
    }

    fn exp(mut self, n: usize) -> GF2p16 {
        self.pow(n)
    }

    //     // // reduces from some polynomial with degree <= 2.
    //     // #[inline]
    //     // fn reduce_from(mut x: [u8; 3]) -> Self {
    //     //     if x[0] != 0 {
    //     //         // divide x by EXT_POLY and use remainder.
    //     //         // i = 0 here.
    //     //         // c*x^(i+j)  = a*x^i*b*x^j
    //     //         x[1] ^= galois_8::mul(1] x[0]);
    //     //         x[2] ^= galois_8::mul(2, x[0]);
    //     //     }

    //     //     Element([x[1], x[2]])
    //     // }

    //     // fn degree(&self) -> usize {
    //     //     if self.0[0] != 0 {
    //     //         1
    //     //     } else {
    //     //         0
    //     //     }
    //     // }
}

// impl From<GF2p16> for Element {
//     fn from(c: GF2p16) -> Self {
//         Element(c)
//     }
// }

// impl Default for Element {
//     fn default() -> Self {
//         Element::zero()
//     }
// }

// impl Add for Element {
//     type Output = Element;

//     fn add(self, other: Self) -> Element {
//         self + other
//     }
// }

// impl Sub for Element {
//     type Output = Element;

//     fn sub(self, other: Self) -> Element {
//         self - other
//     }
// }

// impl Mul for GF2p16 {
//     type Output = GF2p16;

//     fn mul(self, rhs: Self) -> GF2p16 {
//         self * rhs
//     }
// }

// impl Mul<u8> for GF2p16 {
//     type Output = GF2p16;

//     fn mul(self, rhs: u8) -> GF2p16 {
//         GF2p16(self.0 * rhs)
//     }
// }

// impl Div for GF2p16 {
//     type Output = GF2p16;

//     fn div(self, rhs: Self) -> GF2p16 {
//         self / rhs
//     }
// }

// helpers for division.
#[cfg(test)]
mod tests {
    use super::*;
    use quickcheck::Arbitrary;

    impl Arbitrary for Element {
        fn arbitrary<G: quickcheck::Gen>(gen: &mut G) -> Self {
            GF2p16(u16::arbitrary(gen))
        }
    }

    quickcheck! {
        fn qc_add_associativity(a: Element, b: Element, c: Element) -> bool {
            a + (b + c) == (a + b) + c
        }

        fn qc_mul_associativity(a: Element, b: Element, c: Element) -> bool {
            a * (b * c) == (a * b) * c
        }

        fn qc_additive_identity(a: Element) -> bool {
            let zero = Element::zero();
            a - (zero - a) == zero
        }

        fn qc_multiplicative_identity(a: Element) -> bool {
            a.is_zero() || {
                let one = Element::ONE;
                (one / a) * a == one
            }
        }

        fn qc_add_commutativity(a: Element, b: Element) -> bool {
            a + b == b + a
        }

        fn qc_mul_commutativity(a: Element, b: Element) -> bool {
            a * b == b * a
        }

        fn qc_add_distributivity(a: Element, b: Element, c: Element) -> bool {
            a * (b + c) == (a * b) + (a * c)
        }

        fn qc_inverse(a: Element) -> bool {
            a.is_zero() || {
                let inv = GF2p16(1)/a;
                a * inv == Element::constant(1)
            }
        }

        fn qc_exponent_1(a: Element, n: u8) -> bool {
            a.is_zero() || n == 0 || {
                let mut b = a.exp(n as usize);
                for _ in 1..n {
                    b = b / a;
                }

                a == b
            }
        }

        fn qc_exponent_2(a: Element, n: u8) -> bool {
            a.is_zero() || {
                let mut res = true;
                let mut b = Element::constant(1);

                for i in 0..n {
                    res = res && b == a.exp(i as usize);
                    b = b * a;
                }

                res
            }
        }

        fn qc_exp_zero_is_one(a: Element) -> bool {
            a.exp(0) == Element::constant(1)
        }
    }

    #[test]
    #[should_panic]
    fn test_div_b_is_0() {
        let _ = Element::ONE / Element::zero();
    }

    #[test]
    fn zero_to_zero_is_one() {
        assert_eq!(Element::zero().exp(0), Element::constant(1))
    }
}
