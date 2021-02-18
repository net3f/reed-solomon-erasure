//! GF(2^16) interface to its gf-complete library implementation

#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

use std::mem;
use std::ops::{Add, Div, Mul, Sub};
use std::convert::TryInto;

#[cfg(feature = "simd-accel")]
use libc::c_void;

use crate::Field as FieldTrait;

/// An element of `GF(2^16)`.
type Element = u16;
const EXTENSION_DEGREE: i32 = 16;

/// The field GF(2^16).
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq)]
pub struct Field;

static mut GF2_to_16 : Option<gf_t> = None;
//                                         divide: None,
//                                         inverse: None,
//                                         multiply_region: None,
//                                         extract_word: None,
//                                         scratch: None
// };

#[cfg(feature = "simd-accel")]
pub fn init_gf_c_field() {

    unsafe{
        let mut gf: gf_t = { mem::zeroed() };
        let null_ptr: *mut c_void = ::std::ptr::null::<c_void>() as *mut c_void;
        let null_gf_ptr: *mut gf_t = ::std::ptr::null::<gf_t>() as *mut gf_t;
        gf_init_hard(&mut gf, 16, gf_mult_type_t_GF_MULT_SPLIT_TABLE as i32, GF_REGION_ALTMAP as i32, gf_division_type_t_GF_DIVIDE_DEFAULT as i32,  0, 16, 4, null_gf_ptr, null_ptr );
        //gf_init_hard(&mut gf, 16, gf_mult_type_t_GF_MULT_CARRY_FREE as i32, GF_REGION_DEFAULT as i32, gf_division_type_t_GF_DIVIDE_DEFAULT as i32,  0, 0, 0, null_gf_ptr, null_ptr );
        //gf_init_easy(&mut gf, EXTENSION_DEGREE);

        GF2_to_16 = Some(gf);
    }
}

impl crate::Field for Field {
    const ORDER: usize = 65536;

    type Elem = u16;

    fn add(a: u16, b: u16) -> u16 {
        add(a, b)
    }

    fn mul(a: u16, b: u16) -> u16 {
       mul(a, b)
    }

    fn div(a: u16, b: u16) -> u16 {
        div(a, b)
    }

    fn zero() -> u16 {
        0
    }

    fn one() -> u16 {
        1
    }

    fn exp(elem: u16, n: usize) -> u16 {
        exp(elem, n)
    }

    fn nth_internal(n: usize) -> u16 {
        n.try_into().unwrap() //TODO: this should be mod the field poly
    }

    // #[cfg(feature = "simd-accel")]    
    // fn mul_slice(c: u16, input: &[u16], out: &mut [u16]) {
    //     mul_slice(c, input, out)
    // }

    // #[cfg(feature = "simd-accel")]    
    // fn mul_slice_add(c: u16, input: &[u16], out: &mut [u16]) {
    //     mul_slice_xor(c, input, out)
    // }

}

#[cfg(feature = "simd-accel")]
pub fn mul_slice(c: u16, input: &[u16], out: &mut [u16]) {
    unsafe {
        let input_ptr : *mut c_void = &input[0] as *const _ as *const c_void as *mut c_void;
        //let input_ptr : *const c_void = &input[0] as *const _ as *const c_void;
        let out_ptr : *mut c_void = &mut out[0] as *mut _ as *mut c_void;
 
        GF2_to_16.unwrap().multiply_region.w32.unwrap()(&mut GF2_to_16.unwrap(), input_ptr.into(), out_ptr.into(), c.into(), (input.len() * 2) as i32, 0)            
    }

    // gf.multiply_region.w32(&gf, r1, r2, a, 16, 0);
    
    // let low: *const u8 = &MUL_TABLE_LOW[c as usize][0];
    // let high: *const u8 = &MUL_TABLE_HIGH[c as usize][0];

    // assert_eq!(input.len(), out.len());

    // let input_ptr: *const u8 = &input[0];
    // let out_ptr: *mut u8 = &mut out[0];
    // let size: libc::size_t = input.len();

    // let bytes_done: usize =
    //     unsafe { reedsolomon_gal_mul(low, high, input_ptr, out_ptr, size) as usize };

    // mul_slice_pure_rust(c, &input[bytes_done..], &mut out[bytes_done..]);
}

#[cfg(feature = "simd-accel")]
pub fn mul_slice_xor(c: u16, input: &[u16], out: &mut [u16]) {
    unsafe {
        let input_ptr : *mut c_void = &input[0] as *const _ as *const c_void as *mut c_void;
        //let input_ptr : *const c_void = &input[0] as *const _ as *const c_void;
        let out_ptr : *mut c_void = &mut out[0] as *mut _ as *mut c_void;
 
        GF2_to_16.unwrap().multiply_region.w32.unwrap()(&mut GF2_to_16.unwrap(), input_ptr.into(), out_ptr.into(), c.into(), (input.len() * 2) as i32, 1)            
    }
}
            
/// Type alias of ReedSolomon over GF(2^16).
pub type ReedSolomon = crate::ReedSolomon<Field>;

/// Type alias of ShardByShard over GF(2^16).
pub type ShardByShard<'a> = crate::ShardByShard<'a, Field>;

/// Add two elements.
pub fn add(a: u16, b: u16) -> u16 {
    a ^ b
}

/// Subtract `b` from `a`.
#[cfg(test)]
pub fn sub(a: u16, b: u16) -> u16 {
    a ^ b
}

/// Multiply two elements.
pub fn mul(a: u16, b: u16) -> u16 {
    unsafe {
        GF2_to_16.unwrap().multiply.w32.unwrap()(&mut GF2_to_16.unwrap(), a.into(), b.into()).try_into().unwrap()
    }
}

/// Divide one element by another. `b`, the divisor, may not be 0.
pub fn div(a: u16, b: u16) -> u16 {
        unsafe {
            GF2_to_16.unwrap().divide.w32.unwrap()(&mut GF2_to_16.unwrap(), a.into(), b.into()).try_into().unwrap()
        }
}

/// Compute a^n.
pub fn exp(mut elem: u16, n: usize) -> u16 {
    if n == 0 {
        1
    } else if elem == 0 {
        0
    } else {
        let x = elem;
        let mut res: u16 = 0;
        for _ in 1..n {
            elem = mul(elem, x);
        }        
        elem
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use quickcheck::Arbitrary;
    quickcheck! {
        fn qc_add_associativity(a: Element, b: Element, c: Element) -> bool {
            init_gf_c_field();
            add(a , add(b , c)) == add(add(a, b), c)
        }

        fn qc_mul_associativity(a: Element, b: Element, c: Element) -> bool {
            init_gf_c_field();
            mul(a, mul(b, c)) == mul( mul(a, b), c)
        }

        fn qc_additive_identity(a: Element) -> bool {
            init_gf_c_field();
            let zero = 0;
            sub(a, sub(zero, a)) == zero
        }

        fn qc_multiplicative_identity(a: Element) -> bool {
            init_gf_c_field();
            a == 0 || {
                let one = 1;
                mul(div(one, a), a) == one
            }
        }

        fn qc_add_commutativity(a: Element, b: Element) -> bool {
            init_gf_c_field();
            add(a,b) == add(b, a)
        }

        fn qc_mul_commutativity(a: Element, b: Element) -> bool {
            init_gf_c_field();
            mul(a, b) == mul(b, a)
        }

        fn qc_add_distributivity(a: Element, b: Element, c: Element) -> bool {
            init_gf_c_field();
            mul(a ,add(b, c)) == add(mul(a,b), mul (a, c))               
        }

        fn qc_inverse(a: Element) -> bool {
            init_gf_c_field();
            a == 0 || {
                let inv : u16 = div(1,a);
                mul(a, inv) == 1
            }
        }

        fn qc_exponent_1(a: Element, n: u8) -> bool {
            init_gf_c_field();
            a == 0 || n == 0 || {
                let mut b = exp(a, n as usize);
                for _ in 1..n {
                    b = div(b, a);
                }

                a == b
            }
        }

        fn qc_exponent_2(a: Element, n: u8) -> bool {
            init_gf_c_field();
            a == 0 || {
                let mut res = true;
                let mut b = 1;

                for i in 0..n {
                    res = res && b == exp(a, i as usize);
                    b = mul(b, a);
                }

                res
            }
        }

        fn qc_exp_zero_is_one(a: Element) -> bool {
            init_gf_c_field();
            exp(a,0) == 1
        }
    }

    #[test]
    fn gf_complete_inti() {
        init_gf_c_field();
    }

    #[test]
    #[should_panic]
    fn test_div_b_is_0() {
        init_gf_c_field();
        panic!();
        let result : u16 =  div(1 as u16, 0 as u16) as u16;
    }
    
    #[test]
    fn zero_to_zero_is_one() {
        init_gf_c_field();
        assert_eq!(exp(0,0), 1)
    }
}
