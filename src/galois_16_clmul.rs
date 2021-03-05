//! GF(2^16) implementation using clmul cpu instruction

use std::arch::x86_64::*;
use std::convert::TryInto;
use std::ops::{Add, Div, Mul, Sub};

use crate::Field as FieldTrait;

include!(concat!(env!("OUT_DIR"), "/table_g2p16.rs"));

/// An element of `GF(2^16)`.
type Element = u16;
const EXTENSION_DEGREE: i32 = 16;
const prim_poly : u32 = 0x1002d;
//const reducing_poly : u32 = (prim_poly as u64) & 0x1ffff as u64;

/// The field GF(2^16).
#[derive(Debug, Default, Copy, Clone, PartialEq, Eq)]
pub struct Field;

impl crate::Field for Field {
    const ORDER: usize = 65536;

    type Elem = u16;

    fn add(a: u16, b: u16) -> u16 {
        add(a, b)
    }

    fn mul(a: u16, b: u16) -> u16 {
        unsafe {
            mul(a, b)
        }
    }

    fn div(a: u16, b: u16) -> u16 {
        //div(a, b)
        a
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

    // fn mul_slice(c: u16, input: &[u16], out: &mut [u16]) {
    //     mul_slice(c, input, out)
    // }

    // fn mul_slice_add(c: u16, input: &[u16], out: &mut [u16]) {
    //     mul_slice_xor(c, input, out)
    // }

}

// #[cfg(feature = "simd-accel")]
// pub fn mul_slice(c: u16, input: &[u16], out: &mut [u16]) {
//     unsafe {
//         let input_ptr : *mut c_void = &input[0] as *const _ as *const c_void as *mut c_void;
//         //let input_ptr : *const c_void = &input[0] as *const _ as *const c_void;
//         let out_ptr : *mut c_void = &mut out[0] as *mut _ as *mut c_void;
 
//         GF2_to_16.unwrap().multiply_region.w32.unwrap()(&mut GF2_to_16.unwrap(), input_ptr.into(), out_ptr.into(), c.into(), (input.len() * 2) as i32, 0)            
//     }

//     // gf.multiply_region.w32(&gf, r1, r2, a, 16, 0);
    
//     // let low: *const u8 = &MUL_TABLE_LOW[c as usize][0];
//     // let high: *const u8 = &MUL_TABLE_HIGH[c as usize][0];

//     // assert_eq!(input.len(), out.len());

//     // let input_ptr: *const u8 = &input[0];
//     // let out_ptr: *mut u8 = &mut out[0];
//     // let size: libc::size_t = input.len();

//     // let bytes_done: usize =
//     //     unsafe { reedsolomon_gal_mul(low, high, input_ptr, out_ptr, size) as usize };

//     // mul_slice_pure_rust(c, &input[bytes_done..], &mut out[bytes_done..]);
// }

// #[cfg(feature = "simd-accel")]
// pub fn mul_slice_xor(c: u16, input: &[u16], out: &mut [u16]) {
//     unsafe {
//         let input_ptr : *mut c_void = &input[0] as *const _ as *const c_void as *mut c_void;
//         //let input_ptr : *const c_void = &input[0] as *const _ as *const c_void;
//         let out_ptr : *mut c_void = &mut out[0] as *mut _ as *mut c_void;
 
//         GF2_to_16.unwrap().multiply_region.w32.unwrap()(&mut GF2_to_16.unwrap(), input_ptr.into(), out_ptr.into(), c.into(), (input.len() * 2) as i32, 1)            
//     }
// }
            
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
#[target_feature(enable = "sse2", enable = "sse4.1", enable = "pclmulqdq")]
unsafe fn mul(a: u16, b: u16) -> u16 {

        let a_m = _mm_insert_epi32 (_mm_setzero_si128(), a as i32, 0);
        let b_m = _mm_insert_epi32 (a_m, b as i32, 0);

        let prim_poly_m = _mm_set_epi32(0, 0, 0, prim_poly as i32);

        // /* Do the initial multiply */
  
        let mut result = _mm_clmulepi64_si128 (a_m, b_m, 0);

        let mut w = _mm_clmulepi64_si128 (prim_poly_m, _mm_srli_si128 (result, 2), 0);
        result = _mm_xor_si128 (result, w);
        w = _mm_clmulepi64_si128 (prim_poly_m, _mm_srli_si128 (result, 2), 0);
        result = _mm_xor_si128 (result, w);

        /* Extracts 32 bit value from result. */
        return _mm_extract_epi32(result, 0) as u16

}

/// Divide one element by another. `b`, the divisor, may not be 0.
pub fn div(a: u16, b: u16) -> u16 {
    a
//         unsafe {
//             GF2_to_16.unwrap().divide.w32.unwrap()(&mut GF2_to_16.unwrap(), a.into(), b.into()).try_into().unwrap()
//         }
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
            unsafe {
                elem = mul(elem, x);
            }
        }        
        elem
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use quickcheck::Arbitrary;
    use std::arch::x86_64::*;

    // quickcheck! {
    //     fn qc_add_associativity(a: Element, b: Element, c: Element) -> bool {
    //         add(a , add(b , c)) == add(add(a, b), c)
    //     }

    //     fn qc_mul_associativity(a: Element, b: Element, c: Element) -> bool {
    //         mul(a, mul(b, c)) == mul( mul(a, b), c)
    //     }

    //     fn qc_additive_identity(a: Element) -> bool {
    //         let zero = 0;
    //         sub(a, sub(zero, a)) == zero
    //     }

    //     fn qc_multiplicative_identity(a: Element) -> bool {
    //         a == 0 || {
    //             let one = 1;
    //             //mul(div(one, a), a) == one
    //             1==1
    //         }
    //     }

    //     fn qc_add_commutativity(a: Element, b: Element) -> bool {
    //         add(a,b) == add(b, a)
    //     }

    //     fn qc_mul_commutativity(a: Element, b: Element) -> bool {
    //         mul(a, b) == mul(b, a)
    //     }

    //     fn qc_add_distributivity(a: Element, b: Element, c: Element) -> bool {
    //         mul(a ,add(b, c)) == add(mul(a,b), mul (a, c))               
    //     }

    //     fn qc_inverse(a: Element) -> bool {
    //         a == 0 || {
    //             let inv : u16 = div(1,a);
    //             mul(a, inv) == 1
    //         }
    //     }

    //     fn qc_exponent_1(a: Element, n: u8) -> bool {
    //         a == 0 || n == 0 || {
    //             let mut b = exp(a, n as usize);
    //             for _ in 1..n {
    //                 b = div(b, a);
    //             }

    //             a == b
    //         }
    //     }

    //     fn qc_exponent_2(a: Element, n: u8) -> bool {
    //         a == 0 || {
    //             let mut res = true;
    //             let mut b = 1;

    //             for i in 0..n {
    //                 res = res && b == exp(a, i as usize);
    //                 b = mul(b, a);
    //             }

    //             res
    //         }
    //     }

    //     fn qc_exp_zero_is_one(a: Element) -> bool {
    //         exp(a,0) == 1
    //     }
    // }

   #[test]
   fn lots_of_mul() {
        use rand::Rng;

        let mut rng = rand::thread_rng();

        let mut a: u16 = rng.gen();
        let mut b: u16 = rng.gen();
        let mut c: u16 = 0;
        let mut d: u16 = 0;

        const number_of_mul: u32 = 1000000000;

       unsafe {
        for x in 0..number_of_mul {
            c = mul(a, b);
            a = b;
            b = c;
        }
           println!("{}", c);
       }

         // unsafe {
         //     for x in 0..number_of_mul {
             
         //         let a_m = _mm_insert_epi32 (_mm_setzero_si128(), a as i32, 0);
         //         let b_m = _mm_insert_epi32 (a_m, b as i32, 0);

         //         let prim_poly_m = _mm_set_epi32(0, 0, 0, prim_poly as i32);
                 
         //         // /* Do the initial multiply */
  
         //         let mut result = _mm_clmulepi64_si128 (a_m, b_m, 0);


         //         let mut w = _mm_clmulepi64_si128 (prim_poly_m, _mm_srli_si128 (result, 2), 0);
         //         result = _mm_xor_si128 (result, w);
         //         w = _mm_clmulepi64_si128 (prim_poly_m, _mm_srli_si128 (result, 2), 0);
         //         result = _mm_xor_si128 (result, w);

         //         c =  _mm_extract_epi32(result, 0) as u16;
         //         a = b;
         //         b = c;

         //     }
         // }

     }
    
    #[test]
    #[should_panic]
    fn test_div_b_is_0() {
        panic!();
        let result : u16 =  div(1 as u16, 0 as u16) as u16;
    }
    
    #[test]
    fn zero_to_zero_is_one() {
        assert_eq!(exp(0,0), 1)
    }
}
