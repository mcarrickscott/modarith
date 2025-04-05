#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(non_upper_case_globals)]
#![allow(unused_imports)]
// For example: cargo run --release --bin ed448
// python curve_rust.py 64 ED448
// This completes ed448.rs and edwards.rs for this curve
// EdDSA Implementation for curve ED448
// see https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-5.pdf


// to include in library..
/*
use crate::edwards;
use crate::edwards::*;
use crate::edwards::ECP;

use crate::hash;
use crate::hash::*;
use crate::hash::SHA3;
*/

mod edwards;
use edwards::*;
use edwards::ECP;

mod hash;
use hash::*;

/*** Insert code automatically generated in group.rs here ***/
/* Note that much of this code is not needed and can be deleted */

@group@

/*** End of automatically generated code ***/

// number of limbs and bytes in representation
pub const BYTES:usize = NBYTES;
const LIMBS:usize= NLIMBS;
type GEL = [SPINT; LIMBS];



// reduce 114 byte array h to integer r modulo group order q, in constant time
// Consider h as 2^440.(2^440 + y) + z, where x,y and z < q (z is bottom 55 bytes, y is next 55 bytes, x is top 4 bytes)
// Important that x,y and z < q, 55 bytes = 440 bits, q is 446 bits
fn reduce(h:&[u8],r:&mut [SPINT]) {
    let mut buff:[u8;BYTES]=[0;BYTES];    
    let mut x:GEL=[0;LIMBS];
    let mut y:GEL=[0;LIMBS];
    let mut z:GEL=[0;LIMBS];
    let mut c:GEL=[0;LIMBS];

    mod2r(8*(BYTES-1),&mut c);

    for i in 0..55 {  // bottom 55 bytes
        buff[i]=h[i];
    }
    buff[55]=0;
    buff.reverse();
    modimp(&buff,&mut z);

    for i in 0..55 { // middle 55 bytes
        buff[i]=h[i+55];
    }
    buff[55]=0;
    buff.reverse();
    modimp(&buff,&mut y);


    for i in 0..4 {
        buff[i]=h[110+i];
    }
    for i in 4..56 {
        buff[i]=0;
    }
    buff.reverse();
    modimp(&buff,&mut x);

    modmul(&c,&mut x);
    modadd(&y,&mut x);
    modmul(&c,&mut x);
    modadd(&z,&mut x);
    modcpy(&x,r);
}

// general purpose SHAKE256 hash function
// Input ilen bytes, output olen bytes
fn H(ilen:usize,olen:usize,s:&[u8],digest: &mut [u8])
{
    let mut sha3=SHA3::new(SHAKE256);
    for i in 0..ilen { 
        sha3.process(s[i]);
    }
    sha3.shake(digest,olen); 
}

// Input private key - 57 random bytes
// Output public key - 57 bytes
pub fn KEY_PAIR(prv: &[u8],public: &mut [u8]) {
    let mut P=ECP::new();
    ecngen(&mut P);
    let mut s:[u8;BYTES]=[0;BYTES]; 

    H(BYTES+1,BYTES,&prv,&mut s);
// clamp s
    s[0]&=0xFC;
    s[55]|=0x80;

    s.reverse(); // little endian to big endian
 
    ecnmul(&s,&mut P); 

    let sign=ecnget(&mut P,None,Some(&mut s));  // get y coordinate and sign
    s.reverse(); // big endian to little endian
    for i in 0..56 {
        public[i]=s[i];
    }
             
    public[56]=(sign<<7) as u8;
}

const dom4:[u8;10]=[b'S',b'i',b'g',b'E',b'd',b'4',b'4',b'8',0,0];

// input private key, public key, message to be signed. Output signature
pub fn SIGN(prv:&[u8],public: Option<&[u8]>,m:&[u8],sig:&mut [u8]) {
    let mut h:[u8;2*BYTES+2]=[0;2*BYTES+2];  
    let mut ipub:[u8;BYTES+1]=[0;BYTES+1];
    let mut sh:[u8;BYTES]=[0;BYTES];
    let mut s:GEL=[0;LIMBS];
    let mut r:GEL=[0;LIMBS];
    let mut d:GEL=[0;LIMBS];

    let mut sha3=SHA3::new(SHAKE256);

    let mut R=ECP::new();
    ecngen(&mut R);

    if let Some(pb) = public {
        for i in 0..BYTES+1 {
            ipub[i]=pb[i];
        }
    } else {
        KEY_PAIR(prv,&mut ipub);
    }

    H(BYTES+1,2*BYTES+2,&prv,&mut h);

    for i in 0..BYTES {
        sh[i]=h[i];
    }

// derive and clamp s
    sh[0]&=0xFC;
    sh[55]|=0x80;
    sh.reverse();
    modimp(&sh,&mut s);
    
    for i in 0..10 {
        sha3.process(dom4[i]);
    }
    for i in BYTES+1..2*BYTES+2 {
        sha3.process(h[i]);
    }

    for i in 0..m.len() {
        sha3.process(m[i]);
    }
    sha3.shake(&mut h,2*BYTES+2); 

    reduce(&h,&mut r);
    modexp(&r,&mut sh);  // convert to big endian array
    ecnmul(&sh,&mut R);

    let sign=ecnget(&mut R,None,Some(&mut sh));  // get y coordinate and sign
    sh.reverse();              // big endian to little endian

    for i in 0..BYTES {
        sig[i]=sh[i];
    }
    sig[BYTES]=(sign<<7) as u8; // first part of signature

    sha3=SHA3::new(SHAKE256);

    for i in 0..10 {
        sha3.process(dom4[i]);
    }
    for i in 0..BYTES+1 {
        sha3.process(sig[i]);  // R
    }
    for i in 0..BYTES+1 {
        sha3.process(ipub[i]);  // Q
    }
    for i in 0..m.len() {
        sha3.process(m[i]);
    }
    sha3.shake(&mut h,2*BYTES+2);

    reduce(&h,&mut d);
    modmul(&s,&mut d);
    modadd(&r,&mut d);

    modexp(&d,&mut sh);
    sh.reverse();
    for i in 0..BYTES {
        sig[BYTES+i+1]=sh[i];
    }
    sig[2*BYTES+1]=0;           // second part of signature
}

// input public key, message and signature
// NOTE signatures that are of the wrong length should be rejected prior to calling this function
pub fn VERIFY(public: &[u8],m:&[u8],sig:&[u8]) -> bool {
    let mut buff:[u8;BYTES]=[0;BYTES]; 
    let mut sh:[u8;BYTES]=[0;BYTES];
    let mut G=ECP::new();
    let mut R=ECP::new();
    let mut Q=ECP::new();
   
    let mut u:GEL=[0;LIMBS];
    let mut h:[u8;2*BYTES+2]=[0;2*BYTES+2];

    ecngen(&mut G);

// reconstruct point R
    for i in 0..BYTES {
        buff[i]=sig[i];
    }
    buff.reverse();
    let mut sign=(sig[BYTES]>>7) as usize;
    ecnset(sign,None,Some(&buff),&mut R);
    if ecnisinf(&R) {
        return false;
    }

// reconstruct point Q 
    for i in 0..BYTES {
        buff[i]=public[i];
    }
    buff.reverse();
    sign=((public[BYTES]>>7)&1) as usize;
    ecnset(sign,None,Some(&buff),&mut Q);
    if ecnisinf(&Q) {
        return false;
    }


    for i in 0..BYTES {
        buff[i]=sig[i+BYTES+1];
    }
    buff.reverse();

    let mut sha3=SHA3::new(SHAKE256);

    for i in 0..10 {
        sha3.process(dom4[i]);
    }
    for i in 0..BYTES+1 {
        sha3.process(sig[i]);  // R
    }
    for i in 0..BYTES+1 {
        sha3.process(public[i]);  // Q
    }
    for i in 0..m.len() {
        sha3.process(m[i]);
    }
    sha3.shake(&mut h,2*BYTES+2);

    reduce(&h,&mut u); modneg(&mut u); modexp(&u,&mut sh);

    if !modimp(&buff,&mut u) {
        return false;  // out of range
    }
    ecncof(&mut G); ecncof(&mut R); ecncof(&mut Q);
    Q=ecnmul2(&buff,&G,&sh,&Q);

    if ecncmp(&R,&Q) {
        return true;
    }
    return false;
}



// Some utility functions for I/O and debugging

fn char2int(inp: u8) -> u8 {
    if inp>='0' as u8 && inp <='9' as u8 {
        return inp-'0' as u8;
    }
    if inp>='A' as u8 && inp <='F' as u8 {
        return inp-('A' as u8) +10;
    }
    if inp>='a' as u8 && inp <='f' as u8 {
        return inp-('a' as u8) +10;
    }
    return 0;
}

// string s better have even number of characters!
fn from_hex(ilen:usize,s: &str,x: &mut[u8]) {
    let mut pad:[u8;128]=[0;128];
    let c=s.as_bytes();
    let len=c.len();
    let mut lz=2*ilen-len;
    if 2*ilen<len {lz=0;}
    for i in 0..lz {
        pad[i]='0' as u8;
    }
    for i in lz..2*ilen {
        pad[i]=c[i-lz];
    }

    for i in 0..ilen {
        x[i]=char2int(pad[2*i])*16+char2int(pad[2*i+1]);
    }
}

fn printhex(len:usize,array: &[u8]) {
    for i in 0..len {
        print!("{:02X}", array[i])
    }
    println!("")
}

fn main() {
    const SK:&str="c4eab05d357007c632f3dbb48489924d552b08fe0c353a0d4a1f00acda2c463afbea67c5e8d2877c5e3bc397a659949ef8021e954e0a12274e";

    let mut prv:[u8;BYTES+1]=[0;BYTES+1];
    let mut public:[u8;BYTES+1]=[0;BYTES+1];
    let mut sig:[u8;2*BYTES+2]=[0;2*BYTES+2];
    let mut m:[u8;2]=[0;2];

    println!("Run RFC8032 test vector");

    from_hex(BYTES+1,&SK,&mut prv); 
    print!("private key= "); printhex(BYTES+1,&prv);
    KEY_PAIR(&prv,&mut public);
    print!("Public key= "); printhex(BYTES+1,&public);

    m[0]=0x03; // message to be signed
    SIGN(&prv,Some(&public),&m[0..1],&mut sig);
    print!("signature=  "); printhex(2*BYTES+2,&sig); 

    let res=VERIFY(&public,&m[0..1],&sig);
    if res {
        println!("Signature is valid");
    } else {
        println!("Signature is NOT valid");
    }
}
