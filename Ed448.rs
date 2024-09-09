#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(non_upper_case_globals)]
#![allow(unused_imports)]
// For example: cargo run --release --features ED448 --bin Ed448
// python curve_rust.py 64 ED448
// This completes edwards.rs for this curve
// EdDSA Implementation for curve ED448
// see https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-5.pdf

#[cfg(feature="WEIERSTRASS")]
mod weierstrass;
#[cfg(feature="WEIERSTRASS")]
use weierstrass::*;
#[cfg(feature="WEIERSTRASS")]
use weierstrass::ECP;

#[cfg(feature="EDWARDS")]
mod edwards;
#[cfg(feature="EDWARDS")]
use edwards::*;
#[cfg(feature="EDWARDS")]
use edwards::ECP;

mod hash;
use hash::*;
use hash::SHA3;

/*** Insert code automatically generated in group.rs here ***/
/* Note that much of this code is not needed and can be deleted */



/*** End of automatically generated code ***/

// number of limbs and bytes in representation
const BYTES:usize = NBYTES;
const LIMBS:usize= NLIMBS;

// Some utility functions for I/O and debugging

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
// reduce 114 byte array h to integer r modulo group order q, in constant time
// Consider h as 2^472.x + 2^440.y + z, where x,y and z < q (z is first 55 bytes, y is next 4 bytes, x is last 55 bytes)
// precalculate c1=nres(2^472 mod q) and c2=nres(2^440 mod q)
// using utility ed448_order.py 
fn reduce(h:&[u8],r:&mut [SPINT]) {
    let mut buff:[u8;BYTES]=[0;BYTES];    
    let mut x:[SPINT;LIMBS]=[0;LIMBS];
    let mut y:[SPINT;LIMBS]=[0;LIMBS];
    let mut z:[SPINT;LIMBS]=[0;LIMBS];

    for i in 0..55 {
        buff[i]=h[i];
    }
    buff[55]=0;
    buff.reverse();
    modimp(&buff,&mut z);

    for i in 0..4 {
        buff[i]=h[55+i];
    }
    for i in 4..56 {
        buff[i]=0;
    }
    buff.reverse();
    modimp(&buff,&mut y);

    for i in 0..55 {
        buff[i]=h[59+i];
    }
    buff[55]=0;
    buff.reverse();
    modimp(&buff,&mut x);

    if LIMBS==8 {
        let constant_c1: [SPINT; 8] = [0xe3033c23525654,0x7624da8d5b86ce,0x3a503352aa569,0x3337c35f209580,0xae17cf72c9860f,0x9cc14ba3c47c44,0xbcb7e4d070af1a,0x39f823b7292052];
        let constant_c2: [SPINT; 8] = [0xecfdd2acf1d7f0,0x1d5087ce8f0058,0xbe9c8459d676ba,0xff8fc444c3d266,0xa3c47c44ae17ce,0xd070af1a9cc14b,0xb7292052bcb7e4,0x383402a939f823];
        modmul(&constant_c1,&mut x);
        modmul(&constant_c2,&mut y);        
    } else {
        let constant_c1: [SPINT;16] = [0x3525654,0xe3033c2,0xd5b86ce,0x7624da8,0x52aa569,0x3a5033,0xf209580,0x3337c35,0x2c9860f,0xae17cf7,0x3c47c44,0x9cc14ba,0x70af1a,0xbcb7e4d,0x7292052,0x39f823b];
        let constant_c2: [SPINT;16] = [0xcf1d7f0,0xecfdd2a,0xe8f0058,0x1d5087c,0x9d676ba,0xbe9c845,0x4c3d266,0xff8fc44,0x4ae17ce,0xa3c47c4,0xa9cc14b,0xd070af1,0x2bcb7e4,0xb729205,0x939f823,0x383402a];
        modmul(&constant_c1,&mut x);
        modmul(&constant_c2,&mut y);
    }

    modcpy(&x,r); modadd(&y,r);
    modadd(&z,r);
}

// Input private key - 57 random bytes
// Output public key - 57 bytes
pub fn ED448_KEY_GEN(prv: &[u8],public: &mut [u8]) {
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
pub fn ED448_SIGN(prv:&[u8],public:&[u8],m:&[u8],sig:&mut [u8]) {
    let mut h:[u8;2*BYTES+2]=[0;2*BYTES+2];  
    let mut sh:[u8;BYTES]=[0;BYTES];
    let mut s:[SPINT;LIMBS]=[0;LIMBS];
    let mut r:[SPINT;LIMBS]=[0;LIMBS];
    let mut d:[SPINT;LIMBS]=[0;LIMBS];

    let mut sha3=SHA3::new(SHAKE256);

    let mut R=ECP::new();
    ecngen(&mut R);

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

    let mut i=0;
    while m[i]!=0 {
        sha3.process(m[i]);
        i+=1;
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
        sha3.process(public[i]);  // Q
    }
    i=0;
    while m[i]!=0 {
        sha3.process(m[i]);
        i+=1;
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

pub fn ED448_VERIFY(public: &[u8],m:&[u8],sig:&[u8]) -> bool {
    let mut buff:[u8;BYTES]=[0;BYTES]; 
    let mut sh:[u8;BYTES]=[0;BYTES];
    let mut G=ECP::new();
    let mut R=ECP::new();
    let mut Q=ECP::new();
   
    let mut u:[SPINT;LIMBS]=[0;LIMBS];
    let mut h:[u8;2*BYTES+2]=[0;2*BYTES+2];

    ecngen(&mut G);

// reconstruct point R
    for i in 0..BYTES {
        buff[i]=sig[i];
    }
    buff.reverse();
    let mut sign=(sig[BYTES]>>7) as usize;
    ecnset(sign,None,Some(&buff),&mut R);

// reconstruct point Q 
    for i in 0..BYTES {
        buff[i]=public[i];
    }
    buff.reverse();
    sign=((public[BYTES]>>7)&1) as usize;
    ecnset(sign,None,Some(&buff),&mut Q);

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
    let mut i=0;
    while m[i]!=0 {
        sha3.process(m[i]);
        i+=1;
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


fn main() {
    const sk:&str="c4eab05d357007c632f3dbb48489924d552b08fe0c353a0d4a1f00acda2c463afbea67c5e8d2877c5e3bc397a659949ef8021e954e0a12274e";

    let mut prv:[u8;BYTES+1]=[0;BYTES+1];
    let mut public:[u8;BYTES+1]=[0;BYTES+1];
    let mut sig:[u8;2*BYTES+2]=[0;2*BYTES+2];
    let mut m:[u8;2]=[0;2];

    println!("Run RFC8032 test vector");

    from_hex(BYTES+1,&sk,&mut prv); 
    println!("private key= "); printhex(BYTES+1,&prv);
    ED448_KEY_GEN(&prv,&mut public);
    print!("Public key= "); printhex(BYTES+1,&public);

    m[0]=0x03; // message to be signed
    m[1]=0;    // null terminated string
    ED448_SIGN(&prv,&mut public,&m,&mut sig);
    print!("signature=  "); printhex(2*BYTES+2,&sig); 

    let res=ED448_VERIFY(&public,&m,&sig);
    if res {
        println!("Signature is valid");
    } else {
        println!("Signature is NOT valid");
    }
}
