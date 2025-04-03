#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(non_upper_case_globals)]
#![allow(unused_imports)]
// For example: cargo run --release --bin EC256
// python curve_rust.py 64 NIST256
// This completes weierstrass.rs for this curve
// ECDSA Implementation for curve P-256
// see https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-5.pdf


// to include in library..
/*
use crate::weierstrass;
use crate::weierstrass::*;
use crate::weierstrass::ECP;

use crate::hash;
use crate::hash::*;

*/


mod weierstrass;
use weierstrass::*;
use weierstrass::ECP;

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

const COMPRESS:bool = false;

// reduce 40 byte array h to integer r modulo group order q, in constant time
// Consider h as 2^248.x + y, where x and y < q (x is top 9 bytes, y is bottom 31 bytes)
// Important that x and y < q
fn reduce(h:&[u8],r:&mut [SPINT]) {
    let mut buff:[u8;BYTES]=[0;BYTES];    
    let mut x:GEL=[0;LIMBS];
    let mut y:GEL=[0;LIMBS];
    let mut c:GEL=[0;LIMBS];

    mod2r(8*(BYTES-1),&mut c);

    for i in 0..BYTES-1 {
        buff[i]=h[i];
    }
    buff[BYTES-1]=0;
    buff.reverse();
    modimp(&buff,&mut y);

    for i in 0..9 {
        buff[i]=h[BYTES-1+i];
    }
    for i in 9..BYTES {
        buff[i]=0;
    }
    buff.reverse();
    modimp(&buff,&mut x);

    modmul(&c,&mut x); 
    modcpy(&x,r); modadd(&y,r);
}


// Input private key - 32 random bytes
// Output public key - 65 bytes (0x04<x>|<y>), or 33 if compressed (0x02<x>.. or 0x03<x>)
pub fn KEY_PAIR(compress: bool,prv: &[u8],public: &mut [u8]) {
    let mut P=ECP::new();
    let mut x:[u8;BYTES]=[0;BYTES];
    let mut y:[u8;BYTES]=[0;BYTES];
    let fb:u8;
    ecngen(&mut P);

    ecnmul(prv,&mut P); 

    if compress {
        fb=0x02+ecnget(&mut P,&mut x,None) as u8; // 0x02 or 0x03
        for i in 0..BYTES {
            public[1+i]=x[i];
        }
    } else {
        fb=0x04; // no compression
        ecnget(&mut P,&mut x,Some(&mut y));  // get x and y
        for i in 0..BYTES {
            public[1+i]=x[i];
            public[1+i+BYTES]=y[i];
        }
    }
    public[0]=fb;
}

pub fn PREHASH(sha: usize,m: &[u8]) -> [u8; BYTES] {
    let mut thm: [u8; BYTES] = [0; BYTES];
    if sha==32 {
        let mut sh256 = SHA256::new();
        for i in 0..m.len() {
            sh256.process(m[i]);
        }
        let h=sh256.hash();
        for i in 0..BYTES {
            thm[i]=h[i];
        } 
    }
    if sha==48 {
        let mut sh384 = SHA384::new();
        for i in 0..m.len() {
            sh384.process(m[i]);
        }
        let h=sh384.hash();
        for i in 0..BYTES {
            thm[i]=h[i];
        }
    }
    return thm;
}

pub fn SIGN(prv: &[u8],ran: &[u8],thm:&[u8],sig: &mut [u8]) {
    let mut rb:[u8;BYTES]=[0;BYTES];
    let mut sb:[u8;BYTES]=[0;BYTES];
    let mut R=ECP::new();
    let mut e:GEL=[0;LIMBS];
    let mut r:GEL=[0;LIMBS];
    let mut s:GEL=[0;LIMBS];
    let mut k:GEL=[0;LIMBS];

    modimp(thm,&mut e);

    ecngen(&mut R);
    modimp(prv,&mut s);

    reduce(ran,&mut k);
    modexp(&k,&mut rb);
    ecnmul(&rb,&mut R);
    modinv(None,&mut k);

    ecnget(&mut R,&mut rb,None);
    modimp(&rb,&mut r);

    modmul(&r,&mut s);
    modadd(&e,&mut s);
    modmul(&k,&mut s);
    modzer(&mut k);

    modexp(&r,&mut rb);
    modexp(&s,&mut sb);

    for i in 0..BYTES {
        sig[i]=rb[i];
        sig[BYTES+i]=sb[i];
    }
}

// input public key, hashed message and signature
// NOTE signatures that are of the wrong length should be rejected prior to calling this function
pub fn VERIFY(public: &[u8],thm:&[u8],sig:&[u8]) -> bool {
    let mut G=ECP::new();
    let mut Q=ECP::new();

    let mut rb:[u8;BYTES]=[0;BYTES];
    let mut sb:[u8;BYTES]=[0;BYTES];
    let mut u:[u8;BYTES]=[0;BYTES];
    let mut v:[u8;BYTES]=[0;BYTES];
    
    let mut e:GEL=[0;LIMBS];
    let mut r:GEL=[0;LIMBS];
    let mut s:GEL=[0;LIMBS];
    let mut rds:GEL=[0;LIMBS];
    
    modimp(thm,&mut e);

    ecngen(&mut G);

// import from signature
    for i in 0..BYTES {
        rb[i]=sig[i];
        sb[i]=sig[BYTES+i];
    }    

    if !modimp(&rb,&mut r) {
        return false; // if not in range
    }
    if !modimp(&sb,&mut s) {
        return false;
    }

    if modis0(&r) || modis0(&s) {
        return false;
    }
    modinv(None,&mut s); modcpy(&r,&mut rds);
    modmul(&s,&mut rds); modexp(&rds,&mut v);  // export to byte array
    modmul(&e,&mut s); modexp(&s,&mut u); 

    for i in 0..BYTES {
        rb[i]=public[1+i];
        sb[i]=public[1+i+BYTES];
    }

    if public[0]==0x04 {
        ecnset(0,&rb,Some(&sb),&mut Q);
    } else {
        ecnset((public[0]&1) as usize,&rb,None,&mut Q);
    }

    let mut R=ecnmul2(&u,&G,&v,&Q);
    if ecnisinf(&R) {
        return false;
    }
    ecnget(&mut R,&mut rb,None);

    modimp(&rb,&mut e);
    if modcmp(&r,&e) {
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

// test vector for FIPS 186-3 ECDSA Signature Generation
fn main() {
    const SK:&str="519b423d715f8b581f4fa8ee59f4771a5b44c8130b4e3eacca54a56dda72b464";
    const RAN:&str="94a1bbb14b906a61a280f245f9e93c7f3b4a6247824f5d33b9670787642a68deb9670787642a68de";
    const MSG:&str="44acf6b7e36c1342c2c5897204fe09504e1e2efb1a900377dbc4e7a6a133ec56";
    let mut k:[u8;BYTES+8]=[0;BYTES+8];
    let mut m:[u8;BYTES]=[0;BYTES];
    let mut prv:[u8;BYTES]=[0;BYTES];
    let mut public:[u8;2*BYTES+1]=[0;2*BYTES+1];
    let mut sig:[u8;2*BYTES]=[0;2*BYTES];
    let compress:bool=true;
    println!("Run test vector");
    from_hex(BYTES,&SK,&mut prv);
    print!("private key= "); printhex(BYTES,&prv);
    from_hex(BYTES+8,&RAN,&mut k);
    from_hex(BYTES,&MSG,&mut m);
    KEY_PAIR(compress,&prv,&mut public);  
    print!("Public key= ");
    if compress {
        printhex(BYTES+1,&public);
    } else {
        printhex(2*BYTES+1,&public);
    }
    SIGN(&prv,&k,&m[0..32],&mut sig);
    print!("signature= "); printhex(2*BYTES,&sig);

    let res=VERIFY(&public,&m[0..32],&sig);
    if res {
        println!("Signature is valid");
    } else {
        println!("Signature is NOT valid");
    }
}

