#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(non_upper_case_globals)]
// Weierstrass curve support 
// Use python scripts to generate code for NIST256, or your own curve
//
// Mike Scott 4th September 2024
// TII
//
// code for 32/64-bit processor for NIST256 curve can be generated  by 
//
// python curve_rust.py 32/64 NIST256
//

// make sure decoration and generic are both set to False in monty_rust.py or pseudo_rust.py

/*** Insert automatically generated code for modulus field.rs here ***/

@field@

/*** Insert automatically generated definition for curve curve.rs here ***/

@curve@

/*** End of automatically generated code ***/


/*** Insert automatically generated point definition point.rs here ***/

@point@

/*** End of automatically generated code ***/

fn teq(b: usize,c: usize) -> usize {
    let mut x=b^c;
    x-=1;
    return (x>>(usize::BITS-1))&1;
}

pub fn ecncpy(Q: &ECP,P: &mut ECP) {
    modcpy(&Q.x,&mut P.x);
    modcpy(&Q.y,&mut P.y);
    modcpy(&Q.z,&mut P.z);
}

pub fn ecnneg(P: &mut ECP) {
    modneg(&mut P.y);
}

// add Q to P
// standard projective method from EFD - https://www.hyperelliptic.org/EFD/
pub fn ecnadd(Q: &ECP,P: &mut ECP) {
    let mut b:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut t0:[SPINT;NLIMBS]=[0;NLIMBS];     
    let mut t1:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut t2:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut t3:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut t4:[SPINT;NLIMBS]=[0;NLIMBS];   

    modcpy(&P.x,&mut t0); modmul(&Q.x,&mut t0);
    modcpy(&P.y,&mut t1); modmul(&Q.y,&mut t1);
    modcpy(&P.z,&mut t2); modmul(&Q.z,&mut t2);
    
    modcpy(&P.x,&mut t3); modadd(&P.y,&mut t3);
    modcpy(&Q.x,&mut t4); modadd(&Q.y,&mut t4);
    modmul(&t4,&mut t3);

    modcpy(&t0,&mut t4); modadd(&t1,&mut t4);
    modsub(&mut t4,&mut t3);
    modcpy(&P.y,&mut t4); modadd(&P.z,&mut t4);

    modcpy(&Q.y,&mut b); modadd(&Q.z,&mut b);
    modmul(&b,&mut t4);
    modcpy(&t1,&mut b); modadd(&t2,&mut b);
    
    modsub(&b, &mut t4);
    modadd(&P.z,&mut P.x);
    modcpy(&Q.z,&mut P.y); modadd(&Q.x,&mut P.y);

    modmul(&P.y,&mut P.x);
    modcpy(&t0,&mut P.y); modadd(&t2,&mut P.y);
    modsub(&P.x,&mut P.y); modneg(&mut P.y);// P.y is negative!

    if CONSTANT_A==0 {
        modcpy(&t0,&mut P.x); modadd(&t0,&mut P.x);
        modadd(&P.x,&mut t0);
#[cfg(feature="MLI")]
        if CONSTANT_B>0 {
            modmli(3*CONSTANT_B as usize,&mut t2);
            modmli(3*CONSTANT_B as usize,&mut P.y);
        }
#[cfg(feature="MLI")]
        if CONSTANT_B<0 {
            modmli((-3*CONSTANT_B) as usize,&mut t2); modneg(&mut t2);
            modmli((-3*CONSTANT_B) as usize,&mut P.y); modneg(&mut P.y);
        }
#[cfg(not(feature="MLI"))]
        if CONSTANT_B==0 {
            modcpy(&constant_b3,&mut b);
            modmul(&b, &mut t2);
            modmul(&b,&mut P.y); 
        }
        modcpy(&t1,&mut P.z); modadd(&t2,&mut P.z);
        modsub(&t2,&mut t1);

        modcpy(&P.y,&mut P.x); modmul(&t4,&mut P.x);
        modcpy(&t3,&mut t2); modmul(&t1,&mut t2);
        modsub(&t2,&mut P.x); modneg(&mut P.x);
        modmul(&t0,&mut P.y);
        modmul(&P.z,&mut t1);
        modadd(&t1,&mut P.y);
        modmul(&t3,&mut t0);
        modmul(&t4,&mut P.z);
        modadd(&t0,&mut P.z);
    } else {
#[cfg(feature="MLI")]
        if CONSTANT_B>0 {
            modcpy(&t2,&mut P.z); modmli(CONSTANT_B as usize,&mut P.z);
            modcpy(&P.y,&mut P.x); modsub(&P.z,&mut P.x);
            modmli(CONSTANT_B as usize,&mut P.y);
        }
#[cfg(feature="MLI")]
        if CONSTANT_B<0 {
            modcpy(&t2,&mut P.z); modmli((-CONSTANT_B) as usize,&mut P.z);
            modcpy(&P.y,&mut P.x); modadd(&P.z,&mut P.x);
            modmli((-CONSTANT_B) as usize,&mut P.y); modneg(&mut P.y);
        }
#[cfg(not(feature="MLI"))]
        if CONSTANT_B==0 {
            modcpy(&constant_b,&mut b);
            modcpy(&b,&mut P.z); modmul(&t2,&mut P.z);
            modcpy(&P.y,&mut P.x); modsub(&P.z,&mut P.x);
            modmul(&b,&mut P.y);
        }
        modcpy(&P.x,&mut P.z); modadd(&P.x,&mut P.z);

        modadd(&P.z,&mut P.x);
        modcpy(&t1,&mut P.z); modsub(&P.x,&mut P.z);
        modadd(&t1,&mut P.x);

        modcpy(&t2,&mut t1);  modadd(&t2,&mut t1);
        modadd(&t1,&mut t2);

        modsub(&t2,&mut P.y);
        modsub(&t0,&mut P.y);
        modcpy(&P.y,&mut t1); modadd(&P.y,&mut t1); 

        modadd(&t1,&mut P.y);
        modcpy(&t0,&mut t1); modadd(&t0,&mut t1);
        modadd(&t1,&mut t0);

        modsub(&t2,&mut t0);
        modcpy(&t4,&mut t1); modmul(&P.y,&mut t1);
        modcpy(&t0,&mut t2); modmul(&P.y,&mut t2);

        modcpy(&P.x,&mut P.y); modmul(&P.z,&mut P.y);
        modadd(&t2,&mut P.y);
        modmul(&t3,&mut P.x);

        modsub(&t1,&mut P.x);
        modmul(&t4,&mut P.z);
        modcpy(&t3,&mut t1); modmul(&t0,&mut t1);

        modadd(&t1,&mut P.z);
    }
}

pub fn ecnsub(Q: &ECP,P: &mut ECP) {
    let mut W = ECP::new();  
    ecncpy(Q,&mut W); ecnneg(&mut W);
    ecnadd(&W,P);
}

// double P
// standard projective method from EFD - https://www.hyperelliptic.org/EFD/
pub fn ecndbl(P: &mut ECP) {
    let mut b:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut t0:[SPINT;NLIMBS]=[0;NLIMBS];     
    let mut t1:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut t2:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut t3:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut t4:[SPINT;NLIMBS]=[0;NLIMBS]; 

    if CONSTANT_A==0 {
        modcpy(&P.y,&mut t0); modsqr(&mut t0);
        modcpy(&t0,&mut t3); modadd(&t0,&mut t3);
        modcpy(&t3,&mut t1); modadd(&t3,&mut t1);
        modcpy(&t1,&mut t3); modadd(&t1,&mut t3);
        modcpy(&P.x,&mut t4); modmul(&P.y,&mut t4);
        modcpy(&P.y,&mut t1); modmul(&P.z,&mut t1);
        modcpy(&P.z,&mut t2); modsqr(&mut t2);
#[cfg(feature="MLI")]
        if CONSTANT_B>0 {
            modmli(3*CONSTANT_B as usize,&mut t2);
        }
#[cfg(feature="MLI")]
        if CONSTANT_B<0 {
            modmli((-3*CONSTANT_B) as usize,&mut t2); modneg(&mut t2);
        }
#[cfg(not(feature="MLI"))]
        if CONSTANT_B==0 {
            modcpy(&constant_b3,&mut b);
            modmul(&b,&mut t2);
        }
        modcpy(&t2,&mut P.x); modmul(&t3,&mut P.x);
        modcpy(&t0,&mut P.y); modadd(&t2,&mut P.y);
        modcpy(&t3,&mut P.z); modmul(&t1,&mut P.z);

        modcpy(&t2,&mut t1); modadd(&t2,&mut t1);
        modadd(&t1,&mut t2);
        modsub(&t2,&mut t0);
        modmul(&t0,&mut P.y);
        modadd(&P.x,&mut P.y);
        modcpy(&t0,&mut P.x); modmul(&t4,&mut P.x);
        modcpy(&P.x,&mut t0); modadd(&t0,&mut P.x);
    } else {
        modcpy(&P.x,&mut t0); modsqr(&mut t0);
        modcpy(&P.y,&mut t1); modsqr(&mut t1);
        modcpy(&P.z,&mut t2); modsqr(&mut t2);

        modcpy(&P.x,&mut t3); modmul(&P.y,&mut t3);
        modcpy(&P.y,&mut t4); modmul(&P.z,&mut t4);
        modcpy(&t3,&mut b); modadd(&b,&mut t3);
        modmul(&P.x,&mut P.z);
        modcpy(&P.z,&mut b); modadd(&b,&mut P.z);
#[cfg(feature="MLI")]
        if CONSTANT_B>0 {
            modcpy(&t2,&mut P.y); modmli(CONSTANT_B as usize,&mut P.y);
            modsub(&P.z,&mut P.y);
            modmli(CONSTANT_B as usize,&mut P.z);
        }
#[cfg(feature="MLI")]
        if CONSTANT_B<0 {
            modcpy(&t2,&mut P.y); modmli((-CONSTANT_B) as usize,&mut P.y);  modneg(&mut P.y);
            modsub(&P.z,&mut P.y);
            modmli((-CONSTANT_B) as usize,&mut P.z); modneg(&mut P.z);
        }
#[cfg(not(feature="MLI"))]
        if CONSTANT_B==0 {
            modcpy(&constant_b,&mut b);
            modcpy(&t2,&mut P.y); modmul(&b,&mut P.y);
            modsub(&P.z,&mut P.y);
            modmul(&b,&mut P.z);
        }
        modcpy(&P.y,&mut P.x); modadd(&P.y,&mut P.x);
        modadd(&P.x,&mut P.y);
        modcpy(&t1,&mut P.x); modsub(&P.y,&mut P.x);

        modadd(&t1,&mut P.y);
        modmul(&P.x,&mut P.y);
        modmul(&t3,&mut P.x);

        modcpy(&t2,&mut t3); modadd(&t2,&mut t3);
        modadd(&t3,&mut t2);

        modsub(&t2,&mut P.z);
        modsub(&t0,&mut P.z);
        modcpy(&P.z,&mut t3); modadd(&P.z,&mut t3);

        modadd(&t3,&mut P.z);
        modcpy(&t0,&mut t3); modadd(&t0,&mut t3);
        modadd(&t3,&mut t0);

        modsub(&t2,&mut t0);
        modmul(&P.z,&mut t0);
        modadd(&t0,&mut P.y);

        modcpy(&t4,&mut b); modadd(&b,&mut t4);
        modmul(&t4,&mut P.z);

        modsub(&P.z,&mut P.x);
        modcpy(&t4,&mut P.z); modmul(&t1,&mut P.z);
        modcpy(&P.z,&mut b); modadd(&b,&mut P.z);
        modcpy(&P.z,&mut b); modadd(&b,&mut P.z);
    }
}

// set to infinity
pub fn ecninf(P: &mut ECP) {
    modzer(&mut P.x);
    modone(&mut P.y);
    modzer(&mut P.z);
}

// test for infinity
pub fn ecnisinf(P: &ECP) -> bool {
    return modis0(&P.x) && modis0(&P.z);
}

// set to affine
pub fn ecnaffine(P: &mut ECP) {
    let mut i:[SPINT;NLIMBS]=[0;NLIMBS]; 
    modcpy(&P.z,&mut i); modinv(None,&mut i);
    modone(&mut P.z);
    modmul(&i,&mut P.x);
    modmul(&i,&mut P.y);
}

// move Q to P if d=1
fn ecncmv(d: usize,Q: &ECP,P: &mut ECP) {
    modcmv(d,&Q.x,&mut P.x);
    modcmv(d,&Q.y,&mut P.y);
    modcmv(d,&Q.z,&mut P.z);
}

// return true if equal, else false
pub fn ecncmp(Q: &ECP,P: &ECP) -> bool {
    let mut a:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut b:[SPINT;NLIMBS]=[0;NLIMBS];  
    modcpy(&P.x,&mut a); modmul(&Q.z,&mut a);
    modcpy(&Q.x,&mut b); modmul(&P.z,&mut b);
    if !modcmp(&a,&b) {
        return false;
    }
    modcpy(&P.y,&mut a); modmul(&Q.z,&mut a);
    modcpy(&Q.y,&mut b); modmul(&P.z,&mut b);
    if !modcmp(&a,&b) {
        return false;
    }
    return true;
}

// extract (x,y) from point, if y is NULL compress and just return x and sign of y
pub fn ecnget(P: &mut ECP,x: &mut [u8],y: Option<&mut [u8]>) -> usize {
    let mut sx:[SPINT;NLIMBS]=[0;NLIMBS];  
    let mut sy:[SPINT;NLIMBS]=[0;NLIMBS]; 
    ecnaffine(P);
    modcpy(&P.x,&mut sx);
    modexp(&sx,x);
    if let Some(ry) = y {
        modcpy(&P.y,&mut sy);
        modexp(&sy,ry);
        return 0;
    } else {
        return modsign(&P.y); 
    }
}

// weierstrass set point function
// sets P=O if point not on curve
// if y!=NULL tries to set (x,y)
// if y==NULL calculates y (decompresses x) and selects sign from s=0/1

fn setxy(s: usize,x: &[SPINT],y: Option<&[SPINT]>,P: &mut ECP) {
    let mut t:[SPINT;NLIMBS]=[0;NLIMBS];     
    let mut v:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut h:[SPINT;NLIMBS]=[0;NLIMBS];    
    modcpy(x,&mut P.x);
    modcpy(x,&mut v); modsqr(&mut v);
    modmul(x,&mut v); // x^3
    if CONSTANT_A==-3 {
        modsub(x,&mut v);
        modsub(x,&mut v);
        modsub(x,&mut v); // x^3-3x
    }  
#[cfg(feature="MLI")]
    if CONSTANT_B>0 {
        modint(CONSTANT_B as usize,&mut t);
        modadd(&t,&mut v); // V=1+dx^2
    }
#[cfg(feature="MLI")]
    if CONSTANT_B<0 {
        modint((-CONSTANT_B) as usize,&mut t);
        modsub(&t,&mut v);
    }
#[cfg(not(feature="MLI"))]
    if CONSTANT_B==0 {
        modadd(&constant_b,&mut v);
    }
    if let Some(ry)=y {
        modcpy(&ry,&mut t); modsqr(&mut t);
        if modcmp(&t,&v) {
            modcpy(&ry,&mut P.y);
            modone(&mut P.z);
        } else {
            ecninf(P);
        }
        return;
    }
    modpro(&v,&mut h);
    if !modqr(Some(&h),&v)
    { // point not on curve
        ecninf(P);
        return;
    }
    modsqrt(&v,Some(&h),&mut P.y);
    let d=(modsign(&P.y)-s)&1;
    modcpy(&P.y,&mut t); modneg(&mut t);
    modcmv(d,&t,&mut P.y);
    modone(&mut P.z);
}

// multiply point by small curve cofactor (here assumed to be 4 or 8)
pub fn ecncof(_P: &mut ECP) {
}

// Is (x,y) of the right order? Must be checked by calling program!
// api visible version, x and y are big endian byte arrays
pub fn ecnset(s: usize,x: &[u8],y: Option<&[u8]>,P: &mut ECP) {
    let mut sx:[SPINT;NLIMBS]=[0;NLIMBS];     
    let mut sy:[SPINT;NLIMBS]=[0;NLIMBS];   
    modimp(&x,&mut sx);
    if let Some(ry)=y {
        modimp(&ry,&mut sy);
        setxy(s,&sx,Some(&sy),P);
        return;
    }        
    setxy(s,&sx,None,P);    
}

// set generator
pub fn ecngen(P: &mut ECP) {
#[cfg(feature="SMX")]
    if CONSTANT_X!=0 {
        let mut sx:[SPINT;NLIMBS]=[0;NLIMBS];  
        modint(CONSTANT_X,&mut sx);
        setxy(0,&sx,None,P);
    } 
#[cfg(not(feature="SMX"))]
    if CONSTANT_X==0 {
        setxy(0,&constant_x,Some(&constant_y),P);
    }
}

// select point from precomputed array in constant time
fn select(b:isize,W: &[ECP],P: &mut ECP) {
    let mut MP=ECP::new();
    let m=b>>(isize::BITS-1);
    let babs=((b^m)-m) as usize;

    ecncmv(teq(babs,0),&W[0],P);
    ecncmv(teq(babs,1),&W[1],P);
    ecncmv(teq(babs,2),&W[2],P);
    ecncmv(teq(babs,3),&W[3],P);
    ecncmv(teq(babs,4),&W[4],P);
    ecncmv(teq(babs,5),&W[5],P);
    ecncmv(teq(babs,6),&W[6],P);
    ecncmv(teq(babs,7),&W[7],P);
    ecncmv(teq(babs,8),&W[8],P);
    
    ecncpy(P,&mut MP);
    ecnneg(&mut MP);
    ecncmv((m&1) as usize,&MP,P);
}

// convert to double naf form
fn dnaf(e: &[u8],f: &[u8],w: &mut [i8]) {
    let mut ce:u8=0;
    let mut cf:u8=0;
    for i in 0..NBYTES {
        let mut m=e[NBYTES-i-1];
        let mut n=m;
        let mut t=3*(n as usize)+(ce as usize);
        ce=(t>>8) as u8;
        n=(t&0xff) as u8;
        let mut p=f[NBYTES-i-1];
        let mut q=p;
        t=3*(q as usize)+(cf as usize);
        cf=(t>>8) as u8;
        q=(t&0xff) as u8;
        for j in 0..8 {
            w[8*i+j]=((n&1)-(m&1)+3*((q&1)-(p&1))) as i8;
            n>>=1; m>>=1; p>>=1; q>>=1;
        }
    }
    for j in 0..8 {
        w[8*NBYTES+j]=((ce&1)+3*(cf&1)) as i8;
        ce>>=1; cf>>=1;
    }
}

pub fn ecnmul(e: &[u8],P: &mut ECP) {
    let mut Q=ECP::new();
    let mut W: [ECP; 9] = [
        ECP::new(),
        ECP::new(),
        ECP::new(),
        ECP::new(),
        ECP::new(),
        ECP::new(),
        ECP::new(),
        ECP::new(),
        ECP::new(),
    ];  
    let mut w: [i8; 2*NBYTES+1] = [0; 2*NBYTES+1];

    ecninf(&mut W[0]);                      // O
    ecncpy(P,&mut W[1]);                    // P
    ecncpy(P,&mut W[2]); ecndbl(&mut W[2]);         // 2P
    ecncpy(&W[2],&mut Q); ecncpy(&Q,&mut W[3]); ecnadd(P,&mut W[3]);  // 3P
    ecncpy(&Q,&mut W[4]); ecndbl(&mut W[4]); // 4P
    ecncpy(&W[4],&mut Q); ecncpy(&Q,&mut W[5]); ecnadd(P,&mut W[5]); // 5P
    ecncpy(&Q,&mut W[8]); ecndbl(&mut W[8]); // 8P
    ecncpy(&W[3],&mut Q); ecndbl(&mut Q); ecncpy(&Q,&mut W[6]); // 6P
    ecncpy(&Q,&mut W[7]); ecnadd(P,&mut W[7]); // 7P

    let mut j=0;
    for i in 0..NBYTES {
        let c=e[NBYTES-i-1] as i8;
        w[j]=c&0x0f;
        w[j+1]=(c>>4)&0xf;
        j+=2;
    }
    w[2*NBYTES]=0;
    for j in 0..2*NBYTES {
        let t=7-w[j];
        let m=(t>>4)&1;
        w[j]-=m<<4;
        w[j+1]+=m;
    }
    select(w[2*NBYTES] as isize,&W,P);
    for i in (0..2*NBYTES).rev() {
        select(w[i] as isize,&W,&mut Q);
        ecndbl(P);
        ecndbl(P);
        ecndbl(P);
        ecndbl(P);
        ecnadd(&Q,P);
    }
}

// double point multiplication R=eP+fQ
// not constant time
pub fn ecnmul2(e: &[u8],P: &ECP,f: &[u8],Q: &ECP) -> ECP {
    let mut R=ECP::new();
    let mut W: [ECP; 5] = [
        ECP::new(),
        ECP::new(),
        ECP::new(),
        ECP::new(),
        ECP::new(),
    ]; 
    let mut w: [i8; 8*NBYTES+8] = [0; 8*NBYTES+8]; 
    ecninf(&mut W[0]);      // O
    ecncpy(P,&mut W[1]);    // P
    ecncpy(Q,&mut W[3]);    // Q
    ecncpy(Q,&mut W[2]); ecnsub(P,&mut W[2]);    // Q-P
    ecncpy(Q,&mut W[4]); ecnadd(P,&mut W[4]);    // Q+P

    dnaf(&e,&f,&mut w);
    let mut i=8*NBYTES+7;
    while w[i]==0 {  // ignore leading zeros
        i-=1;
    }
    ecninf(&mut R);
    while i>=1 {
        ecndbl(&mut R);
        let j=w[i];
        if j>0 {
            ecnadd(&W[j as usize],&mut R);
        }
        if j<0 {
            ecnsub(&W[(-j) as usize],&mut R);
        }
        i-=1;
    }
    return R;
}


