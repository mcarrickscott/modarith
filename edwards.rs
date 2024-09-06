#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(non_upper_case_globals)]
// Edwards curve support 
// Use python scripts to generate code for ED25519 or ED448, or your own curve
//
// Mike Scott 4th September 2024
// TII
//
// code for 16/32/64-bit processor for ED25519 curve can be generated  by 
//
// python curve_rust.py 16/32/64 ED25519
//
// code for 16/32/64-bit processor for ED448 curve can be generated  by
//
// python curve_rust.py 16/32/64 ED448

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
    modneg(&mut P.x);
}

// add Q to P
// standard projective method from EFD - https://www.hyperelliptic.org/EFD/
pub fn ecnadd(Q: &ECP,P: &mut ECP) {
    let mut a:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut b:[SPINT;NLIMBS]=[0;NLIMBS];     
    let mut c:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut d:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut e:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut f:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut g:[SPINT;NLIMBS]=[0;NLIMBS]; 

    modcpy(&Q.z,&mut a); modmul(&P.z,&mut a);
    modcpy(&a,&mut b); modsqr(&mut b);
    modcpy(&Q.x,&mut c); modmul(&P.x,&mut c);
    modcpy(&Q.y,&mut d); modmul(&P.y,&mut d);
    modcpy(&c,&mut e); modmul(&d,&mut e);
#[cfg(feature="MLI")]
    if CONSTANT_B>0 {
        modmli(CONSTANT_B as usize,&mut e);
        modcpy(&b, &mut f); modsub(&e,&mut f);
        modcpy(&b, &mut g); modadd(&e,&mut g);
    }
#[cfg(feature="MLI")]
    if CONSTANT_B<0 {
        modmli(-CONSTANT_B as usize,&mut e);
        modcpy(&b, &mut f); modadd(&e,&mut f);
        modcpy(&b, &mut g); modsub(&e,&mut g);
    }
#[cfg(not(feature="MLI"))]
    if CONSTANT_B==0 {
        modmul(&constant_b,&mut e);
        modcpy(&b, &mut f); modsub(&e,&mut f);
        modcpy(&b, &mut g); modadd(&e,&mut g);
    }
    modcpy(&P.x,&mut b); modadd(&P.y,&mut b);
    modcpy(&Q.x,&mut e); modadd(&Q.y,&mut e);
    modcpy(&b,&mut P.x); modmul(&e,&mut P.x);
    modsub(&c,&mut P.x);
    modsub(&d,&mut P.x);
    modmul(&f,&mut P.x);
    modmul(&a,&mut P.x);
    if CONSTANT_A == -1 {
        modcpy(&d,&mut P.y); modadd(&c,&mut P.y);
    } else {
        modcpy(&d,&mut P.y); modsub(&c,&mut P.y);
    }
    modmul(&a,&mut P.y);
    modmul(&g,&mut P.y);
    modcpy(&f,&mut P.z); modmul(&g,&mut P.z);
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
    let mut c:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut d:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut e:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut f:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut h:[SPINT;NLIMBS]=[0;NLIMBS];  
    let mut j:[SPINT;NLIMBS]=[0;NLIMBS];    

    modcpy(&P.x,&mut b); modadd(&P.y,&mut b);
    modsqr(&mut b);
    modcpy(&P.x,&mut c); modsqr(&mut c);
    modcpy(&P.y,&mut d); modsqr(&mut d);
    modcpy(&P.z,&mut h); modsqr(&mut h);
    modcpy(&h,&mut j); modadd(&j,&mut h);
    modcpy(&c,&mut e);

    if CONSTANT_A == -1 {
        modneg(&mut e);
    } 
    modcpy(&e,&mut f); modadd(&d,&mut f);
    modcpy(&f,&mut j); modsub(&h,&mut j);
    modcpy(&b,&mut P.x); modsub(&c,&mut P.x);
    modsub(&d,&mut P.x);
    modmul(&j,&mut P.x);
    modcpy(&e,&mut P.y); modsub(&d,&mut P.y);
    modmul(&f,&mut P.y);
    modcpy(&f,&mut P.z); modmul(&j,&mut P.z);
}

// set to infinity
pub fn ecninf(P: &mut ECP) {
    modzer(&mut P.x);
    modone(&mut P.y);
    modone(&mut P.z);
}

// test for infinity
pub fn ecnisinf(P: &ECP) -> bool {
    return modis0(&P.x) && modcmp(&P.y,&P.z);
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

// return 1 if equal, else 0
pub fn ecncmp(Q: &ECP,P: &ECP) -> usize {
    let mut a:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut b:[SPINT;NLIMBS]=[0;NLIMBS];  
    modcpy(&P.x,&mut a); modmul(&Q.z,&mut a);
    modcpy(&Q.x,&mut b); modmul(&P.z,&mut b);
    if !modcmp(&a,&b) {
        return 0;
    }
    modcpy(&P.y,&mut a); modmul(&Q.z,&mut a);
    modcpy(&Q.y,&mut b); modmul(&P.z,&mut b);
    if !modcmp(&a,&b) {
        return 0;
    }
    return 1;
}

// extract (x,y) from point, if y is NULL compress and just return x and sign of y, if x is NULL compress and just return y and sign of x
pub fn ecnget(P: &mut ECP,x: Option<&mut [u8]>,y: Option<&mut [u8]>) -> usize {
    let mut sx:[SPINT;NLIMBS]=[0;NLIMBS];  
    let mut sy:[SPINT;NLIMBS]=[0;NLIMBS]; 
    let mut ynull:bool=false;
    let mut xnull:bool=false;
    ecnaffine(P);
    if let Some(rx) = x {
        modcpy(&P.x,&mut sx);
        modexp(&sx,rx);
    } else {
        xnull=true;
    }
    if let Some(ry) = y {
        modcpy(&P.y,&mut sy);
        modexp(&sy,ry);
    } else {
        ynull=true;
    }
    if ynull {
        return modsign(&P.y);
    }
    if xnull {
        return modsign(&P.x);
    }
    return 0;
}

// general purpose set point function
// sets P=O if point not on curve
// if x and y are not NULL tries to set (x,y)
// if y==NULL tries to set from x and sign s of y (decompression)
// if x==NULL tries to set from y and sign s of x

fn setxy(s: usize,x: Option<&[SPINT]>,y: Option<&[SPINT]>,P: &mut ECP) {
    let mut sx:[SPINT;NLIMBS]=[0;NLIMBS];     
    let mut sy:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut o:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut u:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut v:[SPINT;NLIMBS]=[0;NLIMBS];   
    let mut h:[SPINT;NLIMBS]=[0;NLIMBS];  
    modone(&mut o);
    if let Some(ry)=y {
        if let Some(rx)=x {
            modcpy(&rx,&mut sx); modsqr(&mut sx);
            modcpy(&ry,&mut sy); modsqr(&mut sy);
            if CONSTANT_A == -1 {
                modcpy(&sy,&mut u); modsub(&sx,&mut u);
            } else {
                modcpy(&sy,&mut u); modadd(&sx,&mut u);
            }
            modcpy(&sx,&mut v); modmul(&sy,&mut v);
#[cfg(feature="MLI")]
            if CONSTANT_B>0 {
                modmli(CONSTANT_B as usize,&mut v);
                modadd(&o,&mut v);
            }
#[cfg(feature="MLI")]
            if CONSTANT_B<0 {
                modmli(-CONSTANT_B as usize,&mut v);
                modsub(&o,&mut v); modneg(&mut v);
            }
#[cfg(not(feature="MLI"))]
            if CONSTANT_B==0 {
                modmul(&constant_b,&mut v);
                modadd(&o,&mut v);
            }
            if modcmp(&u,&v) {
                modcpy(&rx,&mut P.x);
                modcpy(&ry,&mut P.y);
                modone(&mut P.z);
            } else {
                ecninf(P);
            }
            return;
        }
    }  

    if let Some(rx)=x {
        modcpy(&rx,&mut sx); modsqr(&mut sx);
        if CONSTANT_A == -1 {
            modcpy(&o,&mut u); modadd(&sx, &mut u);
        } else {
            modcpy(&o,&mut u); modsub(&sx, &mut u); 
        }
        modcpy(&sx,&mut v);
    }
    if let Some(ry)=y {        
        modcpy(&ry,&mut sy); modsqr(&mut sy);
        modcpy(&o,&mut u); modsub(&sy,&mut u);
        if CONSTANT_A == -1 {
            modneg(&mut o);
        } else {
            modcpy(&sy,&mut v);
        }
    }
#[cfg(feature="MLI")]
    if CONSTANT_B>0 {
        modmli(CONSTANT_B as usize,&mut v);
        modsub(&o,&mut v); modneg(&mut v);
    }
#[cfg(feature="MLI")]
    if CONSTANT_B<0 {
        modmli(-CONSTANT_B as usize,&mut v);
        modadd(&o,&mut v);
    }
#[cfg(not(feature="MLI"))]
    if CONSTANT_B==0 {
        modmul(&constant_b,&mut v);
        modsub(&o,&mut v); modneg(&mut v);
    }
    modcpy(&u,&mut o); modsqr(&mut o);
    modmul(&o,&mut u);
    modmul(&v,&mut u);
    modpro(&u,&mut h);
    if !modqr(Some(&h),&u) {
        ecninf(P);
        return;
    }
    modsqrt(&u,Some(&h),&mut v);
    modinv(Some(&h),&mut u);
    modmul(&v, &mut u);
    modmul(&o, &mut u);
    let d=(modsign(&u)-s)&1;
    modcpy(&u, &mut v); modneg(&mut v);
    modcmv(d,&v, &mut u);

    if let Some(rx)=x {
        modcpy(&u, &mut P.y);
        modcpy(&rx,&mut P.x);
    }
    if let Some(ry)=y {
        modcpy(&u, &mut P.x);
        modcpy(&ry,&mut P.y);
    }
    modone(&mut P.z);
}

// multiply point by small curve cofactor (here assumed to be 4 or 8)
pub fn ecncof(P: &mut ECP) {
    for _ in 0..COF {
        ecndbl(P);
    }
}

// Is (x,y) of the right order? Must be checked by calling program!
// api visible version, x and y are big endian byte arrays
pub fn ecnset(s: usize,x: Option<&[u8]>,y: Option<&[u8]>,P: &mut ECP) {
    let mut sx:[SPINT;NLIMBS]=[0;NLIMBS];     
    let mut sy:[SPINT;NLIMBS]=[0;NLIMBS];  
    if let Some(ry)=y {
        if let Some(rx)=x {
            modimp(&rx,&mut sx);
            modimp(&ry,&mut sy);
            setxy(s,Some(&sx),Some(&sy),P);
            return;
        }
    }
    if let Some(rx)=x {
        modimp(rx,&mut sx);
        setxy(s,Some(&sx),None,P);
    }
    if let Some(ry)=y {
        modimp(ry,&mut sy);
        setxy(s,None,Some(&sy),P);
    }
}

// set generator
pub fn ecngen(P: &mut ECP) {
#[cfg(feature="SMX")]
    if CONSTANT_X!=0 {
        let mut sx:[SPINT;NLIMBS]=[0;NLIMBS];  
        modint(CONSTANT_X,&mut sx);
        setxy(0,Some(&sx),None,P);
    } 
#[cfg(not(feature="SMX"))]
    if CONSTANT_X==0 {
        setxy(0,Some(&constant_x),Some(&constant_y),P);
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

