#![allow(dead_code)]
//
// Rust Fp2 functions
// M.Scott 2026
//
// To compile:-
// rustc -O fp2.rs

// run monty_rust.py or pseudo_rust.py script for chosen modulus
// This code assumes modulus is prime = 3 mod 4

/*** Insert automatically generated code for modulus field.rs here ***/



/*** End of automatically generated code ***/

#[derive(Copy, Clone)]
pub struct Fp2 {
    a: [SPINT; NLIMBS],
    b: [SPINT; NLIMBS],
}

impl Fp2 {
// constructors
    pub const fn new() -> Fp2 {
        Fp2 {
            a: [0;NLIMBS],
            b: [0;NLIMBS],
        }
    }

    pub fn zero(&mut self) {
        modzer(&mut self.a);
	modzer(&mut self.b);
    }

    pub fn new_ints(x: usize,y: usize) -> Fp2 {
	let mut f = Fp2::new();
	modint(x,&mut f.a);
	modint(y,&mut f.b);
	f
    }

    pub fn new_fps(x: &[SPINT] ,y: &[SPINT]) -> Fp2 {
	let mut f = Fp2::new();
	modcpy(x,&mut f.a);
	modcpy(y,&mut f.b);
	f
    }

    pub fn new_bigs(x: &[SPINT] ,y: &[SPINT]) -> Fp2 {
	let mut f = Fp2::new();
	modcpy(x,&mut f.a); nres(&mut f.a);
	modcpy(y,&mut f.b); nres(&mut f.b);
	f
    }

    pub fn new_fp(x: &[SPINT]) -> Fp2 {
	let mut f = Fp2::new();
	modcpy(x,&mut f.a);
	modzer(&mut f.b);
	f
    }

    pub fn new_big(x: &[SPINT]) -> Fp2 {
	let mut f = Fp2::new();
	modcpy(x,&mut f.a); nres(&mut f.a);
	modzer(&mut f.b);
	f
    }

    pub fn new_copy(x: &Fp2) -> Fp2 {
        let mut f = Fp2::new();
        modcpy(&x.a,&mut f.a);
        modcpy(&x.b,&mut f.b);
        f
    }

// test self=0 ? 
    pub fn is0(&self) -> bool {
	modis0(&self.a) && modis0(&self.b)
    }

    pub fn is1(&self) -> bool {
	modis1(&self.a) && modis0(&self.b)
    }

    pub fn export(&self,x: &mut [u8],y: &mut[u8]) {
	modexp(&self.a,x);
	modexp(&self.b,y);
    }

    pub fn import(&mut self,x: &[u8],y: &[u8]) {
	modimp(x,&mut self.a);
	modimp(y,&mut self.b);
    }

    pub fn cmp(&self,x: &Fp2) -> bool {
	modcmp(&self.a,&x.a) && modcmp(&self.b,&x.b)
    }

    pub fn cmv(&mut self,d: usize,x:&Fp2) {
	modcmv(d,&x.a,&mut self.a);
	modcmv(d,&x.b,&mut self.b);
    }

    pub fn csw(&mut self,d: usize,x:&mut Fp2) {
	modcsw(d,&mut x.a,&mut self.a);
	modcsw(d,&mut x.b,&mut self.b);
    }

    pub fn split(&self,r: &mut [SPINT],i: &mut [SPINT]) {
	modcpy(&self.a,r); modcpy(&self.b,i);
    }

    pub fn copy(&mut self,x:&Fp2) {
	modcpy(&x.a,&mut self.a);
	modcpy(&x.b,&mut self.b);
    }

    pub fn sign(&self) -> usize {
	let p1=modsign(&self.a);
	let p2=modsign(&self.b);
	let u=modis0(&self.a);
	return p1^((p1^p2)&(u as usize));
    }

    pub fn neg(&mut self) {
	modneg(&mut self.a);
	modneg(&mut self.b);
    }

    pub fn conj(&mut self) {
	modneg(&mut self.b);
    }

    pub fn half(&mut self) {
	modhaf(&mut self.a);
	modhaf(&mut self.b);
    }

    pub fn times_i(&mut self) {
        let mut t: [SPINT; NLIMBS] = [0; NLIMBS];
	modcpy(&self.b,&mut t);
	modneg(&mut t);
	modcpy(&self.a,&mut self.b);
	modcpy(&t,&mut self.a);
    }

    pub fn size(&self) -> [SPINT;NLIMBS] {
        let mut w1: [SPINT; NLIMBS] = [0; NLIMBS];
        let mut w2: [SPINT; NLIMBS] = [0; NLIMBS];	

	modcpy(&self.a,&mut w1); modsqr(&mut w1);
	modcpy(&self.b,&mut w2); modsqr(&mut w2);
	modadd(&w1,&mut w2);
	return w2;
    }

    pub fn add(&mut self,x: &Fp2) {
        modadd(&x.a,&mut self.a);
        modadd(&x.b,&mut self.b);
    } 
    
    pub fn dbl(&mut self) {
	let x=Fp2::new_copy(&self);
	self.add(&x);
    }    
    
    pub fn sub(&mut self,x: &Fp2) {
        modsub(&x.a,&mut self.a);
        modsub(&x.b,&mut self.b);
    }    
 
    pub fn rsub(&mut self,x: &Fp2) {
        self.neg();
        self.add(x);
    }  

    pub fn imul(&mut self,x: usize) {
	modmli(x,&mut self.a);
	modmli(x,&mut self.b);
    }

    pub fn pmul(&mut self,x: &[SPINT]) {
	modmul(x,&mut self.a);
	modmul(x,&mut self.b);
    }

    pub fn sqr(&mut self) {
	let mut c: [SPINT; NLIMBS] = [0; NLIMBS];	
	let mut d: [SPINT; NLIMBS] = [0; NLIMBS];
	modcpy(&self.a,&mut c); modaddf(&self.a,&mut c);
	modcpy(&self.a,&mut d); modaddf(&self.b,&mut d);
	modsubf(&self.b,&mut self.a);
	modmul(&d,&mut self.a);
	modmul(&c,&mut self.b);
    }

// NOTE: faster "lazy" method may apply if prime modulus is not "shaped"
    pub fn mul(&mut self,y: &Fp2) {
	let mut c: [SPINT; NLIMBS] = [0; NLIMBS];	
	let mut d: [SPINT; NLIMBS] = [0; NLIMBS];
	let mut f: [SPINT; NLIMBS] = [0; NLIMBS];
	modcpy(&self.a,&mut c); modmul(&y.a,&mut c);
	modcpy(&self.b,&mut d); modmul(&y.b,&mut d);
	modcpy(&y.a,&mut f); modaddf(&y.b,&mut f);
	modaddf(&self.a,&mut self.b);
	modcpy(&c,&mut self.a); modsub(&d,&mut self.a);
	modadd(&d,&mut c);
	modmul(&f,&mut self.b);
	modsub(&c,&mut self.b);
    }

    pub fn inv(&mut self,h:Option<&[SPINT]>) {
	let mut w1: [SPINT; NLIMBS] = [0; NLIMBS];	
	let mut w2: [SPINT; NLIMBS] = [0; NLIMBS];
	modcpy(&self.a,&mut w1); modsqr(&mut w1);
	modcpy(&self.b,&mut w2); modsqr(&mut w2);
	modadd(&w2,&mut w1);
	modinv(h,&mut w1);
	modmul(&w1,&mut self.a);
	modneg(&mut w1);
	modmul(&w1,&mut self.b);
    }

    pub fn qr(&self,h:Option<&[SPINT]>) -> bool {
	let mut w1: [SPINT; NLIMBS] = [0; NLIMBS];	
	let mut w2: [SPINT; NLIMBS] = [0; NLIMBS];
	modcpy(&self.a,&mut w1); modsqr(&mut w1);
	modcpy(&self.b,&mut w2); modsqr(&mut w2);	
	modadd(&w2,&mut w1);
	return modqr(h,&w1);
    }

    pub fn sqrt(&mut self,h:Option<&[SPINT]>) {
	let mut w1: [SPINT; NLIMBS] = [0; NLIMBS];	
	let mut w2: [SPINT; NLIMBS] = [0; NLIMBS];
	let mut w3: [SPINT; NLIMBS] = [0; NLIMBS];	
	let mut w4: [SPINT; NLIMBS] = [0; NLIMBS];
	let mut hint: [SPINT; NLIMBS] = [0; NLIMBS];
	if self.is0() || self.is1() {
	    return;
	}
	modcpy(&self.a,&mut w1); modsqr(&mut w1);
	modcpy(&self.b,&mut w2); modsqr(&mut w2);	
	modadd(&w2,&mut w1);
	modsqrt(h,&mut w1);         // sqrt(a^2+b^2)  - could make use of an input hint to avoid exp!
	modcpy(&self.a,&mut w2); modadd(&w1,&mut w2);
	modhaf(&mut w2);

	modcpy(&self.b,&mut w1);
	modhaf(&mut w1);

	modpro(&w2,&mut hint);
	let qr=modqr(Some(&hint),&w2) as usize;

	modcpy(&hint,&mut w3); modneg(&mut w3);
	modcpy(&w2,&mut w4); modneg(&mut w4);

	modcmv(1-qr,&w4,&mut w2);
	modcmv(1-qr,&w3,&mut hint);

	modcpy(&w2,&mut self.a); modsqrt(Some(&hint),&mut self.a);
	modcpy(&w2,&mut w3); modinv(Some(&hint),&mut w3);
	modmul(&self.a,&mut w3);
	modcpy(&w3,&mut self.b); modmul(&w1,&mut self.b);
	modcpy(&self.a,&mut w4);

	modcmv(1-qr,&self.b,&mut self.a);
	modcmv(1-qr,&w4,&mut self.b);

	let sgn=self.sign();
	let mut nz=Fp2::new_copy(&self);
	nz.neg();
	self.cmv(sgn,&nz);
    }
}

// ------------------------------------

// Give it a good work out

fn printhex(array: &[u8]) {
    for i in 0..array.len() {
        print!("{:02X}", array[i])
    }
    println!("")
}

fn main() {
    if PMOD8 == 1 || PMOD8 == 5 {
	println!("The prime modulus must be 3 mod 4");
	return;
    }
    let mut buffr:[u8;NBYTES]=[0;NBYTES];
    let mut buffi:[u8;NBYTES]=[0;NBYTES]; 
    
    let mut x = Fp2::new_ints(0x18,0x34);
    let mut y = Fp2::new_ints(0x87,0x22);
    let mut z = Fp2::new();
    let mut w1 = Fp2::new();
    let mut w2 = Fp2::new();

    let mut hint: [SPINT; NLIMBS] = [0; NLIMBS];
    for _ in 0..10000 {
        w1.copy(&x); w1.add(&y);
        w2.copy(&x); w2.sub(&y);
        z.copy(&w1); z.mul(&w2); // z=(x+y)(x-y)
        let sgn1=z.sign();
        z.sqr();                 // z=[(x+y)(x-y)]^2

        let t=z.size();          // t=|x|
        modpro(&t,&mut hint);    // create a hint for z
        z.inv(Some(&hint));      // z=1/[(x+y)(x-y)]^2
        z.half();
        z.dbl();
    
        z.sqrt(None);            // z=1/[(x+y)(x-y)]
        z.inv(None);             // z=(x+y)(x-y)
        let sgn2=z.sign();
        if sgn1 != sgn2 {
            z.neg();             // wrong root
        }
        
        w1.copy(&y); w1.sqr();   // w1=y^2
        z.add(&w1);              // z=x^2
        w2.copy(&x); w2.sqr();
        if z.cmp(&w2) {
            println!("passed");
        } else {
            println!("failed");
            return;
        }
	x.copy(&w1);
	y.copy(&z);
    }
    z.export(&mut buffr,&mut buffi);
    printhex(&buffr);
    printhex(&buffi);
}
