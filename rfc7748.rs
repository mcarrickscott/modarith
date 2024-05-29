// Program to implement RFC7748 - https://datatracker.ietf.org/doc/html/rfc7748
// Montgomery curve key exchange code, as used by TLS
// Use associated python scripts to generate code for C25519 or C448, but easily modified for other Montgomery curves
//
// A good Montgomery curve can be found by running the sagemath script bowe.sage
//
// Mike Scott 23rd November 2023
// TII
// code for 16/32/64-bit processor for C25519 curve can be generated  by 
//
// python pseudo_rust.py 16/32/64 C25519
// or
// python monty_rust.py 16/32/64 C25519
//
// code for 16/32/64-bit processor for C448 curve can be generated  by
//
// python monty_rust.py 16/32/64 C448

// rustc -O rfc7748.rs

// NOTE: Code is set up for C25519. See // ****** for places where changes are needed for X448

/*** Insert automatically generated code for modulus code.rs here ***/



/*** End of automatically generated code ***/

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
fn from_hex(s: &str,x: &mut[u8]) -> usize {
    let c=s.as_bytes();
    let mut i=0;
    let mut j=0;

    while j<c.len() {
        x[i]=char2int(c[j])*16+char2int(c[j+1]);
        i+=1;
        j+=2;
    }
    return i;
}

fn printhex(array: &[u8]) {
    for i in 0..array.len() {
        print!("{:02X}", array[i])
    }
    println!("")
}

// output a modulo number in hex
#[allow(dead_code)]
fn output(x: &[SPINT]) {
    let mut b:[u8;NBYTES]=[0;NBYTES];
    let mut c:[SPINT;NLIMBS]=[0;NLIMBS];
    modcpy(x,&mut c);
    modexp(&c,&mut b);
    printhex(&b);
}

// Describe Montgomery Curve parameters
// ******
//X25519
const A24: usize = 121665;
const COF: usize = 3;
const TWIST_SECURE: bool = true;

//X448
//const A24: usize = 39081;
//const COF: usize = 2;
//const TWIST_SECURE: bool = true;

// clamp input - see RFC7748
fn clamp(bk: &mut [u8]) {
    let s=(8-(NBITS%8))%8;
    bk[0]&=(-(1<<COF)) as u8;
    let mask=(0xff as u8)>>s;
    bk[NBYTES-1]&=mask;
    bk[NBYTES-1]|=(0x80 as u8)>>s;
}

// extract a bit from a byte array
fn bit(n: usize,a: &[u8]) -> isize {
    return ((a[n / 8] & (1 << (n % 8))) >> (n%8)) as isize;
}

// RFC7748 - Montgomery curve
// bv=bk*bu, bu,bv are x coordinates on elliptic curve
fn rfc7748(bk: &[u8],bu: &[u8],bv: &mut [u8]) {
    let mut swap = 0 as isize;
    let mut ck:[u8;NBYTES]=[0;NBYTES];
    let mut u:[SPINT;NLIMBS]=[0;NLIMBS];
    let mut x1:[SPINT;NLIMBS]=[0;NLIMBS];
    let mut x2:[SPINT;NLIMBS]=[0;NLIMBS];
    let mut x3:[SPINT;NLIMBS]=[0;NLIMBS];
    let mut z2:[SPINT;NLIMBS]=[0;NLIMBS];
    let mut z3:[SPINT;NLIMBS]=[0;NLIMBS];
    let mut a:[SPINT;NLIMBS]=[0;NLIMBS];
    let mut b:[SPINT;NLIMBS]=[0;NLIMBS];
    let mut aa:[SPINT;NLIMBS]=[0;NLIMBS];
    let mut bb:[SPINT;NLIMBS]=[0;NLIMBS];
    let mut c:[SPINT;NLIMBS]=[0;NLIMBS];
    let mut d:[SPINT;NLIMBS]=[0;NLIMBS];
    let mut e:[SPINT;NLIMBS]=[0;NLIMBS];

    for i in 0..NBYTES {
        ck[i]=bk[i];
    }
// clamp input
    clamp(&mut ck);

// import into internal representation
    modimp(&bu,&mut u);

    modcpy(&u,&mut x1);  // x_1=u
    modone(&mut x2);     // x_2=1
    modzer(&mut z2);     // z_2=0
    modcpy(&u,&mut x3);  // x_3=u
    modone(&mut z3);     // z_3=1

    for i in (0..NBITS).rev() {
        let kt=bit(i,&ck);
        swap^=kt;
        modcsw(swap,&mut x2,&mut x3);
        modcsw(swap,&mut z2,&mut z3);
        swap=kt;
            
        modcpy(&x2,&mut a); modadd(&z2,&mut a);   // A = x_2 + z_2
        modcpy(&a,&mut aa); modsqr(&mut aa);      // AA = A^2
        modcpy(&x2,&mut b); modsub(&z2,&mut b);   // B = x_2 - z_2
        modcpy(&b,&mut bb); modsqr(&mut bb);      // BB = B^2

        modcpy(&aa,&mut e); modsub(&bb,&mut e);   // E = AA - BB
        modcpy(&x3,&mut c); modadd(&z3,&mut c);   // C = x_3 + z_3

        modcpy(&x3,&mut d); modsub(&z3,&mut d);   // D = x_3 - z_3
        modmul(&a,&mut d);  // DA = D * A
        modmul(&b,&mut c);  // CB = C * B
 
        modcpy(&d,&mut x3); modadd(&c,&mut x3); modsqr(&mut x3); // x_3 = (DA + CB)^2

        modcpy(&d,&mut z3); modsub(&c,&mut z3); modsqr(&mut z3); modmul(&x1,&mut z3); // z_3 = x_1 * (DA - CB)^2
        modcpy(&aa,&mut x2); modmul(&bb,&mut x2); // x_2 = AA * BB
       
        modcpy(&e,&mut z2);
        modmli(A24,&mut z2); 
              
        modadd(&aa,&mut z2); modmul(&e,&mut z2);   // z_2 = E * (AA + a24 * E)
    }
    modcsw(swap,&mut x2,&mut x3);
    modcsw(swap,&mut z2,&mut z3);

    if TWIST_SECURE {
        modpro(&z2,&mut a);       
        modinv(Some(&a),&mut z2);    // sufficient for twist secure curves like C25519 and C448 
    } else { // Do cheap point validation here - see https://eprint.iacr.org/2020/1497
        modcpy(&u,&mut b); modmul(&z2,&mut b);
        modcpy(&b,&mut a); modmul(&z2,&mut a);
        modpro(&a,&mut e);
        modcpy(&a,&mut c);
        modcpy(&e,&mut d); modmul(&z2,&mut d);
        modsqr(&mut d);
        modmul(&u,&mut d);

        for _ in 0..COF-2 {  // COF is 2 or 3
            modsqr(&mut c);
            modmul(&a,&mut c);
        }
        for _ in 0..COF {
            modsqr(&mut e);
        }
        modmul(&e,&mut c);
        modcpy(&c,&mut z2); modmul(&b,&mut z2);

        for _ in 0..COF-2 {
            modsqr(&mut d);
        }
        modone(&mut a); modadd(&a,&mut d); modfsb(&mut d); modshr(1,&mut d); // 1 for QR, else 0
        modmul(&d,&mut x2); // set to zero for bad input point
    }
    modmul(&z2,&mut x2);   

    modexp(&x2,bv);
    bv.reverse(); // little endian
}

// a test vector for x25519 or x448 from RFC7748
fn main() {
// ******
//C25519
    const SK:&str="77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a";
    const SU:&str="0000000000000000000000000000000000000000000000000000000000000009";
//C448
//    const SK:&str="9a8f4925d1519f5775cf46b04b5800d4ee9ee8bae8bc5565d498c28dd9c9baf574a9419744897391006382a6f127ab1d9ac2d8c0a598726b";
//    const SU:&str="0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000005";

    let mut bk:[u8;NBYTES]=[0;NBYTES];
    let mut bu:[u8;NBYTES]=[0;NBYTES];
    let mut bv:[u8;NBYTES]=[0;NBYTES];

// convert to byte arrays
    from_hex(&SU,&mut bu);
    from_hex(&SK,&mut bk);

    rfc7748(&bk,&bu,&mut bv);

    println!("{}",&SK);
    printhex(&bv);

	unsafe {
	    let pre = std::arch::x86_64::_rdtsc();
	    for _i in 0..5000 {
            rfc7748(&bk,&bu,&mut bv);
            rfc7748(&bk,&bv,&mut bu);
	    }
	    let post = std::arch::x86_64::_rdtsc();
	    println!("Clock cycles= {}",((post-pre)/10_000));
        printhex(&bu);
	}
}
