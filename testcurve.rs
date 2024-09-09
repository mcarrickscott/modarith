#![allow(non_snake_case)]
#![allow(dead_code)]
#![allow(non_upper_case_globals)]
// For example: cargo run --release --features ED25519 --bin testcurve
// test and timing program for Edwards and Weierstrass curves
//
// run curve_rust.py script to complete edwards.rs or weierstrass.rs for ED25519
// For example python curve_rust.py 64 ED25519

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

use std::time::Instant;

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
fn from_hex(s: &str,x: &mut[u8]) {
    let mut pad:[u8;2*NBYTES]=[0;2*NBYTES];
    let c=s.as_bytes();
    let len=c.len();
    let mut lz=2*NBYTES-len;
    if 2*NBYTES<len {lz=0;}
    for i in 0..lz {
        pad[i]='0' as u8;
    }
    for i in lz..2*NBYTES {
        pad[i]=c[i-lz];
    }

    for i in 0..NBYTES {
        x[i]=char2int(pad[2*i])*16+char2int(pad[2*i+1]);
    }
}

fn printhex(array: &[u8]) {
    for i in 0..array.len() {
        print!("{:02X}", array[i])
    }
    println!("")
}

// output a point (x,y)
fn outputxy(P: &mut ECP) {
    if ecnisinf(P) {
        println!("P= O");
    } else {
        let mut x: [u8; NBYTES] = [0; NBYTES]; 
        let mut y: [u8; NBYTES] = [0; NBYTES]; 
#[cfg(feature="WEIERSTRASS")]
        ecnget(P,&mut x, Some(&mut y));
#[cfg(feature="EDWARDS")]
        ecnget(P,Some(&mut x), Some(&mut y));
        print!("Px= "); printhex(&x); 
        print!("Py= "); printhex(&y); 
    }
}
#[cfg(feature="ED25519")]
const ORDER:&str="1000000000000000000000000000000014DEF9DEA2F79CD65812631A5CF5D3ED";  
#[cfg(feature="ED25519")] 
const r1:&str="66876CB6C86C76660666789A376F6790956A0D6A507657196D75D610E0C9D7B";
#[cfg(feature="ED25519")]
const r2:&str="9978934937938999F9998765C8909870B885907FDF03764C13B05B94EE93672";
#[cfg(feature="ED25519")]
const n1:&str="20347457078878f77b707c070707077a07707b7b07070707223252357134272";
#[cfg(feature="ED25519")]
const n2:&str="35279279432f249b298a876788d86294e02842092769136c086038b1812383a";

#[cfg(feature="ED448")]
const ORDER:&str="3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF7CCA23E9C44EDB49AED63690216CC2728DC58F552378C292AB5844F3";
#[cfg(feature="ED448")]
const r1:&str=   "1F6868465567867838578786787978654735836726356262066876CB6C86C76660666789A376F6790956A0D6A507657196D75D610E0C9D7B";
#[cfg(feature="ED448")]
const r2:&str=   "209797B9AA987987C7A878798786879AB8CA7C98D9CA9D9DF997893410435C8363E873C00B5F40171816219BE8BE29E38CA165319D4BA778";
#[cfg(feature="ED448")]
const n1:&str=   "21347457078878f77b707c070707077a07707b7b0707070722325235713427294582948924358f9898529852985298085690868659860866";
#[cfg(feature="ED448")]
const n2:&str=   "35279279432f249b298a876788d86294e02842092769136c086038b1812383a5875793576874589355835357878467862937894039158588";

#[cfg(feature="NIST256")]
const ORDER:&str="FFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551"; // order
#[cfg(feature="NIST256")]
const r1:&str=  "166876CB6C86C76660666789A376F6790956A0D6A507657196D75D610E0C9D7B"; // random
#[cfg(feature="NIST256")]
const r2:&str=   "E99789339379389A9F9998765C890986B39059D7021039135CE26D61EE5687D6"; // order-r1
#[cfg(feature="NIST256")]
const n1:&str=   "C20347457078878f77b707c070707077a07707b7b07070707223252357134272"; // random
#[cfg(feature="NIST256")]
const n2:&str=  "D35279279432f249b298a876788d86294e02842092769136c086038b1812383a"; // random

#[cfg(feature="NUMS256W")]
const ORDER:&str="FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE43C8275EA265C6020AB20294751A825";
#[cfg(feature="NUMS256W")]    
const r1:&str=  "166876CB6C86C76660666789A376F6790956A0D6A507657196D75D610E0C9D7B";
#[cfg(feature="NUMS256W")]
const r2:&str=   "E9978934937938999F9998765C890986DAE5E19F451EF6EE89D3C2C839450AAA";
#[cfg(feature="NUMS256W")]
const n1:&str=   "120347457078878f77b707c070707077a07707b7b07070707223252357134272";
#[cfg(feature="NUMS256W")]
const n2:&str=  "235279279432f249b298a876788d86294e02842092769136c086038b1812383a";    

#[cfg(feature="NUMS256E")]
const ORDER:&str="4000000000000000000000000000000041955AA52F59439B1A47B190EEDD4AF5";
#[cfg(feature="NUMS256E")]
const r1:&str=  "166876CB6C86C76660666789A376F6790956A0D6A507657196D75D610E0C9D7B";
#[cfg(feature="NUMS256E")]
const r2:&str=   "29978934937938999F9998765C890987383EB9CE8A51DE298370542FE0D0AD7A";
#[cfg(feature="NUMS256E")]
const n1:&str=   "21347457078878f77b707c070707077a07707b7b070707072232523571342729";
#[cfg(feature="NUMS256E")]
const n2:&str=  "35279279432f249b298a876788d86294e02842092769136c086038b1812383a5";


#[cfg(feature="SECP256K1")]
const ORDER:&str="FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141";
#[cfg(feature="SECP256K1")]    
const r1:&str=  "166876CB6C86C76660666789A376F6790956A0D6A507657196D75D610E0C9D7B";
#[cfg(feature="SECP256K1")]
const r2:&str=   "E9978934937938999F9998765C890985B1583C100A413ACA28FB012BC229A3C6";
#[cfg(feature="SECP256K1")]
const n1:&str=   "120347457078878f77b707c070707077a07707b7b07070707223252357134272";
#[cfg(feature="SECP256K1")]
const n2:&str=  "235279279432f249b298a876788d86294e02842092769136c086038b1812383a";   


const NBYTES: usize=ORDER.len()/2;

fn main() {
 
    let mut r:[u8;NBYTES]=[0;NBYTES];
    let mut a:[u8;NBYTES]=[0;NBYTES];    
    let mut b:[u8;NBYTES]=[0;NBYTES];       
    let mut P=ECP::new();
    let mut Q=ECP::new();
// convert to byte array
    from_hex(&ORDER,&mut r);
    ecngen(&mut P);
    outputxy(&mut P);
    ecncpy(&P,&mut Q); 
    
    ecnmul(&r,&mut P);
    if ecnisinf(&P) {
        println!("MUL test passed OK");
    } else {
        println!("MUL test FAILED");
    }
    from_hex(&r1,&mut a);
    from_hex(&r2,&mut b);
    P=ecnmul2(&a,&Q,&b,&Q);
    if ecnisinf(&P) {
        println!("MUL2 test passed OK");
    } else {
        println!("MUL2 test FAILED");
    }    
    
    from_hex(&n1,&mut a); // random
    from_hex(&n2,&mut b);

    ecncpy(&Q,&mut P);    
    
    println!("Timing point multiplication");
    unsafe {
        let pre = std::arch::x86_64::_rdtsc();
        let begin = Instant::now();
    	for _ in 0..10000 {
        	ecnmul(&a,&mut P);
    	}        
	let post = std::arch::x86_64::_rdtsc();
	let elapsed = begin.elapsed();
	let dur = (elapsed.as_secs() * 1_000_000_000) + (elapsed.subsec_nanos()) as u64;      
	println!("Clock cycles= {} Microsecs= {}",((post-pre)/10000),dur/10000000);  
   	outputxy(&mut P);	
    }	
    println!("Timing double point multiplication");
    unsafe {
        let pre = std::arch::x86_64::_rdtsc();
        let begin = Instant::now();
    	for _ in 0..10000 {
            P=ecnmul2(&a,&P,&b,&Q);
    	}        
	let post = std::arch::x86_64::_rdtsc();
	let elapsed = begin.elapsed();
	let dur = (elapsed.as_secs() * 1_000_000_000) + (elapsed.subsec_nanos()) as u64;   
	println!("Clock cycles= {} Microsecs= {}",((post-pre)/10000),dur/10000000);  
   	outputxy(&mut P);	   
    }
}
