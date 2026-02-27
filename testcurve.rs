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
    let mut pad:[u8;2*FBYTES]=[0;2*FBYTES];
    let c=s.as_bytes();
    let len=c.len();
    let mut lz=2*FBYTES-len;
    if 2*FBYTES<len {lz=0;}
    for i in 0..lz {
        pad[i]='0' as u8;
    }
    for i in lz..2*FBYTES {
        pad[i]=c[i-lz];
    }

    for i in 0..FBYTES {
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
        let mut x: [u8; FBYTES] = [0; FBYTES]; 
        let mut y: [u8; FBYTES] = [0; FBYTES]; 
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

#[cfg(feature="NIST384")]
const ORDER:&str="ffffffffffffffffffffffffffffffffffffffffffffffffc7634d81f4372ddf581a0db248b0a77aecec196accc52973"; // order
#[cfg(feature="NIST384")]
const r1:&str=  "bd9c66b3ad3c2d6d1a3d1fa7bc8960a923b8c1e9392456de3eb13b9046685257bdd640fb06671ad11c80317fa3b1799d"; // random
#[cfg(feature="NIST384")]
const r2:&str=   "4263994c52c3d292e5c2e05843769f56dc473e16c6dba92188b211f1adcedb879a43ccb742498ca9d06be7eb2913afd6"; // order-r1
#[cfg(feature="NIST384")]
const n1:&str=   "9a1de644815ef6d13b8faa1837f8a88b17fc695a07a0ca6e0822e8f36c031199972a846916419f828b9d2434e465e150"; // random
#[cfg(feature="NIST384")]
const n2:&str=  "4737819096da1dac72ff5d2a386ecbe06b65a6a48b8148f6b38a088ca65ed389b74d0fb132e706298fadc1a606cb0fb3"; // random

#[cfg(feature="NIST521")]
const ORDER:&str="1fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffa51868783bf2f966b7fcc0148f709a5d03bb5c9b8899c47aebb6fb71e91386409"; // order
#[cfg(feature="NIST521")]
const r1:&str=  "d8972a846916419f828b9d2434e465e150bd9c66b3ad3c2d6d1a3d1fa7bc8960a923b8c1e9392456de3eb13b9046685257bdd640fb06671ad11c80317fa3b1799d"; // random
#[cfg(feature="NIST521")]
const r2:&str=   "12768d57b96e9be607d7462dbcb1b9a1eaf4263994c52c3d292e5c2e05843769f512dcdc59a860b3f8d411ac5b8b0a153787ddf88bd83352cdd9eef859eed86ea6c"; // order-r1
#[cfg(feature="NIST521")]
const n1:&str=   "e5386ecbe06b65a6a48b8148f6b38a088ca65ed389b74d0fb132e706298fadc1a606cb0fb39a1de644815ef6d13b8faa1837f8a88b17fc695a07a0ca6e0822e8f3"; // random
#[cfg(feature="NIST521")]
const n2:&str=  "acc37459eef50bea63371ecd7b27cd813047229389571aa8766c307511b2b9437a28df6ec4ce4a2bbdc241330b01a9e71fde8a774bcf36d58b4737819096da1dac"; // random

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

#[cfg(feature="SQISIGN_1")]
const ORDER:&str="13FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF098677E8D0D856DA332BA970DCFDEA1";
#[cfg(feature="SQISIGN_1")]
const r1:&str=  "10876CB6C86C76660666789A376F6790956A0D6A507657196D75D610E0C9D7B";  
#[cfg(feature="SQISIGN_1")]
const r2:&str=   "378934937938999f9998765c890986e741c6a7e8061ffc0c5b5d35ffc34126"; 
#[cfg(feature="SQISIGN_1")]
const n1:&str=   "7347457078878f77b707c070707077a07707b7b07070707223252357134272"; 
#[cfg(feature="SQISIGN_1")]
const n2:&str=  "3279279432f249b298a876788d86294e02842092769136c086038b1812383a";   

#[cfg(feature="SQISIGN_3")]
const ORDER:&str="104000000000000000000000000000000000000000000000303A69B3514879CD109A98F29F0D04F09F855D4F3C6A7037";
#[cfg(feature="SQISIGN_3")]
const r1:&str=  "6076CB6C86C76660666789A376F6790956A0D6A507657196D75D610E0C9D7B87508750F98765C890986DAE5E19F451E";  
#[cfg(feature="SQISIGN_3")]
const r2:&str=   "a38934937938999f9998765c890986f6a95f295af89a8e6c2c493a2707ea2149b9223e30696a86795fe82695acb2b19";    
#[cfg(feature="SQISIGN_3")]
const n1:&str=   "93347457078878f77b707c070707077a07707b7b070707072232523571342720986DAE5E19F451EF6EE89D3C2C0986D";    
#[cfg(feature="SQISIGN_3")]
const n2:&str=  "235279279432f249b298a876788d86294e02842092769136c086038b1812383a13427294582948924358f98985956A0";

#[cfg(feature="SQISIGN_5")]
const ORDER:&str="6C00000000000000000000000000000000000000000000000000000000000002C8858DA0CB07C5ABCADABC1BEE86F8C9101174D8A115AD57E5F0228C2D0871";
#[cfg(feature="SQISIGN_5")]
const r1:&str=  "47235279279432f249b298a876788d86294e02842092769136c086038b1812383a8765C890985B1583C100A413ACA28FB0654735836726356262066876CB6C";  
#[cfg(feature="SQISIGN_5")]
const r2:&str=   "24dcad86d86bcd0db64d675789877279d6b1fd7bdf6d896ec93f79fc74e7edca8dfe27d83a6f6a964719bb77dada56395fac2da31dae8722838e1c23b63d05";    
#[cfg(feature="SQISIGN_5")]
const n1:&str=   "32120347457078878f77b707c070707077a07707b7b0707070722325235713427270707077a07707b7b07070707223252307707b7b070707072232523688A3";   
#[cfg(feature="SQISIGN_5")] 
const n2:&str=  "19235279279432f249b298a876788d86294e02842092769136c086038b1812383a13427294582948924358f98985956A077b707c070707077a07707b7b0707";

fn main() {
 
    let mut r:[u8;FBYTES]=[0;FBYTES];
    let mut a:[u8;FBYTES]=[0;FBYTES];    
    let mut b:[u8;FBYTES]=[0;FBYTES];       
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
