// Some useful hash functions

// SHA2

const SHA256_H0: u32 = 0x6A09E667;
const SHA256_H1: u32 = 0xBB67AE85;
const SHA256_H2: u32 = 0x3C6EF372;
const SHA256_H3: u32 = 0xA54FF53A;
const SHA256_H4: u32 = 0x510E527F;
const SHA256_H5: u32 = 0x9B05688C;
const SHA256_H6: u32 = 0x1F83D9AB;
const SHA256_H7: u32 = 0x5BE0CD19;

const SHA256_K: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

pub struct SHA256 {
    length: [u32; 2],
    h: [u32; 8],
    w: [u32; 64],
}

impl SHA256 {
    fn s(n: u32, x: u32) -> u32 {
        ((x) >> n) | ((x) << (32 - n))
    }
    fn r(n: u32, x: u32) -> u32 {
        (x) >> n
    }

    fn ch(x: u32, y: u32, z: u32) -> u32 {
        (x & y) ^ (!(x) & z)
    }

    fn maj(x: u32, y: u32, z: u32) -> u32 {
        (x & y) ^ (x & z) ^ (y & z)
    }
    fn sig0(x: u32) -> u32 {
        SHA256::s(2, x) ^ SHA256::s(13, x) ^ SHA256::s(22, x)
    }

    fn sig1(x: u32) -> u32 {
        SHA256::s(6, x) ^ SHA256::s(11, x) ^ SHA256::s(25, x)
    }

    fn theta0(x: u32) -> u32 {
        SHA256::s(7, x) ^ SHA256::s(18, x) ^ SHA256::r(3, x)
    }

    fn theta1(x: u32) -> u32 {
        SHA256::s(17, x) ^ SHA256::s(19, x) ^ SHA256::r(10, x)
    }

    pub fn as_bytes(&self, array: &mut [u8]) {
        let mut ptr = 0;
        for i in 0..2 {
            let mut t = self.length[i];
            array[ptr] = (t % 256) as u8;
            t /= 256;
            ptr += 1;
            array[ptr] = (t % 256) as u8;
            t /= 256;
            ptr += 1;
            array[ptr] = (t % 256) as u8;
            t /= 256;
            ptr += 1;
            array[ptr] = t as u8;
            ptr += 1;
        }
        for i in 0..8 {
            let mut t = self.h[i];
            array[ptr] = (t % 256) as u8;
            t /= 256;
            ptr += 1;
            array[ptr] = (t % 256) as u8;
            t /= 256;
            ptr += 1;
            array[ptr] = (t % 256) as u8;
            t /= 256;
            ptr += 1;
            array[ptr] = t as u8;
            ptr += 1;
        }
        for i in 0..64 {
            let mut t = self.w[i];
            array[ptr] = (t % 256) as u8;
            t /= 256;
            ptr += 1;
            array[ptr] = (t % 256) as u8;
            t /= 256;
            ptr += 1;
            array[ptr] = (t % 256) as u8;
            t /= 256;
            ptr += 1;
            array[ptr] = t as u8;
            ptr += 1;
        }
    }

    pub fn from_bytes(&mut self, array: &[u8]) {
        let mut ptr = 0;
        for i in 0..2 {
            let mut t = array[ptr + 3] as u32;
            t = 256 * t + (array[ptr + 2] as u32);
            t = 256 * t + (array[ptr + 1] as u32);
            t = 256 * t + (array[ptr] as u32);
            self.length[i] = t;
            ptr += 4;
        }
        for i in 0..8 {
            let mut t = array[ptr + 3] as u32;
            t = 256 * t + (array[ptr + 2] as u32);
            t = 256 * t + (array[ptr + 1] as u32);
            t = 256 * t + (array[ptr] as u32);
            self.h[i] = t;
            ptr += 4;
        }
        for i in 0..64 {
            let mut t = array[ptr + 3] as u32;
            t = 256 * t + (array[ptr + 2] as u32);
            t = 256 * t + (array[ptr + 1] as u32);
            t = 256 * t + (array[ptr] as u32);
            self.w[i] = t;
            ptr += 4;
        }
    }

    fn transform(&mut self) {
        /* basic transformation step */
        for j in 16..64 {
            self.w[j] = SHA256::theta1(self.w[j - 2])
                .wrapping_add(self.w[j - 7])
                .wrapping_add(SHA256::theta0(self.w[j - 15]))
                .wrapping_add(self.w[j - 16]);
        }
        let mut a = self.h[0];
        let mut b = self.h[1];
        let mut c = self.h[2];
        let mut d = self.h[3];
        let mut e = self.h[4];
        let mut f = self.h[5];
        let mut g = self.h[6];
        let mut hh = self.h[7];
        for j in 0..64 {
            /* 64 times - mush it up */
            let t1 = hh
                .wrapping_add(SHA256::sig1(e))
                .wrapping_add(SHA256::ch(e, f, g))
                .wrapping_add(SHA256_K[j])
                .wrapping_add(self.w[j]);
            let t2 = SHA256::sig0(a).wrapping_add(SHA256::maj(a, b, c));
            hh = g;
            g = f;
            f = e;
            e = d.wrapping_add(t1);
            d = c;
            c = b;
            b = a;
            a = t1.wrapping_add(t2);
        }
        self.h[0] = self.h[0].wrapping_add(a);
        self.h[1] = self.h[1].wrapping_add(b);
        self.h[2] = self.h[2].wrapping_add(c);
        self.h[3] = self.h[3].wrapping_add(d);
        self.h[4] = self.h[4].wrapping_add(e);
        self.h[5] = self.h[5].wrapping_add(f);
        self.h[6] = self.h[6].wrapping_add(g);
        self.h[7] = self.h[7].wrapping_add(hh);
    }

    /* Initialise Hash function */
    pub fn init(&mut self) {
        /* initialise */
        for i in 0..64 {
            self.w[i] = 0
        }
        self.length[0] = 0;
        self.length[1] = 0;
        self.h[0] = SHA256_H0;
        self.h[1] = SHA256_H1;
        self.h[2] = SHA256_H2;
        self.h[3] = SHA256_H3;
        self.h[4] = SHA256_H4;
        self.h[5] = SHA256_H5;
        self.h[6] = SHA256_H6;
        self.h[7] = SHA256_H7;
    }

    pub fn new() -> SHA256 {
        let mut nh = SHA256 {
            length: [0; 2],
            h: [0; 8],
            w: [0; 64],
        };
        nh.init();
        nh
    }

    pub fn new_copy(hh: &SHA256) -> SHA256 {
        let mut nh = SHA256 {
            length: [0; 2],
            h: [0; 8],
            w: [0; 64],
        };
        nh.length[0] = hh.length[0];
        nh.length[1] = hh.length[1];
        for i in 0..64 {
            nh.w[i] = hh.w[i];
        }
        for i in 0..8 {
            nh.h[i] = hh.h[i];
        }
        nh
    }

    /* process a single byte */
    pub fn process(&mut self, byt: u8) {
        /* process the next message byte */
        let cnt = ((self.length[0] / 32) % 16) as usize;
        self.w[cnt] <<= 8;
        self.w[cnt] |= byt as u32;
        self.length[0] += 8;
        if self.length[0] == 0 {
            self.length[1] += 1;
            self.length[0] = 0
        }
        if (self.length[0] % 512) == 0 {
            self.transform()
        }
    }

    /* process an array of bytes */

    pub fn process_array(&mut self, b: &[u8]) {
        for i in 0..b.len() {
            self.process(b[i])
        }
    }

    /* process a 32-bit integer */
    pub fn process_num(&mut self, n: i32) {
        self.process(((n >> 24) & 0xff) as u8);
        self.process(((n >> 16) & 0xff) as u8);
        self.process(((n >> 8) & 0xff) as u8);
        self.process((n & 0xff) as u8);
    }

    /* Generate 32-byte Hash */
    pub fn hash(&mut self) -> [u8; 32] {
        /* pad message and finish - supply digest */
        let mut digest: [u8; 32] = [0; 32];
        let len0 = self.length[0];
        let len1 = self.length[1];
        self.process(0x80);
        while (self.length[0] % 512) != 448 {
            self.process(0)
        }
        self.w[14] = len1;
        self.w[15] = len0;
        self.transform();
        for i in 0..32 {
            /* convert to bytes */
            digest[i] = ((self.h[i / 4] >> (8 * (3 - i % 4))) & 0xff) as u8;
        }
        self.init();
        digest
    }

    pub fn continuing_hash(&self) -> [u8; 32] {
        let mut sh = SHA256::new_copy(self);
        sh.hash()
    }
}

// SHA3

pub const HASH224: usize = 28;
pub const HASH256: usize = 32;
pub const HASH384: usize = 48;
pub const HASH512: usize = 64;
pub const SHAKE128: usize = 16;
pub const SHAKE256: usize = 32;

const ROUNDS: usize = 24;

const RC: [u64; 24] = [
    0x0000000000000001,
    0x0000000000008082,
    0x800000000000808A,
    0x8000000080008000,
    0x000000000000808B,
    0x0000000080000001,
    0x8000000080008081,
    0x8000000000008009,
    0x000000000000008A,
    0x0000000000000088,
    0x0000000080008009,
    0x000000008000000A,
    0x000000008000808B,
    0x800000000000008B,
    0x8000000000008089,
    0x8000000000008003,
    0x8000000000008002,
    0x8000000000000080,
    0x000000000000800A,
    0x800000008000000A,
    0x8000000080008081,
    0x8000000000008080,
    0x0000000080000001,
    0x8000000080008008,
];

pub struct SHA3 {
    length: usize,
    rate: usize,
    len: usize,
    //s: [[u64; 5]; 5],
    s: [u64; 25],
}

impl SHA3 {
    fn rotl(x: u64, n: u64) -> u64 {
        ((x) << n) | ((x) >> (64 - n))
    }

    fn transform(&mut self) {
        for k in 0..ROUNDS {
            let c0 = self.s[0] ^ self.s[5] ^ self.s[10] ^ self.s[15] ^ self.s[20];
            let c1 = self.s[1] ^ self.s[6] ^ self.s[11] ^ self.s[16] ^ self.s[21];
            let c2 = self.s[2] ^ self.s[7] ^ self.s[12] ^ self.s[17] ^ self.s[22];
            let c3 = self.s[3] ^ self.s[8] ^ self.s[13] ^ self.s[18] ^ self.s[23];
            let c4 = self.s[4] ^ self.s[9] ^ self.s[14] ^ self.s[19] ^ self.s[24];

            let d0 = c4 ^ SHA3::rotl(c1, 1);
            let d1 = c0 ^ SHA3::rotl(c2, 1);
            let d2 = c1 ^ SHA3::rotl(c3, 1);
            let d3 = c2 ^ SHA3::rotl(c4, 1);
            let d4 = c3 ^ SHA3::rotl(c0, 1);

            let b00 = self.s[0] ^ d0;
            let b02 = SHA3::rotl(self.s[1] ^ d1, 1);
            let b04 = SHA3::rotl(self.s[2] ^ d2, 62);
            let b01 = SHA3::rotl(self.s[3] ^ d3, 28);
            let b03 = SHA3::rotl(self.s[4] ^ d4, 27);

            let b13 = SHA3::rotl(self.s[5] ^ d0, 36);
            let b10 = SHA3::rotl(self.s[6] ^ d1, 44);
            let b12 = SHA3::rotl(self.s[7] ^ d2, 6);
            let b14 = SHA3::rotl(self.s[8] ^ d3, 55);
            let b11 = SHA3::rotl(self.s[9] ^ d4, 20);

            let b21 = SHA3::rotl(self.s[10] ^ d0, 3);
            let b23 = SHA3::rotl(self.s[11] ^ d1, 10);
            let b20 = SHA3::rotl(self.s[12] ^ d2, 43);
            let b22 = SHA3::rotl(self.s[13] ^ d3, 25);
            let b24 = SHA3::rotl(self.s[14] ^ d4, 39);

            let b34 = SHA3::rotl(self.s[15] ^ d0, 41);
            let b31 = SHA3::rotl(self.s[16] ^ d1, 45);
            let b33 = SHA3::rotl(self.s[17] ^ d2, 15);
            let b30 = SHA3::rotl(self.s[18] ^ d3, 21);
            let b32 = SHA3::rotl(self.s[19] ^ d4, 8);

            let b42 = SHA3::rotl(self.s[20] ^ d0, 18);
            let b44 = SHA3::rotl(self.s[21] ^ d1, 2);
            let b41 = SHA3::rotl(self.s[22] ^ d2, 61);
            let b43 = SHA3::rotl(self.s[23] ^ d3, 56);
            let b40 = SHA3::rotl(self.s[24] ^ d4, 14);

            self.s[0] = b00 ^ (!b10 & b20);
            self.s[1] = b10 ^ (!b20 & b30);
            self.s[2] = b20 ^ (!b30 & b40);
            self.s[3] = b30 ^ (!b40 & b00);
            self.s[4] = b40 ^ (!b00 & b10);

            self.s[5] = b01 ^ (!b11 & b21);
            self.s[6] = b11 ^ (!b21 & b31);
            self.s[7] = b21 ^ (!b31 & b41);
            self.s[8] = b31 ^ (!b41 & b01);
            self.s[9] = b41 ^ (!b01 & b11);

            self.s[10] = b02 ^ (!b12 & b22);
            self.s[11] = b12 ^ (!b22 & b32);
            self.s[12] = b22 ^ (!b32 & b42);
            self.s[13] = b32 ^ (!b42 & b02);
            self.s[14] = b42 ^ (!b02 & b12);

            self.s[15] = b03 ^ (!b13 & b23);
            self.s[16] = b13 ^ (!b23 & b33);
            self.s[17] = b23 ^ (!b33 & b43);
            self.s[18] = b33 ^ (!b43 & b03);
            self.s[19] = b43 ^ (!b03 & b13);

            self.s[20] = b04 ^ (!b14 & b24);
            self.s[21] = b14 ^ (!b24 & b34);
            self.s[22] = b24 ^ (!b34 & b44);
            self.s[23] = b34 ^ (!b44 & b04);
            self.s[24] = b44 ^ (!b04 & b14);

            self.s[0] ^= RC[k];
        }
    }

    /* Initialise Hash function */
    pub fn init(&mut self, olen: usize) {
        /* initialise */
        for i in 0..25 {
            self.s[i] = 0;
        }
        self.length = 0;
        self.len = olen;
        self.rate = 200 - 2 * olen;
    }

    pub fn new(olen: usize) -> SHA3 {
        let mut nh = SHA3 {
            length: 0,
            rate: 0,
            len: 0,
            s: [0; 25],
        };
        nh.init(olen);
        nh
    }

    pub fn new_copy(hh: &SHA3) -> SHA3 {
        let mut nh = SHA3 {
            length: 0,
            rate: 0,
            len: 0,
            s: [0; 25],
        };
        nh.length = hh.length;
        nh.len = hh.len;
        nh.rate = hh.rate;
        for i in 0..25 {
            nh.s[i] = hh.s[i];
        }
        nh
    }

    /* process a single byte */
    pub fn process(&mut self, byt: u8) {
        /* process the next message byte */
        let cnt = self.length;
        let b = cnt % 8;
        let ind = cnt / 8;
        self.s[ind] ^= (byt as u64) << (8 * b);
        self.length += 1;
        if self.length == self.rate {
            self.length = 0;
            self.transform();
        }
    }

    /* process an array of bytes */

    pub fn process_array(&mut self, b: &[u8]) {
        for i in 0..b.len() {
            self.process(b[i])
        }
    }

    /* process a 32-bit integer */
    pub fn process_num(&mut self, n: i32) {
        self.process(((n >> 24) & 0xff) as u8);
        self.process(((n >> 16) & 0xff) as u8);
        self.process(((n >> 8) & 0xff) as u8);
        self.process((n & 0xff) as u8);
    }

    pub fn squeeze(&mut self, buff: &mut [u8], olen: usize) {
        let mut m = 0;
        let nb = olen / self.rate;

        for _ in 0..nb {
            for i in 0..self.rate / 8 {
                let mut el = self.s[i];
                for _ in 0..8 {
                    buff[m] = (el & 0xff) as u8;
                    m += 1;
                    el >>= 8;
                }
            }
            self.transform();
        }

        let mut i = 0;
        while m < olen {
            let mut el = self.s[i];
            i += 1;
            for _ in 0..8 {
                buff[m] = (el & 0xff) as u8;
                m += 1;
                if m >= olen {
                    break;
                }
                el >>= 8;
            }
        }
    }

    /* Generate 32-byte Hash */
    pub fn hash(&mut self, digest: &mut [u8]) {
        /* pad message and finish - supply digest */
        let q = self.rate - self.length;
        if q == 1 {
            self.process(0x86);
        } else {
            self.process(0x06);
            while self.length != self.rate - 1 {
                self.process(0x00)
            }
            self.process(0x80);
        }
        let hlen = self.len;
        self.squeeze(digest, hlen);
    }

    pub fn continuing_hash(&mut self, digest: &mut [u8]) {
        let mut sh = SHA3::new_copy(self);
        sh.hash(digest)
    }

    pub fn shake(&mut self, digest: &mut [u8], olen: usize) {
        let q = self.rate - self.length;
        if q == 1 {
            self.process(0x9f);
        } else {
            self.process(0x1f);
            while self.length != self.rate - 1 {
                self.process(0x00)
            }
            self.process(0x80);
        }
        self.squeeze(digest, olen);
    }

    pub fn continuing_shake(&mut self, digest: &mut [u8], olen: usize) {
        let mut sh = SHA3::new_copy(self);
        sh.shake(digest, olen);
    }
}

