// Define required intrinsics depending on chosen engine, in case AVX-512/IFMA
// Using n 64-bit lanes in a n*64 bit wide SIMD extension, where n is 2, 4 or 8

#ifndef INTRINSICS_H
#define INTRINSICS_H

// Key
// S - Signed
// U - Unsigned
// SU - Signed to Unsigned
// Signed addition/subtraction are in fact the same operation. But parameter types will be different.

#if SIMD_ENGINE == AVX128_IFMA

#define MR_ZERO() _mm_setzero_si128()
#define MR_AND(X,Y) _mm_and_si128(X,Y)
#define MR_ANDNOT(X,Y) _mm_andnot_si128(X,Y)
#define MR_OR(X,Y) _mm_or_si128(X,Y)
#define MR_XOR(X,Y) _mm_xor_si128(X,Y)
#define MR_SET_ALL_LANES_TO_CONSTANT(C) _mm_set1_epi64x(C)
#define MR_SET_EACH_LANE_TO_CONSTANT(C) _mm_set_epi64x(C[1],C[0])
#define MR_ADD64U(X,Y) _mm_add_epi64(X,Y)
#define MR_MUL_DOUBLE_ADD(A,B,C,D) _mm_madd52lo_epu64(_mm_madd52lo_epu64(_mm_setzero_si128(),A,B),C,D)
#define MR_MULADD_HI(T,X,Y) _mm_madd52hi_epu64(T,X,Y)
#define MR_MULADD_LO(T,X,Y) _mm_madd52lo_epu64(T,X,Y)
#define MR_MULADD_LO_64(T,X,Y) _mm_add_epi64(T,_mm_mullo_epi64(X,Y))
#define MR_ADDSHL(T,X,S) _mm_add_epi64(T,_mm_and_si128(_mm_slli_epi64(X,S),_mm_set1_epi64x(((int64_t)1<<52)-1)) )
#define MR_ADDSHR(T,X,S) _mm_add_epi64(T,_mm_srli_epi64(X,52-S))
#define MR_SHR52(LO,HI) _mm_add_epi64(_mm_srai_epi64(LO,52),HI)
#define MR_MULHI52(X,Y) _mm_madd52hi_epu64(_mm_setzero_si128(),X,Y)    
#define MR_MULLO52(X,Y) _mm_madd52lo_epu64(_mm_setzero_si128(),X,Y)
#define MR_MULLO(X,Y) _mm_mullo_epi64(X,Y)
#define MR_SHL(L,H,S) _mm_add_epi64(_mm_srai_epi64(L,S),_mm_slli_epi64(H,(52-S)))
#define MR_SUB64U(X,Y) _mm_sub_epi64(X,Y)
#define MR_SHR64U(X,C) _mm_srli_epi64(X,C)
#define MR_SHL64U(X,C) _mm_slli_epi64(X,C)
#define MR_SRA64S(X,C) _mm_srai_epi64(X,C)
#define MR_EXTRACT(X,L) ((unsigned char)(_mm_extract_epi16(X,L*4)&0xFF))

#endif


#if SIMD_ENGINE == AVX256_IFMA

#define MR_ZERO() _mm256_setzero_si256()
#define MR_AND(X,Y) _mm256_and_si256(X,Y)
#define MR_ANDNOT(X,Y) _mm256_andnot_si256(X,Y)
#define MR_OR(X,Y) _mm256_or_si256(X,Y)
#define MR_XOR(X,Y) _mm256_xor_si256(X,Y)
#define MR_SET_ALL_LANES_TO_CONSTANT(C) _mm256_set1_epi64x(C)
#define MR_SET_EACH_LANE_TO_CONSTANT(C) _mm256_set_epi64x(C[3],C[2],C[1],C[0])
#define MR_ADD64U(X,Y) _mm256_add_epi64(X,Y)
#define MR_MUL_DOUBLE_ADD(A,B,C,D) _mm256_madd52lo_epu64(_mm256_madd52lo_epu64(_mm256_setzero_si128(),A,B),C,D)
#define MR_MULADD_HI(T,X,Y) _mm256_madd52hi_epu64(T,X,Y)
#define MR_MULADD_LO(T,X,Y) _mm256_madd52lo_epu64(T,X,Y)
#define MR_MULADD_LO_64(T,X,Y) _mm256_add_epi64(T,_mm256_mullo_epi64(X,Y))
#define MR_ADDSHL(T,X,S) _mm256_add_epi64(T,_mm256_and_si256(_mm256_slli_epi64(X,S),_mm256_set1_epi64x(((int64_t)1<<52)-1)) )
#define MR_ADDSHR(T,X,S) _mm256_add_epi64(T,_mm256_srli_epi64(X,52-S))
#define MR_SHR52(LO,HI) _mm256_add_epi64(_mm256_srai_epi64(LO,52),HI)
#define MR_MULHI52(X,Y) _mm256_madd52hi_epu64(_mm256_setzero_si256(),X,Y)    
#define MR_MULLO52(X,Y) _mm256_madd52lo_epu64(_mm256_setzero_si256(),X,Y)
#define MR_MULLO(X,Y) _mm256_mullo_epi64(X,Y)
#define MR_SHL(L,H,S) _mm256_add_epi64(_mm256_srai_epi64(L,S),_mm256_slli_epi64(H,(52-S)))
#define MR_SUB64U(X,Y) _mm256_sub_epi64(X,Y)
#define MR_SHR64U(X,C) _mm256_srli_epi64(X,C)
#define MR_SHL64U(X,C) _mm256_slli_epi64(X,C)
#define MR_SRA64S(X,C) _mm256_srai_epi64(X,C)
#define MR_EXTRACT(X,L) ((unsigned char)(_mm256_extract_epi16(X,L*4)&0xFF))

#endif

#if SIMD_ENGINE == AVX512_IFMA

#define MR_ZERO() _mm512_setzero_si512()
#define MR_AND(X,Y) _mm512_and_si512(X,Y)
#define MR_ANDNOT(X,Y) _mm512_andnot_si512(X,Y)
#define MR_OR(X,Y) _mm512_or_si512(X,Y)
#define MR_XOR(X,Y) _mm512_xor_si512(X,Y)
#define MR_SET_ALL_LANES_TO_CONSTANT(C) _mm512_set1_epi64(C)
#define MR_SET_EACH_LANE_TO_CONSTANT(C) _mm512_set_epi64(C[7],C[6],C[5],C[4],C[3],C[2],C[1],C[0])
#define MR_ADD64U(X,Y) _mm512_add_epi64(X,Y)
#define MR_MUL_DOUBLE_ADD(A,B,C,D) _mm512_madd52lo_epu64(_mm512_madd52lo_epu64(_mm512_setzero_si128(),A,B),C,D)
#define MR_MULADD_HI(T,X,Y) _mm512_madd52hi_epu64(T,X,Y)
#define MR_MULADD_LO(T,X,Y) _mm512_madd52lo_epu64(T,X,Y)
#define MR_MULADD_LO_64(T,X,Y) _mm512_add_epi64(T,_mm512_mullo_epi64(X,Y))
#define MR_ADDSHL(T,X,S) _mm512_add_epi64(T,_mm512_and_si512(_mm512_slli_epi64(X,S),_mm512_set1_epi64(((int64_t)1<<52)-1)) )
#define MR_ADDSHR(T,X,S) _mm512_add_epi64(T,_mm512_srli_epi64(X,52-S))
#define MR_SHR52(LO,HI) _mm512_add_epi64(_mm512_srai_epi64(LO,52),HI)
#define MR_MULHI52(X,Y) _mm512_madd52hi_epu64(_mm512_setzero_si512(),X,Y)    
#define MR_MULLO52(X,Y) _mm512_madd52lo_epu64(_mm512_setzero_si512(),X,Y)
#define MR_MULLO(X,Y) _mm512_mullo_epi64(X,Y)
#define MR_SHL(L,H,S) _mm512_add_epi64(_mm512_srai_epi64(L,S),_mm512_slli_epi64(H,(52-S)))
#define MR_SUB64U(X,Y) _mm512_sub_epi64(X,Y)
#define MR_SHR64U(X,C) _mm512_srli_epi64(X,C)
#define MR_SHL64U(X,C) _mm512_slli_epi64(X,C)
#define MR_SRA64S(X,C) _mm512_srai_epi64(X,C)
#define MR_EXTRACT(X,L) ((unsigned char)(_mm256_extract_epi16(_mm512_extracti64x4_epi64(X,L/4),(L%4)*4)&0xFF))

#endif


#endif