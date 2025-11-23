// Define required intrinsics depending on chosen engine
// Using n 64-bit lanes in a n*64 bit wide SIMD extension, where n is 2, 4 or 8

#ifndef INTRINSICS_H
#define INTRINSICS_H

// Key
// S - Signed
// U - Unsigned
// SU - Signed to Unsigned
// Signed addition/subtraction are in fact the same operation. But parameter types will be different.

#if SIMD_ENGINE == SSE4

#define MR_ZERO() _mm_setzero_si128()
#define MR_DZERO() _mm_setzero_si128()
#define MR_AND(X,Y) _mm_and_si128(X,Y)
#define MR_ANDNOT(X,Y) _mm_andnot_si128(X,Y)
#define MR_OR(X,Y) _mm_or_si128(X,Y)
#define MR_XOR(X,Y) _mm_xor_si128(X,Y)
#define MR_SET_ALL_LANES_TO_CONSTANT(C) _mm_set1_epi64x(C)
#define MR_SET_EACH_LANE_TO_CONSTANT(C) _mm_set_epi32(0,C[1],0,C[0])
#define MR_ADD64U(X,Y) _mm_add_epi64(X,Y)
#define MR_ADD32S(X,Y) _mm_add_epi32(X,Y)
#define MR_ADD32U(X,Y) _mm_add_epi32(X,Y)
#define MR_MULADDU(T,X,Y) ({__m128i PP=_mm_mul_epu32(X,Y); _mm_add_epi64(T,PP);})
#define MR_MULADDS(T,X,Y) ({__m128i PP=_mm_mul_epi32(X,Y); _mm_add_epi64(T,PP);})
#define MR_MUL64_CONSTANT(X,C) ({__m128i S=MR_SET_ALL_LANES_TO_CONSTANT(C); __m128i PP1=_mm_mul_epu32(X,S); __m128i PP2=_mm_mul_epu32(_mm_srli_epi64(X,32),S); _mm_add_epi64(PP1,_mm_slli_epi64(PP2,32));})
#define MR_MUL32_CONSTANT(X,C) ({__m128i S=MR_SET_ALL_LANES_TO_CONSTANT(C); _mm_mul_epu32(X,S);})
#define MR_MUL32_W_CONSTANT(X,C) ({__m128i S=MR_SET_ALL_LANES_TO_CONSTANT(C); _mm_mul_epu32(X,S);})
#define MR_MULADD32_CONSTANT(T,X,C) ({__m128i S=MR_SET_ALL_LANES_TO_CONSTANT(C); _mm_add_epi64(T,_mm_mul_epu32(X,S));})
#define MR_MUL32U(X,Y) _mm_mul_epu32(X,Y)
#define MR_MUL32U_N(X,Y) _mm_mul_epu32(X,Y)
#define MR_SUB32U(X,Y) _mm_sub_epi32(X,Y)
#define MR_SUB32S(X,Y) _mm_sub_epi32(X,Y)
#define MR_SUB64U(X,Y) _mm_sub_epi64(X,Y)
#define MR_ADD6432U(X,Y) _mm_add_epi64(X,Y)
#define MR_SUB6432U(X,Y) _mm_sub_epi64(X,Y)
#define MR_SRA32S(X,C) _mm_srai_epi32(X,C)
#define MR_SHR32U(X,C) _mm_srli_epi32(X,C)
#define MR_SHR64U(X,C) _mm_srli_epi64(X,C)
#define MR_SHL64U(X,C) _mm_slli_epi64(X,C)
#define MR_SHL32U(X,C) _mm_slli_epi32(X,C)
#define MR_CAST32_SU(X) (X)
#define MR_CAST32_US(X) (X)
#define MR_CAST64_SU(X) (X)
#define MR_CAST64_US(X) (X)
#define MR_CAST6432(X) (X)
#define MR_CAST3264(X) (X)
#define MR_EXTRACT(X,L) ((unsigned char)(_mm_extract_epi16(X,L*4)&0xFF))

#endif


#if SIMD_ENGINE == AVX2

#define MR_ZERO() _mm256_setzero_si256()
#define MR_DZERO() _mm256_setzero_si256()
#define MR_AND(X,Y) _mm256_and_si256(X,Y)
#define MR_ANDNOT(X,Y) _mm256_andnot_si256(X,Y)
#define MR_OR(X,Y) _mm256_or_si256(X,Y)
#define MR_XOR(X,Y) _mm256_xor_si256(X,Y)
#define MR_SET_ALL_LANES_TO_CONSTANT(C) _mm256_set1_epi64x(C)
#define MR_SET_EACH_LANE_TO_CONSTANT(C) _mm256_set_epi32(0,C[3],0,C[2],0,C[1],0,C[0])
#define MR_ADD64U(X,Y) _mm256_add_epi64(X,Y)
#define MR_ADD32S(X,Y) _mm256_add_epi32(X,Y)
#define MR_ADD32U(X,Y) _mm256_add_epi32(X,Y)
#define MR_MULADDU(T,X,Y) ({__m256i PP=_mm256_mul_epu32(X,Y); _mm256_add_epi64(T,PP);})
#define MR_MULADDS(T,X,Y) ({__m256i PP=_mm256_mul_epi32(X,Y); _mm256_add_epi64(T,PP);})
#define MR_MUL64_CONSTANT(X,C) ({__m256i S=MR_SET_ALL_LANES_TO_CONSTANT(C); __m256i PP1=_mm256_mul_epu32(X,S); __m256i PP2=_mm256_mul_epu32(_mm256_srli_epi64(X,32),S); _mm256_add_epi64(PP1,_mm256_slli_epi64(PP2,32));})
#define MR_MUL32_CONSTANT(X,C) ({__m256i S=MR_SET_ALL_LANES_TO_CONSTANT(C); _mm256_mul_epu32(X,S);})
#define MR_MUL32_W_CONSTANT(X,C) ({__m256i S=MR_SET_ALL_LANES_TO_CONSTANT(C); _mm256_mul_epu32(X,S);})
#define MR_MULADD32_CONSTANT(T,X,C) ({__m256i S=MR_SET_ALL_LANES_TO_CONSTANT(C); _mm256_add_epi64(T,_mm256_mul_epu32(X,S));})
#define MR_MUL32U(X,Y) _mm256_mul_epu32(X,Y)
#define MR_MUL32U_N(X,Y) _mm256_mul_epu32(X,Y)
#define MR_SUB32U(X,Y) _mm256_sub_epi32(X,Y)
#define MR_SUB32S(X,Y) _mm256_sub_epi32(X,Y)
#define MR_SUB64U(X,Y) _mm256_sub_epi64(X,Y)
#define MR_ADD6432U(X,Y) _mm256_add_epi64(X,Y)
#define MR_SUB6432U(X,Y) _mm256_sub_epi64(X,Y)
#define MR_SRA32S(X,C) _mm256_srai_epi32(X,C)
#define MR_SHR32U(X,C) _mm256_srli_epi32(X,C)
#define MR_SHR64U(X,C) _mm256_srli_epi64(X,C)
#define MR_SHL64U(X,C) _mm256_slli_epi64(X,C)
#define MR_SHL32U(X,C) _mm256_slli_epi32(X,C)
#define MR_CAST32_SU(X) (X)
#define MR_CAST32_US(X) (X)
#define MR_CAST64_SU(X) (X)
#define MR_CAST64_US(X) (X)
#define MR_CAST6432(X) (X)
#define MR_CAST3264(X) (X)
#define MR_EXTRACT(X,L) ((unsigned char)(_mm256_extract_epi16(X,L*4)&0xFF))

#endif

// Here we assume the property AVX-512DQ is available. If so we get a faster MR_MUL64_CONSTANT(X,C). See below
// Note that if this property (which requires AVX-512 support) is available a similar optimization could be appled to AVX2 and SSE (see above)

#if SIMD_ENGINE == AVX512

#define MR_ZERO() _mm512_setzero_si512()
#define MR_DZERO() _mm512_setzero_si512()
#define MR_AND(X,Y) _mm512_and_si512(X,Y)
#define MR_ANDNOT(X,Y) _mm512_andnot_si512(X,Y)
#define MR_OR(X,Y) _mm512_or_si512(X,Y)
#define MR_XOR(X,Y) _mm512_xor_si512(X,Y)
#define MR_SET_ALL_LANES_TO_CONSTANT(C) _mm512_set1_epi64(C)
#define MR_SET_EACH_LANE_TO_CONSTANT(C) _mm512_set_epi32(0,C[7],0,C[6],0,C[5],0,C[4],0,C[3],0,C[2],0,C[1],0,C[0])
#define MR_ADD64U(X,Y) _mm512_add_epi64(X,Y)
#define MR_ADD32S(X,Y) _mm512_add_epi32(X,Y)
#define MR_ADD32U(X,Y) _mm512_add_epi32(X,Y)
#define MR_MULADDU(T,X,Y) ({__m512i PP=_mm512_mul_epu32(X,Y); _mm512_add_epi64(T,PP);})
#define MR_MULADDS(T,X,Y) ({__m512i PP=_mm512_mul_epi32(X,Y); _mm512_add_epi64(T,PP);})
#define MR_MUL64_CONSTANT(X,C) ({__m512i S=MR_SET_ALL_LANES_TO_CONSTANT(C); _mm512_mullo_epi64(X,S);})
//#define MR_MUL64_CONSTANT(X,C) ({__m512i S=MR_SET_ALL_LANES_TO_CONSTANT(C); __m512i PP1=_mm512_mul_epu32(X,S); __m512i PP2=_mm512_mul_epu32(_mm512_srli_epi64(X,32),S); _mm512_add_epi64(PP1,_mm512_slli_epi64(PP2,32));})
#define MR_MUL32_CONSTANT(X,C) ({__m512i S=MR_SET_ALL_LANES_TO_CONSTANT(C); _mm512_mul_epu32(X,S);})
#define MR_MUL32_W_CONSTANT(X,C) ({__m512i S=MR_SET_ALL_LANES_TO_CONSTANT(C); _mm512_mul_epu32(X,S);})
#define MR_MULADD32_CONSTANT(T,X,C) ({__m512i S=MR_SET_ALL_LANES_TO_CONSTANT(C); _mm512_add_epi64(T,_mm512_mul_epu32(X,S));})
#define MR_MUL32U(X,Y) _mm512_mul_epu32(X,Y)
#define MR_MUL32U_N(X,Y) _mm512_mul_epu32(X,Y)
#define MR_SUB32U(X,Y) _mm512_sub_epi32(X,Y)
#define MR_SUB32S(X,Y) _mm512_sub_epi32(X,Y)
#define MR_SUB64U(X,Y) _mm512_sub_epi64(X,Y)
#define MR_ADD6432U(X,Y) _mm512_add_epi64(X,Y)
#define MR_SUB6432U(X,Y) _mm512_sub_epi64(X,Y)
#define MR_SRA32S(X,C) _mm512_srai_epi32(X,C)
#define MR_SHR32U(X,C) _mm512_srli_epi32(X,C)
#define MR_SHR64U(X,C) _mm512_srli_epi64(X,C)
#define MR_SHL64U(X,C) _mm512_slli_epi64(X,C)
#define MR_SHL32U(X,C) _mm512_slli_epi32(X,C)
#define MR_CAST32_SU(X) (X)
#define MR_CAST32_US(X) (X)
#define MR_CAST64_SU(X) (X)
#define MR_CAST64_US(X) (X)
#define MR_CAST6432(X) (X)
#define MR_CAST3264(X) (X)
#define MR_EXTRACT(X,L) ((unsigned char)(_mm256_extract_epi16(_mm512_extracti64x4_epi64(X,L/4),(L%4)*4)&0xFF))

#endif



#if SIMD_ENGINE == NEON

#define MR_ZERO() vdup_n_u32(0)
#define MR_DZERO() vdupq_n_u64(0)
#define MR_AND(X,Y) vand_u32(X,Y)
#define MR_ANDNOT(X,Y) vand_u32(vmvn_u32(X),Y)
#define MR_OR(X,Y) vorr_u32(X,Y)
#define MR_XOR(X,Y) veor_u32(X,Y)
#define MR_SET_ALL_LANES_TO_CONSTANT(C) vdup_n_u32(C)
#define MR_SET_EACH_LANE_TO_CONSTANT(C) vld1_u32(C)
#define MR_ADD64U(X,Y) vaddq_u64(X,Y)
#define MR_ADD32S(X,Y) vadd_s32(X,Y)
#define MR_ADD32U(X,Y) vadd_u32(X,Y)
#define MR_MULADDU(T,X,Y) vmlal_u32(T,X,Y)
#define MR_MULADDS(T,X,Y) vmlal_s32(T,X,Y) 
#define MR_MUL64_CONSTANT(X,C) ({uint64x2_t PP1=vmull_n_u32(vmovn_u64(X),C); uint64x2_t PP2=vmull_n_u32(vmovn_u64(vshrq_n_u64(X,32)),C);vaddq_u64(PP1,vshlq_n_u64(PP2,32));})
#define MR_MUL32_CONSTANT(X,C) vmul_n_u32(X,C)
#define MR_MUL32_W_CONSTANT(X,C) vmull_n_u32(X,C)
#define MR_MULADD32_CONSTANT(T,X,C) vmlal_n_u32(T,X,C)
#define MR_MUL32U(X,Y) vmull_u32(X,Y)
#define MR_MUL32U_N(X,Y) vmul_u32(X,Y)
#define MR_SUB32U(X,Y) vsub_u32(X,Y)
#define MR_SUB32S(X,Y) vsub_s32(X,Y)
#define MR_SUB64U(X,Y) vsubq_u64(X,Y)
#define MR_ADD6432U(X,Y) vaddw_u32(X,Y)
#define MR_SUB6432U(X,Y) vsubw_u32(X,Y)
#define MR_SRA32S(X,C) vshr_n_s32(X,C)
#define MR_SHR32U(X,C) vshr_n_u32(X,C)
#define MR_SHR64U(X,C) vshrq_n_u64(X,C)
#define MR_SHL64U(X,C) vshlq_n_u64(X,C)
#define MR_SHL32U(X,C) vshl_n_u32(X,C)
#define MR_CAST32_SU(X) vreinterpret_u32_s32(X)
#define MR_CAST32_US(X) vreinterpret_s32_u32(X)
#define MR_CAST64_SU(X) vreinterpretq_u64_s64(X)
#define MR_CAST64_US(X) vreinterpretq_s64_u64(X)
#define MR_CAST6432(X) vmovn_u64(X)
#define MR_CAST3264(X) vmovl_u32(X)
#define MR_EXTRACT(X,L) ({uint32_t b; vst1_lane_u32(&b,X,L); (unsigned char)(b&0xff);})

#endif


#endif