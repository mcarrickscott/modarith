// Some useful hash functions

#include <stdint.h>

#ifndef HASH_H
#define HASH_H

/* SHA2 */

/**
 * @brief SHA256 hash function instance */
typedef struct
{
    uint32_t length[2]; /**< 64-bit input length */
    uint32_t h[8];      /**< Internal state */
    uint32_t w[64];	/**< Internal state */
    int hlen;		/**< Hash length in bytes */
} hash256;


/**
 * @brief SHA384-512 hash function instance */
typedef struct
{
    uint64_t length[2]; /**< 64-bit input length */
    uint64_t h[8];      /**< Internal state */
    uint64_t w[80];	/**< Internal state */
    int hlen;           /**< Hash length in bytes */
} hash512;

/**
 * @brief SHA384 hash function instance */
typedef hash512 hash384;



/* SHA3 */

/**
 * @brief SHA3 hash function instance */
typedef struct
{
    int length;   /**< 64-bit input length */
    uint64_t S[25];  /**< Internal state */
    int rate;          /**< TODO */
    int len;           /**< Hash length in bytes */
} sha3;

#define SHA3_HASH224 28 /**< SHA3 224 bit hash */
#define SHA3_HASH256 32 /**< SHA3 256 bit hash */
#define SHA3_HASH384 48 /**< SHA3 384 bit hash */
#define SHA3_HASH512 64 /**< SHA3 512 bit hash */

#define SHAKE128 16 /**< SHAKE128   hash */
#define SHAKE256 32 /**< SHAKE256 hash */

/* Hash function */
/**	@brief Initialise an instance of SHA256
 *
	@param H an instance SHA256
 */
extern void HASH256_init(hash256 *H);
/**	@brief Add a byte to the hash
 *
	@param H an instance SHA256
	@param b byte to be included in hash
 */
extern void HASH256_process(hash256 *H, int b);
/**	@brief Generate 32-byte final hash
 *
	@param H an instance SHA256
	@param h is the output 32-byte hash
 */
extern void HASH256_hash(hash256 *H, char *h);

/**	@brief Generate 32-byte intermediate hash
 *
	@param H an instance SHA256
	@param h is the output 32-byte hash
 */
extern void HASH256_continuing_hash(hash256 *H, char *h);


/**	@brief Initialise an instance of SHA384
 *
	@param H an instance SHA384
 */
extern void HASH384_init(hash384 *H);
/**	@brief Add a byte to the hash
 *
	@param H an instance SHA384
	@param b byte to be included in hash
 */
extern void HASH384_process(hash384 *H, int b);
/**	@brief Generate 48-byte final hash
 *
	@param H an instance SHA384
	@param h is the output 48-byte hash
 */
extern void HASH384_hash(hash384 *H, char *h);

/**	@brief Generate 48-byte intermediate hash
 *
	@param H an instance SHA384
	@param h is the output 48-byte hash
 */
extern void HASH384_continuing_hash(hash384 *H, char *h);

/**	@brief Initialise an instance of SHA512
 *
	@param H an instance SHA512
 */
extern void HASH512_init(hash512 *H);
/**	@brief Add a byte to the hash
 *
	@param H an instance SHA512
	@param b byte to be included in hash
 */
extern void HASH512_process(hash512 *H, int b);
/**	@brief Generate 64-byte final hash
 *
	@param H an instance SHA512
	@param h is the output 64-byte hash
 */
extern void HASH512_hash(hash512 *H, char *h);

/**	@brief Generate 64-byte intermediate hash
 *
	@param H an instance SHA512
	@param h is the output 64-byte hash
 */
extern void HASH512_continuing_hash(hash512 *H, char *h);


/**	@brief Initialise an instance of SHA3
 *
	@param H an instance SHA3
	@param t the instance type
 */
extern void  SHA3_init(sha3 *H, int t);
/**	@brief process a byte for SHA3
 *
	@param H an instance SHA3
	@param b a byte of date to be processed
 */
extern void  SHA3_process(sha3 *H, int b);
/**	@brief create fixed length final hash output of SHA3
 *
	@param H an instance SHA3
	@param h a byte array to take hash
 */
extern void  SHA3_hash(sha3 *H, char *h);

/**	@brief create fixed length intermediate hash output of SHA3
 *
	@param H an instance SHA3
	@param h a byte array to take hash
 */
extern void  SHA3_continuing_hash(sha3 *H, char *h);

/**	@brief create variable length final hash output of SHA3
 *
	@param H an instance SHA3
	@param h a byte array to take hash
	@param len is the length of the hash
 */
extern void  SHA3_shake(sha3 *H, char *h, int len);

/**	@brief create variable length intermediate hash output of SHA3
 *
	@param H an instance SHA3
	@param h a byte array to take hash
	@param len is the length of the hash
 */
extern void  SHA3_continuing_shake(sha3 *H, char *h, int len);

/**	@brief generate further hash output of SHA3
 *
	@param H an instance SHA3
	@param h a byte array to take hash
	@param len is the length of the hash
 */
extern void  SHA3_squeeze(sha3 *H, char *h, int len);

#endif