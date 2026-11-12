#include "fixed.h"

/* a + b with logic only: the sum without carries is a ^ b, the carries are
 * (a & b) << 1, and they are added in again until none is left.
 */
uint32_t op1(uint32_t a, uint32_t b, int *n)
{
    int passes = 0;
    while (b) {
        uint32_t carry = (a & b) << 1;
        a ^= b;
        b = carry;
        passes++;
    }
    *n = passes;
    return a;
}

/* a - b modulo 2^32, which is what a + (~b + 1) computes. */
uint32_t op2(uint32_t a, uint32_t b)
{
    return a - b;
}

/* floor((a + b) / 2) without overflow: the bits both operands share count in
 * full, the bits only one of them has count half.
 */
uint32_t op3(uint32_t a, uint32_t b)
{
    return (a & b) + ((a ^ b) >> 1);
}

/* Q11.20 multiply. The 64-bit product is built from 16-bit halves, so each
 * partial product fits in 32 bits and needs only __mulsi3, never __muldi3.
 * The unsigned product of the magnitudes is negated when the signs differ;
 * then half a unit is added and the sum is shifted right by 20, which is
 * the arithmetic shift of the signed product. Only the low 32 bits are kept.
 */
int32_t q_mul(int32_t a, int32_t b)
{
    uint32_t ua = a < 0 ? 0u - (uint32_t) a : (uint32_t) a;
    uint32_t ub = b < 0 ? 0u - (uint32_t) b : (uint32_t) b;
    uint32_t al = ua & 0xffffu, ah = ua >> 16;
    uint32_t bl = ub & 0xffffu, bh = ub >> 16;

    uint32_t ll = al * bl, lh = al * bh, hl = ah * bl, hh = ah * bh;

    /* lo + 2^32 hi = ll + 2^16 (lh + hl) + 2^32 hh */
    uint32_t mid = lh + hl;
    uint32_t mid_carry = mid < lh; /* worth 2^48 */
    uint32_t lo = ll + (mid << 16);
    uint32_t hi = hh + (mid >> 16) + (mid_carry << 16) + (lo < ll);

    if ((a < 0) != (b < 0)) {
        hi = ~hi + (lo == 0);
        lo = 0u - lo;
    }

    uint32_t half = 1u << (Q - 1);
    uint32_t sum = lo + half;
    hi += sum < lo;
    return (int32_t) ((sum >> Q) | (hi << (32 - Q)));
}
