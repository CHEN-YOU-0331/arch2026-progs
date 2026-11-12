#include "bf16.h"

#define BF16_SIGN 0x8000u
#define BF16_INF 0x7F80u
#define BF16_NAN 0x7FC0u
#define BF16_BIAS 127

uint16_t bf16_from_f32(uint32_t f32bits)
{
    /* Inf and NaN keep their top 16 bits as they are. */
    if (((f32bits >> 23) & 0xFFu) == 0xFFu)
        return (uint16_t) (f32bits >> 16);

    /* Round to nearest, ties to even: 0x7FFF rounds everything above the
     * halfway point up, and bit 16 adds the last unit for a tie on odd.
     */
    uint32_t lsb = (f32bits >> 16) & 1u;
    return (uint16_t) ((f32bits + 0x7FFFu + lsb) >> 16);
}

uint16_t bf16_mul(uint16_t a, uint16_t b)
{
    uint32_t sign = (uint32_t) (a ^ b) & BF16_SIGN;
    int32_t exp_a = (a >> 7) & 0xFF, exp_b = (b >> 7) & 0xFF;
    uint32_t mant_a = a & 0x7Fu, mant_b = b & 0x7Fu;

    /* a first, then b; see the contract for the order. */
    if (exp_a == 0xFF) {
        if (mant_a)
            return a;
        if (!exp_b && !mant_b)
            return BF16_NAN;
        return (uint16_t) (sign | BF16_INF);
    }
    if (exp_b == 0xFF) {
        if (mant_b)
            return b;
        if (!exp_a && !mant_a)
            return BF16_NAN;
        return (uint16_t) (sign | BF16_INF);
    }
    if ((!exp_a && !mant_a) || (!exp_b && !mant_b))
        return (uint16_t) sign;

    /* A denormal is normalized by shifting its leading one up to bit 7; it
     * then counts as exponent 1 minus the shift.
     */
    if (!exp_a) {
        exp_a = 1;
        while (!(mant_a & 0x80u)) {
            mant_a <<= 1;
            exp_a--;
        }
    } else {
        mant_a |= 0x80u;
    }
    if (!exp_b) {
        exp_b = 1;
        while (!(mant_b & 0x80u)) {
            mant_b <<= 1;
            exp_b--;
        }
    } else {
        mant_b |= 0x80u;
    }

    /* Both significands lie in [128, 256), so the product lies in
     * [2^14, 2^16): one renormalization step at most, then truncation.
     */
    uint32_t prod = mant_a * mant_b;
    int32_t exp = exp_a + exp_b - BF16_BIAS;
    uint32_t frac;
    if (prod & 0x8000u) {
        frac = (prod >> 8) & 0x7Fu;
        exp++;
    } else {
        frac = (prod >> 7) & 0x7Fu;
    }

    if (exp >= 0xFF)
        return (uint16_t) (sign | BF16_INF);
    if (exp <= 0) {
        if (exp < -6)
            return (uint16_t) sign;
        frac >>= 1 - exp;
        exp = 0;
    }
    return (uint16_t) (sign | ((uint32_t) exp << 7) | frac);
}

uint16_t bf16_sqrt(uint16_t a)
{
    uint32_t exp = (a >> 7) & 0xFFu, mant = a & 0x7Fu;

    if (exp == 0xFF) {
        if (mant)
            return a; /* NaN */
        if (a & BF16_SIGN)
            return BF16_NAN; /* -Inf */
        return a;            /* +Inf */
    }
    if (!exp && !mant)
        return 0; /* either zero gives +0 */
    if (a & BF16_SIGN)
        return BF16_NAN;
    if (!exp)
        return 0; /* denormal */

    /* sqrt(2^e * m) = 2^(e/2) * sqrt(m). With e = exp - 127 odd the
     * significand is doubled so the exponent can be halved exactly; in both
     * cases the new biased exponent is floor((exp + 127) / 2).
     */
    int32_t e = (int32_t) exp - BF16_BIAS;
    int32_t new_exp = (int32_t) ((exp + BF16_BIAS) >> 1);
    uint32_t m = 0x80u | mant;
    if (e & 1)
        m <<= 1;

    /* The largest mid in [90, 256] with mid * mid / 128 <= m. */
    uint32_t low = 90, high = 256, result = 128;
    while (low <= high) {
        uint32_t mid = (low + high) >> 1;
        if ((mid * mid) / 128 <= m) {
            result = mid;
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }

    if (result >= 256) {
        result >>= 1;
        new_exp++;
    } else {
        while (result < 128 && new_exp > 1) {
            result <<= 1;
            new_exp--;
        }
    }
    if (new_exp >= 0xFF)
        return BF16_INF;
    if (new_exp <= 0)
        return 0;
    return (uint16_t) (((uint32_t) new_exp << 7) | (result & 0x7Fu));
}
