#include "gcd.h"

/* Binary GCD (Stein's algorithm). gcd(2a, 2b) = 2 gcd(a, b), gcd(2a, b) =
 * gcd(a, b) for odd b, and gcd(a, b) = gcd(a, b - a) for odd a < b, where
 * b - a is even. The loop only shifts and subtracts, and each subtraction is
 * followed by at least one shift, so the larger operand loses a bit per pass.
 */
uint32_t gcd_u32(uint32_t a, uint32_t b)
{
    if (a == 0)
        return b;
    if (b == 0)
        return a;

    /* The common factor of two. */
    unsigned shift = 0;
    while (((a | b) & 1u) == 0) {
        a >>= 1;
        b >>= 1;
        shift++;
    }

    while ((a & 1u) == 0)
        a >>= 1;

    /* a is odd from here on. */
    for (;;) {
        while ((b & 1u) == 0)
            b >>= 1;
        if (a == b)
            break;
        if (a > b) {
            uint32_t t = a;
            a = b;
            b = t;
        }
        b -= a;
    }
    return a << shift;
}
