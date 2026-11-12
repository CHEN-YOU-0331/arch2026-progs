# BFloat16 routines: hand-written RV32I.
#
# uint16_t bf16_from_f32_asm(uint32_t f32bits)   a0 -> a0
# uint16_t bf16_mul_asm(uint16_t a, uint16_t b)  a0, a1 -> a0
# uint16_t bf16_sqrt_asm(uint16_t a)             a0 -> a0
#
# RV32I and Zicsr only: no mul, div, rem, no floating point, no compressed
# instructions. Tests: tests/tbf16.c.
#
# The uint16_t arguments arrive zero-extended (the calling convention extends
# them to 32 bits), so a sign test is an unsigned compare against 0x8000.

    .text

# bf16_from_f32: Inf and NaN keep their top half; anything else rounds to
# nearest, ties to even, by adding 0x7FFF plus bit 16 before the shift.
    .globl  bf16_from_f32_asm
    .type   bf16_from_f32_asm, @function
    .align  2
bf16_from_f32_asm:
    lui     t1, 0x7f800             # exponent field
    and     t0, a0, t1
    beq     t0, t1, 1f
    srli    t0, a0, 16
    andi    t0, t0, 1
    addi    t0, t0, -1              # bit 16 - 1
    add     a0, a0, t0
    lui     t1, 0x8                 # + 0x8000: 0x7FFF + bit 16 in all
    add     a0, a0, t1
1:  srli    a0, a0, 16
    ret
    .size   bf16_from_f32_asm, .-bf16_from_f32_asm

# bf16_mul
#
#   a2, a3   exponent fields of a and b (later: effective exponents)
#   a4, a6   significands of a and b, with the leading one at bit 7
#   a5       product sign (bit 15)
#   a7       significand product
    .globl  bf16_mul_asm
    .type   bf16_mul_asm, @function
    .align  2
bf16_mul_asm:
    xor     a5, a0, a1
    lui     t6, 0x8
    and     a5, a5, t6              # sign
    srli    a2, a0, 7
    andi    a2, a2, 0xff
    srli    a3, a1, 7
    andi    a3, a3, 0xff
    li      t0, 0xff
    beq     a2, t0, .Lmul_a_special
    beq     a3, t0, .Lmul_b_special
    slli    t1, a0, 17              # a is a zero
    beqz    t1, .Lmul_zero
    slli    t1, a1, 17              # b is a zero
    beqz    t1, .Lmul_zero
    andi    a4, a0, 0x7f
    andi    a6, a1, 0x7f
    beqz    a2, .Lmul_a_denorm
    ori     a4, a4, 0x80
.Lmul_a_ready:
    beqz    a3, .Lmul_b_denorm
    ori     a6, a6, 0x80
.Lmul_b_ready:
    # 8 x 8 product by shift and add. Bit 7 of b is always set, so the
    # product starts at a << 7 and bits 0..6 of b are tested one by one.
    slli    a7, a4, 7
    andi    t0, a6, 0x01
    beqz    t0, 1f
    add     a7, a7, a4
1:  andi    t0, a6, 0x02
    beqz    t0, 1f
    slli    t1, a4, 1
    add     a7, a7, t1
1:  andi    t0, a6, 0x04
    beqz    t0, 1f
    slli    t1, a4, 2
    add     a7, a7, t1
1:  andi    t0, a6, 0x08
    beqz    t0, 1f
    slli    t1, a4, 3
    add     a7, a7, t1
1:  andi    t0, a6, 0x10
    beqz    t0, 1f
    slli    t1, a4, 4
    add     a7, a7, t1
1:  andi    t0, a6, 0x20
    beqz    t0, 1f
    slli    t1, a4, 5
    add     a7, a7, t1
1:  andi    t0, a6, 0x40
    beqz    t0, 1f
    slli    t1, a4, 6
    add     a7, a7, t1
1:  # exp = exp_a + exp_b - 127, plus one when the product reaches bit 15;
    # the fraction is bits 7..13 or 8..14 of the product, truncated.
    add     a2, a2, a3
    srli    t0, a7, 15
    add     a2, a2, t0
    addi    a2, a2, -127
    addi    t1, t0, 7
    srl     a7, a7, t1
    andi    a7, a7, 0x7f
    addi    t1, a2, -1              # exponent outside [1, 254]?
    li      t2, 254
    bgeu    t1, t2, .Lmul_range
    slli    a2, a2, 7
    or      a0, a5, a2
    or      a0, a0, a7
    ret

.Lmul_range:
    blez    a2, .Lmul_small
    li      t0, 0x7f80              # overflow: Inf with the product sign
    or      a0, a5, t0
    ret
.Lmul_small:
    li      t0, -6
    blt     a2, t0, .Lmul_zero
    li      t0, 1                   # denormal: fraction >> (1 - exp)
    sub     t0, t0, a2
    srl     a7, a7, t0
    or      a0, a5, a7
    ret

.Lmul_zero:
    mv      a0, a5
    ret

.Lmul_a_denorm:                     # leading one up to bit 7, exp 1 - shift
    li      a2, 1
1:  slli    a4, a4, 1
    addi    a2, a2, -1
    andi    t0, a4, 0x80
    beqz    t0, 1b
    j       .Lmul_a_ready

.Lmul_b_denorm:
    li      a3, 1
1:  slli    a6, a6, 1
    addi    a3, a3, -1
    andi    t0, a6, 0x80
    beqz    t0, 1b
    j       .Lmul_b_ready

.Lmul_a_special:                    # a is Inf or NaN
    andi    t0, a0, 0x7f
    bnez    t0, .Lmul_ret           # a NaN a is returned unchanged
    slli    t1, a1, 17
    beqz    t1, .Lmul_nan           # Inf times zero
.Lmul_inf:
    li      t0, 0x7f80
    or      a0, a5, t0
.Lmul_ret:
    ret

.Lmul_b_special:                    # b is Inf or NaN, a is not
    andi    t0, a1, 0x7f
    bnez    t0, .Lmul_ret_b
    slli    t1, a0, 17
    bnez    t1, .Lmul_inf
.Lmul_nan:
    li      a0, 0x7fc0
    ret
.Lmul_ret_b:
    mv      a0, a1
    ret
    .size   bf16_mul_asm, .-bf16_mul_asm

# bf16_sqrt
#
# For a positive normal a with exponent field E and significand m in
# [128, 256), the contract asks for the largest r with floor(r * r / 128) <= m'
# where m' is m, doubled when E is even (an odd unbiased exponent). That is
# r * r <= 128 m' + 127, i.e. r = isqrt(128 m' + 127), which lies in
# [128, 255] for every m' in [128, 512): the search bounds never bind and no
# renormalization is needed. The square root is taken digit by digit (two
# bits of the radicand per result bit), so no multiply is needed. The result
# exponent is floor((E + 127) / 2), always in [64, 190].
    .globl  bf16_sqrt_asm
    .type   bf16_sqrt_asm, @function
    .align  2
bf16_sqrt_asm:
    srli    a1, a0, 7
    andi    a1, a1, 0xff            # E
    li      t0, 0xff
    beq     a1, t0, .Lsqrt_special
    slli    t1, a0, 17
    beqz    t1, .Lsqrt_zero         # +0 and -0 give +0
    lui     t2, 0x8
    bgeu    a0, t2, .Lsqrt_nan      # negative nonzero
    beqz    a1, .Lsqrt_zero         # denormal
    andi    a2, a0, 0x7f
    ori     a2, a2, 0x80            # m
    andi    t0, a1, 1
    xori    t0, t0, 1
    sll     a2, a2, t0              # m' (E even: odd unbiased exponent)
    slli    a2, a2, 7
    ori     a2, a2, 127             # x = 128 m' + 127, in [2^14, 2^16)

    # Digit-by-digit square root. The top digit is always 1 because
    # x >= 2^14, so it starts at: x -= 2^14, res = 2^14, bit = 2^12.
    lui     t0, 0x4                 # 0x4000
    sub     a2, a2, t0
    mv      a3, t0                  # res
    lui     t1, 0x1                 # bit = 0x1000
    add     t2, a3, t1
    srli    a3, a3, 1
    bltu    a2, t2, 1f
    sub     a2, a2, t2
    add     a3, a3, t1
1:  addi    t2, a3, 0x400
    srli    a3, a3, 1
    bltu    a2, t2, 1f
    sub     a2, a2, t2
    addi    a3, a3, 0x400
1:  addi    t2, a3, 0x100
    srli    a3, a3, 1
    bltu    a2, t2, 1f
    sub     a2, a2, t2
    addi    a3, a3, 0x100
1:  addi    t2, a3, 0x40
    srli    a3, a3, 1
    bltu    a2, t2, 1f
    sub     a2, a2, t2
    addi    a3, a3, 0x40
1:  addi    t2, a3, 0x10
    srli    a3, a3, 1
    bltu    a2, t2, 1f
    sub     a2, a2, t2
    addi    a3, a3, 0x10
1:  addi    t2, a3, 0x4
    srli    a3, a3, 1
    bltu    a2, t2, 1f
    sub     a2, a2, t2
    addi    a3, a3, 0x4
1:  addi    t2, a3, 0x1
    srli    a3, a3, 1
    bltu    a2, t2, 1f
    addi    a3, a3, 0x1
1:  # r = a3 in [128, 255]: the fraction is r - 128, so the exponent goes in
    # one lower and the leading one of r carries it back up.
    addi    a1, a1, 125
    srli    a1, a1, 1
    slli    a1, a1, 7
    add     a0, a1, a3
    ret

.Lsqrt_special:                     # Inf or NaN
    andi    t0, a0, 0x7f
    bnez    t0, .Lsqrt_ret          # NaN unchanged
    lui     t2, 0x8
    bltu    a0, t2, .Lsqrt_ret      # +Inf unchanged
.Lsqrt_nan:
    li      a0, 0x7fc0
.Lsqrt_ret:
    ret
.Lsqrt_zero:
    li      a0, 0
    ret
    .size   bf16_sqrt_asm, .-bf16_sqrt_asm
