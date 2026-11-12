# op1, op2, op3 and q_mul: hand-written RV32I.
#
# uint32_t op1_asm(uint32_t a, uint32_t b, int *n)   a0, a1, a2 -> a0
# uint32_t op2_asm(uint32_t a, uint32_t b)           a0, a1 -> a0
# uint32_t op3_asm(uint32_t a, uint32_t b)           a0, a1 -> a0
# int32_t  q_mul_asm(int32_t a, int32_t b)           a0, a1 -> a0
#
# RV32I and Zicsr only: no mul, div, rem, no floating point, no compressed
# instructions. q_mul forms a 64-bit product in a register pair. Tests:
# tests/tfixed.c.

    .text

# op1: add with logic only. Each pass: carry = (a & b) << 1, a ^= b, b = carry.
# The loop is tested at the bottom, so a pass costs and, xor, slli, addi, bnez.
    .globl  op1_asm
    .type   op1_asm, @function
    .align  2
op1_asm:
    li      t1, 0                   # passes
    beqz    a1, 2f
1:  and     t0, a0, a1
    xor     a0, a0, a1
    slli    a1, t0, 1
    addi    t1, t1, 1
    bnez    a1, 1b
2:  sw      t1, 0(a2)
    ret
    .size   op1_asm, .-op1_asm

# op2: a - b modulo 2^32, which is what a + (~b + 1) computes.
    .globl  op2_asm
    .type   op2_asm, @function
    .align  2
op2_asm:
    sub     a0, a0, a1
    ret
    .size   op2_asm, .-op2_asm

# op3: floor((a + b) / 2) = (a & b) + ((a ^ b) >> 1), no intermediate overflow.
    .globl  op3_asm
    .type   op3_asm, @function
    .align  2
op3_asm:
    and     t0, a0, a1
    xor     a0, a0, a1
    srli    a0, a0, 1
    add     a0, a0, t0
    ret
    .size   op3_asm, .-op3_asm

# q_mul: Q11.20 product, (a * b + 2^19) >> 20 arithmetically, low 32 bits.
#
# The unsigned product of the magnitudes is built by shift and add over the
# bits of the smaller magnitude, which ends the loop as early as possible; the
# larger one is the 64-bit multiplicand. The product is then negated when the
# signs differ, half a unit is added, and bits 20..51 are returned.
#
#   a0       multiplicand, low word (the larger magnitude)
#   a4       multiplicand, high word
#   a1       multiplier (the smaller magnitude), consumed from bit 0
#   a2, a3   product, low and high word
#   a5       a ^ b: its sign is the sign of the product
    .globl  q_mul_asm
    .type   q_mul_asm, @function
    .align  2
q_mul_asm:
    xor     a5, a0, a1
    srai    t0, a0, 31              # |a|
    xor     a0, a0, t0
    sub     a0, a0, t0
    srai    t0, a1, 31              # |b|
    xor     a1, a1, t0
    sub     a1, a1, t0
    bgeu    a0, a1, 1f              # multiplier = the smaller magnitude
    mv      t0, a0
    mv      a0, a1
    mv      a1, t0
1:  li      a2, 0
    li      a3, 0
    li      a4, 0
    beqz    a1, 4f
2:  andi    t0, a1, 1
    beqz    t0, 3f
    add     a2, a2, a0              # product += multiplicand
    sltu    t0, a2, a0
    add     a3, a3, a4
    add     a3, a3, t0
3:  srli    t0, a0, 31              # multiplicand <<= 1
    slli    a0, a0, 1
    slli    a4, a4, 1
    or      a4, a4, t0
    srli    a1, a1, 1
    bnez    a1, 2b
4:  bgez    a5, 5f                  # negate when the signs differ
    seqz    t0, a2
    not     a3, a3
    add     a3, a3, t0
    neg     a2, a2
5:  lui     t0, 0x80                # + 2^19, half a unit
    add     a2, a2, t0
    sltu    t0, a2, t0
    add     a3, a3, t0
    srli    a2, a2, 20              # bits 20..51
    slli    a3, a3, 12
    or      a0, a2, a3
    ret
    .size   q_mul_asm, .-q_mul_asm
