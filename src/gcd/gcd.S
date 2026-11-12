# gcd_u32_asm: binary GCD, hand-written RV32I.
#
# uint32_t gcd_u32_asm(uint32_t a, uint32_t b)   a in a0, b in a1, result in a0
#
# RV32I and Zicsr only: no mul, div, rem, no floating point, no compressed
# instructions. Tests: tests/tgcd.c.
#
# Stein's algorithm. After the common power of two is removed and both
# operands are odd, each pass subtracts the smaller from the larger and strips
# the zeros of the difference, which is even, so it shifts at least once.
#
# Instruction count is the cost (rv32emu has no timing model), so the code is
# laid out to avoid jumps rather than branches:
#   * the two subtract-and-strip paths each end with the comparison of the next
#     pass, and stay on their own side with a backward branch;
#   * a strip loop costs srli, andi, beqz per zero bit.
#
#   a0, a1   the operands (a, b)
#   a2       the common power of two
#   t0       scratch

    .text

    .globl  gcd_u32_asm
    .type   gcd_u32_asm, @function
    .align  2
gcd_u32_asm:
    beqz    a0, .Lret_b             # gcd(0, b) = b
    beqz    a1, .Lret               # gcd(a, 0) = a

    # Common power of two: shift both while both are even.
    or      t0, a0, a1
    andi    t0, t0, 1
    li      a2, 0
    bnez    t0, .Lodd_a
.Lcommon:
    srli    a0, a0, 1
    srli    a1, a1, 1
    addi    a2, a2, 1
    or      t0, a0, a1
    andi    t0, t0, 1
    beqz    t0, .Lcommon

.Lodd_a:                            # strip a until it is odd
    andi    t0, a0, 1
    bnez    t0, .Lodd_b
1:  srli    a0, a0, 1
    andi    t0, a0, 1
    beqz    t0, 1b
.Lodd_b:                            # strip b until it is odd
    andi    t0, a1, 1
    bnez    t0, .Lcmp
1:  srli    a1, a1, 1
    andi    t0, a1, 1
    beqz    t0, 1b

.Lcmp:                              # both odd
    beq     a0, a1, .Ldone
    bltu    a1, a0, .Lagt
.Lbgt:                              # b > a: b -= a, strip b
    sub     a1, a1, a0
1:  srli    a1, a1, 1
    andi    t0, a1, 1
    beqz    t0, 1b
    beq     a0, a1, .Ldone
    bltu    a0, a1, .Lbgt
.Lagt:                              # a > b: a -= b, strip a
    sub     a0, a0, a1
1:  srli    a0, a0, 1
    andi    t0, a0, 1
    beqz    t0, 1b
    beq     a0, a1, .Ldone
    bltu    a1, a0, .Lagt
    j       .Lbgt

.Ldone:
    sll     a0, a0, a2
.Lret:
    ret
.Lret_b:
    mv      a0, a1
    ret
    .size   gcd_u32_asm, .-gcd_u32_asm
