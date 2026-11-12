# Optimal 2x2x2 solver: hand-written RV32I.
#
# int rubik_solve_asm(const char *state, char *out, size_t cap)
#                                                  a0, a1, a2 -> a0
#
# RV32I and Zicsr only: no mul, div, rem, no floating point, no compressed
# instructions. Tests: tests/trubik.c.
#
# The algorithm is the one of rubik.c: iterative deepening for limit = 6..10
# over move sequences of length limit - 5, each sequence ending in a lookup in
# the table of the states within 5 moves, which also supplies the last moves;
# a state at distance 11 is marked in that table and solved through its first
# unmarked neighbor with the limit at 10. The tables are rubik_tab.S, the same
# data as rubik_tab.c (layout in rubik_tab.h).
#
# A state is carried as two row addresses: P = rubik_asm_perm + 8 * p and
# O = rubik_asm_orient + 8 * o. A row is four halfwords; halfword f < 3 holds
# (next rank << 3) | lower bound, so masking its low three bits gives 8 * next
# rank, the byte offset of the next row, and the low three bits are the
# bound. Halfword 3 of a permutation row is the byte offset of its group in
# rubik_asm_near; the next row's halfword 3 ends the group. A table entry is
# (o << 4) | code, so (O - rubik_asm_orient) << 1 is the key it is compared
# with.
#
# Registers held through the whole call (callee-saved):
#   s0   rubik_asm_perm
#   s1   rubik_asm_orient
#   s2   rubik_asm_near
#   s3   move counter of the distance-11 loop
#   s5   P of the state the deepening starts from
#   s6   O of the state the deepening starts from
#   s7   limit of the deepening
#   s8   where the next move is written
#   s9   moves written so far
#   s10  out
#   s11  cap
# The move buffer is 12 bytes at 0(sp).

    .text

# A child in .Lsearch: a turn of face f, the (t + 1)-th quarter turn from the
# node. s9, s10 hold P, O of the previous turn (the node itself for t = 0),
# s6 the moves left for the child. A child whose bound exceeds them is
# skipped; otherwise its move is written and the search goes one level down.
.macro SCHILD f, t
    lhu     t1, 2*\f(s9)
    lhu     t2, 2*\f(s10)
    andi    t3, t1, -8
    add     s9, s0, t3              # P of the child
    andi    t3, t2, -8
    add     s10, s1, t3             # O of the child
    andi    t1, t1, 7               # dp
    bltu    s6, t1, .Lsc\@
    andi    t2, t2, 7               # dor
    bltu    s6, t2, .Lsc\@
    li      t0, 3*\f+\t
    sb      t0, 0(s8)
    mv      a0, s9
    mv      a1, s10
    mv      a2, s6
    li      a3, \f
    addi    a4, s8, 1
    jal     .Lsearch
    bnez    a0, .Lsearch_ret
.Lsc\@:
.endm

.macro SFACE f
    li      t0, \f
    beq     s7, t0, .Lsf\@          # never the face turned last
    mv      s9, s4
    mv      s10, s5
    SCHILD  \f, 0
    SCHILD  \f, 1
    SCHILD  \f, 2
.Lsf\@:
.endm

# A child in .Lfrontier: as SCHILD, but the child is at the end of the
# sequence, so it is looked up in the table instead of searched. a5, a6 hold
# P, O of the previous turn; a2 is 5, the moves left for the child. The group
# of the child's permutation is scanned in increasing order up to the first
# entry not below the key; it is a hit when that entry has the key and a code
# below 10 (code 15 marks distance 11, which is no hit).
.macro FCHILD f, t
    lhu     t1, 2*\f(a5)
    lhu     t2, 2*\f(a6)
    andi    t3, t1, -8
    andi    t4, t2, -8
    add     a5, s0, t3              # P of the child
    add     a6, s1, t4              # O of the child
    andi    t1, t1, 7               # dp
    bltu    a2, t1, .Lfc\@
    andi    t2, t2, 7               # dor
    bltu    a2, t2, .Lfc\@
    lhu     t5, 6(a5)               # group of the permutation
    lhu     t6, 14(a5)
    beq     t5, t6, .Lfc\@
    add     t5, t5, s2
    add     t6, t6, s2
    slli    t4, t4, 1               # key: o << 4
.Lfl\@:
    lhu     t0, 0(t5)
    bgeu    t0, t4, .Lfk\@
    addi    t5, t5, 2
    bne     t5, t6, .Lfl\@
    j       .Lfc\@
.Lfk\@:
    sub     t0, t0, t4              # the code if the keys match
    sltiu   t1, t0, 10
    beqz    t1, .Lfc\@
    li      t1, 3*\f+\t
    j       .Lfrontier_hit
.Lfc\@:
.endm

.macro FFACE f
    li      t0, \f
    beq     a3, t0, .Lff\@
    mv      a5, a0
    mv      a6, a1
    FCHILD  \f, 0
    FCHILD  \f, 1
    FCHILD  \f, 2
.Lff\@:
.endm

# Inner products of the Lehmer code: t4 counts the cubies after position i
# that are smaller than the one at i (in t1). Recursion over j = i+1..6.
.macro LESS j
    .if \j < 7
    lbu     t2, \j(a0)
    sltu    t3, t2, t1
    add     t4, t4, t3
    LESS    "(\j+1)"
    .endif
.endm

.macro SMALLER i
    lbu     t1, \i(a0)
    li      t4, 0
    LESS    "(\i+1)"
.endm

    .globl  rubik_solve_asm
    .type   rubik_solve_asm, @function
    .align  2
rubik_solve_asm:
    addi    sp, sp, -64
    sw      ra, 60(sp)
    sw      s0, 56(sp)
    sw      s1, 52(sp)
    sw      s2, 48(sp)
    sw      s3, 44(sp)
    sw      s4, 40(sp)
    sw      s5, 36(sp)
    sw      s6, 32(sp)
    sw      s7, 28(sp)
    sw      s8, 24(sp)
    sw      s9, 20(sp)
    sw      s10, 16(sp)
    sw      s11, 12(sp)
    mv      s10, a1
    mv      s11, a2

    # Validate. Each character is checked before the next is read, so a
    # short string stops at its NUL (which is out of range).
    mv      t0, a0
    addi    t4, a0, 7
    li      t5, 7
    li      t6, 0                   # cubies seen, one bit each
1:  lbu     t1, 0(t0)
    addi    t1, t1, -'1'
    bgeu    t1, t5, .Linvalid
    li      t2, 1
    sll     t2, t2, t1
    and     t3, t6, t2
    bnez    t3, .Linvalid           # repeated cubie
    or      t6, t6, t2
    addi    t0, t0, 1
    bne     t0, t4, 1b

    addi    t4, a0, 13
    li      t5, 3
    li      a3, 0                   # orientation rank
    li      a4, 0                   # orientation sum
2:  lbu     t1, 0(t0)
    addi    t1, t1, -'1'
    bgeu    t1, t5, .Linvalid
    add     a4, a4, t1
    slli    t2, a3, 1               # rank = 3 * rank + digit
    add     a3, a3, t2
    add     a3, a3, t1
    addi    t0, t0, 1
    bne     t0, t4, 2b
    lbu     t1, 0(t0)               # the seventh orientation
    addi    t1, t1, -'1'
    bgeu    t1, t5, .Linvalid
    add     a4, a4, t1
    lbu     t1, 1(t0)
    bnez    t1, .Linvalid           # longer than 14 characters
    li      t2, 0x1249              # bits 0, 3, 6, 9, 12: sums divisible by 3
    srl     t2, t2, a4
    andi    t2, t2, 1
    beqz    t2, .Linvalid

    # Lehmer code p = ((((c0 * 6 + c1) * 5 + c2) * 4 + c3) * 3 + c4) * 2 + c5,
    # where ci counts the later cubies smaller than cubie i.
    SMALLER 0
    mv      a5, t4
    SMALLER 1
    slli    t2, a5, 1               # * 6
    add     a5, a5, t2
    slli    a5, a5, 1
    add     a5, a5, t4
    SMALLER 2
    slli    t2, a5, 2               # * 5
    add     a5, a5, t2
    add     a5, a5, t4
    SMALLER 3
    slli    a5, a5, 2               # * 4
    add     a5, a5, t4
    SMALLER 4
    slli    t2, a5, 1               # * 3
    add     a5, a5, t2
    add     a5, a5, t4
    SMALLER 5
    slli    a5, a5, 1               # * 2
    add     a5, a5, t4

    lla     s0, rubik_asm_perm
    lla     s1, rubik_asm_orient
    lla     s2, rubik_asm_near
    slli    a5, a5, 3
    add     s5, s0, a5              # P of the cube
    slli    a3, a3, 3
    add     s6, s1, a3              # O of the cube
    mv      s8, sp
    li      s9, 0

    mv      a0, s5
    mv      a1, s6
    jal     .Lnear
    li      t0, 15
    beq     a2, t0, .Lfar
    bgez    a2, .Lwalk              # within 5 moves: the table has it all
    li      s7, 6

    # Iterative deepening from (s5, s6), the limit in s7; the moves go to s8.
.Ldeepen:
    mv      a0, s5
    mv      a1, s6
    mv      a2, s7
    li      a3, 3                   # no face turned yet
    mv      a4, s8
    jal     .Lsearch
    bnez    a0, .Lfound
    addi    s7, s7, 1
    li      t0, 10
    bgeu    t0, s7, .Ldeepen
    j       .Linvalid               # not reached: every state is within 11

.Lfound:                            # limit - 5 moves written, then the table
    addi    t0, s7, -5
    add     s8, s8, t0
    add     s9, s9, t0
    mv      a0, a1
    mv      a1, a2
    mv      a2, a3

    # Follow the table: a0, a1 a state in it, a2 its code.
.Lwalk:
    sltiu   t0, a2, 9
    beqz    t0, .Lemit              # 9: solved
    sb      a2, 0(s8)
    addi    s8, s8, 1
    addi    s9, s9, 1
    jal     .Lapply
    jal     .Lnear
    j       .Lwalk

    # Distance 11: the first neighbor not marked is at distance 10.
.Lfar:
    li      s3, 0
1:  mv      a0, s5
    mv      a1, s6
    mv      a2, s3
    jal     .Lapply
    jal     .Lnear
    li      t0, 15
    bne     a2, t0, 2f
    addi    s3, s3, 1
    j       1b
2:  sb      s3, 0(s8)
    addi    s8, s8, 1
    li      s9, 1
    mv      s5, a0
    mv      s6, a1
    li      s7, 10
    j       .Ldeepen

    # s9 moves at 0(sp). Each takes its letter and a space or the NUL, plus
    # a second character for a half or counterclockwise turn.
.Lemit:
    lla     t5, .Lmove_tab
    li      t1, 0                   # second characters
    li      t0, 0
1:  beq     t0, s9, 2f
    add     t2, sp, t0
    lbu     t2, 0(t2)
    add     t2, t2, t5
    lbu     t2, 27(t2)
    snez    t2, t2
    add     t1, t1, t2
    addi    t0, t0, 1
    j       1b
2:  slli    t0, s9, 1               # 2n, or 1 for the solved cube
    add     t0, t0, t1
    bnez    s9, 3f
    li      t0, 1
3:  bgeu    s11, t0, 4f
    li      a0, -2                  # does not fit
    beqz    s11, .Lret
    sb      zero, 0(s10)
    j       .Lret
4:  mv      t3, s10
    li      t0, 0
    li      t4, ' '
5:  beq     t0, s9, 7f
    beqz    t0, 6f
    sb      t4, 0(t3)
    addi    t3, t3, 1
6:  add     t2, sp, t0
    lbu     t2, 0(t2)
    add     t2, t2, t5
    lbu     t1, 18(t2)              # letter
    sb      t1, 0(t3)
    addi    t3, t3, 1
    lbu     t1, 27(t2)              # '2', '\'' or nothing
    beqz    t1, 8f
    sb      t1, 0(t3)
    addi    t3, t3, 1
8:  addi    t0, t0, 1
    j       5b
7:  sb      zero, 0(t3)
    mv      a0, s9
    j       .Lret

.Linvalid:
    li      a0, -1
.Lret:
    lw      ra, 60(sp)
    lw      s0, 56(sp)
    lw      s1, 52(sp)
    lw      s2, 48(sp)
    lw      s3, 44(sp)
    lw      s4, 40(sp)
    lw      s5, 36(sp)
    lw      s6, 32(sp)
    lw      s7, 28(sp)
    lw      s8, 24(sp)
    lw      s9, 20(sp)
    lw      s10, 16(sp)
    lw      s11, 12(sp)
    addi    sp, sp, 64
    ret
    .size   rubik_solve_asm, .-rubik_solve_asm

# .Lsearch: depth-first search below a node.
#   in   a0 P, a1 O, a2 moves left (6 or more), a3 face turned last (3: none),
#        a4 where the move to a child goes
#   out  a0 1 and a1 P, a2 O, a3 code of the table state reached, or a0 0
# With 6 moves left the children end the sequence: .Lfrontier. Otherwise the
# node keeps P, O in s4, s5, the child's moves left in s6, the last face in
# s7 and the move address in s8 across the calls; s9, s10 follow the turns.
.Lsearch:
    li      t0, 6
    beq     a2, t0, .Lfrontier
    addi    sp, sp, -32
    sw      ra, 28(sp)
    sw      s4, 24(sp)
    sw      s5, 20(sp)
    sw      s6, 16(sp)
    sw      s7, 12(sp)
    sw      s8, 8(sp)
    sw      s9, 4(sp)
    sw      s10, 0(sp)
    mv      s4, a0
    mv      s5, a1
    addi    s6, a2, -1
    mv      s7, a3
    mv      s8, a4
    SFACE   0
    SFACE   1
    SFACE   2
    li      a0, 0
.Lsearch_ret:
    lw      ra, 28(sp)
    lw      s4, 24(sp)
    lw      s5, 20(sp)
    lw      s6, 16(sp)
    lw      s7, 12(sp)
    lw      s8, 8(sp)
    lw      s9, 4(sp)
    lw      s10, 0(sp)
    addi    sp, sp, 32
    ret

# .Lfrontier: the last level, a leaf routine. Same interface as .Lsearch with
# 6 moves left; a0, a1 keep the node, a5, a6 follow the turns.
.Lfrontier:
    li      a2, 5
    FFACE   0
    FFACE   1
    FFACE   2
    li      a0, 0
    ret
.Lfrontier_hit:                     # t1 move, t0 code, a5, a6 the child
    sb      t1, 0(a4)
    li      a0, 1
    mv      a1, a5
    mv      a2, a6
    mv      a3, t0
    ret

# .Lnear: look a state up in the table.
#   in   a0 P, a1 O        out  a2 code (0..15), or -1 when absent
# a0, a1 are kept; t0..t6 are used.
.Lnear:
    lhu     t5, 6(a0)
    lhu     t6, 14(a0)
    sub     t4, a1, s1
    slli    t4, t4, 1               # key: o << 4
    li      a2, -1
    beq     t5, t6, 2f
    add     t5, t5, s2
    add     t6, t6, s2
1:  lhu     t0, 0(t5)
    bgeu    t0, t4, 3f
    addi    t5, t5, 2
    bne     t5, t6, 1b
2:  ret
3:  sub     t0, t0, t4
    li      t1, 16
    bgeu    t0, t1, 2b
    mv      a2, t0
    ret

# .Lapply: apply move a2 (0..8) to a0 P, a1 O in place; t0..t3 are used.
.Lapply:
    lla     t0, .Lmove_tab
    add     t0, t0, a2
    lbu     t1, 0(t0)               # 2 * face: the column
    lbu     t2, 9(t0)               # quarter turns
1:  add     t3, a0, t1
    lhu     t3, 0(t3)
    andi    t3, t3, -8
    add     a0, s0, t3
    add     t3, a1, t1
    lhu     t3, 0(t3)
    andi    t3, t3, -8
    add     a1, s1, t3
    addi    t2, t2, -1
    bnez    t2, 1b
    ret

    .section .rodata
# Per move R R2 R' B B2 B' D D2 D': column (2 * face), quarter turns, letter,
# second character (0 for none).
.Lmove_tab:
    .byte   0, 0, 0, 2, 2, 2, 4, 4, 4
    .byte   1, 2, 3, 1, 2, 3, 1, 2, 3
    .ascii  "RRRBBBDDD"
    .byte   0, '2', '\'', 0, '2', '\'', 0, '2', '\''
