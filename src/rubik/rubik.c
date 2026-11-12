#include "rubik.h"
#include "rubik_tab.h"

/* Optimal 2x2x2 solver: iterative deepening down to a table of every state
 * within NEAR moves of solved.
 *
 * rubik_near (see rubik_tab.h) holds the 12,224 states within NEAR = 5 moves,
 * each with a move that brings it one move closer, and marks the 2,644 states
 * at distance 11. Anything else is at distance 6..10 and is found by a
 * depth-first search for a state of the table:
 *
 *   for limit = 6, 7, ..., 10:
 *       search every move sequence of length limit - NEAR from the start
 *       whose states keep max(dp, dor) within the moves left; at the end of
 *       each sequence, look the state up in rubik_near
 *
 * A hit at the end of a sequence of length limit - NEAR is a state within NEAR
 * moves, so it gives a solution of at most limit moves. No shorter limit
 * succeeded, and every state within NEAR moves is in the table, so the
 * solution has exactly limit moves and limit is the distance: the first hit
 * is optimal. The table then supplies the last NEAR moves one at a time.
 *
 * A state at distance 11 is recognized in the table. Each of its neighbors is
 * at distance 10 or 11; the first one that is not marked is solved with the
 * limit set to 10 at once, after the move that reaches it.
 *
 * dp and dor are exact distances of the permutation and of the orientation
 * alone, so max(dp, dor) never overestimates and the pruning never cuts an
 * optimal sequence. Only one face is excluded after a move: turning the same
 * face again is never shorter than a single turn of that face.
 */

#define NEAR 5   /* rubik_near holds every state within NEAR moves */
#define LIMIT 10 /* the longest search; distance 11 comes from the table */
#define MAX_MOVES 11
#define NO_FACE 3

/* Move m turns face m / 3 by m % 3 + 1 quarter turns; tables avoid the
 * division, which RV32I has no instruction for.
 */
static const uint8_t move_face[9] = {0, 0, 0, 1, 1, 1, 2, 2, 2};
static const uint8_t move_turns[9] = {1, 2, 3, 1, 2, 3, 1, 2, 3};
static const char move_name[9][3] = {"R",  "R2", "R'", "B", "B2",
                                     "B'", "D",  "D2", "D'"};

/* The code of state (p, o) in rubik_near, or -1 when it is not there. A group
 * holds at most 17 entries, so a linear scan that stops at the first key not
 * below the wanted one is enough.
 */
static int near_code(unsigned p, unsigned o)
{
    unsigned key = o << 4;
    for (unsigned i = rubik_perm[p][3]; i < rubik_perm[p + 1][3]; i++) {
        unsigned e = rubik_near[i];
        if (e >= key)
            return e < key + 16 ? (int) (e & 15) : -1;
    }
    return -1;
}

static void apply(unsigned *p, unsigned *o, int move)
{
    int face = move_face[move];
    for (int t = 0; t < move_turns[move]; t++) {
        *p = rubik_perm[*p][face] >> 3;
        *o = rubik_orient[*o][face] >> 3;
    }
}

/* Follow the table from a state of it to solved, writing the moves. */
static int walk(unsigned p, unsigned o, int code, uint8_t *moves)
{
    int n = 0;
    while (code >= 0 && code < RUBIK_CODE_SOLVED && n < NEAR) {
        moves[n++] = (uint8_t) code;
        apply(&p, &o, code);
        code = near_code(p, o);
    }
    return n;
}

/* Depth-first search below (p, o), which was reached in `depth` moves and has
 * `left` moves left; `last` is the face turned last. Returns the length of
 * the solution written to moves[], or 0 when there is none within the limit.
 */
static int search(unsigned p,
                  unsigned o,
                  int left,
                  int last,
                  uint8_t *moves,
                  int depth)
{
    for (int face = 0; face < 3; face++) {
        if (face == last)
            continue;
        unsigned q = p, r = o;
        /* One, two and three quarter turns of the face, in turn. */
        for (int t = 0; t < 3; t++) {
            unsigned vq = rubik_perm[q][face], vr = rubik_orient[r][face];
            q = vq >> 3;
            r = vr >> 3;
            /* The child has left - 1 moves; prune when a projection needs
             * more.
             */
            if ((int) (vq & 7) >= left || (int) (vr & 7) >= left)
                continue;
            moves[depth] = (uint8_t) (face * 3 + t);
            if (left - 1 == NEAR) {
                int code = near_code(q, r);
                if (code >= 0 && code != RUBIK_CODE_FAR)
                    return depth + 1 + walk(q, r, code, moves + depth + 1);
            } else {
                int n = search(q, r, left - 1, face, moves, depth + 1);
                if (n)
                    return n;
            }
        }
    }
    return 0;
}

static int deepen(unsigned p, unsigned o, int first, uint8_t *moves)
{
    for (int limit = first; limit <= LIMIT; limit++) {
        int n = search(p, o, limit, NO_FACE, moves, 0);
        if (n)
            return n;
    }
    return 0; /* not reached: every state is within 11 moves */
}

static int solve(unsigned p, unsigned o, uint8_t *moves)
{
    int code = near_code(p, o);
    if (code == RUBIK_CODE_FAR) {
        for (int m = 0; m < 9; m++) {
            unsigned q = p, r = o;
            apply(&q, &r, m);
            if (near_code(q, r) == RUBIK_CODE_FAR)
                continue;
            moves[0] = (uint8_t) m;
            return 1 + deepen(q, r, LIMIT, moves + 1);
        }
    }
    if (code >= 0)
        return walk(p, o, code, moves);
    return deepen(p, o, NEAR + 1, moves);
}

/* Validate the 14 characters and rank the state. Every character is checked
 * before the next one is read, so a short string ends the scan at its NUL.
 */
static int parse(const char *state, unsigned *pp, unsigned *po)
{
    uint8_t cubie[7];
    unsigned seen = 0, sum = 0, o = 0, p = 0;

    for (int i = 0; i < 7; i++) {
        unsigned d = (unsigned) (unsigned char) state[i] - '1';
        if (d >= 7 || ((seen >> d) & 1u))
            return 0;
        seen |= 1u << d;
        cubie[i] = (uint8_t) d;
    }
    for (int i = 0; i < 7; i++) {
        unsigned d = (unsigned) (unsigned char) state[7 + i] - '1';
        if (d >= 3)
            return 0;
        sum += d;
        if (i < 6)
            o = o * 3 + d;
    }
    /* The orientations must sum to 0, 3, 6, 9 or 12: bits of 0x1249. */
    if (state[14] != '\0' || !((0x1249u >> sum) & 1u))
        return 0;

    /* Lehmer code: for each cubie, how many later cubies are smaller. */
    for (int i = 0; i < 7; i++) {
        unsigned smaller = 0;
        for (int j = i + 1; j < 7; j++)
            smaller += cubie[j] < cubie[i];
        p = p * (unsigned) (7 - i) + smaller;
    }
    *pp = p;
    *po = o;
    return 1;
}

int rubik_solve(const char *state, char *out, size_t cap)
{
    unsigned p, o;
    uint8_t moves[MAX_MOVES];

    if (!parse(state, &p, &o))
        return -1;
    int n = solve(p, o, moves);

    /* Names of one or two characters, single spaces between them, a NUL. */
    size_t need = 1;
    for (int i = 0; i < n; i++)
        need += (size_t) (move_name[moves[i]][1] ? 2 : 1) + (i > 0);
    if (need > cap) {
        if (cap)
            out[0] = '\0';
        return -2;
    }
    char *w = out;
    for (int i = 0; i < n; i++) {
        if (i)
            *w++ = ' ';
        *w++ = move_name[moves[i]][0];
        if (move_name[moves[i]][1])
            *w++ = move_name[moves[i]][1];
    }
    *w = '\0';
    return n;
}
