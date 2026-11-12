/*
 * rubiktab: generate the lookup tables of the 2x2x2 solver in src/rubik/.
 *
 *   cc -O2 -o rubiktab tools/rubiktab.c
 *   ./rubiktab c > src/rubik/rubik_tab.c     tables for rubik_solve
 *   ./rubiktab asm > src/rubik/rubik_tab.S   the same tables for the assembly,
 *                                            with column 3 of perm[] doubled
 *                                            into a byte offset
 *
 * The cube model is the one of src/rubik/rubik.h (sysprog21/minirubik): seven
 * cubies, quarter turns of R, B and D. A state is ranked as a permutation
 * (Lehmer code, 0..5039) and an orientation (base 3 over the first six
 * cubies, 0..728); the seventh orientation follows from the other six.
 *
 * Three tables come out, all of 16-bit words:
 *
 *   perm[5041][4]    perm[p][f], f < 3: the permutation after a quarter turn
 *                    of face f, as (rank << 3) | dp, where dp is the fewest
 *                    moves that solve that permutation alone (at most 7).
 *                    perm[p][3]: the index in near[] of the first entry with
 *                    permutation p. The extra row 5040 closes the last range.
 *   orient[729][4]   orient[o][f], f < 3: the orientation after a quarter turn
 *                    of face f, as (rank << 3) | dor, where dor is the fewest
 *                    moves that solve the orientation alone (at most 6).
 *                    orient[o][3] is 0; it only makes a row 8 bytes long.
 *   near[14868]      every state within 5 moves of solved (12,224) and every
 *                    state at distance 11 (2,644), grouped by permutation and
 *                    sorted by orientation inside a group, as (o << 4) | code.
 *                    code 0..8 is a move that brings the state one move closer
 *                    to solved (R R2 R' B B2 B' D D2 D'), 9 marks the solved
 *                    state, 15 marks a state at distance 11.
 *
 * The full breadth-first search over 3,674,160 states runs here, on the host;
 * the solver only carries the 14,868 states above and searches the rest.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
    CUBIES = 7,
    NPERM = 5040,
    NORIENT = 729,
    NSTATES = NPERM * NORIENT,
    CODE_SOLVED = 9,
    CODE_FAR = 15,
};

typedef struct {
    uint8_t p[CUBIES], o[CUBIES];
} cube_t;

/* Each destination takes a cubie from source[face][destination]. */
static const uint8_t source[3][CUBIES] = {
    {1, 4, 2, 0, 3, 5, 6},
    {0, 1, 2, 4, 5, 6, 3},
    {0, 2, 5, 3, 1, 4, 6},
};
static const uint8_t twist[3][CUBIES] = {
    {1, 2, 0, 2, 1, 0, 0},
    {0, 0, 0, 1, 2, 1, 2},
    {0, 0, 0, 0, 0, 0, 0},
};

static cube_t quarter_turn(const cube_t *c, int face)
{
    cube_t r;
    for (int i = 0; i < CUBIES; i++) {
        int from = source[face][i];
        r.p[i] = c->p[from];
        r.o[i] = (uint8_t) ((c->o[from] + twist[face][i]) % 3);
    }
    return r;
}

static int rank_perm(const cube_t *c)
{
    int p = 0;
    for (int i = 0; i < CUBIES; i++) {
        int smaller = 0;
        for (int j = i + 1; j < CUBIES; j++)
            smaller += c->p[j] < c->p[i];
        p = p * (CUBIES - i) + smaller;
    }
    return p;
}

static int rank_orient(const cube_t *c)
{
    int o = 0;
    for (int i = 0; i < CUBIES - 1; i++)
        o = o * 3 + c->o[i];
    return o;
}

static void unrank_perm(int p, cube_t *c)
{
    uint8_t avail[CUBIES] = {0, 1, 2, 3, 4, 5, 6};
    int f = 720;
    for (int i = 0; i < CUBIES; i++) {
        int q = p / f;
        p %= f;
        c->p[i] = avail[q];
        for (int j = q; j + 1 < CUBIES - i; j++)
            avail[j] = avail[j + 1];
        if (i < CUBIES - 2)
            f /= CUBIES - 1 - i;
    }
}

static void unrank_orient(int o, cube_t *c)
{
    int sum = 0;
    for (int i = CUBIES - 1; i-- > 0;) {
        c->o[i] = (uint8_t) (o % 3);
        sum += c->o[i];
        o /= 3;
    }
    c->o[CUBIES - 1] = (uint8_t) ((3 - sum % 3) % 3);
}

static uint16_t next_perm[3][NPERM], next_orient[3][NORIENT];
static uint8_t dist_perm[NPERM], dist_orient[NORIENT];
static uint16_t perm[NPERM + 1][4], orient[NORIENT][4];
static uint16_t near[16384];
static int nnear;

/* Breadth-first search of a projection (permutation or orientation alone). */
static void project_bfs(int n, uint16_t next[3][NPERM], uint8_t *d)
{
    int *queue = malloc(sizeof *queue * (size_t) n), head = 0, tail = 1;
    if (!queue)
        exit(1);
    memset(d, 0xFF, (size_t) n);
    queue[0] = 0;
    d[0] = 0;
    while (head < tail) {
        int x = queue[head++];
        for (int f = 0; f < 3; f++) {
            int y = x;
            for (int t = 0; t < 3; t++) {
                y = next[f][y];
                if (d[y] == 0xFF) {
                    d[y] = (uint8_t) (d[x] + 1);
                    queue[tail++] = y;
                }
            }
        }
    }
    free(queue);
}

static uint32_t apply(uint32_t s, int move)
{
    int p = (int) (s / NORIENT), o = (int) (s % NORIENT), f = move / 3;
    for (int t = 0; t <= move % 3; t++) {
        p = next_perm[f][p];
        o = next_orient[f][o];
    }
    return (uint32_t) p * NORIENT + (uint32_t) o;
}

static void build(void)
{
    cube_t c = {{0}, {0}};
    for (int r = 0; r < NPERM; r++) {
        unrank_perm(r, &c);
        for (int f = 0; f < 3; f++) {
            cube_t n = quarter_turn(&c, f);
            next_perm[f][r] = (uint16_t) rank_perm(&n);
        }
    }
    for (int i = 0; i < CUBIES; i++)
        c.p[i] = (uint8_t) i;
    for (int r = 0; r < NORIENT; r++) {
        unrank_orient(r, &c);
        for (int f = 0; f < 3; f++) {
            cube_t n = quarter_turn(&c, f);
            next_orient[f][r] = (uint16_t) rank_orient(&n);
        }
    }
    project_bfs(NPERM, next_perm, dist_perm);
    /* The orientation table is narrower; widen it for the shared routine. */
    static uint16_t wide[3][NPERM];
    for (int f = 0; f < 3; f++)
        memcpy(wide[f], next_orient[f], sizeof next_orient[f]);
    project_bfs(NORIENT, wide, dist_orient);

    uint8_t *dist = malloc(NSTATES);
    uint32_t *queue = malloc(sizeof *queue * NSTATES), head = 0, tail = 1;
    if (!dist || !queue)
        exit(1);
    memset(dist, 0xFF, NSTATES);
    queue[0] = 0;
    dist[0] = 0;
    while (head < tail) {
        uint32_t x = queue[head++];
        for (int m = 0; m < 9; m++) {
            uint32_t y = apply(x, m);
            if (dist[y] == 0xFF) {
                dist[y] = (uint8_t) (dist[x] + 1);
                queue[tail++] = y;
            }
        }
    }
    if (tail != NSTATES || dist[queue[NSTATES - 1]] != 11) {
        fputs("rubiktab: unexpected state graph\n", stderr);
        exit(1);
    }
    free(queue);

    for (int p = 0; p < NPERM; p++) {
        perm[p][3] = (uint16_t) nnear;
        for (int o = 0; o < NORIENT; o++) {
            uint32_t s = (uint32_t) p * NORIENT + (uint32_t) o;
            int d = dist[s], code = -1;
            if (d == 0) {
                code = CODE_SOLVED;
            } else if (d == 11) {
                code = CODE_FAR;
            } else if (d <= 5) {
                for (int m = 0; m < 9 && code < 0; m++)
                    if (dist[apply(s, m)] == d - 1)
                        code = m;
            }
            if (code >= 0)
                near[nnear++] = (uint16_t) (o << 4 | code);
        }
    }
    perm[NPERM][3] = (uint16_t) nnear;
    free(dist);

    for (int p = 0; p < NPERM; p++)
        for (int f = 0; f < 3; f++) {
            int q = next_perm[f][p];
            perm[p][f] = (uint16_t) (q << 3 | dist_perm[q]);
        }
    for (int o = 0; o < NORIENT; o++)
        for (int f = 0; f < 3; f++) {
            int q = next_orient[f][o];
            orient[o][f] = (uint16_t) (q << 3 | dist_orient[q]);
        }
}

/* Rows of four words are braced, two rows to a line; a flat array is
 * written eight words to a line.
 */
static void emit_c_array(const char *decl, const uint16_t *v, int n, int rows)
{
    printf("%s = {\n", decl);
    for (int i = 0; i < n; i += 8) {
        printf("   ");
        for (int j = i; j < n && j < i + 8; j++) {
            const char *open = rows && j % 4 == 0 ? "{" : "";
            const char *close = rows && j % 4 == 3 ? "}" : "";
            printf("%s%s0x%04X%s,", rows && j % 4 ? "" : " ", open, v[j],
                   close);
            if (rows && j % 4 != 3)
                printf(" ");
        }
        printf("\n");
    }
    printf("};\n\n");
}

static void emit_asm_array(const char *name, const uint16_t *v, int n)
{
    printf("    .globl  %s\n", name);
    printf("    .type   %s, @object\n", name);
    printf("    .balign 8\n%s:\n", name);
    for (int i = 0; i < n; i += 8) {
        printf("    .half   ");
        for (int j = i; j < n && j < i + 8; j++)
            printf("0x%04X%s", v[j], j + 1 < n && j + 1 < i + 8 ? ", " : "");
        printf("\n");
    }
    printf("    .size   %s, .-%s\n\n", name, name);
}

int main(int argc, char **argv)
{
    if (argc != 2 || (strcmp(argv[1], "c") && strcmp(argv[1], "asm"))) {
        fputs("usage: rubiktab c|asm\n", stderr);
        return 2;
    }
    build();
    if (!strcmp(argv[1], "c")) {
        printf(
            "/* Generated by tools/rubiktab.c. Do not edit. The layout is "
            "described\n * in rubik_tab.h.\n */\n"
            "#include \"rubik_tab.h\"\n\n/* clang-format off */\n");
        emit_c_array("const uint16_t rubik_perm[RUBIK_NPERM + 1][4]",
                     &perm[0][0], (NPERM + 1) * 4, 1);
        emit_c_array("const uint16_t rubik_orient[RUBIK_NORIENT][4]",
                     &orient[0][0], NORIENT * 4, 1);
        emit_c_array("const uint16_t rubik_near[RUBIK_NNEAR]", near, nnear, 0);
        printf("/* clang-format on */\n");
    } else {
        /* The assembly indexes rubik_near by byte offset. */
        for (int p = 0; p <= NPERM; p++)
            perm[p][3] = (uint16_t) (perm[p][3] * 2);
        printf(
            "# Generated by tools/rubiktab.c. Do not edit. The layout is "
            "described in\n# rubik_tab.h; these are the same tables under the "
            "names of the assembly,\n# except that column 3 of "
            "rubik_asm_perm is a byte offset into rubik_asm_near.\n\n    "
            ".section .rodata\n\n");
        emit_asm_array("rubik_asm_perm", &perm[0][0], (NPERM + 1) * 4);
        emit_asm_array("rubik_asm_orient", &orient[0][0], NORIENT * 4);
        emit_asm_array("rubik_asm_near", near, nnear);
    }
    fprintf(stderr, "rubiktab: %d near entries\n", nnear);
    return 0;
}
