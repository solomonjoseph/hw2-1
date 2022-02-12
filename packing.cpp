#include <immintrin.h>
#include <string.h>


#define bool int
#define true 1
#define false 0

#define max(a, b) (((a) < (b)) ? (b) : (a))
#define min(a, b) (((a) < (b)) ? (a) : (b))

const int levels = 3;
const int ms[] = { 1, 4, 32, 128 };
const int ns[] = { 1, 4, 32, 128 };

static int calc_blocked_size(int lda, int block_size) {
    int newblocks = lda / block_size;
    if (lda % block_size != 0) {
        newblocks += 1;
    }
    return newblocks * block_size;
}

static inline int offset_from_coords(int r, int c, int H, int W, int dr, int dc, bool row_wise) {
    if (row_wise) {
        return r * W + c * min(dr, H-r);
    } else {
        return c * H + r * min(dc, W-c);
    }
}

static int recursive_pack_helper(double* dest, double* src, int lda, int M, int N, int level, const int* ms, const int* ns, bool col_wise) {
    int dr = ms[level];
    int dc = ns[level];

    int r_max = calc_blocked_size(M, ms[max(level, 1)]);
    int c_max = calc_blocked_size(N, ns[max(level, 1)]);

    int i = 0;
    for (int r = 0; r < r_max; r += dr) {
        for (int c = 0; c < c_max; c += dc) {
            int target = col_wise ? (c*lda + r) : (c + r*lda);
            if (level == 0) {
                if (r < M && c < N) {
                    dest[i++] = src[target];
                    src[target] = -1;
                } else {
                    dest[i++] = 0; // maybe memset before to get less branching?
                }
            } else {
                int M_next = min(dr, M-r);
                int N_next = min(dc, N-c);
                i += recursive_pack_helper(dest + i, src + target, lda, M_next, N_next, level-1, ms, ns, col_wise);
            }
        }
    }

    return i;
}

static double* recursive_pack(double* src, int lda, int levels, const int* ms, const int* ns, bool col_wise) {
    int M = calc_blocked_size(lda, ms[1]);
    int N = calc_blocked_size(lda, ns[1]);
    double* dest = _mm_malloc(M * N * sizeof(double), 64); // align to 64-byte boundary
    recursive_pack_helper(dest, src, lda, lda, lda, levels, ms, ns, col_wise);
    return dest;
}

static int recursive_unpack_helper(double* dest, double* src, int lda, int M, int N, int level, const int* ms, const int* ns, bool col_wise) {
    int dr = ms[level];
    int dc = ns[level];

    int r_max = calc_blocked_size(M, ms[max(level, 1)]);
    int c_max = calc_blocked_size(N, ns[max(level, 1)]);

    int i = 0;
    for (int r = 0; r < r_max; r += dr) {
        for (int c = 0; c < c_max; c += dc) {
            int target = col_wise ? (c*lda + r) : (c + r*lda);
            if (level == 0) {
                if (r < M && c < N) {
                    dest[target] = src[i++];
                } else {
                    i++;
                }
            } else {
                int M_next = min(dr, M-r);
                int N_next = min(dc, N-c);
                i += recursive_unpack_helper(dest + target, src + i, lda, M_next, N_next, level-1, ms, ns, col_wise);
            }
        }
    }

    return i;
}

static void recursive_unpack(double* src, double* dest, int lda, int levels, const int* ms, const int* ns, bool col_wise) {
    recursive_unpack_helper(dest, src, lda, lda, lda, levels, ms, ns, col_wise);
}
