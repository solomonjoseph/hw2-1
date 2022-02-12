#include "common.h"
#include <cmath>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <cstdint>

/*
We will preduce B bins for N particles.
The naive runtime is O(N^2). When split into bins,
we end up with O(B * (N/B)^2) = O(N^2/B) work.
Thus, we set B = aN to get a total of O(N/a) work.
We want each bin to roughly fill up the cache,
so the general goal is to maximize a by increasing it until
it starts to make things work.

Ideally, B = aN, but in reality we need a power of 2 bins.
Can try rounding either up or down, but I think rounding down the bin size
makes the most sense. Instead of sqrt(aN) bins on each side, we need
the nearest power of two *above* sqrt(aN).

Even more in reality, bin capacity fitting in cache is probably the sole most 
important factor for performance, so that will dictate things. Luckily density
is constant, meaning that size and N increase at the same rate, leading 
bin size to directly correspond to a.
*/


//https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
inline int round_up_pow2(const unsigned int n) {
    unsigned int v = n;
    --v;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    return ++v;
}

/*
Bin capacity is what determines actual performance due to cache limitations...
So I think actually fixing bin capacity is best, calculate others from there.
Number of bins is then N / bin_capacity -> bins per side is that rounded 
*/
template <unsigned int bin_capacity = 400, class particle = particle_t>
struct bin_store {
    const unsigned int N;
    const unsigned int num_bins_per_side;
    const double bin_width;
    const double size;
    const unsigned int num_bins;

    //Actual memory backing for the bins.
    //Will change this data type if we rearrange the struct or anything.
    //Will need to be bin_capacity * num_bins big.
    particle *bins;

    static inline unsigned int compute_bins_per_side(const unsigned int N) {
        unsigned int num_bins = N / bin_capacity;
        unsigned int bps = ceil(sqrt(num_bins));
        return round_up_pow2(bps);
    }

    //Access the ith vector in 
    static inline unsigned int index(const unsigned int x, const unsigned int y, const unsigned int i) {
        //TODO: implement Z curve mapping here
        return 0;
    }

    bin_store(const unsigned int N, const double size) : N(N), size(size), num_bins_per_side(compute_bins_per_side(N)),
              num_bins = num_bins_per_side * num_bins_per_side, bin_width(size / num_bins_per_side) {
        bins = align(64) new particle[num_bins * bin_capacity];
    }

    ~bin_store() {
        delete bins;
    }

    //TODO: implement operator[] for get/set
};

unsigned int num_bins;
unsigned int bin_width;

using std::vector;
using std::set;
using std::max;

vector<vector<set<int>>> bins;

template<int bin_exp>
class BinStore {

};

// Apply the force from neighbor to particle
void apply_force(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

void apply_force_bin(particle_t* parts, int bin_x, int bin_y) {
    set<int>::iterator itr_i, itr_j;
    set<int> bin = bins[bin_x][bin_y];
    for (itr_i = bin.begin(); itr_i != bin.end(); itr_i++) {
        for (itr_j = bin.begin(); itr_j != bin.end(); itr_j++) {
            apply_force(parts[*itr_i], parts[*itr_j]);
        }
    }
}

void apply_force_between_bins(particle_t* parts, int bin_one_x, int bin_one_y, int bin_two_x, int bin_two_y) {
    set<int>::iterator itr_i, itr_j;
    set<int> bin_one = bins[bin_one_x][bin_one_y];
    set<int> bin_two = bins[bin_two_x][bin_two_y];
    for (itr_i = bin_one.begin(); itr_i != bin_one.end(); itr_i++) {
        for (itr_j = bin_two.begin(); itr_j != bin_two.end(); itr_j++) {
            apply_force(parts[*itr_i], parts[*itr_j]);
        }
    }
}

// Integrate the ODE
void move(particle_t& p, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }
}


void init_simulation(particle_t* parts, int num_parts, double size) {
	// You can use this space to initialize static, global data objects
    // that you may need. This function will be called once before the
    // algorithm begins. Do not do any particle simulation here
    bin_size = max(min_r, sqrt(cutoff / density));
    num_bins = size / bin_size;
    bins = vector<vector<set<int>>>(num_bins, vector<set<int>>(num_bins, set<int>()));
    for (int i = 0; i < num_parts; i++) {
        int x_bin = std::min(num_bins - 1, (int) (parts[i].x / bin_size));
        int y_bin = std::min(num_bins - 1, (int) (parts[i].y / bin_size));
        bins[x_bin][y_bin].insert(i);
    }
}


void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // Compute Forces
    for (int i = 0; i < num_bins; i++) {
        for (int j = 0; j < num_bins; j++) {
            apply_force_bin(parts, i, j);
            if (j > 0) {
                apply_force_between_bins(parts, i, j, i, j-1);
                if (i > 0) {
                    apply_force_between_bins(parts, i, j, i-1, j);
                    apply_force_between_bins(parts, i, j, i-1, j-1);
                }
                if (i < num_bins - 1) {
                    apply_force_between_bins(parts, i, j, i+1, j-1);
                }
            }
        }
    }

    // Move Particles
    for (int i = 0; i < num_parts; ++i) {
        double x_old = parts[i].x;
        double y_old = parts[i].y;
        move(parts[i], size);
        if (x_old != parts[i].x || y_old != parts[i].y) {
            int x_bin_old = std::min(num_bins - 1, (int) (x_old / bin_size));
            int y_bin_old = std::min(num_bins - 1, (int) (y_old / bin_size));
            int x_bin_new = std::min(num_bins - 1, (int) (parts[i].x / bin_size));
            int y_bin_new = std::min(num_bins - 1, (int) (parts[i].y / bin_size));
            bins[x_bin_old][y_bin_old].erase(i);
            bins[x_bin_new][y_bin_new].insert(i);
        }
    }
}
