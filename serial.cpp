#include "common.h"
#include "immintrin.h"
#include "smmintrin.h"
#include <cmath>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <cstdint>
#include <iostream>

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

using std::max;
using std::min;
constexpr double cutoff_squared = cutoff * cutoff;

// Apply the force from neighbor to particle
inline void apply_force(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact - too far or are same particle, then no
    if (r2 > cutoff * cutoff || (dx == 0 && dy == 0))
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

// Apply the force from neighbor to particle
inline void apply_force_symmetric(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact - too far or are same particle, then no
    if (r2 > cutoff * cutoff || (dx == 0 && dy == 0))
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
    neighbor.ax -= coef * dx;
    neighbor.ay -= coef * dy;
}

// Integrate the ODE
inline void move(particle_t& p, double size) {
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

    p.ax = p.ay = 0.0;
}

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

struct improved_particle_t {
    particle_t part;
    double last_t;
    unsigned int id;
    bool valid;
};

struct bin {
    improved_particle_t *particles = nullptr;
    unsigned int capacity = 0, count = 0;
    bool own_mem = false;

    bin() {}

    void push_back(improved_particle_t new_part) {
        if (count >= capacity) {
            improved_particle_t *old_parts = particles;
            particles = new improved_particle_t[capacity * 2];
            std::copy(old_parts, old_parts + capacity, particles);
            capacity *= 2;
            if (own_mem) {
                delete[] old_parts;
            } else {
                own_mem = true;
            }
        }
        particles[count++] = new_part;
    }

    //order doesn't matter, so keep things packed
    //IMPORTANT: do NOT remove _during_ simulation. remove particles BETWEEN simulation steps
    void remove(unsigned int index) {
        if (particles[index].valid) {
            particles[index] = particles[--count];
            particles[count].valid = false;
        }
    }

    improved_particle_t &operator[](int i) {
        return particles[i];
    }

    void flush(particle_t *orig_backing) {
        for (unsigned int i = 0; i < count; ++i) {
            improved_particle_t this_part = particles[i];
            orig_backing[this_part.id] = this_part.part;
        }
    }

    ~bin() {
        if (own_mem) {
            delete particles;
        }
    }
};

/*
1. initialize our data structure, pack in the improved_particle_t's from initialization data
2. simulation is as before, but now we traverse bins in Z order so that neighbors are near
3. unpack back into the original values based on the ID's
*/

/*
Questions for Ben:
    - why is class particle templated
    - why struct over class
*/

/*
Bin capacity is what determines actual performance due to cache limitations...
So I think actually fixing bin capacity is best, calculate others from there.
Number of bins is then N / bin_capacity -> bins per side is that rounded 
*/
struct bin_store {
    const unsigned int N;
    unsigned int num_bins_per_side;
    double bin_width;
    const double size;
    unsigned int num_bins;
    unsigned int bin_capacity;
    unsigned int steps;

    //Actual memory backing for the bins.
    //Will change this data type if we rearrange the struct or anything.
    //Will need to be bin_capacity * num_bins big.
    std::vector<bin> bins;
    improved_particle_t *particle_mem;
    particle_t *base_array;

    inline unsigned int compute_bins_per_side(const unsigned int N) {
        unsigned int num_bins = N / bin_capacity;
        unsigned int bps = ceil(sqrt(num_bins));
        return round_up_pow2(bps);
    }
    
    static inline unsigned int calc_Z_order(uint16_t i, uint16_t j) {
        //Attribution: https://stackoverflow.com/a/14853492
        static constexpr unsigned int MASKS[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF};
        static constexpr unsigned int SHIFTS[] = {1, 2, 4, 8};

        unsigned int x = i;  // Interleave lower 16 bits of x and y, so the bits of x
        unsigned int y = j;  // are in the even positions and bits from y in the odd;

        //TODO: microopt could be AVX-ifying this
        x = (x | (x << SHIFTS[3])) & MASKS[3];
        x = (x | (x << SHIFTS[2])) & MASKS[2];
        x = (x | (x << SHIFTS[1])) & MASKS[1];
        x = (x | (x << SHIFTS[0])) & MASKS[0];

        y = (y | (y << SHIFTS[3])) & MASKS[3];
        y = (y | (y << SHIFTS[2])) & MASKS[2];
        y = (y | (y << SHIFTS[1])) & MASKS[1];
        y = (y | (y << SHIFTS[0])) & MASKS[0];

        const unsigned int result = x | (y << 1);
        return result;
    }

    //Access the ith vector in the bin at 
    inline bin &get_bin(const unsigned int x, const unsigned int y) {
        // return &bins[calc_Z_order(x, y) * bin_capacity];
        return bins[calc_Z_order(x, y)];
    }

    inline bin &get_bin_from_coord(const double x, const double y) {
        unsigned int x_int = x / bin_width;
        unsigned int y_int = y / bin_width;
        return get_bin(x_int, y_int);
    }

    void write_back() {
        for (int i = 0; i < num_bins; ++i) {
            bins[i].flush(base_array);
        }
    }

    bin_store(const unsigned int N, const double size, const unsigned int bin_capacity, particle_t* parts) : N(N), size(size), bin_capacity(bin_capacity) {
        steps = 0;
        num_bins_per_side = compute_bins_per_side(N);
        bin_width = size / num_bins_per_side;
        num_bins = num_bins_per_side * num_bins_per_side;
        std::cout << N << ", " << num_bins << ", " << bin_capacity << ", " << size << ", " << bin_width << std::endl;
        base_array = parts;
        particle_mem = new improved_particle_t[bin_capacity * num_bins];
        //std::cout << num_bins << std::endl;
        bins = std::vector<bin>(num_bins, bin());
        for (uint16_t i = 0; i < num_bins_per_side; ++i) {
            for (uint16_t j = 0; j < num_bins_per_side; ++j) {
                unsigned int ind = calc_Z_order(i, j);
                bins[ind].particles = &particle_mem[bin_capacity * ind];
                bins[ind].capacity = bin_capacity;
                bins[ind].count = 0;
                bins[ind].own_mem = false;
                //std::cout << "Assigning " << i << ", " << j << " to " << ind << std::endl;
            }
        }
        
        for (unsigned int i = 0; i < N; i++) {
            parts[i].ax = parts[i].ay = 0.0;
            bin &curr_bin = get_bin_from_coord(parts[i].x, parts[i].y);
            curr_bin.push_back((improved_particle_t) {parts[i], 0.0, i, true});
        }
    }

    static inline double d2(double x1, double y1, double x2, double y2) {
        double dx = x1 - x2;
        double dy = y1 - y2;
        return dx * dx + dy * dy;
    }

    void simulate_one_step() {
        double left_edge, right_edge, top_edge, bottom_edge;
        std::vector<improved_particle_t *> buf;
        for (int i = 0; i < num_bins_per_side; ++i) {
            for (int j = 0; j < num_bins_per_side; ++j) {
                buf.clear();
                bin &my_bin = get_bin(i, j);
                left_edge = i * bin_width;
                right_edge = (i + 1) * bin_width;
                top_edge = j * bin_width;
                bottom_edge = (j + 1) * bin_width; 
                if (i - 1 >= 0 && j + 1 < num_bins_per_side) {//bottom left
                    bin &bottom_left = get_bin(i - 1, j + 1);
                    for (int k = 0; k < bottom_left.count; ++k) {
                        improved_particle_t *other = &bottom_left[k];
                        if (d2(left_edge, bottom_edge, other->part.x, other->part.y) <= cutoff_squared) {
                            buf.push_back(other);
                        }
                    }
                }
                if (j + 1 < num_bins_per_side) {//bottom
                    bin &bottom = get_bin(i, j + 1);
                    for (int k = 0; k < bottom.count; ++k) {
                        improved_particle_t *other = &bottom[k];
                        if (other->part.y - bottom_edge <= cutoff) {
                            buf.push_back(other);
                        }
                    }
                }
                if (i + 1 < num_bins_per_side && j + 1 < num_bins_per_side) {//bottom right
                    bin &bottom_right = get_bin(i + 1, j + 1);
                    for (int k = 0; k < bottom_right.count; ++k) {
                        improved_particle_t *other = &bottom_right[k];
                        if (d2(right_edge, bottom_edge, other->part.x, other->part.y) <= cutoff_squared) {
                            buf.push_back(other);
                        }
                    }
                }
                if (i + 1 < num_bins_per_side) {//right
                    bin &right = get_bin(i + 1, j);
                    for (int k = 0; k < right.count; ++k) {
                        improved_particle_t *other = &right[k];
                        if (other->part.x - right_edge <= cutoff) {
                            buf.push_back(other);
                        }
                    }
                }

                //Now have collected all the particles we might be updating, so do the update
                for (int p = 0; p < my_bin.count; ++p) {
                    int s = buf.size();
                    particle_t &left = my_bin[p].part;
                    for (int q = 0; q < s; ++q) {
                        particle_t &right = buf[q]->part;
                        apply_force_symmetric(left, right);
                    }
                    for (int q = p + 1; q < my_bin.count; ++q) {
                        particle_t &right = my_bin[q].part;
                        apply_force_symmetric(left, right);
                    }
                }
            }
        }

        double x_offset, y_offset;
        // Move Particles
        for (int i = 0; i < num_bins_per_side; ++i) {
            for (int j = 0; j < num_bins_per_side; ++j) {
                bin &current = get_bin(i, j);
                for (int k = 0; k < current.count; ++k) {
                    particle_t &part = current[k].part;
                    move(part, size);
                }
            }
        }

        int new_i, new_j, current_index, new_index;
        double x_t, y_t;
        for (int i = 0; i < num_bins_per_side; ++i) {
            for (int j = 0; j < num_bins_per_side; ++j) {
                bin &current = get_bin(i, j);
                left_edge = i * bin_width;
                right_edge = (i + 1) * bin_width;
                top_edge = j * bin_width;
                bottom_edge = (j + 1) * bin_width; 
                for (int k = 0; k < current.count; ++k) {
                    particle_t &part = current[k].part;
                    x_t = part.x;
                    y_t = part.y;
                    new_i = i;
                    new_j = j;

                    while(x_t < left_edge) {
                        x_t += bin_width;
                        --new_i;
                    }
                    while(x_t > right_edge) {
                        x_t -= bin_width;
                        ++new_i;
                    }
                    while(y_t < top_edge) {
                        y_t += bin_width;
                        --new_j;
                    }
                    while(y_t > bottom_edge) {
                        y_t -= bin_width;
                        ++new_j;
                    }
                    
                    if (new_i != i || new_j != j) {
                        improved_particle_t tmp = current[k];
                        current.remove(k);
                        bin &other = get_bin(new_i, new_j);
                        other.push_back(tmp);
                    }
                }
            }
        }
        write_back();
    }

    ~bin_store() {
        delete particle_mem;
    }
};

bin_store *bins;


void init_simulation(particle_t* parts, int num_parts, double size) {
	// You can use this space to initialize static, global data objects
    // that you may need. This function will be called once before the
    // algorithm begins. Do not do any particle simulation here
    bins = new bin_store(num_parts, size, 32, parts);
}


void simulate_one_step(particle_t* parts, int num_parts, double size) {
    bins->simulate_one_step();
}
