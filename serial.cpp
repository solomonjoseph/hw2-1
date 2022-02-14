#include "common.h"
#include <cmath>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <cstdio>

int num_bins;
double BIN_SIZE;

#define GET_INDEX(X, Y, NUM_BINS) (X * NUM_BINS + Y)

using std::vector;
using std::unordered_set;
using std::max;

// vector<vector<set<int>>> bins;

unordered_set<int>* bins;

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
    unordered_set<int>::iterator itr_i, itr_j;
    unordered_set<int> bin = bins[GET_INDEX(bin_x, bin_y, num_bins)];
    for (itr_i = bin.begin(); itr_i != bin.end(); itr_i++) {
        for (itr_j = bin.begin(); itr_j != bin.end(); itr_j++) {
            apply_force(parts[*itr_i], parts[*itr_j]);
        }
    }
}

void apply_force_between_bins(particle_t* parts, int bin_one_x, int bin_one_y, int bin_two_x, int bin_two_y) {
    unordered_set<int>::iterator itr_i, itr_j;
    unordered_set<int> bin_one = bins[GET_INDEX(bin_one_x, bin_one_y, num_bins)];
    unordered_set<int> bin_two = bins[GET_INDEX(bin_two_x, bin_two_y, num_bins)];
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
    BIN_SIZE = std::min(size, 5*cutoff);
    num_bins = size / BIN_SIZE;
    bins = new unordered_set<int>[num_bins * num_bins];
    // bins = vector<vector<set<int>>>(num_bins, vector<set<int>>(num_bins, set<int>()));
    for (int i = 0; i < num_parts; i++) {
        int x_bin = std::min(num_bins - 1, (int) (parts[i].x / BIN_SIZE));
        int y_bin = std::min(num_bins - 1, (int) (parts[i].y / BIN_SIZE));
        bins[GET_INDEX(x_bin, y_bin, num_bins)].insert(i);
    }
    for (int i = 0; i < num_parts; i++) {
        printf("%d\n", bins[i].size());
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
            int x_bin_old = std::min(num_bins - 1, (int) (x_old / BIN_SIZE));
            int y_bin_old = std::min(num_bins - 1, (int) (y_old / BIN_SIZE));
            int x_bin_new = std::min(num_bins - 1, (int) (parts[i].x / BIN_SIZE));
            int y_bin_new = std::min(num_bins - 1, (int) (parts[i].y / BIN_SIZE));
            bins[GET_INDEX(x_bin_old, y_bin_old, num_bins)].erase(i);
            bins[GET_INDEX(x_bin_new, y_bin_new, num_bins)].insert(i);
        }
    }
}
