
// MODIFIED: Sequential LNS (Large Neighborhood Search) for TRP
// FIX: Converted from 5-Way Swarm to purely Sequential execution (Worker 0 only).
// FIX: Epochs are now tracked by K full iterations rather than time intervals.
// FIX: Removed multithreading locks; greedy acceptance naturally tracks the global best.
//
// [CHANGE v2] LOCAL SEARCH RESTRUCTURED:
//   OLD: a RESTRICTED 2-opt was woven INSIDE the reinsertion loop. After every single
//        node insertion, only 4 sweeps were evaluated, all ANCHORED at the insertion
//        point (best_idx-1 / best_idx), i.e. the anchor edge was FIXED.
//   NEW: no 2-opt during insertion at all. Once ALL removed nodes have been reinserted,
//        a FULL 2-opt is run over the completed tour. The inner sweep algebra is kept
//        EXACTLY the same as the restricted version; the only difference is that the
//        anchored edge (a,b) = (tour[i], tour[i+1]) is now ITERATED over every i,
//        instead of being pinned to the insertion point.
//   Everything else (datasets, parameters, parsing, metrics, logging, epochs, threading,
//   acceptance rule) is untouched.

#include <iostream>      // For standard console input/output
#include <vector>        // For dynamic array data structures
#include <cmath>         // For mathematical operations (sqrt, cos, acos)
#include <fstream>       // For file reading and writing (TSPLIB parsing, logging)
#include <string>        // For string manipulation
#include <sstream>       // For string stream processing (parsing lines)
#include <algorithm>     // For standard algorithms (shuffle, reverse, find)
#include <chrono>        // For time tracking and deadlines
#include <random>        // For the Mersenne Twister random number generator
#include <iomanip>       // For formatted output (setprecision, setw)
#include <limits>        // For maximum double values (INF_VAL)
#include <cstdlib>       // For standard library utilities
#include <functional>    // For std::function (if needed for callbacks)
#include <numeric>       // For numeric operations
#include <queue>         // For queue structures
#include <unordered_map> // For hash map structures
#include <deque>         // For double-ended queues
#include <thread>        // [MODIFIED] Added for the separated logger/monitor thread
#include <atomic>        // [MODIFIED] Added for thread-safe cross-thread atomic flags

using namespace std;     // Using standard namespace for convenience

// Define a large constant to represent infinity for cost comparisons
const double INF_VAL = numeric_limits<double>::max() / 2.0;
// Threshold for O(1) matrix usage (kept intact from original)
const int O1_MATRIX_THRESHOLD = 0;

// Structure to hold 2D coordinates of a node
struct Point { double x, y; };

// Function to write messages to the log file (removed mutex as it's sequential now)
void write_to_log(const string& message) {
    // Open the log file in append mode
    ofstream log_file("log_report_ROAR_regularDataset_modified_trivialCombination.txt", ios_base::app);
    if (log_file.is_open()) { // Check if file stream successfully opened
        log_file << message << "\n"; // Write the message and a newline
    } else {
        // Output error to standard error stream if logging fails
        cerr << "[ERROR] Unable to open log_report_ROAR_regularDataset_modified_trivialCombination.txt for writing!" << endl;
    }
}

// Structure representing the metadata of a TSP dataset
struct TSPInstance {
    string name;            // Name of the dataset (e.g., "berlin52")
    int nodes;              // Number of nodes/cities
    long long optimum;      // Known optimum (if applicable)
    bool depot_inclusive;   // Whether the tour must return to the depot
    bool floor_rounding;    // Whether distances should be floored
};

// Structure holding the parsed points and distance metric type
struct ParsedData {
    vector<Point> coords;     // The actual coordinate data
    string edge_weight_type;  // The type of distance metric to use
};

// --- TSPLIB DISTANCE METRICS (Kept exactly intact) ---

// Standard Euclidean distance, rounded to nearest integer
long long get_dist_euc_2d_nearest(Point p1, Point p2) {
    double dx = p1.x - p2.x; double dy = p1.y - p2.y; // Calculate differences
    return static_cast<long long>(std::sqrt(dx * dx + dy * dy) + 0.5); // Euclidean formula with rounding
}
// Euclidean distance, floored to integer
long long get_dist_euc_2d_floor(Point p1, Point p2) {
    double dx = p1.x - p2.x; double dy = p1.y - p2.y; // Calculate differences
    return static_cast<long long>(std::floor(std::sqrt(dx * dx + dy * dy))); // Euclidean formula with floor
}
// Euclidean distance, ceiling to integer
long long get_dist_ceil_2d(Point p1, Point p2) {
    double dx = p1.x - p2.x; double dy = p1.y - p2.y; // Calculate differences
    return static_cast<long long>(std::ceil(std::sqrt(dx * dx + dy * dy))); // Euclidean formula with ceil
}
// Pseudo-Euclidean distance (ATT metric from TSPLIB)
long long get_dist_att(Point p1, Point p2) {
    double dx = p1.x - p2.x; double dy = p1.y - p2.y; // Calculate differences
    double rij = std::sqrt((dx * dx + dy * dy) / 10.0); // Pseudo-Euclidean formula
    long long tij = static_cast<long long>(std::round(rij)); // Round to nearest long long
    return (tij < rij) ? tij + 1 : tij; // Specific ATT rounding rule
}
// Helper to convert geographic coordinates to radians
double convert_to_rad(double x) {
    double PI = 3.141592; int deg = static_cast<int>(x); // Extract degrees
    double min = x - deg; // Extract minutes
    return PI * (deg + 5.0 * min / 3.0) / 180.0; // Conversion formula
}
// Geographical distance metric (Great Circle distance)
long long get_dist_geo(Point p1, Point p2) {
    double RRR = 6378.388; // Radius of the Earth in TSPLIB standard
    double rad_x1 = convert_to_rad(p1.x); double rad_y1 = convert_to_rad(p1.y); // Convert p1
    double rad_x2 = convert_to_rad(p2.x); double rad_y2 = convert_to_rad(p2.y); // Convert p2
    double q1 = std::cos(rad_y1 - rad_y2); // Cosine of longitudinal difference
    double q2 = std::cos(rad_x1 - rad_x2); double q3 = std::cos(rad_x1 + rad_x2); // Latitudinal terms
    double val = 0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3); // Great circle formula term
    val = std::max(-1.0, std::min(1.0, val)); // Clamp value to [-1, 1] to avoid acos domain errors
    return static_cast<long long>(RRR * std::acos(val) + 1.0); // Final Geo distance
}

// --- ON-THE-FLY DISTANCE EVALUATOR ---
class DistanceEvaluator {
    const vector<Point>& coords; // Reference to the coordinate array
    string edge_type;            // Metric string
    bool floor_rounding;         // Floor preference
public:
    // Constructor initializing the evaluator constraints
    DistanceEvaluator(const vector<Point>& c, const string& et, bool fr)
        : coords(c), edge_type(et), floor_rounding(fr) {}

    // Operator overload to calculate distance between index i and j on the fly
    double operator()(int i, int j) const {
        if (i == j) return 0.0; // Distance to itself is 0
        // Route to the appropriate mathematical metric based on string tag
        if (edge_type == "GEO") return (double)get_dist_geo(coords[i], coords[j]);
        if (edge_type == "ATT") return (double)get_dist_att(coords[i], coords[j]);
        if (edge_type == "CEIL_2D") return (double)get_dist_ceil_2d(coords[i], coords[j]);
        // Default to EUC_2D, checking floor condition
        return floor_rounding ? (double)get_dist_euc_2d_floor(coords[i], coords[j]) : (double)get_dist_euc_2d_nearest(coords[i], coords[j]);
    }
    // Returns total number of nodes
    size_t size() const { return coords.size(); }
};

// --- CORE UTILS ---

// Verifies if a sequence of nodes is a valid Hamiltonian path starting at 0
bool is_hamiltonian_path(const vector<int>& tour, int expected_nodes) {
    if ((int)tour.size() != expected_nodes) return false; // Size mismatch
    vector<bool> visited(expected_nodes, false); // Tracking array
    for (int city : tour) {
        // Check bounds and duplicates
        if (city < 0 || city >= expected_nodes || visited[city]) return false;
        visited[city] = true; // Mark as visited
    }
    return tour[0] == 0; // Ensure it starts at the depot (node 0)
}

// Calculates the pure cumulative Traveling Repairman Problem (TRP) cost
double calculate_trp_cost(const DistanceEvaluator& dm, const vector<int>& tour, bool depot_inclusive) {
    if (tour.empty()) return INF_VAL; // Guard against empty tours
    double total_arrival_time = 0.0; // Accumulator for all arrivals (TRP objective)
    double current_time = 0.0;       // Tracking the current time along the route (TSP distance)
    for (size_t i = 1; i < tour.size(); ++i) { // Loop through the sequence
        current_time += dm(tour[i-1], tour[i]); // Add edge to current time
        total_arrival_time += current_time;     // Add current time to the cumulative TRP objective
    }
    if (depot_inclusive) { // If we must return to depot
        current_time += dm(tour.back(), tour[0]); // Add return edge
        total_arrival_time += current_time;       // Include return time in objective
    }
    return total_arrival_time; // Return final latency sum
}

// Struct for O(1) subsequence concatenation (Kept intact as requested)
struct Subseq {
    long long n; double d; double l; int first; int last;
};

// Function to combine two subsequences algebraically (Kept intact)
Subseq combine_subseq(const Subseq& a, const Subseq& b, const DistanceEvaluator& dm) {
    if (a.n == 0) return b; // Empty left side
    if (b.n == 0) return a; // Empty right side
    Subseq res;
    res.n = a.n + b.n; // Combined node count
    double edge = dm(a.last, b.first); // Edge connecting the two blocks
    res.d = a.d + edge + b.d; // Combined pure distance
    res.l = a.l + b.l + b.n * (a.d + edge); // Combined TRP latency using algebraic shift
    res.first = a.first; res.last = b.last; // Set new boundaries
    return res;
}

// Parses standard TSPLIB format files
ParsedData fetch_and_parse_tsplib(const string& filename) {
    ifstream file(filename); ParsedData data;
    data.edge_weight_type = "EUC_2D"; // Default metric
    string line; bool reading_nodes = false;
    if (!file.is_open()) return data; // Return empty if file not found
    while (getline(file, line)) { // Read line by line
        if (line.find("EOF") != string::npos) break; // End of file marker
        if (line.find("EDGE_WEIGHT_TYPE") != string::npos) { // Extract distance metric
            if (line.find("GEO") != string::npos) data.edge_weight_type = "GEO";
            else if (line.find("ATT") != string::npos) data.edge_weight_type = "ATT";
            else if (line.find("CEIL_2D") != string::npos) data.edge_weight_type = "CEIL_2D";
            else if (line.find("EUC_2D") != string::npos) data.edge_weight_type = "EUC_2D";
        }
        if (line.find("NODE_COORD_SECTION") != string::npos) { reading_nodes = true; continue; } // Trigger point reading
        if (reading_nodes) {
            stringstream ss(line); int id; double x, y;
            if (ss >> id >> x >> y) data.coords.push_back({x, y}); // Parse ID and coords
        }
    }
    return data;
}

// Determines the time limit allowed based on dataset size
double calculate_dynamic_limit(int n, double base_val) {
    if (n < 200) return base_val; // Small datasets get base time
    if (n >= 1000) return 10.0 * base_val; // Large get 10x
    double progress = static_cast<double>(n - 200) / 800.0; // Interpolation factor
    double multiplier = 1.0 + (progress * 9.0); // Scale between 1x and 10x
    return base_val * multiplier;
}

// [MODIFIED] Re-implemented time interval according to strict rules (1s to 30s linear interpolation)
double get_epoch_interval(int nodes) {
    // if (nodes <= 1000) return 1.0;
    if (nodes >= 100000) return 30.0;
    return 1.0 + 29.0 * ((static_cast<double>(nodes) - 1000.0) / 99000.0);
    // return 1;
}

// Determines hard timeout for the search based on node count
double get_search_timeout(int nodes) {
    if (nodes <= 1000) return 100.0; // 100s for up to 1000 nodes
    if (nodes >= 100000) return 1000.0; // Max cap 1000s
    return 100.0 + 900.0 * ((nodes - 1000.0) / 99000.0); // Linear interpolation
}

// Generates an initial draft route using the Nearest Neighbor heuristic
vector<int> generate_nn_draft_trp(const DistanceEvaluator& dm) {
    int N = dm.size(); if (N <= 1) return {}; // Guard for empty/trivial graphs
    vector<int> tour; tour.reserve(N); vector<bool> visited(N, false); // Allocation
    int curr = 0; tour.push_back(curr); visited[curr] = true; // Start at node 0
    for (int step = 1; step < N; ++step) { // For remaining N-1 steps
        int best_next = -1; double best_dist = INF_VAL; // Track best local move
        for (int i = 0; i < N; ++i) { // Scan all nodes
            if (!visited[i]) { // If unvisited
                double d = dm(curr, i); // Calculate distance
                if (d < best_dist) { best_dist = d; best_next = i; } // Update nearest
            }
        }
        curr = best_next; tour.push_back(curr); visited[curr] = true; // Move to nearest
    }
    return tour; // Return generated draft
}

// ========================================================================
// --- [NEW BLOCK] FULL TRP-AWARE 2-OPT (REPLACES THE RESTRICTED IN-LOOP 2-OPT) ---
// ========================================================================
//
// This function replaces the four restricted sweeps that used to live inside the
// reinsertion loop (i_pre / j_pre / i_suf / j_suf). Those four blocks all shared the
// SAME inner sweep algebra; their only role was to try the four ways of anchoring a
// reversal onto the freshly-inserted node.
//
// Here we keep that inner sweep algebra EXACTLY as it was, and simply promote the
// anchor to a loop variable:
//
//     for i = 0 .. p_size-3          <-- [CHANGED] anchored edge (a,b) is now ITERATED
//         for j = i+2 .. p_size-1    <-- unchanged sweep, incremental accumulators
//
// Correctness of reusing the same algebra: reversing tour[i+1 .. j] breaks the edges
// (a=tour[i], b=tour[i+1]) and (c=tour[j], d=tour[j+1]) and reconnects a-c and b-d.
//   * The (a,b) -> (a,c) length change delays every node from index i+1 onward
//     => multiplier K_b = p_size - i - 1 (+1 for the depot return leg).
//   * The (c,d) -> (b,d) length change delays every node from index j+1 onward
//     => multiplier K_c = p_size - j - 1 (+1 for the depot return leg).
//   * The internal latency of the segment is swapped from its forward accumulation
//     (total_latency_first) to its reversed accumulation (total_latency_reverse).
// Since K_b = (j - i) + K_c, these three terms together give the EXACT TRP delta,
// which is why the identical formula stays valid for an arbitrary anchor i.
//
// Search policy: BEST-IMPROVEMENT, SINGLE PASS. We loop over ALL O(n^2) possible 2-opt
// reversals exactly ONCE, record the single best improving reversal, and — if one exists —
// apply it and return. This is NOT iterated local search: we do NOT sweep again, and we do
// NOT continue until a 2-opt local optimum is reached. At most ONE reversal is applied per
// call. The destroy/recreate loop is what drives repeated improvement across iterations.
// Index 0 (the depot) can never move because the reversal always starts at i+1 >= 1.
//
// [FIX - EPOCH REPORTING] The incumbent cost is published to the monitor thread as soon as
// the reversal is applied, so a long O(n^2) sweep on a large instance does not leave the
// monitor observing a frozen cost (which previously fired a spurious EARLY STOP).
// This is safe (it never advertises a cost that will not be realised): we only publish when
// tour_cost < shared_best_cost <= current_cost == old_cost, and since the reversal strictly
// lowered the cost, the strictly-greedy acceptance test below is GUARANTEED to accept it.
// The 2-opt phase is only ever reached when the rebuilt tour is complete, so the failsafe
// rollback cannot discard an already-published value.
bool apply_full_two_opt_trp(
    const DistanceEvaluator& dm,
    vector<int>& tour,                       // [in/out] tour to be improved in place
    double& tour_cost,                       // [in/out] TRP cost of 'tour', kept in sync
    bool depot_inclusive,
    const std::atomic<bool>& stop_flag,      // honour the hard timeout / early stop
    std::atomic<double>& shared_best_cost)   // live channel to the monitor thread
{
    int p_size = static_cast<int>(tour.size());
    if (p_size < 4) return false;            // No reversal is possible on a tiny tour

    // [CHANGED] SINGLE PASS ONLY. The enclosing 'while (!stop_flag)' loop that used to repeat
    // the sweep until no improving move remained has been REMOVED, per the requirement that a
    // full 2-opt loops over the n^2 moves exactly once and applies only the single best move.

    double original_baseline = tour_cost;    // Absolute cost of the CURRENT tour
    double current_best = tour_cost;         // Best cost seen during this single sweep
    int best_s = -1, best_e = -1;            // Reversal bounds of the best move found
    bool swapped = false;                    // Was any improving move found this sweep?

    // [KEPT INTACT] Same lambda as the restricted version: turns the three algebraic
    // terms into an absolute candidate cost and records it if strictly better.
    auto check_and_record = [&](int i, int j, double first_term, double second_term, double tl_rev, double tl_first) {
        double test_cost = original_baseline + first_term + second_term + tl_rev - tl_first;
        if (test_cost < current_best - 1e-5) { // Strict improvement threshold
            current_best = test_cost; best_s = i + 1; best_e = j; swapped = true;
        }
    };

    // [CHANGED] ANCHOR LOOP. Previously 'i' was fixed to best_idx-1 / best_idx (restricted).
    // Now every edge (tour[i], tour[i+1]) in the tour is used as the anchor in turn, giving
    // the full O(n^2) set of 2-opt moves in one pass.
    for (int i = 0; i + 2 < p_size; ++i) {

        // Bail out mid-sweep if the monitor thread signalled the deadline
        if (stop_flag.load(std::memory_order_relaxed)) break;

        // [KEPT INTACT] Per-anchor accumulators, reset for each new anchor i.
        double cumulative_length_tsp = 0.0;   // Running forward length inside the segment
        double total_latency_reverse = 0.0;   // Internal latency if the segment is reversed
        double total_latency_first   = 0.0;   // Internal latency in the current orientation

        int K_b = p_size - i - 1 + (depot_inclusive ? 1 : 0); // Downstream shift multiplier for edge (a,b)
        int a = tour[i], b = tour[i + 1];                     // The anchored edge being broken
        double d_ab = dm(a, b);                               // Hoisted: constant across the j sweep

        // [KEPT INTACT] Inner sweep: extend the reversed segment one node at a time,
        // updating the latency accumulators incrementally so each (i,j) costs O(1).
        for (int j = i + 2; j < p_size; ++j) {
            double dist = dm(tour[j - 1], tour[j]);           // Edge entering node j
            cumulative_length_tsp += dist;                    // Accumulate path distance
            total_latency_reverse += dist * (j - i - 1);      // Latency of reversed segment
            total_latency_first   += cumulative_length_tsp;   // Latency of forward segment

            int K_c = p_size - j - 1 + (depot_inclusive ? 1 : 0); // Secondary shift multiplier
            int c = tour[j];                                  // Tail of the segment
            int d = (j + 1 < p_size) ? tour[j + 1]            // Successor of the segment tail,
                                     : (depot_inclusive ? tour[0] : -1); // or the depot / nothing

            double f_term = (dm(a, c) - d_ab) * K_b;                       // First edge swap delta
            double s_term = (d != -1) ? (dm(b, d) - dm(c, d)) * K_c : 0.0; // Second edge swap delta
            check_and_record(i, j, f_term, s_term, total_latency_reverse, total_latency_first);
        }
    }

    // [CHANGED] No improving reversal found in the single pass => leave the tour unchanged.
    if (!swapped) return false;

    // Apply the SINGLE best reversal found in this one pass, then return (no re-sweep).
    reverse(tour.begin() + best_s, tour.begin() + best_e + 1);

    // [NOTE] The algebra already yields the exact new cost (current_best), but we recompute
    // once here (O(n), negligible next to the O(n^2) sweep) to avoid any floating-point drift.
    tour_cost = calculate_trp_cost(dm, tour, depot_inclusive);

    // [FIX - EPOCH REPORTING] Publish the improvement immediately so the monitor thread does
    // not see a frozen cost during the long sweep. Single worker => no race on the store.
    double global_best = shared_best_cost.load(std::memory_order_relaxed);
    if (tour_cost < global_best - 1e-5) {
        shared_best_cost.store(tour_cost, std::memory_order_relaxed);
    }

    return true;
}

// ========================================================================
// --- SEQUENTIAL LNS (LARGE NEIGHBORHOOD SEARCH) METAHEURISTIC ---
// ========================================================================

// Replaces perform_parallel_hybrid_search. Strictly executes Worker 0 logic.
vector<int> perform_sequential_lns_search(
    const DistanceEvaluator& dm, const vector<int>& seed_tour, double k_prop,
    chrono::steady_clock::time_point deadline, bool depot_inclusive,
    int epoch_patience, int loops_per_epoch)
{
    int N = dm.size(); // Total nodes
    if (N <= 3 || seed_tour.empty()) return seed_tour; // Return immediately if trivial

    // Calculate how many nodes to rip out per LNS cycle based on k_prop
    int actual_k = max(1, static_cast<int>(k_prop * (N - 1)));
    actual_k = min({actual_k, N - 2, 300}); // Bound the removal limit

    // Random Number Generator seeded for determinism/variance
    mt19937 rng(1337 + chrono::steady_clock::now().time_since_epoch().count());

    // Because it strictly accepts improvements, current_tour IS the global best.
    vector<int> current_tour = seed_tour;
    // Calculate initial cost
    double current_cost = calculate_trp_cost(dm, seed_tour, depot_inclusive);

    // [MODIFIED] Cross-thread synchronization variables to replace inline timing checks
    // stop_flag replaces chrono::now() > deadline in the hot loops.
    // shared_best_cost allows the monitor thread to read the current best cost asynchronously.
    std::atomic<bool> stop_flag{false};
    std::atomic<double> shared_best_cost{current_cost};

    // Epoch tracking variables
    int loop_counter = 0; // Tracks the inner loops to trigger K-based epochs [Kept variable intact per constraint]
    double epoch_time_interval = get_epoch_interval(N);
    auto start_time = chrono::steady_clock::now();

    // [MODIFIED] Generate a separated monitor thread.
    // This thread strictly handles sleeping for the epoch interval, reporting to the log,
    // and managing the hard timeout, eliminating all such timing checks from the main loop.
    std::thread monitor_thread([&]() {
        int shared_epoch = 0;
        int non_improving_epochs = 0;
        double last_epoch_cost = shared_best_cost.load(std::memory_order_relaxed);

        while (!stop_flag.load(std::memory_order_relaxed)) {
            auto now = chrono::steady_clock::now();
            auto next_epoch_time = now + chrono::duration<double>(epoch_time_interval);

            // [MODIFIED] Check if the next epoch is further than the hard search timeout
            if (next_epoch_time >= deadline) {
                // Wake up exactly at the search timeout to trigger the stop of all threads
                std::this_thread::sleep_until(deadline);
                stop_flag.store(true, std::memory_order_relaxed);
                break;
            } else {
                // Sleep for epoch_interval time
                std::this_thread::sleep_until(next_epoch_time);

                // Safety exit if the main worker thread finished unexpectedly
                if (stop_flag.load(std::memory_order_relaxed)) break;

                // Wake up, access and report to log the best cost
                double observed_cost = shared_best_cost.load(std::memory_order_relaxed);
                double elapsed_sec = chrono::duration<double>(chrono::steady_clock::now() - start_time).count();

                // Console output
                cout << "\n    [Sync Epoch " << shared_epoch + 1
                     << " | Time: " << fixed << setprecision(2) << elapsed_sec << "s"
                     << "] Best Cost: " << fixed << setprecision(2) << observed_cost << flush;

                // File output
                write_to_log("Sync Epoch " + to_string(shared_epoch + 1) +
                             " | Best Cost: " + to_string(observed_cost));

                // Check for strict improvement against the last epoch baseline
                if (observed_cost < last_epoch_cost - 1e-5) {
                    last_epoch_cost = observed_cost; // Update baseline
                    non_improving_epochs = 0;        // Reset patience counter
                } else {
                    non_improving_epochs++; // Increment patience counter
                    // Early stopping mechanism
                    if (non_improving_epochs >= epoch_patience) {
                        cout << " [EARLY STOP: " << epoch_patience << " epochs without improvement]" << flush;
                        write_to_log("Terminated due to " + to_string(epoch_patience) + " non-improving epochs.");
                        stop_flag.store(true, std::memory_order_relaxed); // Trigger the stop of the worker thread
                        break;
                    }
                }
                shared_epoch++; // Increment epoch ID
            }
        }
    });

    // ================== WORKER 0: LNS EXECUTION ==================
    // [MODIFIED] The discrete epoch timing block was removed from here and shifted to the thread above.

    // [MODIFIED] Main execution loop: run until the monitor thread triggers stop_flag
    while (!stop_flag.load(std::memory_order_relaxed)) {

        // Save the baseline state before applying mutations
        vector<int> old_tour = current_tour;
        double old_cost = current_cost;

        // Pick a random number of nodes to remove (between 1 and actual_k)
        int curr_k = uniform_int_distribution<int>(1, actual_k)(rng);

        // Generate a list of valid candidate indices to remove (ignoring index 0 / depot)
        vector<int> candidates; candidates.reserve(N - 1);
        for (int i = 1; i < N; ++i) candidates.push_back(i);
        shuffle(candidates.begin(), candidates.end(), rng); // Randomize candidate indices

        // Track which specific node IDs are being ripped out
        vector<bool> is_removed(N, false);
        vector<int> removed_nodes;
        for (int i = 0; i < curr_k; ++i) {
            int node = old_tour[candidates[i]]; // Identify the node
            is_removed[node] = true;            // Mark it as removed
            removed_nodes.push_back(node);      // Add to removed list
        }

        // Build the partial tour consisting only of nodes that were NOT removed
        vector<int> partial_tour; partial_tour.reserve(N);
        for (int city : old_tour) {
            if (!is_removed[city]) partial_tour.push_back(city);
        }

        // ---------------- RECREATE PHASE: PURE GREEDY REINSERTION ----------------
        // [CHANGED] This loop is now ONLY greedy insertion. Everything that used to follow
        // an insertion (the calculate_trp_cost baseline refresh, the check_and_record lambda,
        // and the four restricted 2-opt sweeps anchored at the insertion point) has been
        // REMOVED from here and replaced by a single full 2-opt run after the loop ends.
        // Consequence: reinsertion is now much cheaper per node (no O(N) cost recomputation
        // and no local sweeps per insertion), and the local search sees the finished tour.
        for (int node : removed_nodes) {
            // [MODIFIED] Changed absolute time check to atomic flag check to prevent overrunning
            if (stop_flag.load(std::memory_order_relaxed)) break;

            double best_increase = INF_VAL; // Track the best marginal insertion cost
            int best_idx = -1;              // Track the best index to insert the node
            int p_size = partial_tour.size(); // Current size of the partial tour
            double cum_time = 0.0;            // Tracks distance elapsed up to insertion point

            // Find the best position to insert the current 'node'
            for (int i = 1; i <= p_size; ++i) {
                int prev = partial_tour[i-1]; // Node before insertion
                int curr = (i < p_size) ? partial_tour[i] : (depot_inclusive ? 0 : -1); // Node after insertion

                double increase = 0.0; // The calculated TRP latency increase for this slot
                if (curr != -1) {
                    // Triangle inequality change (dist added to routing)
                    double delta = dm(prev, node) + dm(node, curr) - dm(prev, curr);
                    // TRP Shift: Immediate time + cascading delay applied to all remaining nodes
                    increase = cum_time + dm(prev, node) + (p_size - i + (depot_inclusive ? 1 : 0)) * delta;
                    cum_time += dm(prev, curr); // Update cum_time for the next iteration
                } else {
                    // Inserting at the very end with no depot return
                    increase = cum_time + dm(prev, node);
                }

                // If this is the best slot found so far, record it
                if (increase < best_increase) {
                    best_increase = increase;
                    best_idx = i;
                }
            }
            if (best_idx == -1) break; // Failsafe

            // Execute the insertion at the optimally found index
            partial_tour.insert(partial_tour.begin() + best_idx, node);

            // [REMOVED] p_size++ / local_best_dist / original_baseline / current_best /
            // [REMOVED] swapped / best_s / best_e / check_and_record lambda
            // [REMOVED] restricted 2-opt blocks 1..4 (i_pre, j_pre, i_suf, j_suf)
            // [REMOVED] the conditional reverse() that applied the restricted move
            // They are all superseded by the full 2-opt executed once, below.
        }
        // -------------- END OF RECREATE PHASE --------------

        // After attempting all removals and reinsertions, verify the tour is complete
        if (partial_tour.size() == (size_t)N) {

            double new_cost = calculate_trp_cost(dm, partial_tour, depot_inclusive); // Final cost check

            // [NEW] FULL 2-OPT IMPROVEMENT PHASE — SINGLE PASS.
            // Runs ONCE over the COMPLETE rebuilt tour: it sweeps every one of the O(n^2)
            // 2-opt moves a single time, then applies only the single best improving reversal
            // (or none). It does NOT iterate to a 2-opt local optimum. It refines
            // 'partial_tour' in place and keeps 'new_cost' in sync, so the acceptance test
            // below is unchanged. 'shared_best_cost' is handed in so the applied move is
            // reported to the monitor thread immediately rather than only at acceptance.
            apply_full_two_opt_trp(dm, partial_tour, new_cost, depot_inclusive, stop_flag, shared_best_cost);

            // STRICTLY GREEDY ACCEPTANCE: Only accept if strictly better
            if (new_cost < old_cost - 1e-5) {
                current_tour = partial_tour;
                current_cost = new_cost;
                // [MODIFIED] Store updated best cost safely for the monitor thread to read
                // (Still required: the 2-opt only publishes when IT improves the tour. An LNS
                // iteration whose gain came purely from the ruin-and-recreate step, with no
                // 2-opt reversal applied, would otherwise never be reported.)
                shared_best_cost.store(current_cost, std::memory_order_relaxed);
            } else {
                // Reject mutation, reset to previous best
                current_tour = old_tour;
                current_cost = old_cost;
            }
        } else {
            // Failsafe rollback if tour rebuilding failed
            current_tour = old_tour;
            current_cost = old_cost;
        }
        loop_counter++; // Increment the main loop counter for the epoch system
    }

    // [MODIFIED] Ensure synchronization and complete shutdown of the monitor thread before exiting
    stop_flag.store(true, std::memory_order_relaxed);
    if (monitor_thread.joinable()) {
        monitor_thread.join();
    }
    write_to_log("Total LNS Loops Executed: " + to_string(loop_counter));
    return current_tour; // Return the best sequence found
}

// --- MAIN BENCHMARK LOOP ---

int main() {
    // Pipeline parameters
    double p_sec = 1.0;          // Base timeout scaling factor
    double k_prop = 0.2;         // Default LNS removal proportion
    int epoch_patience = 10;      // Default epochs without improvement to trigger early stop
    int loops_per_epoch = 100;   // Default K loops to represent an epoch

    // UI Greetings
    cout << "=== TRP Scalability Benchmark (ROAR) ===" << endl;

    string in;
    // Parameter Inputs
    cout << "LNS Max Removal Nodes (Prop) [Default 0.2]: "; getline(cin, in); if(!in.empty()) k_prop = stod(in);
    cout << "Epoch Patience for Early Stop [Default 10]: "; getline(cin, in); if(!in.empty()) epoch_patience = stoi(in);
    cout << "Loops per Epoch (K) [Default 100]: "; getline(cin, in); if(!in.empty()) loops_per_epoch = stoi(in);
    cout << "Auto-advance to the next dataset without waiting? (y/n) [Default n]: ";

    // Write configuration to log
    write_to_log("=== TRP Scalability Benchmark (ROAR) ===");
    write_to_log("LNS Max Removal Nodes (Prop): " + to_string(k_prop));
    write_to_log("Epoch Patience for Early Stop: " + to_string(epoch_patience));
    write_to_log("Loops per Epoch (K): " + to_string(loops_per_epoch));

    // Read auto-advance choice
    string auto_input; getline(cin, auto_input);
    bool auto_advance = (auto_input == "y" || auto_input == "Y");

    // Array of TSPLIB datasets to process
    vector<TSPInstance> datasets = {
        {"berlin52", 52, 1000, true, false},
        {"st70", 70, 1000, true, false},
        {"kroA100", 100, 1000, true, false},
        {"pr107", 107, 1000, true, false},
        {"ch150", 150, 1000, true, false},
        {"d198", 198, 1000, true, false},
        {"ts225", 225, 1000, true, false},
        {"a280", 280, 1000, true, false},
        {"pcb442", 442, 1000, true, false},
        {"att532", 532, 1000, true, false},
        {"rat783", 783, 1000, true, false},
        {"dsj1000", 1000, 1000, true, false},
        {"u1432", 1432, 1000, true, false},
        {"d1655", 1655, 1000, true, false},
        {"pr2392", 2392, 1000, true, false},
        // All datasets with > 3000 nodes included below
        {"pcb3038", 3038, 1000, true, false},
        {"fnl4461", 4461, 1000, true, false},
        {"rl5934", 5934, 1000, true, false},
        {"pla7397", 7397, 1000, true, false},
        {"rl11849", 11849, 1000, true, false}
    };

    int number_execution = 10; // Number of times to loop the entire benchmark suite
    for (int i = 0; i < number_execution; i++) {
        write_to_log("=== EXECUTION " + to_string(i+1) + " ===");

        // Open CSV to store execution results
        ofstream csv_file("trp_results_ROAR_regularDataset_modified_trivialCombination.csv", ios_base::app);
        csv_file << "Dataset,Nodes,Method,RawCost,FinalCost,Time(s)\n";

        const int num_methods = 1;
        string methods[] = { "ROAR_trivialCombination" };

        // Loop through each dataset
        for (const auto& ds : datasets) {
            // Verify file exists
            if (!ifstream(ds.name + ".tsp").good()) {
                cout << "[WARNING] Dataset " << ds.name << ".tsp not found. Skipping." << endl;
                continue;
            }
            // Parse dataset
            ParsedData parsed = fetch_and_parse_tsplib(ds.name + ".tsp");
            if (parsed.coords.empty() || parsed.coords.size() < 2) {
                cout << "[WARNING] Parsing failed for: " << ds.name << ". Skipping." << endl;
                continue;
            }

            // Establish dynamic time limit based on node count
            double dynamic_no_improve = calculate_dynamic_limit(ds.nodes, p_sec);
            DistanceEvaluator dm(parsed.coords, parsed.edge_weight_type, ds.floor_rounding);

            // Output vectors for results
            vector<long long> raw_results(num_methods, -1);
            vector<long long> final_results(num_methods, -1);
            vector<double> times(num_methods, 0.0);

            cout << string(105, '=') << "\nPROCESSING: " << ds.name << " (" << ds.nodes << " nodes)" << endl;

            for (int m = 0; m < num_methods; ++m) {
                cout << "  [Running] " << left << setw(15) << methods[m] << "... " << flush;
                write_to_log("=== Processing Dataset: " + ds.name + " | Method: " + methods[m] + " ===");

                // Track execution time
                auto s1 = chrono::steady_clock::now();
                double search_timeout = get_search_timeout(ds.nodes);

                // [MODIFIED] final_dl now solely relies on get_search_timeout, get_dl lambda is completely removed.
                auto final_dl = chrono::steady_clock::now() + chrono::milliseconds(static_cast<long long>(search_timeout * 1000.0));

                // Generate Initial draft tour using Nearest Neighbor
                vector<int> seed_tour = generate_nn_draft_trp(dm);

                // Verify the initial draft
                if (!seed_tour.empty() && is_hamiltonian_path(seed_tour, ds.nodes)) {
                    raw_results[m] = static_cast<long long>(calculate_trp_cost(dm, seed_tour, ds.depot_inclusive));
                    cout << "[Draft Cost: " << raw_results[m] << "] -> " << flush;
                } else {
                    cout << "[Draft Cost: INVALID] -> " << flush;
                }

                // If draft generation succeeded, run the Sequential LNS
                if (!seed_tour.empty()) {
                    auto final_tour = perform_sequential_lns_search(dm, seed_tour, k_prop, final_dl, ds.depot_inclusive, epoch_patience, loops_per_epoch);

                    // Verify the final refined tour
                    if (is_hamiltonian_path(final_tour, ds.nodes)) {
                        final_results[m] = static_cast<long long>(calculate_trp_cost(dm, final_tour, ds.depot_inclusive));
                        cout << "\n    -> [Refined Cost: " << final_results[m] << "] ";

                        // Output final route to log
                        stringstream route_ss;
                        route_ss << "[ROUTE FOUND] Dataset: " << ds.name << " | Method: " << methods[m] << "\n";
                        for (size_t j = 0; j < final_tour.size(); ++j) {
                            route_ss << final_tour[j] << (j == final_tour.size() - 1 ? "" : " ");
                        }
                        write_to_log(route_ss.str());
                    } else {
                        // Debug for failed structural integrity
                        cout << "\n    [DEBUG] Tour is invalid: Size=" << final_tour.size() << " Expected=" << ds.nodes << endl;
                        cout << "    -> [Refined Cost: INVALID] ";
                    }
                }

                // Finalize time taken
                times[m] = chrono::duration<double>(chrono::steady_clock::now() - s1).count();
                cout << "Done in " << fixed << setprecision(2) << times[m] << "s\n";
            }

            // Print neatly formatted table for the dataset result
            cout << "\n--- RESULTS FOR " << ds.name  << endl;
            cout << left << setw(18) << "Method" << setw(30) << "Raw Cost" << setw(30) << "Final Cost" << setw(15) << "Time(s)" << endl;
            cout << string(105, '-') << endl;

            for(int j = 0; j < num_methods; j++) {
                stringstream row_ss;
                if (final_results[j] == -1 || raw_results[j] == -1) { // If failure
                    row_ss << left << setw(18) << methods[j] << "INVALID RUN\n";
                    csv_file << ds.name << "," << ds.nodes << "," << methods[j] << ",INVALID,INVALID," << times[j] << "\n";
                } else { // If successful
                    row_ss << left << setw(18) << methods[j]
                        << setw(30) << raw_results[j]
                        << setw(30) << final_results[j]
                        << setprecision(3) << times[j] << "s\n";
                    // Append data to CSV
                    csv_file << ds.name << "," << ds.nodes << "," << methods[j]
                            << "," << raw_results[j] << "," << final_results[j] << "," << times[j] << "\n";
                }
                cout << row_ss.str();
                write_to_log(row_ss.str()); // Log the table row
            }
            // Pause if auto advance is disabled
            if (!auto_advance) {
                cout << "\nPress [Enter] to continue..." << flush;
                string dummy; getline(cin, dummy);
            }
        }
        csv_file.close(); // Close CSV properly
    }
    return 0; // End program cleanly
}
