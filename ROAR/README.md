# ROAR — Restricted 2-Opt And Ruin-and-recreate
Reference implementation of **ROAR**, the CPU-native ruin-and-recreate heuristic for the
Travelling Repairman Problem (TRP / Minimum Latency Problem) introduced in:



This folder holds **original code** written for the paper — it is not adapted from any third-party
codebase. The only external material is the instance data (see [Instances](#instances)).
IMPORTANT: The programs contain multiple artefacts and comment from the initial draft versions. Ignore dead code (not executed in main) and strange comment.
---

## 1. Contents

### Source files

| File | Paper artefact | Instance set driven |
|---|---|---|
| `main_roar.cpp` | **ROAR** (Algorithms 1 & 2) | *Regular* set, 52 – 11,849 nodes (Table I) |
| `main_ROAR_large.cpp` | **ROAR** | *Large* set, 13,509 – 33,810 nodes (Table II) |
| `main_roar_extreme.cpp` | **ROAR** | *Extreme* set, 32,892 – 200,000 nodes (Table III) |
| `main_trivialCombination.cpp` | **LNS+2-opt**, the post-repair-local-search version of Section V-F | *Regular* set |
| `main_trivialCombination_large.cpp` | **LNS+2-opt**, the post-repair-local-search version of Section V-F | *Large* set |

> `main_roar.cpp`, `main_ROAR_large.cpp` and `main_roar_extreme.cpp` are **byte-for-byte identical
> except for three things**: the log filename, the results-CSV filename, and the `datasets` array in
> `main()`. There is intentionally no shared header — each is a self-contained translation unit so a
> run can be launched with a single `g++` invocation.

### Build / run helpers

| File | Purpose |
|---|---|
| `run_all.sh` | Compiles and runs the three ROAR binaries back-to-back with default parameters |
| `.vscode/tasks.json` | VS Code `g++ -g -O3` build task for the active file |

### Instances

37 TSPLIB-format `.tsp` files live flat in this directory because every binary opens
`<dataset>.tsp` relative to its own working directory.

* **Regular (20)** — `berlin52 st70 kroA100 pr107 ch150 d198 ts225 a280 pcb442 att532 rat783
  dsj1000 u1432 d1655 pr2392 pcb3038 fnl4461 rl5934 pla7397 rl11849`
* **Large (6)** — `usa13509 fma21553 lsb22777 bbz25234 boa28924 pla33810`
* **Extreme (11)** — `xib32892 pba38478 ics39603 dan59296 pla85900 mona-lisa100K vangogh120K
  venus140K pareja160K courbet180K earring200K`

Sources: TSPLIB (Reinelt), the Waterloo VLSI collection (Rohe), and the Waterloo TSP-art
collection (Bosch / Kaplan & Bosch). See the paper's references [16]–[20].

### Result artefacts (produced by a run, committed for reproducibility)

| File | Content |
|---|---|
| `trp_results_ROAR_<set>Dataset_modified.csv` | One row per (dataset, run): `Dataset,Nodes,Method,RawCost,FinalCost,Time(s)` |
| `log_report_ROAR_<set>Dataset_modified.txt` | Per-epoch incumbent trace, final routes, `Total LNS Loops Executed` |
| `*_trivialCombination.*` | Same two artefacts for the ROAR-PLS ablation |

`RawCost` is the latency of the nearest-neighbour seed (the *Draft (NN)* column of Table III);
`FinalCost` is the latency ROAR returns.

### Binaries

`main_roar`, `main_ROAR_large`, `main_roar_extreme`, `main_trivialCombination` (VS Code task output)
and `run_regular`, `run_large`, `run_extreme` (`run_all.sh` output) are checked-in Linux x86-64
executables. Rebuild them rather than trusting them.

---

## 2. Building and running

### Requirements

C++11 compiler with `<thread>` (GCC 7+ / Clang 6+) and pthreads. No external libraries.
Peak RSS is `O(N)` — a 200,000-node run stays well under 100 MB.

### One-shot: all three instance sets

```bash
cd ROAR
chmod +x run_all.sh
./run_all.sh
```

`run_all.sh` compiles with `g++ -O3 -pthread` and pipes `"\n\n\ny\n"` into each binary, i.e. accepts
all three default hyperparameters and enables auto-advance.

### Manual

```bash
g++ -O3 -pthread main_roar.cpp             -o run_regular
g++ -O3 -pthread main_ROAR_large.cpp       -o run_large
g++ -O3 -pthread main_roar_extreme.cpp     -o run_extreme
g++ -O3 -pthread main_trivialCombination.cpp -o run_pls
g++ -O3 -pthread main_trivialCombination_large.cpp -o run_pls_large


**The outer loop runs the whole instance list 10 times** (`number_execution = 10` in `main()`), which
is where the paper's "mean over 10 independent runs" comes from. Results are *appended* to the CSV,
so delete or move the CSV before a fresh benchmark.

Wall-clock cost of one full pass: roughly `sum over instances of T(N)` — about 30 min for the regular
set, 32 min for the large set, 2.2 h for the extreme set. Multiply by 10 for the full benchmark.

---

## 3. Algorithm walkthrough (`main_roar.cpp`)

Line references are to `main_roar.cpp`; the other three files are structurally identical.

### 3.1 Objective — depot-inclusive latency

`calculate_trp_cost()` (L145–158) accumulates the arrival time at every node **and** the arrival back
at the depot:

```
Z = sum_{k=1}^{N-1} A_k  +  (A_{N-1} + c(pi_{N-1}, pi_0))
```

This is the **Depot-Inclusive** variant of Section III, the same metric used by Silva et al. and by
the GILS-RVND baseline in `../GILS_RVND`. All three algorithms in this repository are scored under it.

### 3.2 Distance metrics (L66–128)

`DistanceEvaluator` computes edge costs **on the fly** — there is deliberately no `O(N^2)` matrix,
which is what lets ROAR run at 200,000 nodes where GILS-RVND goes OOM.

| `EDGE_WEIGHT_TYPE` | Rule | Instances affected |
|---|---|---|
| `EUC_2D` | `(long long)(sqrt(dx^2+dy^2) + 0.5)` — nearest integer | most |
| `CEIL_2D` | `ceil(...)` | `dsj1000` (reported as `dsj1000ceil` in Table I) |
| `ATT` | TSPLIB pseudo-Euclidean | `att532` |
| `GEO` | TSPLIB great-circle, with the `acos` argument clamped to `[-1,1]` | none in these sets |

The metric is read from the file header by `fetch_and_parse_tsplib()` (L179–199). These rules match
`GILS_RVND/MLP-master/src/readData.cpp` exactly, so costs are directly comparable.

### 3.3 Initialization (L226–241)

`generate_nn_draft_trp()` — classical nearest neighbour from the depot, `O(N^2)` scan.
Deterministic: ties break to the lowest index. This is Section IV-A, and it is also the
*Draft (NN)* baseline of Table III.

### 3.4 The ROAR cycle (L340–534)

**Removal budget** (L257–258) — Eq. (2):

```cpp
int actual_k = max(1, (int)(k_prop * (N - 1)));
actual_k = min({actual_k, N - 2, 300});          // alpha_destroy, N-2, C_limit
```

`C_limit = 300` is hard-coded, matching Section V-D.

**Ruin** (L347–367) — draw `N_remove ~ U{1, ..., actual_k}`, shuffle positions `1..N-1` (the depot at
position 0 is never removable) and excise the first `N_remove`.

**Recreate** (L370–401) — cheapest insertion, a literal transcription of **Algorithm 1**:

```cpp
for (int i = 1; i <= p_size; ++i) {
    prev = partial_tour[i-1];
    curr = (i < p_size) ? partial_tour[i] : 0;              // b <- pi'_{i mod m}: depot when i == m
    delta = dm(prev,node) + dm(node,curr) - dm(prev,curr);  // Eq. (3)
    increase = cum_time + dm(prev,node) + (p_size - i + 1) * delta;  // Eq. (4)
    cum_time += dm(prev,curr);                              // T <- T + c(a,b)
}
```

`cum_time` is `T_{i-1}`, carried incrementally, so each candidate slot is `O(1)` and a full scan is
`O(|pi'|)`.

**Restricted 2-opt** (L429–508) — runs immediately after *every* insertion (the ILO scheme of
Section IV-C.3). Four contiguous families, each anchoring one fracture edge on an edge incident to
the freshly inserted node at index `best_idx`:

| Block | Anchor | Broken edge held fixed | Sweep direction |
|---|---|---|---|
| 1 (`i_pre = best_idx-1`) | entry edge as `e1` | `(pi[best_idx-1], pi[best_idx])` | `j` ascends |
| 2 (`j_pre = best_idx-1`) | entry edge as `e2` | `(pi[best_idx-1], pi[best_idx])` | `i` descends |
| 3 (`i_suf = best_idx`) | exit edge as `e1` | `(pi[best_idx], pi[best_idx+1])` | `j` ascends |
| 4 (`j_suf = best_idx`) | exit edge as `e2` | `(pi[best_idx], pi[best_idx+1])` | `i` descends |

All four evaluate Eq. (5):

```
Delta_Z = delta_entry * K_b  +  delta_exit * K_c  +  TRP_rev  -  TRP_first
```

with `K_b = p_size - i - 1 + 1` and `K_c = p_size - j - 1 + 1` (the trailing `+1` is the depot-return
leg). The forward blocks use the ascending recurrence of Eqs. (6)–(8) verbatim:

```cpp
cumulative_length_tsp += dist;                       // L_tsp   (Eq. 6)
total_latency_reverse += dist * (j - i_pre - 1);     // TRP_rev (Eq. 7)
total_latency_first   += cumulative_length_tsp;      // TRP_1st (Eq. 8)
```

The backward blocks use the mirror recurrence noted in the paper as derivable "with minor
modifications" — as `i` decrements, `cum_rev_bw` accumulates `sum_{t=i+2}^{j} e_t` and feeds
`TRP_rev`, while `TRP_first` gains `e_{i+2} * (j - i - 1)`. Every candidate is `O(1)` time with
`O(1)` auxiliary state, which is the paper's Section IV-D contribution.

The single best improving reversal across all four families is applied (L511); if none improves, the
tour is left alone.

**Acceptance** (L515–532) — strictly greedy: `new_cost < old_cost - 1e-5` or roll back to
`old_tour`. Because acceptance is monotone, `current_tour` *is* the incumbent best; there is no
separate best-so-far copy.

### 3.5 Termination (L282–334, L211–223)

A dedicated **monitor thread** owns all timing; the worker loop only polls an atomic `stop_flag`, so
no clock is read in the hot path.

* **Epoch length** `E(N)` — Eq. (9):
  `E(N) = 1 + 29*(N-1000)/99000` for `N <= 100,000`, else `30` seconds.
* **Wall-clock budget** `T(N)` — Eq. (10):
  `T(N) = 100` for `N <= 1000`; `100 + 900*(N-1000)/99000` (= `100 + (N-1000)/110`) up to
  `N = 100,000`; `1000` beyond.
* **Patience** — the run halts after `epoch_patience = 10` consecutive epochs without a strict
  improvement, or when `T(N)` elapses, whichever comes first.

`get_search_timeout()` is identical in `GILS_RVND/MLP-master/src/main.cpp` and in
`PDLSH/source/benchmark.py`, so all three methods share the same budget.

---

## 4. LNS+2-opt (`main_trivialCombination.cpp`)

The Section V-F ablation. Identical to `main_roar.cpp` except:

* the four restricted 2-opt blocks are **removed from the reinsertion loop** — recreate is now pure
  greedy insertion;
* after the tour is complete, `apply_full_two_opt_trp()` (L296–379) sweeps the **full `O(n^2)`**
  2-opt neighbourhood once, using the *same* delta algebra with the anchor `i` promoted to a loop
  variable, and applies the single best improving reversal.

This isolates the effect of *where* the local search sits (interleaved ILO vs. post-repair) while
holding the ruin-and-recreate frame constant.


## 5. Reproducing the paper tables
Note: due to the large size, the log report for large and extreme dataset is omitted. The result files are still available.
| Table | Command | Then |
|---|---|---|
| I (regular) | `printf "\n\n\ny\n" \| ./run_regular` | copy `trp_results_ROAR_regularDataset_modified.csv` to `../processingData/` |
| II (large) | `printf "\n\n\ny\n" \| ./run_large` | copy `trp_results_ROAR_largeDataset_modified.csv` to `../processingData/` |
| III (extreme) | `printf "\n\n\ny\n" \| ./run_extreme` | copy `trp_results_ROAR_extremeDataset_modified.csv` to `../processingData/` |
| V-F ablation | `printf "\n\n\ny\n" \| ./run_pls` | copy `trp_results_ROAR_regularDataset_modified_trivialCombination.csv` to `../processingData/` |

Aggregation, gap computation and the loop-count summary are handled by the scripts in
[`../processingData`](../processingData/README.md).

---

## 6. Attribution

* **Algorithm and code** — original work by the ROAR authors. No third-party solver code is
  incorporated in this folder.
* **Instances** — TSPLIB (G. Reinelt), Waterloo VLSI collection (A. Rohe), Waterloo TSP-art
  collection (R. Bosch; C. S. Kaplan & R. Bosch). Redistributed here unmodified for reproducibility;
  original terms apply.
* **Baselines** — GILS-RVND is in [`../GILS_RVND`](../GILS_RVND/README.md) (adapted from
  `franciscunha/MLP`); PDLSH is in [`../PDLSH`](../PDLSH/README.md) (Yelmewad & Talawar's official
  release).
