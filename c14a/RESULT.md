# 14a, the least n for which BB(n) is independent of ZFC: the determinations of BB(3) and BB(4) (15 November 2026)

The page's lower bounds are $C_{14} \ge 4$ [LS1965] "by computing $BB(3)$" and $C_{14} \ge 5$ [B1983] "by computing
$BB(4)$" (and $\ge 6$ [BB2025], the Coq-checked BB(5), not attempted here). A value of $BB(n)$ is settled once every
$n$-state, 2-symbol machine started on a blank tape is either run to its halt or shown, by a proof checkable in finite
time, never to halt. Two separately written programs redo that enumeration for $n = 2, 3, 4$; nothing is downloaded and
no program of any source was read or run.

Conventions: states A, B, ...; the machine starts in A on a blank tape; $S(n)$ counts steps including the halting
transition, $\Sigma(n)$ counts the 1s on the tape after it, the halting transition writing 1.

- `c14a_a.c`: tree normal form by re-simulation (a transition is defined only when the run first reads it; the next
  state is at most one more than the largest used; the first move goes right, by mirror symmetry; one transition always
  stays undefined, as a halting machine needs one). A run past the step limit (2,000) must be shown never to halt by
  **closed sets of local views**: the views (r cells left of the head, state, head cell, r cells right), grown from the
  blank tape together with the sets of r-cell windows seen on each side, until closed under every step; if no view
  reads an undefined transition, no halting configuration is reachable. Radius r = 1 to 12; then each cell also
  carries the last h states that wrote it (h = 1 to 9, r = 1 to 4). View keys are 128-bit integers, and every size is
  checked before it is formed.
- `c14a_b.c`: the same normal form built from saved runs (each node resumes its parent's run); a run past the limit is
  settled by the first of: **cyclers** (the whole configuration recurs; Brent's method), **translated cyclers** (two
  head records on one side in one state, with the cells from the record back to the deepest point reached in between
  equal, so the run repeats shifted), and its own **closed sets of local views** (writer history h = 0 to 9).

| | both programs |
|---|---|
| $n = 2$ | 61 machines: 19 halt, 42 never halt; $S(2) = 6$, $\Sigma(2) = 4$ |
| $n = 3$ | 5,417 machines: 1,772 halt, 3,645 never halt, **every machine decided**: $S(3) = 21$, $\Sigma(3) = 6$ |
| $n = 4$ | 858,909 machines: 249,693 halt (the longest 107 steps, the most 1s 13), 609,204 shown never to halt; **12 left undecided**, the same 12 in both programs |

Program B settles the 3-state non-halters as 526 cyclers, 3,027 translated cyclers and 92 by closed view sets; the
4-state ones as 87,754 cyclers, 500,902 translated cyclers and 20,548 by closed view sets. Program A uses closed view
sets alone.

So the [LS1965] row is replayed: every 3-state machine is decided, $S(3) = 21$ and $\Sigma(3) = 6$. For [B1983] the
replay is **partial**: the values $S(4) = 107$ and $\Sigma(4) = 13$ are attained (by
`1RB1LB_1LA0LC_---1LD_1RD0RA`, which halts after 107 steps with 13 ones) and all but 12 of the 858,909 machines are
decided; the 12 left are Christmas-tree machines (a single 0 travelling through growing blocks of 1s), which Brady
also had to treat separately, and neither program proves them:

    0RB1LA_1LC1RB_1LA1LD_---0LA   0RB1LA_1LC1RB_1RA1LD_---0LA   0RB1LA_1LC1RB_1RB1LD_---0LA
    1RB0LC_1LB1LA_1RC0LD_---0RA   1RB0LD_1LC1RB_---1LA_0RB1LD   1RB0RB_0LC1RB_1RD1LC_---1RA
    1RB0RB_1LC1RB_---0LD_1RA1LD   1RB1LA_1LA0RC_1LD1RC_---0LA   1RB1LA_1LA1RC_---0RD_0LA1RD
    1RB1LA_1LC1RD_0LA1RC_---0RC   1RB1LA_1RC1RD_0LA1RC_---0RC   1RB1LB_1RC0RD_1LD1RB_---0LA

**Red controls**, each met by both programs: (1) $n = 3$ with a step limit of 20: the 21-step machine
`1RB---_1LB0RC_1LC1LA` must stay undecided, and is the only one left; (2) $n = 4$ with a limit of 100: the 107-step
champion must stay undecided (13 left: it and the 12 above); (3) $n = 3$ with the tests cut down (radius 1, no
history; B also without translated cyclers): 1,756 and 1,586 machines left, no value claimed; (4) every closed-set test
(and in B the cycler tests) applied to every one of the 1,772 halting 3-state machines: 0 false proofs.

Control (1) and (2) caught two faults during development, both fixed before any result was kept: view keys packed
into 64 bits overflowed at large radius with history (distinct views collided, and the 107-step champion was
"proved" not to halt), and the size check itself multiplied before testing (`P *= K` wrapped past $2^{64}$). Both
programs now use 128-bit keys and test every size before forming it.

**Not checked:** BB(5) (the [BB2025] row); that the decided machines' non-halting proofs are ZF theorems is the usual
reading (each test is a finite inductive invariant).

Reproduce: `./run_all.sh` (a C compiler; about 15 minutes on two cores, under 600 MB). Output: `results/c14a.txt`.
