# Computing infinity vertices

This guide explains the two infinity-vertex workflows implemented in this
repository:

1. the validated one-point stationary reconstruction, which computes
   aggregate primitive infinity cumulants; and
2. the general CJR graph compiler and dynamic-programming solver, which keep
   contact profiles, descendants, and evaluation insertions separate.

These workflows use descendants in different senses. The aggregate
Bloch--Okounkov workflow uses the ordinary stable-map psi class. CJR's local
formula uses the cotangent line on the log-GLSM curve before stabilization.
The compiler supports both: when stabilization contracts a marked zero tail,
the stable-map descendant is recorded at the adjacent infinity contact.
Primary probes are unaffected.

Run every command from the repository root. The code is tested with SageMath
10.7, and the required `admcycles` package is included under `vendor/`.

## 1. Infinity-vertex notation

The generic unknown is an `EffectiveVertex`:

```sage
load("log_glsm_infinity_dp.sage")

vertex = EffectiveVertex(
    genus=2,
    ambient_degree=0,
    contacts=(-2, -2, -1),
    psi_min=0,
    insertions=(0, 0, 1),
    contact_psi=(0, 0, 2),
)

assert vertex.is_balanced()
assert vertex.is_dimension_zero()
print(vertex)
```

Its fields have the following meanings:

| Field | Meaning |
|---|---|
| `genus` | genus of the connected infinity component |
| `ambient_degree` | degree measured by the hyperplane class of `P^2` |
| `contacts` | negative CJR contact orders, sorted increasingly |
| `psi_min` | power of the minimum-degeneracy cotangent class |
| `insertions` | ambient basis powers at the punctures: `0=1`, `1=H`, `2=H^2` |
| `contact_psi` | ordinary cotangent powers at the same punctures |

For the plane cubic, balance is

```text
3*ambient_degree - (2*genus - 2) = sum(contact + 1).
```

The constructor rejects nonnegative contacts, and `is_balanced()` checks this
equation. A balanced plane-cubic infinity vertex has reduced virtual
dimension

```text
len(contacts) + 3*ambient_degree.
```

Thus a numerical invariant must satisfy

```text
psi_min + sum(insertions) + sum(contact_psi)
    = len(contacts) + 3*ambient_degree.
```

`dimension_defect()` returns the difference between the two sides and
`is_dimension_zero()` checks equality. Be aware that `label` is part of the
object's identity; omit it for contact-resolved vertices that must match
factors produced by the compiler.

The variable `q` in the stationary generating series records intrinsic
degree on the elliptic curve. Ambient plane degree is three times the
intrinsic degree.

## 2. Fast validated genus-one and genus-two reconstruction

The quickest complete calculation is:

```bash
sage log_glsm_infinity_vertices.sage --max-degree 8
```

This compares the known Bloch--Okounkov stationary series with the resummed
positive-degree `O(3)`-twisted zero-level theory. It returns

```text
genus 1, basic contact -1                 = -1/24
genus 2, primitive profile combination   =  1/2880
two genus-1 gluing                        =  1/1152
assembled degree-zero infinity           =  7/5760
```

It also checks that all positive-degree residual coefficients vanish and that
the reconstructed series equals the known elliptic-curve series. Machine
readable output is available with

```bash
sage log_glsm_infinity_vertices.sage --max-degree 8 --json
```

Programmatically:

```sage
load("log_glsm_infinity_vertices.sage")

result = reconstruct_genus_two_infinity_vertices(8)

print(result["genus_one_basic_contact_minus_one"])
print(result["genus_two_primitive_profile_combination"])
print(result["assembled_degree_zero_infinity"])
assert all(result["checks"].values())
```

The output is respectively `-1/24`, `1/2880`, and `7/5760`.

### What the genus-two value means

The one-point equation determines an aggregate primitive cumulant. It does
not, by itself, separate the individual genus-two punctured invariants with
profiles `(-3,)` and `(-2,-2)`, or profiles obtained by adding contact `-1`
legs. The contact-resolved compiler keeps these variables distinct and now
has two probe paths.  Generic stabilized stationary descendants can produce
contact-descendant vertices.  The hypersurface reconstruction instead uses
`effective_basis_only`, where those vertices are forbidden and separation is
tested in the smaller evaluation/`psi_min` basis of Theorem 10.1.

## 3. Aggregate one-point cumulants in arbitrary genus

The triangular one-point provider computes the primitive stationary
cumulants through any chosen genus:

```bash
sage cjr_plane_cubic_equation_provider.sage \
  --max-genus 4 \
  --max-degree 6
```

The first values are

```text
g=1: -1/24
g=2:  1/2880
g=3: -1/181440
g=4:  1/9676800
```

They agree with

```text
-B_(2g) / (2g * (2g)!).
```

The program verifies the reconstructed `q`-series coefficient by coefficient.
The same calculation through the API is:

```sage
load("cjr_plane_cubic_equation_provider.sage")

provider = PlaneCubicOnePointEquationProvider(
    max_genus=4,
    verification_degree=6,
)
solver, values = provider.solve()

for genus in range(1, 5):
    vertex = provider.vertices[genus]
    print(genus, vertex.contacts, values[vertex])

checks = provider.verify(values)
assert all(
    check["diagonal_nonzero"]
    and check["bernoulli_value"]
    and check["all_q_coefficients"]
    for check in checks.values()
)
```

As in genus two, these are aggregate one-point primitive cumulants. The
representative contact `(-(2g-1),)` records the degree-zero balance but does
not assert that the aggregate equals that individual punctured invariant.

## 4. Listing allowed and occurring contact profiles

For degree zero, profiles of fixed genus and valence can be enumerated
directly:

```sage
load("log_glsm_infinity_vertices.sage")

print(infinity_contact_profiles(1, 1))
print(infinity_contact_profiles(2, 1))
print(infinity_contact_profiles(2, 2))
```

This prints

```text
((-1,),)
((-3,),)
((-3, -1), (-2, -2))
```

To see every base profile that actually occurs in the compact-type CJR graph
sum for a one-point probe:

```sage
load("cjr_plane_cubic_equation_provider.sage")

provider = PlaneCubicOnePointEquationProvider(2, 4)
dependencies = provider.combinatorial_infinity_dependencies(
    genus=2,
    ambient_degree=0,
)

for dependency in dependencies:
    print(dependency)
```

The entries are `(genus, ambient_degree, contacts)`. This combinatorial list
does not yet include `psi_min` powers or evaluation insertions; those appear
only after expanding a particular graph contribution.

## 5. Compiling contact-resolved localization equations

`PlaneCubicFullEquationProvider` performs the full calculation:

- enumerate all CJR bipartite graphs;
- evaluate stable zero vertices by base-torus localization on `P^2`;
- integrate their Hodge/psi classes with `admcycles`;
- expand edge, gluing, and infinity-descendant factors;
- collect the result as a polynomial in exact `EffectiveVertex` objects; and
- pull stabilized descendants through `st`, including contracted marked tails;
- insert the convention-compatible compact value.

For the genus-two one-point probe:

```sage
load("cjr_full_equation_provider.sage")

provider = PlaneCubicFullEquationProvider(laurent_precision=8)
gw_probe = ProbeSpec.stationary(
    genus=2,
    ambient_degree=0,
    descendants=(2,),
    label="one-point genus two",
)
assert gw_probe.psi_convention == "stabilized"

compilation = provider.compilation(gw_probe)
assert compilation.is_complete
print("graphs:", len(compilation.contributions))
print("collected monomials:", len(compilation.polynomial.terms))

relation = provider.relation(gw_probe, t_power=0)
print("known reduced invariant:", relation.known_gw)  # -7/5760
print("contact-resolved terms:", len(relation.terms))

for coefficient, factors in relation.terms:
    print("coefficient is nonzero:", coefficient != 0)
    for factor in factors:
        print("  ", factor.to_record())
```

This discrete probe has 11 localization graphs. After the stabilization
pullback and collection, its constant row is

```text
-7/5760 = (-1/6) *
  EffectiveVertex(2,0,(-2,-2,-1),psi_min=0,
                  insertions=(0,0,1),contact_psi=(0,0,2)).
```

The marking lies on a contracted marked zero tail in the surviving term, so
its `psi^2` becomes `contact_psi=2` at the contact `-1`. An explicit
`psi_convention="log"` probe remains available as a different diagnostic,
but its `t^0` side is not supplied by Bloch--Okounkov.

For an auditable per-graph dictionary, use

```sage
report = provider.relation_report(
    gw_probe,
    t_power=0,
    require_complete=True,
)
```

The report records graph signatures, automorphism orders, provenance,
unsupported diagnostics, and the collected effective-vertex factors.

### Splitting `(-3)` from `(-2,-2)`

Do not discard the nonconstant Laurent coefficients. The localization
formula is an identity in the R-weight `t`: its `t^k` coefficient has known
right-hand side zero for `k != 0`, so these coefficients provide additional
relations. After enforcing reduced virtual dimension, the `t^2` relation for
the one-point probe above has just one lower-genus factor,

```text
EffectiveVertex(1, 0, (-1,-1), psi_min=0, insertions=(1,1)),
```

with coefficient `-5/8748`. The previously documented genus-two factor in
this row violated virtual dimension and was spurious. Nonconstant Laurent
coefficients remain useful relations, but dimension pruning is applied before
the rank calculation.

The previously documented two-row string/dilaton witness mixed the two psi
conventions and is invalid. Stabilization pullback is now implemented, but
the corrected rows live in the larger `contact_psi` basis. The provider still
rejects the obsolete two-column matrix because it omitted those variables.
Use the orchestrator's exact support and rank report instead. If that system
is deficient, additional punctured relations involving contact cotangent
lines are required.

## 6. Automatic finite-closure orchestration

`log_glsm_infinity_orchestrator.sage` automates the dependency and rank work
described above. For a requested exact vertex it:

1. applies the CJR punctured unit-removal axioms to vertices without contact
   descendants;
2. compiles stabilization-corrected descendant probes in increasing marking
   count;
3. adds compact unit relatives as further exact rows if earlier stages
   stall;
4. activates only relations incident to the requested dependency closure;
5. substitutes every solved value and detects linear-ready components;
6. selects independent rows by exact matrix rank;
7. solves full-rank blocks, starting with the highest contact/Laurent layer;
8. repeats until the requested roots are solved; and
9. returns an exact kernel and unsupported-level report if the configured
   finite probe family is insufficient.

### Punctured contact-`-1` equations

Do not confuse the two appearances of string/dilaton in this workflow:

- `include_unit_relatives` and `max_unit_insertions` enlarge the **compact
  graph family** (currently only nonconstant homogeneous coefficients are
  usable when those probes contain descendants);
- `include_punctured_axioms` controls CJR relations among the **unknown
  infinity vertices**.

The latter is on by default. For the plane cubic, `cjr_punctured_axioms.sage`
uses `psi_DF=-3H` to implement equations (8.3), (8.4), and Corollary 8.3 of
[Chen--Janda--Ruan, *Punctured logarithmic R-maps*](https://arxiv.org/abs/2208.04519).
If the last marking is the contact-`-1` unit, it generates

```text
string:   V_plus(k; 1)
          = sum_j |c_j| sum_(a=0)^(k-1)
              (-3)^a V_base(k-1-a; H_j^a)

divisor:  V_plus(k; H)
          = D V_base(k)
            + sum_j |c_j| sum_(a=0)^(k-1)
                (-3)^a V_base(k-1-a; H_j^(a+1))

dilaton:  V_plus(k+1; 1) + 3 V_plus(k; H)
          = (2g-2+n) V_base(k)
```

Here `H_j^a` means that the insertion already present at marking `j` is
multiplied by `H^a`; any result above `H^2` is zero. Forgetting is attempted
only in the stable range `2*g-2+n > 0`. At `k=0`, string gives the
fundamental-class vanishing and divisor gives `D*V_base`.

The same provider supplies the CJR closed low-genus inputs. For `O(3)`, genus
zero vanishes; genus one is empty unless `D=0` and all contacts are `-1`; and
the one-marked base values are

```text
EffectiveVertex(1,0,(-1,),psi_min=0,insertions=(1,)) = 0
EffectiveVertex(1,0,(-1,),psi_min=1,insertions=(0,)) = 1/8
```

The second value follows from
`[R_(1,(-1))]^red=3H^2` and `psi_min=lambda-3H`. Higher-valence genus-one
values are then obtained recursively without compiling a compact probe. Every
fully substituted row is checked for zero residual; a disagreement between a
localization convention and a CJR axiom raises
`InconsistentLocalizationRelationError` rather than being ignored.

To inspect the compact-localization system without these inputs, pass
`--no-punctured-axioms` or set `include_punctured_axioms=False`. This is a
diagnostic switch, not the recommended computation. Checkpoint format version
6 records contact cotangent powers. Version-5 checkpoints are intentionally
rejected because their relations used the smaller pre-pullback vertex basis.

The smallest complete command computes the contact-resolved genus-one base
vertex:

```bash
sage log_glsm_infinity_orchestrator.sage \
  --genus 1 \
  --ambient-degree 0 \
  --contacts=-1 \
  --psi-min 0 \
  --insertions 1 \
  --max-markings 1 \
  --t-powers 0 \
  --no-unit-relatives
```

Use `--contact-psi` to target an ordinary cotangent insertion at each
puncture. For example, the vertex emitted by the corrected genus-two
one-point stationary probe is

```bash
sage log_glsm_infinity_orchestrator.sage \
  --genus 2 \
  --ambient-degree 0 \
  --contacts -2 -2 -1 \
  --psi-min 0 \
  --insertions 0 0 1 \
  --contact-psi 0 0 2 \
  --max-markings 1 \
  --t-powers 0
```

For a genus-two primitive target, use a wider Laurent range and preserve a
checkpoint because the probe compilations are much more expensive:

```bash
sage log_glsm_infinity_orchestrator.sage \
  --genus 2 \
  --ambient-degree 0 \
  --contacts=-3 \
  --psi-min 0 \
  --insertions 1 \
  --max-markings 2 \
  --t-powers 2 1 0 -1 -2 -3 -4 \
  --max-unit-insertions 2 \
  --zero-vertex-cache genus2-o3-zero-vertices.sqlite \
  --hodge-cache genus2-o3-hodge.sqlite \
  --checkpoint-out genus2-infinity.json \
  --json
```

If the process is interrupted, repeat the same command with
`--checkpoint-in genus2-infinity.json`. A blocked checkpoint may also be
loaded with a larger `--max-markings`; its extracted relations and solved
values are reused. Keep the same Laurent-power list when resuming.

The equivalent Sage API is:

```sage
load("log_glsm_infinity_orchestrator.sage")

target = EffectiveVertex(
    2, 0, (-3,), psi_min=0, insertions=(1,)
)
config = InfinityOrchestrationConfig(
    max_markings=2,
    t_powers=(2, 1, 0, -1, -2, -3, -4),
    include_unit_relatives=True,
    max_unit_insertions=2,
    include_punctured_axioms=True,
    max_genus=2,
    max_ambient_degree=0,
)
provider = PlaneCubicFullEquationProvider(
    laurent_precision=8,
    zero_vertex_cache="genus2-o3-zero-vertices.json",
    autosave_zero_vertices=True,
)
orchestrator = InfinityVertexOrchestrator(
    provider, (target,), config=config
)
report = orchestrator.run()

if report["status"] == "complete":
    print(orchestrator.solver.values[target])
else:
    print(report["frontier"]["rank_deficiencies"])
    print(report["frontier"]["nonlinear_relation_count"])
    orchestrator.save_checkpoint("genus2-infinity.json")
```

The orchestrator does not replace a singular system by an arbitrary ordering.
`status="blocked"` means that the finite family specified by `max_markings`,
`t_powers`, and the available known-GW backend has not determined the target.
The frontier distinguishes an exact rank deficiency from nonlinear unsolved
dependencies and vertices that occur in no currently linear relation.

If a probe has residual Chow dimension `delta`, its `t^k` coefficient has
ordinary dimension `delta+k`. The numerical backend retains it only when this
number is nonpositive: equality gives a zero-cycle and a negative value gives
a forced-vanishing homogeneous row. Positive-dimensional coefficients require
an additional test class and are rejected. Adding further negative powers
therefore adds exact scalar rows only after this grading check, without a new
graph compilation. Mixed probes are controlled separately:
`max_unit_insertions=1` is the original string/dilaton-shaped family, while
depth two also includes every dimension-zero probe with one or two canonical
unit insertions (`psi^0(1)` or `psi^1(1)`) and at most `max_markings` point
insertions. Their nonconstant Laurent coefficients are exact homogeneous
log-psi relations. Their `t^0` coefficients are skipped until a compatible
compact log-descendant backend exists. These compact graph relations do not
duplicate the punctured unit-removal equations above.

## 6a. Every genus-three vertex at once

`genus_three_infinity_vertices.sage` wraps the orchestrator for a whole
genus. Two things change at genus three.

**Positive infinity degree appears.** Every contact satisfies `c+1 <= 0`, so
the balance equation forces `3*D <= 2g-2`. This is why genus one and two are
purely degree zero, and why genus three splits into two sectors:

```sage
load("genus_three_infinity_vertices.sage")

print(allowed_infinity_degrees(2))   # (0,)
print(allowed_infinity_degrees(3))   # (0, 1)

print(infinity_contact_profiles_in_degree(3, 2, 0))  # (-5,-1), (-4,-2), (-3,-3)
print(infinity_contact_profiles_in_degree(3, 2, 1))  # (-2,-1)
```

`infinity_contact_profiles_in_degree` agrees with
`infinity_contact_profiles` whenever the degree is zero.

**Zero compact sectors provide additional probes.** A curve contained in the
cubic has ambient degree divisible by three, but the ambient log-GLSM moduli
exist in every plane degree. When the degree is not divisible by three,
equation (9.7) has empty hypersurface side. The positive-degree bipartite
assembly is regression-tested in genus one through ambient degree three: after
substituting the CJR infinity values, the complete localization sums are
respectively `0`, `0`, and `-1`, exactly matching the compact side.

Enumerate the truncated family first; it is instant and tells you how much
work a run implies:

```bash
sage genus_three_infinity_vertices.sage --list-only --max-valence 2
```

A contact `-1` is free in the balance equation, so `--max-valence` is what
makes the family finite. Within a fixed valence, insertions range over
`{1,H,H^2}` and virtual dimension fixes
`psi_min = valence + 3*D - sum(insertions)`. Negative powers vanish, while
`--max-psi-min` may discard a required power above the cap. The listed roots
have `contact_psi=0` and are complete for those valence and psi caps.  In
`effective_basis_only` mode the entire dependency closure retains this
property.  Stationary probes in this mode use the log-domain class and obtain
their compact values through the CJR-I stabilization-boundary comparison;
generic direct stabilized-descendant equations may enlarge the state space
with contact-descendant vertices.

Then compute, always with a budget and a checkpoint:

```bash
sage genus_three_infinity_vertices.sage \
  --max-infinity-degree 0 --max-valence 1 --max-psi-min 0 \
  --max-markings 1 --effective-basis-only \
  --max-chow-unit-insertions 0 \
  --t-powers 0 -1 -2 -3 -4 -5 -6 \
  --zero-vertex-cache genus3-o3-zero-vertices.sqlite \
  --hodge-cache genus3-o3-hodge.sqlite \
  --time-budget 600 \
  --checkpoint-out genus3-degree-zero.json
```

The API mirrors the orchestrator:

```sage
load("genus_three_infinity_vertices.sage")

computation = GenusThreeVertexComputation(
    max_valence=1,
    max_psi_min=0,
    max_markings=1,
    t_powers=(0, -1, -2, -3, -4),
    include_unit_relatives=False,
    include_chow_relations=True,
    max_chow_unit_insertions=0,
    chow_primary_only=True,
    effective_basis_only=True,
)
report = computation.run(
    checkpoint_path="genus3-infinity.json", time_budget=3600
)

for record in report["targets"]:
    print(record["description"], record["value"])
```

All targets are roots of a single orchestrator, so the shared genus-one and
genus-two closure and the compiled relations are built once. A target the
configured probe family does not determine is reported with `value=None` and
appears in `report["orchestration"]["frontier"]`; it is never given a value
by a tie-break.

### Profiled feasibility boundary

The optimized CJR graph enumerator itself is no longer the main obstacle in
degree zero. Representative timings are:

| probe data `(g,n,D)` | graph classes | enumeration time |
|---|---:|---:|
| `(3,0,0)` | 32 | 0.15 s |
| `(3,1,0)` | 85 | 0.34 s |
| `(3,2,0)` | 300 | 3.31 s |
| `(3,0,3)` | 478 | 6.17 s |
| `(3,1,3)` | 2,088 | 115.58 s |

The current extractor uses an exact cached recurrence at `t=infinity`, so it
has no Laurent-precision failure. The zero-vertex evaluator uses the
nonresonant path `(0,epsilon,epsilon^2)` and takes `epsilon -> 0` after each
complete fixed-locus sum. A linear ray such as `epsilon*(0,1,3)` is resonant
at a degree-three unstable nodal vertex, where two flag weights cancel. A
fixed numerical choice such as `(0,1,3)` leaves base-equivariant terms
depending on `t`; it is useful for diagnostics but does not define the
required rows.

Historical timings for stationary and compact-unit layers still measure graph
compilation cost, but their old rank figures are invalid because the matrices
omitted contact-descendant columns. Fresh rank figures must use the current
stabilization-corrected support.

The large degree-zero profile recorded on 2026-08-07 predates the exact
infinity virtual-dimension filter. It overcounted the closure by retaining
dimension-incompatible effective vertices, so its old rank and nullity are not
certificates for the corrected system. For `(g,D,c)=(3,0,-5)` and
`psi_min=0`, the corrected target list contains only the `H` insertion; the
unit and `H^2` variants vanish before equation generation. The current
performance and rank figures below come from a fresh checkpoint. Orchestration
checkpoint versions below 8 are intentionally rejected. Version 7 relations
used the twisted virtual dimension to truncate stable zero-side flag
descendants; the correct bound is the ordinary ambient stable-map dimension
because the equivariant Euler class can contribute powers of `t`.
Hodge-integral and zero-vertex caches remain reusable. Nonequivariant
zero-vertex caches carry the quadratic path in their metadata and filename.

### Current status of the `(3,0,-5)` vertex

The target is

```sage
EffectiveVertex(3, 0, (-5,), psi_min=0, insertions=(1,))
```

where insertion power `1` means `H`. The graded infinity backend assigns a
probe coefficient the record

```text
probe dimension = delta
t-power = k
residual Chow dimension = delta+k
```

and admits it to the numerical DP exactly when `delta+k <= 0`. Consequently
the primary one-point probe, whose defect is four, first becomes numerical at
`t^-4`. Its complete row has 50 terms and includes

```text
(-625/72) * V(3,0,(-5); psi_min=0, H).
```

This is a scalar row, but it is not the one-term equation
`(-625/72)*V=0`: all 49 other terms must be retained. Doing so resolves the
former genus-one contradiction. In particular, the defect-two two-unit row
is

```text
V(1,0,(-1); psi_min=1,1)
- 3 V(1,0,(-1); psi_min=0,H)
- V(1,0,(-1,-1); psi_min=2,1,1) = 0,
```

and the known values `1/8`, `0`, `1/8` satisfy it exactly.

The former strongest retry imported 89 additional lower-genus values and left
a 41-dimensional kernel. Its version 7 checkpoint is now diagnostic only:
the corrected flag-descendant bound can add terms to its compiled relations,
so it cannot be used to quote a current rank or value.

A fresh version-8 run through the first one-marking stage, including ambient
degrees zero and one, is saved as
`results/genus3/chow-v8-positive-d1.json`. It compiles 419 relations, activates
341 of them, and proves 90 vertex values. The target belongs to a component
with 161 columns and 202 candidate rows, of rank 108 and nullity 53. Thus this
certified stage is consistent but does **not** determine the `(3,0,-5)` value.
The two-marking stage exposes a performance boundary rather than an algebraic
failure: localization enters a large general Hodge integral in `admcycles`,
while the calibrated Givental backend enters a comparably large symbolic
R-matrix graph product. Neither path produced a new relation checkpoint in an
eight-minute diagnostic run.

Coefficients with positive residual dimension remain class-valued. Extending
the backend to those coefficients requires representing effective
pushforwards as tautological/Chow classes and pairing them with complementary
test classes.

Diagnostic checkpoints named `g3-d0-minus5-class-*` were produced while
testing the invalid scalar projection and must not be used as equations.

Candidate relations are screened by an incremental exact echelon basis. Each
new row is reduced against existing pivots, so probe selection no longer
rebuilds and reranks the entire Sage matrix for every candidate. The final
`ProbeBlockSelection` still constructs an exact matrix for rank and optional
kernel reporting. For a rank-deficient affine block, reduced row-echelon form
also detects pivot coordinates independent of every free variable. Those
coordinates are saved and substituted even when the rest of the component is
singular; no modular heuristic or arbitrary free-variable choice is used.

Positive ambient probe degrees now pass the first three exact genus-one
checks. The earlier residuals `1/24`, `3/4`, and `-3/2` in ambient degrees
one, two, and three came from truncating the stable zero-side flag series by
the twisted virtual dimension. With the ambient stable-map bound, the full
localization sums are respectively `0`, `0`, and `-1`, exactly their compact
hypersurface sides. Higher-degree runs remain guarded by the exact solved-row
consistency check.

### Precompute and reuse the O(3) zero vertices

The persistent table is implemented. First inventory a finite probe family
without doing any tautological integration:

```bash
sage o3_zero_vertex_precompute.sage \
  --genus 3 --ambient-degree 1 --max-point-markings 1 \
  --cache genus3-o3-zero-vertices.sqlite --collect-only --json
```

Then remove `--collect-only` and give the batch a wall-clock budget:

```bash
sage o3_zero_vertex_precompute.sage \
  --genus 3 --ambient-degree 1 --max-point-markings 1 \
  --cache genus3-o3-zero-vertices.sqlite \
  --hodge-cache genus3-o3-hodge.sqlite --time-budget 21600
```

The inventory backend visits every stable zero option but returns zero before
assembling infinity contributions. The evaluation pass skips values already
in the cache. A `.sqlite`, `.sqlite3`, or `.db` cache commits each exact value
independently in WAL mode; a `.json` cache retains the older atomic whole-file
format. It is consequently safe to resume the same command after interruption.
The equation provider and precompute CLI use
`(0,epsilon,epsilon^2)` followed by the limit `epsilon -> 0`, and keep that
table separate: passing `CACHE.sqlite` to the genus-three driver derives
`CACHE.weights-nonequivariant-quadratic.sqlite`. The precompute CLI writes the same
non-equivariant convention to the same derived path, so the base cache name in
the precompute and orchestration commands should match.
A budget is checked between individual zero-vertex requests; it cannot
interrupt one admcycles evaluation already in progress. The request list
itself is stored in `CACHE.inventory.json`; subsequent runs load that manifest
instead of repeating the CJR enumeration. Use `--inventory PATH` to place it
elsewhere.

After the one-time `--collect-only` pass, deterministic shards can evaluate
the same manifest concurrently. Run one command per terminal, changing only
`--shard-index`:

```bash
sage o3_zero_vertex_precompute.sage \
  --genus 3 --ambient-degree 1 --max-point-markings 1 \
  --cache genus3-o3-zero-vertices.sqlite \
  --hodge-cache genus3-o3-hodge.sqlite \
  --shard-count 4 --shard-index 0 --time-budget 21600
```

Indices `0,1,2,3` partition the sorted request list. The zero-vertex and Hodge
stores are safe to share between those processes; `INSERT OR IGNORE` also
makes a repeated shard harmless. A completed Hodge integral is immediately
available to later requests and later runs.

There are four normalization/reuse layers before this batch reaches
admcycles:

1. provenance-only edge labels do not participate in request equality;
2. scalar insertion factors are removed from persistent zero-vertex keys;
3. marked insertions are put in a canonical order before sparse lifts and
   fixed-graph enumeration; and
4. psi powers are sorted in `HodgeIntegralRequest`.

All stable degree-zero vertices are evaluated directly on
`Mbar_(g,n) x P2`, with sparse lifts restricting the fixed-point sum. Hodge
monomials are reduced by Mumford's relation through genus three. Pure psi,
single `lambda_g`, and `lambda_g lambda_(g-1) lambda_(g-2)` descendants use
closed formulas (string/dilaton determines the descendants of the last top
class); only irreducible mixed monomials reach admcycles. Lambda products,
edge factors, stable-vertex values, Hodge integrals, and final zero vertices
are cached at their natural levels.

For example, the unmarked twisted genus-three degree-zero vertex previously
remained inside `admcycles.common_degenerations` after five minutes. The
specialized backend evaluates it in about 0.07 seconds and agrees exactly with
the three-graph localization sum (about 0.10 seconds). A sparse one-divisor
top-descendant request is smaller already: the direct and explicit paths take
roughly 0.38 and 0.22 seconds respectively, so the optimization is not a
universal constant-factor win. Its purpose is to eliminate the pathological
Hodge products while retaining one exact implementation for every stable
degree-zero request.

Historical ambient-degree-three profiling took about 153 seconds for the
one-mark genus-three stationary probe. The positive-degree assembly is now
consistency-checked in genus one through ambient degree three. Completed zero
vertices remain reusable, but compiled version-7 relations are rejected
because they may contain the former flag-descendant truncation. Pass the same
zero-vertex cache file to a fresh version-8 run with `--zero-vertex-cache`.

The degree-zero family with one stationary probe and all mixed probes through
unit depth two is already substantial: 6 probes produced 1,502 distinct
requests and took 328 seconds to inventory. Of these, 920 are genus zero and
582 have higher genus. All now use the direct constant-map path, while the
SQLite tables ensure any remaining admcycles values are paid only once. This
is why unit depths are separate orchestration stages rather than one
monolithic default family.

### Calibrated Givental--Teleman strategy

The repository contains an exact stable-graph R-matrix engine and a calibrated
positive-degree O(3) theory. To use it through either orchestration command,
add

```bash
--zero-vertex-strategy hybrid
```

The hybrid strategy uses Givental--Teleman for every supported stable request
and retains localization as a fallback. `--zero-vertex-strategy givental`
uses only the calibrated reconstruction. `localization` remains the default
until larger positive-degree workloads have been profiled.

`O3TwistedQuantumRing` retains the raw Picard--Fuchs principal symbol for
diagnostics. The calibrated path deliberately does not use its non-flat
residue metric. `O3CalibratedGiventalData` obtains the small quantum product
from exact genus-zero three-point invariants in the flat basis, Hensel-lifts
its eigenvalues, and solves the Dubrovin recursion. Its q=0 integration
constants are the diagonal Mumford-GRR/quantum-Riemann--Roch R-matrix. The
result is checked for metric self-adjointness and symplecticity at every
requested truncation.

The stable-graph action produces ancestors. Before extracting the requested
degree, `O3CalibratedGiventalBackend` applies the S-matrix
ancestor--descendant transformation to every insertion. Exact tests compare
positive-degree primaries and descendants with virtual localization through
degree two.

### Would a faster implementation language help?

For the `(3,0,3)` graph family, cProfile measured 11.20 seconds and assigned
3.02 seconds (27%) to `_zero_vertex_types`. The matching deterministic kernel
benchmark took 1.68 seconds in Python and 0.092 seconds in Go, an 18.2x native
speedup with identical checksums. Amdahl's law therefore predicts only about
1.34x for that whole graph enumeration if this one loop becomes native.

That is a real but secondary improvement. It does not touch the positive-
degree bottleneck, which is exact Sage rational-function algebra and admcycles
tautological multiplication. Rewriting that layer would amount to
reimplementing the mathematical backend, not translating a hot loop. The
reproducible scripts are in `benchmarks/`; the recommended order is persistent
zero-vertex caching, mathematical specialization of recurring requests, and
only then a small native extension for the combinatorial enumerator.

Use staged checkpoints. First run only `--max-infinity-degree 0`, then resume
with a richer probe family:

```bash
sage genus_three_infinity_vertices.sage \
  --max-infinity-degree 0 --max-valence 1 --max-psi-min 0 \
  --max-markings 1 --max-unit-insertions 2 \
  --t-powers 2 1 0 -1 -2 -3 -4 \
  --zero-vertex-cache genus3-o3-zero-vertices.sqlite \
  --hodge-cache genus3-o3-hodge.sqlite \
  --time-budget 1800 \
  --checkpoint-in genus3-degree-zero.json \
  --checkpoint-out genus3-degree-zero.json
```

Large genus-three runs use compact frontier diagnostics by default: every
deficiency records exact `rank`, `columns`, and `kernel_dimension`, while
`kernel` is `null`. Pass `--full-kernel` only when explicit basis vectors are
needed; converting their rational-function entries to strings can take longer
than compiling the relations. `--time-budget` is checked between stages and
solve blocks, not inside one probe compilation, fixed-locus sum, or admcycles
integral. `--checkpoint-in` resumes relations that reached a completed stage.

## 7. Manual solution of a contact-resolved vertex

The DP solver expects one localization equation per target. If all factors
other than the target are already known and strictly lower in `order_key()`,
the scalar pattern is:

```sage
# Choose an exact target appearing in relation.terms.
target = next(
    factor
    for coefficient, factors in relation.terms
    for factor in factors
    if factor.genus == 2
)

provider.assign_probe(target, log_probe, t_power=-1)

# Supply already known lower contact-resolved vertices here. Values must lie
# in provider.rings.base_field. An empty dictionary is sufficient only when
# the equation has no unsolved lower factors.
lower_values = {}
solver = provider.solver(initial_values=lower_values)

# This succeeds once every lower dependency has either an assigned probe or
# an entry in lower_values.
# value = solver.solve(target)
```

Leaving `value = solver.solve(target)` commented is intentional in this
example: the displayed genus-two relation contains lower genus-one vertices
with their own contact, insertion, and descendant data. The repository does
not identify those with the aggregate `-1/24`; they must be computed by
separate probes or supplied as proven initial values.

Use

```sage
solver.dependency_graph((target,))
```

to inspect the recursive dependencies before starting a long calculation.
The solver raises `NonTriangularLocalizationError` rather than silently using
a non-lower unknown.

### Same-rank blocks

When several contact profiles mix, select independent probe rows and solve
them together:

```sage
targets = tuple(contact_resolved_targets)
lower_values = dict(already_solved_lower_vertices)

candidates = provider.candidate_relations(
    genus=2,
    ambient_degree=0,
    max_markings=3,
    t_powers=(2, 1, 0, -1),
    require_complete=True,
    include_unit_relatives=True,
)

factory = ProbeFactory(provider.rings.base_field)
selection = factory.select_relations(
    targets,
    candidates,
    lower_values=lower_values,
)

print("rank:", selection.rank, "/", len(targets))
if not selection.is_full_rank:
    print("kernel:", selection.kernel_basis())
    print("rejected probes:", selection.rejected)
else:
    for target, selected_relation in zip(targets, selection.relations):
        provider.assign_probe(
            target,
            selected_relation.probe,
            selected_relation.t_power,
        )
    solver = provider.solver(initial_values=lower_values)
    values = solver.solve_block(targets)
```

The names `contact_resolved_targets` and `already_solved_lower_vertices` are
placeholders for the finite block being studied. A singular selection is a
mathematical lack of independent probes, not a request for more graph
enumeration. Add different stabilized descendants, markings, or nonconstant
Laurent coefficients until the exact diagonal matrix has full rank.
Contact-descendant vertices are genuine columns and must not be projected out
of the block.

## 8. Checkpoints

Solved exact values can be serialized without converting them to floating
point numbers:

```sage
checkpoint = solver.checkpoint()
print(checkpoint)

new_solver = provider.solver()
new_solver.restore_checkpoint(checkpoint)
```

The checkpoint records the complete vertex signature and stores values as
exact strings in the configured coefficient field.

For an orchestrated calculation, prefer

```sage
orchestrator.save_checkpoint("infinity.json")

resumed = InfinityVertexOrchestrator(
    provider, (target,), config=config
)
resumed.load_checkpoint("infinity.json")
report = resumed.run()
```

This version also stores every extracted localization relation, the expanded
vertex closure, completed probe stages, solved-block matrices, and failure
diagnostics. Restoring it does not re-enumerate graphs already represented in
the checkpoint.

## 9. What is currently automatic

The following are end-to-end and regression tested:

- genus-one and genus-two aggregate stationary reconstruction;
- aggregate one-point primitive cumulants in arbitrary finite genus;
- CJR graph enumeration with exact automorphisms and contacts;
- full `O(3)`-twisted zero-vertex evaluation;
- exact compilation of a chosen probe into contact-resolved terms; and
- finite dependency-closure discovery;
- progressive stabilization-corrected probe expansion;
- exact full-rank block selection and descending back-substitution; and
- versioned checkpoints containing both relations and solved values.

The fast aggregate routine does not identify its stationary cumulant with an
individual contact-resolved punctured invariant. The generic compiler and
orchestrator keep profiles and contact descendants as distinct variables.
Ordinary stationary GW values now give valid equations, but they do not
guarantee that a chosen finite probe family has full rank for every requested
block.
Always inspect the orchestration frontier when a run returns `blocked`.

### Targeted checkpoint extensions

A richer checkpoint can be extended at a new marking depth without repeating
all of its positive-degree probes:

```bash
sage genus_three_infinity_vertices.sage \
  --max-infinity-degree 0 --max-valence 1 --max-psi-min 0 \
  --max-markings 3 --chow-relations --primary-chow-relations-only \
  --additional-probe-ambient-degrees 1 2 \
  --future-probe-ambient-degree-ceiling 0 \
  --checkpoint-in results/genus3/chow-v8-positive-d2-stage2-convolution.json \
  --checkpoint-out results/genus3/chow-v8-d2-plus-degree0-m3.json
```

The ceiling affects only relations compiled after the checkpoint is loaded.
It does not delete saved relations or alter the serialized mathematical
configuration.  Resuming a checkpoint now also finishes any pending exact
block solve before advancing to a new probe stage.

### Current `(3,0,-5)` frontier

Section 9.3 specializes the infinity contribution to

```text
st_*(t [U_infinity]^red / (-t - psi_min)).
```

Accordingly, Theorem 10.1 uses effective invariants containing evaluation
classes and powers of `psi_min`, but no ordinary cotangent descendants at the
contacts.  Run this basis explicitly with

```bash
sage genus_three_infinity_vertices.sage \
  --max-infinity-degree 0 --max-valence 1 --max-psi-min 0 \
  --max-markings 2 --effective-basis-only \
  --max-chow-unit-insertions 0 \
  --additional-probe-ambient-degrees 1 2 \
  --t-powers 0 -1 -2 -3 -4 -5 -6 \
  --zero-vertex-cache results/genus3/zero-vertices.sqlite \
  --hodge-cache results/genus3/hodge.sqlite \
  --checkpoint-out results/genus3/effective-basis-d2-m2.json
```

`--effective-basis-only` implies primary Chow relations, disables compact
descendant unit relatives, and raises an error if graph compilation emits a
nonempty `contact_psi` tuple.  This hard assertion prevents a stabilized
stable-map descendant experiment from being mistaken for the CJR effective
basis.

The exact run above does not yet determine

```text
EffectiveVertex(g=3,D=0,contacts=(-5,),psi_min=0,insertions=(1,)).
```

The historical version-8 checkpoint contains 1259 compiled relations, 1205 active relations, 692
vertices, and 349 solved values.  After partial coordinate recovery, the
target component has rank 278 on 343 columns (kernel dimension 65).  In the
computed kernel basis, only one vector moves the requested coordinate.  Its
support consists of 18
degree-zero and degree-one effective invariants related by the punctured
contact-`-1` equations; every vertex has `contact_psi=[]`.

The newly implemented one-point boundary-comparison row has compact value
`-31/967680`, target coefficient `-625/72`, and no `contact_psi`.  It raises
that historical component's rank from 278 to 279 but does not remove the one
kernel direction moving the target.  Adding all degree-zero, at-most-two-point
boundary-compared stationary rows solves four additional coordinates and
leaves the target component at rank 279 on 339 columns.

Thus the certified answer remains `UNDETERMINED`, not zero, but the reason is
now precise: this finite family leaves one direction through the target.
Theorem 10.1 computes hypersurface GW invariants from twisted invariants and
effective invariants; it does not assert the converse.  Remark 10.7 observes
that negative Laurent coefficients give relations among effective invariants
and explicitly leaves their computational implications for further study.

### Why "basic effective invariants" do not apply here

An earlier version of this note suggested that computing one *basic effective
invariant* would remove the remaining kernel direction.  That is not supported
for the plane cubic, and the reason is worth recording.

CJR II, Definition 9.19, defines basic effective invariants only for `g >= 2`
and only for `beta` satisfying

```text
(3 - dim X + rk E)*(g-1) - int_beta c_1(K_X (x) det E) = 0.
```

For `X = P^2` and `E = O(3)` we have `K_X (x) det E = O(-3) (x) O(3) = O`, so
the integral vanishes for every `beta` and the left side is `2*(g-1)`.  This is
zero only at `g = 1`, which Definition 9.19 excludes.  **No genus-three basic
effective invariant exists in CJR's sense.**

Assumption 9.22, which Corollary 9.24 needs in order to conclude that every
`g >= 2` effective invariant either vanishes or is determined by the basic
ones, requires `2*(g-1) <= 0`.  It fails for every `g >= 2`.

The controlling quantity is `3 - dim X + rk E = 3 - dim Z`.  CJR's finiteness
results are built for Calabi--Yau threefolds, where it is zero (Example 9.23).
Here `Z` is an elliptic *curve*, so it is `+2`, and

```text
red dim R = vir dim Mbar_{g,n}(Z,beta) + sum_i(c_i+1) = (2g-2+n) + sum_i(c_i+1)
```

grows with `g`.  The family of effective invariants really is infinite, which
is what the closure explosion in these runs reflects: five compiled relations
at genus three already introduce 518 distinct vertices.

`basic_minus_two_vertex` in `cjr_plane_cubic_equation_provider.sage` builds the
all-contact-`-2` vertex by analogy with Definition 9.19.  It is a useful
bookkeeping device, and its star coefficient `1/n_2!` is independently
confirmed by the graph compiler, but it is **not** an instance of that
definition and carries none of its determination properties.

### What has been measured and does not help

*More localization probes.* The zero-marking genus-three degree-zero probe is
absent from `stationary_candidates`, which starts at one marking.  Its `t^-4`
coefficient is a legitimate Chow-scalar row containing both the target
(coefficient `125/24`) and the all-`-2` vertex (coefficient `-1/24`).  Adding
it and its one-marking companion to the 1259-relation checkpoint leaves the
target component at rank 278 on 343 columns, kernel dimension 65, and the same
single kernel direction through the target: the new rows lie in the existing
row space.

*The cycle-valued unit axiom, CJR II Theorem 8.4, equation (8.8).* Pairing it
with a test class `Theta` and applying the projection formula gives

```text
int_{Mbar_{g,n+1}} Theta . LHS
    = int_{Mbar_{g,n}} pi_* Theta . A
    + sum_j |c_j| int_{Mbar_{g,n}} (Theta|_delta) . B_j.
```

The variables are dimension-zero effective invariants, so `A` and `B_j` are
point classes and only `deg Theta = 0` contributes, which returns equation
(8.3) -- already implemented.  Any other `Theta` introduces
tautological-weighted effective *cycles* as new unknowns rather than new
equations.  Remark 9.26 uses (8.8) in the opposite direction, propagating known
dimension-zero values up to higher-dimensional cycles.

### The genus-one closed form as a convention test

`PlaneCubicPuncturedAxioms.genus_one_closed_form()` derives the two
dimension-zero genus-one values from CJR II (9.17) and (9.21) instead of
tabulating them:

```sage
load("cjr_punctured_axioms.sage")
closed = PlaneCubicPuncturedAxioms().genus_one_closed_form()
print(closed["reduced_cycle"])    # 3*H^2
print(closed["psi_min"])          # lam - 3*H
print(closed["psi_min_value"])    # 1/8
print(closed["divisor_value"])    # 0
```

This pins the sign of `psi_min`, the `3H^2` normalization, and the ambient
insertion basis simultaneously.  It is a two-sided check: the divisor value
discriminates candidate reduced cycles that the `psi_min` value alone cannot.
An ansatz omitting the `(+) O - N` term of (9.21) reproduces `1/8` but predicts
`-1/8` for the divisor value, and is thereby rejected.

The older `chow-v8-*` checkpoints combine stationary descendant probes with
stabilization pullbacks and contain ordinary `contact_psi` columns.  They are
useful tests of the generic graph compiler, but their larger ranks and kernels
must not be quoted as the Theorem 10.1 hypersurface reconstruction.

All version-8 checkpoints predate the boundary-comparison rows and are
rejected by the version-9 orchestrator.  Their zero-vertex and Hodge caches
remain reusable.
