# Computing infinity vertices

This guide explains the two infinity-vertex workflows implemented in this
repository:

1. the validated one-point stationary reconstruction, which computes
   aggregate primitive infinity cumulants; and
2. the general CJR graph compiler and dynamic-programming solver, which keep
   contact profiles, descendants, and evaluation insertions separate.

Run every command from the repository root. The code is tested with SageMath
10.7, and the required `admcycles` package is included under `vendor/`.

## 1. Infinity-vertex notation

The generic unknown is an `EffectiveVertex`:

```sage
load("log_glsm_infinity_dp.sage")

vertex = EffectiveVertex(
    genus=2,
    ambient_degree=0,
    contacts=(-2, -2),
    psi_min=0,
    insertions=(0, 0),
)

assert vertex.is_balanced()
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

For the plane cubic, balance is

```text
3*ambient_degree - (2*genus - 2) = sum(contact + 1).
```

The constructor rejects nonnegative contacts, and `is_balanced()` checks this
equation. Be aware that `label` is part of the object's identity; omit it for
contact-resolved vertices that must match factors produced by the compiler.

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
legs. The profiles are nevertheless separable: use the contact-resolved
workflow in Section 5, several Laurent coefficients of the R-equivariant
identity, and a full-rank family of stationary, string, and dilaton probes.

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
- insert the known elliptic-curve Gromov--Witten value.

For the genus-two one-point probe:

```sage
load("cjr_full_equation_provider.sage")

provider = PlaneCubicFullEquationProvider(laurent_precision=8)
probe = ProbeSpec.stationary(
    genus=2,
    ambient_degree=0,
    descendants=(2,),
    label="one-point genus two",
)

compilation = provider.compilation(probe)
assert compilation.is_complete
print("graphs:", len(compilation.contributions))
print("collected monomials:", len(compilation.polynomial.terms))

relation = provider.relation(probe, t_power=0)
print("known reduced value:", relation.known_gw)
print("contact-resolved terms:", len(relation.terms))

for coefficient, factors in relation.terms:
    print("coefficient is nonzero:", coefficient != 0)
    for factor in factors:
        print("  ", factor.to_record())
```

At the current conventions this probe has 11 localization graphs, 53
pre-extraction monomials, and 11 terms after taking the `t^0` coefficient.
Seven of those terms contain a genus-two effective vertex.
The coefficients lie in `QQ(lambda0,lambda1,lambda2)` and may be large; print
the factor records first, and inspect coefficients only when needed.

`relation.known_gw` is the signed reduced log-GLSM value used in the CJR
equation, not necessarily the unsigned Gromov--Witten invariant printed by a
stationary-series helper.

For an auditable per-graph dictionary, use

```sage
report = provider.relation_report(
    probe,
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
relations. For the one-point probe above, the `t^2` relation has just one
genus-two factor,

```text
EffectiveVertex(2, 0, (-2,-2,-1), psi_min=0, insertions=(0,0,1)),
```

with coefficient `-1/6`. Lower-genus and zero-level terms still have to be
substituted, but this illustrates the descending-Laurent/back-substitution
structure: start with the highest nonzero powers of `t`, solve the vertices
with extra contact `-1` legs, and descend toward `t^0` and the primitive
profiles.

Additional probe rows make the contact-profile blocks nonsingular. Generate
the string and dilaton relatives of the one-point probe as follows:

```bash
sage cjr_full_equation_provider.sage --profile-split-rank
```

This compiles two 28-graph probes and can take several minutes. It prints the
matrix below and verifies that it has full rank. The equivalent API is:

```sage
load("cjr_full_equation_provider.sage")

witness = genus_two_profile_split_rank_witness()
print(witness["matrix"])
assert witness["full_rank"]
assert witness["determinant"] == 9
```

The rows are `(-27/2, 4)` and `(-9, 2)`. Their nonzero determinant is a
concrete rank witness that the `(-3,-1)` and `(-2,-2)` sectors are not forced
to remain aggregated. After the additional contact `-1` sectors have been
solved, the one-point `t^0` equation contains the primitive `(-3,)` vertex
with nonzero coefficient; back-substitution therefore separates `(-3,)`
from `(-2,-2)`.

This two-row calculation is a rank witness, not an isolated two-variable
equation: both relations contain other effective vertices. In an actual
calculation, put every unresolved factor of the relevant rank into the block,
solve the higher Laurent/contact layers first, and use `select_relations()`
on the complete dependency closure.

## 6. Automatic finite-closure orchestration

`log_glsm_infinity_orchestrator.sage` automates the dependency and rank work
described above. For a requested exact vertex it:

1. compiles stationary probes in increasing marking count;
2. adds string and dilaton relatives if the stationary stages stall;
3. activates only relations incident to the requested dependency closure;
4. substitutes every solved value and detects linear-ready components;
5. selects independent rows by exact matrix rank;
6. solves full-rank blocks, starting with the highest contact/Laurent layer;
7. repeats until the requested roots are solved; and
8. returns an exact kernel and unsupported-level report if the configured
   finite probe family is insufficient.

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
  --t-powers 2 1 0 -1 -2 \
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
    t_powers=(2, 1, 0, -1, -2),
    include_unit_relatives=True,
    max_genus=2,
    max_ambient_degree=0,
)
provider = PlaneCubicFullEquationProvider(laurent_precision=8)
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

provider.assign_probe(target, probe, t_power=0)

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
enumeration. Add different descendants, markings, or Laurent coefficients
until the exact diagonal matrix has full rank.

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
- progressive stationary/string/dilaton probe expansion;
- exact full-rank block selection and descending back-substitution; and
- versioned checkpoints containing both relations and solved values.

The fast aggregate routine does not identify its stationary cumulant with an
individual contact-resolved punctured invariant. The generic compiler and
orchestrator keep—and can distinguish—the profiles. They do not guarantee
that a chosen finite probe family has full rank for every requested block.
Always inspect the orchestration frontier when a run returns `blocked`.
