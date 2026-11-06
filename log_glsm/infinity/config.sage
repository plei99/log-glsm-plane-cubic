r"""infinity / config."""

from sage.all import *



__all__ = ['InconsistentLocalizationRelationError', 'InfinityOrchestrationConfig', 'OrchestrationLimitError', 'ProbeExpansionStage']


class OrchestrationLimitError(RuntimeError):
    """A caller-supplied finite resource bound was exceeded."""

class InconsistentLocalizationRelationError(ArithmeticError):
    """A fully evaluated localization/axiom row has nonzero residual."""

class ProbeExpansionStage(SageObject):
    """One progressively richer family of dimension-compatible probes."""

    def __init__(self, max_markings, include_unit_relatives=False,
                 max_unit_insertions=None):
        self.max_markings = ZZ(max_markings)
        self.include_unit_relatives = bool(include_unit_relatives)
        if max_unit_insertions is None:
            max_unit_insertions = 1 if self.include_unit_relatives else 0
        self.max_unit_insertions = ZZ(max_unit_insertions)
        if self.max_markings < 1:
            raise ValueError("a probe stage needs at least one marking")
        if self.max_unit_insertions < 0:
            raise ValueError("max_unit_insertions must be nonnegative")
        if self.include_unit_relatives != bool(self.max_unit_insertions):
            raise ValueError(
                "unit stages need a positive max_unit_insertions depth"
            )

    def to_record(self):
        return {
            "max_markings": int(self.max_markings),
            "include_unit_relatives": self.include_unit_relatives,
            "max_unit_insertions": int(self.max_unit_insertions),
        }

    @classmethod
    def from_record(cls, record):
        return cls(
            record["max_markings"],
            record.get("include_unit_relatives", False),
            record.get("max_unit_insertions"),
        )

    def __repr__(self):
        suffix = (
            "+unit<=%s" % self.max_unit_insertions
            if self.include_unit_relatives else ""
        )
        return "ProbeExpansionStage(markings<=%s%s)" % (
            self.max_markings, suffix
        )

class InfinityOrchestrationConfig(SageObject):
    """Finite bounds and probe schedule for an orchestration run."""

    def __init__(self, max_markings=2, t_powers=(0, -1, -2),
                 include_unit_relatives=True, max_unit_insertions=None,
                 include_chow_relations=False,
                 max_chow_unit_insertions=1,
                 chow_primary_only=False,
                 effective_basis_only=False,
                 additional_probe_ambient_degrees=(),
                 include_punctured_axioms=True,
                 max_genus=None,
                 max_ambient_degree=None, max_vertices=2000,
                 max_relations=4000, max_solve_rounds=1000,
                 require_complete=True, fail_on_unsupported=False,
                 include_kernel_basis=True):
        self.max_markings = ZZ(max_markings)
        self.t_powers = tuple(ZZ(power) for power in t_powers)
        self.include_unit_relatives = bool(include_unit_relatives)
        if max_unit_insertions is None:
            max_unit_insertions = 1 if self.include_unit_relatives else 0
        self.max_unit_insertions = ZZ(max_unit_insertions)
        self.include_chow_relations = bool(include_chow_relations)
        self.max_chow_unit_insertions = ZZ(max_chow_unit_insertions)
        self.chow_primary_only = bool(chow_primary_only)
        self.effective_basis_only = bool(effective_basis_only)
        self.additional_probe_ambient_degrees = tuple(sorted(set(
            ZZ(degree) for degree in additional_probe_ambient_degrees
        )))
        self.include_punctured_axioms = bool(include_punctured_axioms)
        self.max_genus = None if max_genus is None else ZZ(max_genus)
        self.max_ambient_degree = (
            None if max_ambient_degree is None
            else ZZ(max_ambient_degree)
        )
        self.max_vertices = ZZ(max_vertices)
        self.max_relations = ZZ(max_relations)
        self.max_solve_rounds = ZZ(max_solve_rounds)
        self.require_complete = bool(require_complete)
        self.fail_on_unsupported = bool(fail_on_unsupported)
        self.include_kernel_basis = bool(include_kernel_basis)

        if self.max_markings < 1:
            raise ValueError("max_markings must be positive")
        if self.max_unit_insertions < 0:
            raise ValueError("max_unit_insertions must be nonnegative")
        if self.max_chow_unit_insertions < 0:
            raise ValueError(
                "max_chow_unit_insertions must be nonnegative"
            )
        if any(degree < 0 for degree in self.additional_probe_ambient_degrees):
            raise ValueError("additional probe degrees must be nonnegative")
        if self.include_unit_relatives != bool(self.max_unit_insertions):
            raise ValueError(
                "include_unit_relatives and max_unit_insertions disagree"
            )
        if self.effective_basis_only and (
                not self.include_chow_relations
                or not self.chow_primary_only
                or self.include_unit_relatives):
            raise ValueError(
                "effective-basis-only reconstruction requires primary Chow "
                "relations and excludes descendant unit relatives"
            )
        if not self.t_powers:
            raise ValueError("at least one Laurent coefficient is required")
        if len(set(self.t_powers)) != len(self.t_powers):
            raise ValueError("t_powers cannot contain duplicates")
        if self.max_genus is not None and self.max_genus < 0:
            raise ValueError("max_genus must be nonnegative")
        if self.max_ambient_degree is not None \
                and self.max_ambient_degree < 0:
            raise ValueError("max_ambient_degree must be nonnegative")
        if self.max_vertices < 1 or self.max_relations < 1:
            raise ValueError("resource bounds must be positive")
        if self.max_solve_rounds < 1:
            raise ValueError("max_solve_rounds must be positive")

    def stages(self):
        stages = [
            ProbeExpansionStage(markings, False)
            for markings in range(1, self.max_markings + 1)
        ]
        if self.include_unit_relatives:
            for unit_depth in range(1, self.max_unit_insertions + 1):
                stages.extend(
                    ProbeExpansionStage(markings, True, unit_depth)
                    for markings in range(1, self.max_markings + 1)
                )
        return tuple(stages)

    def to_record(self):
        return {
            "max_markings": int(self.max_markings),
            "t_powers": [int(power) for power in self.t_powers],
            "include_unit_relatives": self.include_unit_relatives,
            "max_unit_insertions": int(self.max_unit_insertions),
            "include_chow_relations": self.include_chow_relations,
            "max_chow_unit_insertions": int(
                self.max_chow_unit_insertions
            ),
            "chow_primary_only": self.chow_primary_only,
            "effective_basis_only": self.effective_basis_only,
            "additional_probe_ambient_degrees": [
                int(degree)
                for degree in self.additional_probe_ambient_degrees
            ],
            "include_punctured_axioms": self.include_punctured_axioms,
            "max_genus": (
                None if self.max_genus is None else int(self.max_genus)
            ),
            "max_ambient_degree": (
                None if self.max_ambient_degree is None
                else int(self.max_ambient_degree)
            ),
            "max_vertices": int(self.max_vertices),
            "max_relations": int(self.max_relations),
            "max_solve_rounds": int(self.max_solve_rounds),
            "require_complete": self.require_complete,
            "fail_on_unsupported": self.fail_on_unsupported,
            "include_kernel_basis": self.include_kernel_basis,
        }

    @classmethod
    def from_record(cls, record):
        return cls(**record)
