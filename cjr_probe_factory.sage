r"""Exact rank-aware probe selection for infinity-vertex blocks."""

load("log_glsm_conventions.sage")
load("log_glsm_infinity_dp.sage")


class CompiledLocalizationRelation(SageObject):
    r"""A target-independent extracted relation ``known=twisted+terms``."""

    def __init__(self, probe, t_power, coefficient_field, known_gw,
                 twisted_zero_level=0, terms=(), compilation=None,
                 is_complete=None):
        self.probe = probe
        self.t_power = ZZ(t_power)
        self.coefficient_field = coefficient_field
        self.known_gw = coefficient_field(known_gw)
        self.twisted_zero_level = coefficient_field(twisted_zero_level)
        self.terms = tuple(
            (coefficient_field(coefficient), tuple(factors))
            for coefficient, factors in terms
        )
        self.compilation = compilation
        self._is_complete_override = (
            None if is_complete is None else bool(is_complete)
        )

    @property
    def is_complete(self):
        if self._is_complete_override is not None:
            return self._is_complete_override
        return self.compilation is None or self.compilation.is_complete

    def equation_for(self, target):
        return LocalizationEquation(
            target,
            self.known_gw,
            twisted_zero_level=self.twisted_zero_level,
            terms=self.terms,
            probe_label="%s [t^%s]" % (self.probe, self.t_power),
            coefficient_field=self.coefficient_field,
        )

    def to_record(self):
        """Serialize an extracted relation without its heavy graph objects."""
        return {
            "probe": self.probe.to_record(),
            "t_power": int(self.t_power),
            "known_gw": str(self.known_gw),
            "twisted_zero_level": str(self.twisted_zero_level),
            "is_complete": self.is_complete,
            "terms": [
                {
                    "coefficient": str(coefficient),
                    "factors": [factor.to_record() for factor in factors],
                }
                for coefficient, factors in self.terms
            ],
        }

    @classmethod
    def from_record(cls, record, coefficient_field=QQ):
        return cls(
            ProbeSpec.from_record(record["probe"]),
            record["t_power"],
            coefficient_field,
            known_gw=coefficient_field(record["known_gw"]),
            twisted_zero_level=coefficient_field(
                record.get("twisted_zero_level", "0")
            ),
            terms=tuple(
                (
                    coefficient_field(term["coefficient"]),
                    tuple(EffectiveVertex.from_record(factor)
                          for factor in term.get("factors", ())),
                )
                for term in record.get("terms", ())
            ),
            compilation=None,
            is_complete=record.get("is_complete", True),
        )


class ProbeBlockSelection(SageObject):
    """Independent relation rows selected for a same-rank target block."""

    def __init__(self, targets, relations, rows, coefficient_field,
                 rejected=()):
        self.targets = tuple(targets)
        self.relations = tuple(relations)
        self.rows = tuple(tuple(coefficient_field(value) for value in row)
                          for row in rows)
        self.coefficient_field = coefficient_field
        self.rejected = tuple(rejected)
        self.matrix = Matrix(coefficient_field, self.rows) if self.rows else Matrix(
            coefficient_field, 0, len(self.targets)
        )

    @property
    def rank(self):
        return self.matrix.rank()

    @property
    def is_full_rank(self):
        return self.rank == len(self.targets)

    def kernel_basis(self):
        return tuple(tuple(vector) for vector in self.matrix.right_kernel().basis())


class ProbeFactory(SageObject):
    """Generate stationary candidates and choose exact independent rows."""

    def __init__(self, coefficient_field=QQ):
        self.coefficient_field = coefficient_field

    def stationary_candidates(self, genus, ambient_degree, max_markings=4):
        genus = ZZ(genus)
        ambient_degree = ZZ(ambient_degree)
        max_markings = ZZ(max_markings)
        if genus < 1 or max_markings < 1:
            raise ValueError("stationary candidate generation needs g,n positive")
        total_psi = 2 * genus - 2
        candidates = []
        for markings in range(1, max_markings + 1):
            profiles = set(
                tuple(sorted(profile, reverse=True))
                for profile in IntegerVectors(total_psi, markings)
            )
            for profile in sorted(profiles, reverse=True):
                candidates.append(ProbeSpec.stationary(
                    genus, ambient_degree, profile,
                    label="stationary candidate",
                ))
        return tuple(candidates)

    def unit_relatives(self, probe):
        r"""Return dimension-zero string and dilaton relatives of a probe.

        For a stationary probe, adding ``tau_1(1)`` gives its dilaton
        relative.  Adding ``tau_0(1)`` and raising one point-descendant by
        one gives a string relative.  Their known elliptic-curve values are
        computed by :class:`EllipticProbeValueBackend`, while their CJR graph
        rows need not be dependent term by term.  They are therefore useful
        for separating contact-profile blocks.
        """
        if not isinstance(probe, ProbeSpec):
            raise TypeError("unit relatives require a ProbeSpec")
        if any(insertion.kind != "point" for insertion in probe.insertions):
            raise ValueError("unit relatives are generated from stationary probes")

        relatives = []
        for index, insertion in enumerate(probe.insertions):
            raised = EllipticInsertion("point", insertion.psi_power + 1)
            insertions = (
                (EllipticInsertion("unit", 0),)
                + probe.insertions[:index]
                + (raised,)
                + probe.insertions[index + 1:]
            )
            relatives.append(ProbeSpec(
                probe.genus, probe.ambient_degree, insertions,
                connected=probe.connected, label="string relative",
            ))
        relatives.append(ProbeSpec(
            probe.genus, probe.ambient_degree,
            (EllipticInsertion("unit", 1),) + probe.insertions,
            connected=probe.connected, label="dilaton relative",
        ))
        return tuple(relatives)

    def diagonal_row(self, relation, targets, lower_values=None):
        targets = tuple(targets)
        target_set = set(targets)
        lower_values = dict(lower_values or {})
        row = [relation.coefficient_field.zero()] * len(targets)
        nonlinear = False
        unresolved = set()
        for coefficient, factors in relation.terms:
            block_factors = [factor for factor in factors if factor in target_set]
            if len(block_factors) > 1:
                nonlinear = True
                continue
            multiplier = coefficient
            skipped = False
            for factor in factors:
                if factor in target_set and not skipped:
                    skipped = True
                    continue
                if factor not in lower_values:
                    unresolved.add(factor)
                    multiplier = None
                    break
                multiplier *= relation.coefficient_field(lower_values[factor])
            if multiplier is not None and block_factors:
                row[targets.index(block_factors[0])] += multiplier
        return tuple(row), nonlinear, tuple(sorted(
            unresolved, key=lambda item: item.order_key()
        ))

    def select_relations(self, targets, relations, lower_values=None):
        targets = tuple(targets)
        if not targets:
            raise ValueError("a probe block needs at least one target")
        selected_relations = []
        selected_rows = []
        rejected = []
        current_rank = 0
        for relation in relations:
            if relation.coefficient_field is not self.coefficient_field:
                raise TypeError("all candidate relations must use the factory field")
            if not relation.is_complete:
                rejected.append((relation, "incomplete graph compilation"))
                continue
            row, nonlinear, unresolved = self.diagonal_row(
                relation, targets, lower_values
            )
            if nonlinear:
                rejected.append((relation, "nonlinear in target block"))
                continue
            if unresolved:
                rejected.append((relation, "unresolved lower factors %s" % (unresolved,)))
                continue
            proposed = selected_rows + [row]
            proposed_rank = Matrix(self.coefficient_field, proposed).rank()
            if proposed_rank > current_rank:
                selected_relations.append(relation)
                selected_rows.append(row)
                current_rank = proposed_rank
            else:
                rejected.append((relation, "row is dependent"))
            if current_rank == len(targets):
                break
        return ProbeBlockSelection(
            targets, selected_relations, selected_rows,
            self.coefficient_field, rejected=rejected,
        )
