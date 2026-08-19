"""Pure-Python baseline for the CJR zero-vertex classification hot loop."""

from dataclasses import dataclass
import argparse
import time


@dataclass(frozen=True)
class Case:
    genera: tuple[int, ...]
    degrees: tuple[int, ...]
    edge_zeros: tuple[int, ...]
    marking_zeros: tuple[int, ...]


def cases(count: int) -> tuple[Case, ...]:
    state = 0xC0FFEE

    def draw() -> int:
        nonlocal state
        state = (1664525 * state + 1013904223) & 0xFFFFFFFF
        return state

    answer = []
    for _ in range(count):
        zero_count = 1 + draw() % 8
        genera = tuple(1 if draw() % 17 == 0 else 0
                       for _ in range(zero_count))
        degrees = tuple(1 if draw() % 19 == 0 else 0
                        for _ in range(zero_count))
        edge_zeros = tuple(draw() % zero_count
                           for _ in range(1 + draw() % 12))
        marking_zeros = tuple(draw() % zero_count
                              for _ in range(draw() % 5))
        answer.append(Case(genera, degrees, edge_zeros, marking_zeros))
    return tuple(answer)


def classify(case: Case) -> int:
    zero_count = len(case.genera)
    edge_valences = [0] * zero_count
    marking_valences = [0] * zero_count
    for zero in case.edge_zeros:
        edge_valences[zero] += 1
    for zero in case.marking_zeros:
        marking_valences[zero] += 1
    checksum = 0
    for zero in range(zero_count):
        edge_valence = edge_valences[zero]
        marking_valence = marking_valences[zero]
        valence = edge_valence + marking_valence
        if case.genera[zero] or case.degrees[zero] or valence > 2:
            kind = 1
        elif edge_valence == 1 and marking_valence == 0:
            kind = 2
        elif edge_valence == 1 and marking_valence == 1:
            kind = 3
        elif edge_valence == 2 and marking_valence == 0:
            kind = 4
        else:
            kind = 5
        checksum += kind * (zero + 1)
    return checksum


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=4096)
    parser.add_argument("--rounds", type=int, default=500)
    arguments = parser.parse_args()
    workload = cases(arguments.cases)
    started = time.perf_counter()
    checksum = 0
    for _ in range(arguments.rounds):
        for case in workload:
            checksum += classify(case)
    elapsed = time.perf_counter() - started
    calls = arguments.cases * arguments.rounds
    print(f"language=python calls={calls} seconds={elapsed:.6f} "
          f"ns_per_call={elapsed * 1e9 / calls:.1f} checksum={checksum}")


if __name__ == "__main__":
    main()
