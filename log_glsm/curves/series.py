r"""curves / series."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable
import math

__all__ = ['Number']


Number = int | Fraction

def _as_fraction(value: Number | str) -> Fraction:
    if isinstance(value, Fraction):
        return value
    return Fraction(value)

def _harmonic(number: int) -> Fraction:
    return sum((Fraction(1, k) for k in range(1, number + 1)), Fraction(0))

def _divisor_power_sum(number: int, power: int) -> int:
    if number <= 0:
        raise ValueError("divisor power sums are defined here for positive integers")
    if power < 0:
        raise ValueError("the divisor power must be nonnegative")
    total = 0
    root = math.isqrt(number)
    for divisor in range(1, root + 1):
        if number % divisor == 0:
            total += divisor**power
            quotient = number // divisor
            if quotient != divisor:
                total += quotient**power
    return total

def _divisor_sum(number: int) -> int:
    return _divisor_power_sum(number, 1)

def _bernoulli_number(index: int) -> Fraction:
    """Return ``B_index`` with the convention ``B_1=-1/2``."""

    if type(index) is not int or index < 0:
        raise ValueError("the Bernoulli index must be a nonnegative integer")
    values = [Fraction(1)]
    for degree in range(1, index + 1):
        values.append(
            -sum(
                Fraction(math.comb(degree + 1, k)) * values[k]
                for k in range(degree)
            )
            / Fraction(degree + 1)
        )
    return values[index]

@dataclass(frozen=True)
class _Series:
    """A truncated power series over Q with a fixed inclusive order."""

    coefficients: tuple[Fraction, ...]

    def __init__(self, coefficients: Iterable[Number]) -> None:
        object.__setattr__(
            self, "coefficients", tuple(_as_fraction(value) for value in coefficients)
        )
        if not self.coefficients:
            raise ValueError("a truncated series needs at least its constant term")

    @property
    def order(self) -> int:
        return len(self.coefficients) - 1

    @classmethod
    def zero(cls, order: int) -> "_Series":
        return cls([0] * (order + 1))

    @classmethod
    def one(cls, order: int) -> "_Series":
        return cls([1] + [0] * order)

    @classmethod
    def variable(cls, order: int) -> "_Series":
        if order < 1:
            return cls.zero(order)
        return cls([0, 1] + [0] * (order - 1))

    def __getitem__(self, degree: int) -> Fraction:
        return self.coefficients[degree]

    def _check_order(self, other: "_Series") -> None:
        if self.order != other.order:
            raise ValueError("power series must have the same truncation order")

    def __add__(self, other: "_Series") -> "_Series":
        self._check_order(other)
        return _Series(a + b for a, b in zip(self.coefficients, other.coefficients))

    def __sub__(self, other: "_Series") -> "_Series":
        self._check_order(other)
        return _Series(a - b for a, b in zip(self.coefficients, other.coefficients))

    def scale(self, scalar: Number) -> "_Series":
        scalar = _as_fraction(scalar)
        return _Series(scalar * value for value in self.coefficients)

    def __mul__(self, other: "_Series") -> "_Series":
        self._check_order(other)
        answer = [Fraction(0)] * (self.order + 1)
        for left_degree, left in enumerate(self.coefficients):
            if not left:
                continue
            for right_degree in range(self.order - left_degree + 1):
                right = other[right_degree]
                if right:
                    answer[left_degree + right_degree] += left * right
        return _Series(answer)

    def reciprocal(self) -> "_Series":
        if not self[0]:
            raise ZeroDivisionError("a series with zero constant term is not invertible")
        answer = [1 / self[0]] + [Fraction(0)] * self.order
        for degree in range(1, self.order + 1):
            convolution = sum(
                self[k] * answer[degree - k] for k in range(1, degree + 1)
            )
            answer[degree] = -convolution / self[0]
        return _Series(answer)

    def __truediv__(self, other: "_Series") -> "_Series":
        return self * other.reciprocal()

    def exp(self) -> "_Series":
        if self[0]:
            raise ValueError("exact formal exp is implemented for zero constant term")
        answer = [Fraction(1)] + [Fraction(0)] * self.order
        for degree in range(1, self.order + 1):
            answer[degree] = sum(
                k * self[k] * answer[degree - k] for k in range(1, degree + 1)
            ) / degree
        return _Series(answer)

    def log(self) -> "_Series":
        if self[0] != 1:
            raise ValueError("exact formal log is implemented for constant term one")
        derivative = _Series(
            [(degree + 1) * self[degree + 1] for degree in range(self.order)] + [0]
        )
        quotient = derivative / self
        answer = [Fraction(0)] * (self.order + 1)
        for degree in range(1, self.order + 1):
            answer[degree] = quotient[degree - 1] / degree
        return _Series(answer)

    def compose(self, inner: "_Series") -> "_Series":
        self._check_order(inner)
        if inner[0]:
            raise ValueError("the inner series in a formal composition must vanish at zero")
        answer = _Series.zero(self.order)
        power = _Series.one(self.order)
        for coefficient in self.coefficients:
            if coefficient:
                answer = answer + power.scale(coefficient)
            power = power * inner
        return answer

    def shifted(self, amount: int = 1) -> "_Series":
        if amount < 0:
            raise ValueError("only nonnegative shifts are supported")
        if amount > self.order:
            return _Series.zero(self.order)
        return _Series([0] * amount + list(self.coefficients[: self.order + 1 - amount]))

    def as_strings(self) -> list[str]:
        return [str(value) for value in self.coefficients]

def _format_series(
    series: _Series,
    variable: str,
    *,
    start_degree: int = 0,
    include_zero: bool = False,
) -> str:
    pieces: list[tuple[str, str]] = []
    for degree in range(start_degree, series.order + 1):
        coefficient = series[degree]
        if not coefficient and not include_zero:
            continue
        sign = "-" if coefficient < 0 else "+"
        magnitude = abs(coefficient)
        if degree == 0:
            monomial = str(magnitude)
        else:
            variable_part = variable if degree == 1 else f"{variable}^{degree}"
            monomial = variable_part if magnitude == 1 else f"{magnitude}*{variable_part}"
        pieces.append((sign, monomial))
    if not pieces:
        return "0"
    first_sign, first = pieces[0]
    answer = ("-" if first_sign == "-" else "") + first
    for sign, monomial in pieces[1:]:
        answer += f" {sign} {monomial}"
    return answer

def _integer_partitions(number: int, maximum: int | None = None):
    """Yield integer partitions in weakly decreasing order."""

    if number == 0:
        yield ()
        return
    if maximum is None or maximum > number:
        maximum = number
    for first in range(maximum, 0, -1):
        for rest in _integer_partitions(number - first, first):
            yield (first,) + rest
