from enum import StrEnum, auto


class Unit(StrEnum):
    Bp = auto()
    Kbp = auto()
    Mbp = auto()

    def convert_value(self, value: int | float, round_to: int = 3) -> float:
        if self == Unit.Bp:
            new_value = value
        elif self == Unit.Kbp:
            new_value = value / 1_000
        else:
            new_value = value / 1_000_000

        return round(new_value, round_to)
