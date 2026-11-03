"""Deterministic scientific unit normalization for chemical concentrations and exposure durations.

Enforces:
1. Exact unit conversion only when dimensions and units are unambiguous and scientifically compatible.
2. Conversion between molar concentration units (M, mM, uM, nM, pM).
3. Conversion between mass concentration units (mg/L, ug/mL, etc.).
4. Strictly PROHIBITS converting between mass concentration (e.g. ug/mL) and molarity (e.g. uM)
   without an explicit, verified molecular weight.
5. Conversion between standard time units (s, min, h, day, week).
6. Preserves original values and units without alteration.
7. Emits structured status codes ('CONVERTED', 'IDENTICAL', 'INCOMPATIBLE_DIMENSIONS', 'UNKNOWN_UNIT').
"""

from dataclasses import dataclass
from enum import StrEnum
import math
import re
from typing import Any


class UnitConversionStatus(StrEnum):
    IDENTICAL = "IDENTICAL"
    CONVERTED = "CONVERTED"
    INCOMPATIBLE_DIMENSIONS = "INCOMPATIBLE_DIMENSIONS"
    UNKNOWN_UNIT = "UNKNOWN_UNIT"
    INVALID_VALUE = "INVALID_VALUE"


@dataclass(frozen=True, slots=True)
class NormalizedDose:
    original_value: float | str | None
    original_unit: str | None
    normalized_value: float | None
    normalized_unit: str | None
    conversion_factor: float | None
    status: UnitConversionStatus
    message: str


@dataclass(frozen=True, slots=True)
class NormalizedDuration:
    original_value: float | str | None
    original_unit: str | None
    normalized_value: float | None
    normalized_unit: str | None
    conversion_factor: float | None
    status: UnitConversionStatus
    message: str


# Canonical reference unit for molar concentration is micromolar (uM)
_MOLAR_TO_UM: dict[str, float] = {
    "m": 1_000_000.0,
    "mol/l": 1_000_000.0,
    "molar": 1_000_000.0,
    "mm": 1_000.0,
    "mmol/l": 1_000.0,
    "millimolar": 1_000.0,
    "um": 1.0,
    "µm": 1.0,
    "\u00b5m": 1.0,  # micro sign
    "\u03bcm": 1.0,  # greek small letter mu
    "umol/l": 1.0,
    "µmol/l": 1.0,
    "micromolar": 1.0,
    "nm": 0.001,
    "nmol/l": 0.001,
    "nanomolar": 0.001,
    "pm": 0.000001,
    "pmol/l": 0.000001,
    "picomolar": 0.000001,
}

# Mass concentration normalized to ug/mL (= mg/L)
_MASS_TO_UG_PER_ML: dict[str, float] = {
    "mg/l": 1.0,
    "ug/ml": 1.0,
    "µg/ml": 1.0,
    "ug/l": 0.001,
    "µg/l": 0.001,
    "ng/ml": 0.001,
    "mg/ml": 1000.0,
    "g/l": 1000.0,
}

# Canonical reference unit for exposure duration is hour (h)
_TIME_TO_HOURS: dict[str, float] = {
    "h": 1.0,
    "hr": 1.0,
    "hrs": 1.0,
    "hour": 1.0,
    "hours": 1.0,
    "d": 24.0,
    "day": 24.0,
    "days": 24.0,
    "w": 168.0,
    "wk": 168.0,
    "wks": 168.0,
    "week": 168.0,
    "weeks": 168.0,
    "m": 1.0 / 60.0,
    "min": 1.0 / 60.0,
    "mins": 1.0 / 60.0,
    "minute": 1.0 / 60.0,
    "minutes": 1.0 / 60.0,
    "s": 1.0 / 3600.0,
    "sec": 1.0 / 3600.0,
    "second": 1.0 / 3600.0,
    "seconds": 1.0 / 3600.0,
}


def _clean_unit_string(u: str | None) -> str:
    if not u:
        return ""
    clean = u.strip().lower()
    # Normalize Greek mu or micro sign to 'u'
    clean = clean.replace("µ", "u").replace("\u00b5", "u").replace("\u03bc", "u")
    return clean


def parse_numeric_value(val: Any) -> float | None:
    """Parse numeric float from float, int, or string, returning None if unparseable."""
    if val is None:
        return None
    if isinstance(val, (int, float)):
        if math.isnan(val) or math.isinf(val):
            return None
        return float(val)
    if isinstance(val, str):
        s = val.strip().replace(",", "")
        # Match leading float
        m = re.match(r"^[-+]?([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)", s)
        if m:
            try:
                f = float(m.group(1))
                return f if not (math.isnan(f) or math.isinf(f)) else None
            except ValueError:
                return None
    return None


def normalize_concentration(
    value: Any,
    unit: str | None,
    target_unit: str = "uM",
) -> NormalizedDose:
    """Normalize chemical concentration/dose safely.

    Default target_unit is 'uM' (micromolar) for molar units.
    Mass units are normalized to 'ug/mL'.
    Cross-dimensional conversion (mass to molar) is blocked and flagged as INCOMPATIBLE_DIMENSIONS.
    """
    orig_val_raw = value
    orig_unit_raw = unit

    if value is None or (isinstance(value, str) and not value.strip()):
        return NormalizedDose(
            original_value=orig_val_raw,
            original_unit=orig_unit_raw,
            normalized_value=None,
            normalized_unit=None,
            conversion_factor=None,
            status=UnitConversionStatus.INVALID_VALUE,
            message="Missing or empty concentration value.",
        )

    parsed_num = parse_numeric_value(value)
    if parsed_num is None:
        return NormalizedDose(
            original_value=orig_val_raw,
            original_unit=orig_unit_raw,
            normalized_value=None,
            normalized_unit=None,
            conversion_factor=None,
            status=UnitConversionStatus.INVALID_VALUE,
            message=f"Cannot parse numeric concentration from '{value}'.",
        )

    if not unit or not str(unit).strip():
        # Unit missing: preserve value, cannot normalize unit
        return NormalizedDose(
            original_value=orig_val_raw,
            original_unit=orig_unit_raw,
            normalized_value=parsed_num,
            normalized_unit=None,
            conversion_factor=1.0,
            status=UnitConversionStatus.UNKNOWN_UNIT,
            message="Concentration unit is missing; normalized unit cannot be established.",
        )

    clean_u = _clean_unit_string(unit)
    target_clean = _clean_unit_string(target_unit)

    # 1. Check if already in target unit
    if clean_u == target_clean:
        return NormalizedDose(
            original_value=orig_val_raw,
            original_unit=orig_unit_raw,
            normalized_value=parsed_num,
            normalized_unit=target_unit,
            conversion_factor=1.0,
            status=UnitConversionStatus.IDENTICAL,
            message=f"Value already in target unit '{target_unit}'.",
        )

    # 2. Check molar concentration
    if clean_u in _MOLAR_TO_UM:
        factor_to_um = _MOLAR_TO_UM[clean_u]
        val_in_um = parsed_num * factor_to_um

        if target_clean in _MOLAR_TO_UM:
            factor_target_to_um = _MOLAR_TO_UM[target_clean]
            converted_val = val_in_um / factor_target_to_um
            overall_factor = factor_to_um / factor_target_to_um
            return NormalizedDose(
                original_value=orig_val_raw,
                original_unit=orig_unit_raw,
                normalized_value=converted_val,
                normalized_unit=target_unit,
                conversion_factor=overall_factor,
                status=UnitConversionStatus.CONVERTED,
                message=f"Converted {parsed_num} {unit} -> {converted_val:.6g} {target_unit}.",
            )
        else:
            # Target is mass or incompatible
            return NormalizedDose(
                original_value=orig_val_raw,
                original_unit=orig_unit_raw,
                normalized_value=val_in_um,
                normalized_unit="uM",
                conversion_factor=factor_to_um,
                status=UnitConversionStatus.INCOMPATIBLE_DIMENSIONS,
                message=f"Source unit '{unit}' is molar, but target '{target_unit}' is incompatible. Normalized to uM.",
            )

    # 3. Check mass concentration
    if clean_u in _MASS_TO_UG_PER_ML:
        factor_to_ug_ml = _MASS_TO_UG_PER_ML[clean_u]
        val_in_ug_ml = parsed_num * factor_to_ug_ml

        if target_clean in _MASS_TO_UG_PER_ML:
            factor_target = _MASS_TO_UG_PER_ML[target_clean]
            converted_val = val_in_ug_ml / factor_target
            overall_factor = factor_to_ug_ml / factor_target
            return NormalizedDose(
                original_value=orig_val_raw,
                original_unit=orig_unit_raw,
                normalized_value=converted_val,
                normalized_unit=target_unit,
                conversion_factor=overall_factor,
                status=UnitConversionStatus.CONVERTED,
                message=f"Converted {parsed_num} {unit} -> {converted_val:.6g} {target_unit}.",
            )
        else:
            # Cannot convert mass concentration to molar without molecular weight!
            return NormalizedDose(
                original_value=orig_val_raw,
                original_unit=orig_unit_raw,
                normalized_value=val_in_ug_ml,
                normalized_unit="ug/mL",
                conversion_factor=factor_to_ug_ml,
                status=UnitConversionStatus.INCOMPATIBLE_DIMENSIONS,
                message=f"Cannot convert mass concentration '{unit}' to molar unit '{target_unit}' without molecular weight. Preserved in ug/mL.",
            )

    # 4. Unknown unit
    return NormalizedDose(
        original_value=orig_val_raw,
        original_unit=orig_unit_raw,
        normalized_value=parsed_num,
        normalized_unit=None,
        conversion_factor=None,
        status=UnitConversionStatus.UNKNOWN_UNIT,
        message=f"Unrecognized concentration unit '{unit}'. Original value preserved without unit conversion.",
    )


def normalize_duration(
    value: Any,
    unit: str | None,
    target_unit: str = "h",
) -> NormalizedDuration:
    """Normalize exposure duration safely.

    Default target_unit is 'h' (hours).
    """
    orig_val_raw = value
    orig_unit_raw = unit

    if value is None or (isinstance(value, str) and not value.strip()):
        return NormalizedDuration(
            original_value=orig_val_raw,
            original_unit=orig_unit_raw,
            normalized_value=None,
            normalized_unit=None,
            conversion_factor=None,
            status=UnitConversionStatus.INVALID_VALUE,
            message="Missing or empty exposure duration value.",
        )

    parsed_num = parse_numeric_value(value)
    if parsed_num is None:
        return NormalizedDuration(
            original_value=orig_val_raw,
            original_unit=orig_unit_raw,
            normalized_value=None,
            normalized_unit=None,
            conversion_factor=None,
            status=UnitConversionStatus.INVALID_VALUE,
            message=f"Cannot parse numeric duration from '{value}'.",
        )

    if not unit or not str(unit).strip():
        return NormalizedDuration(
            original_value=orig_val_raw,
            original_unit=orig_unit_raw,
            normalized_value=parsed_num,
            normalized_unit=None,
            conversion_factor=1.0,
            status=UnitConversionStatus.UNKNOWN_UNIT,
            message="Duration unit is missing; normalized unit cannot be established.",
        )

    clean_u = str(unit).strip().lower()
    target_clean = str(target_unit).strip().lower()

    if clean_u == target_clean:
        return NormalizedDuration(
            original_value=orig_val_raw,
            original_unit=orig_unit_raw,
            normalized_value=parsed_num,
            normalized_unit=target_unit,
            conversion_factor=1.0,
            status=UnitConversionStatus.IDENTICAL,
            message=f"Duration already in target unit '{target_unit}'.",
        )

    if clean_u in _TIME_TO_HOURS:
        factor_to_h = _TIME_TO_HOURS[clean_u]
        val_in_h = parsed_num * factor_to_h

        if target_clean in _TIME_TO_HOURS:
            factor_target = _TIME_TO_HOURS[target_clean]
            converted_val = val_in_h / factor_target
            overall_factor = factor_to_h / factor_target
            return NormalizedDuration(
                original_value=orig_val_raw,
                original_unit=orig_unit_raw,
                normalized_value=converted_val,
                normalized_unit=target_unit,
                conversion_factor=overall_factor,
                status=UnitConversionStatus.CONVERTED,
                message=f"Converted {parsed_num} {unit} -> {converted_val:.6g} {target_unit}.",
            )

    return NormalizedDuration(
        original_value=orig_val_raw,
        original_unit=orig_unit_raw,
        normalized_value=parsed_num,
        normalized_unit=None,
        conversion_factor=None,
        status=UnitConversionStatus.UNKNOWN_UNIT,
        message=f"Unrecognized duration unit '{unit}'. Original value preserved without conversion.",
    )
