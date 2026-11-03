"""Unit tests for scientific unit normalization of concentration and exposure duration."""

import math
import pytest

from src.metadata.units import (
    UnitConversionStatus,
    normalize_concentration,
    normalize_duration,
    parse_numeric_value,
)


class TestNumericParsing:
    def test_parse_valid_numbers(self):
        assert parse_numeric_value(10) == 10.0
        assert parse_numeric_value(10.5) == 10.5
        assert parse_numeric_value("200") == 200.0
        assert parse_numeric_value("  200.5 uM ") == 200.5
        assert parse_numeric_value("1,000.0") == 1000.0
        assert parse_numeric_value("1e-4") == 0.0001

    def test_parse_invalid_or_missing_numbers(self):
        assert parse_numeric_value(None) is None
        assert parse_numeric_value("") is None
        assert parse_numeric_value("not_a_number") is None
        assert parse_numeric_value(float("nan")) is None
        assert parse_numeric_value(float("inf")) is None


class TestConcentrationNormalization:
    def test_identical_molar_unit(self):
        res = normalize_concentration(200.0, "uM", target_unit="uM")
        assert res.status == UnitConversionStatus.IDENTICAL
        assert res.original_value == 200.0
        assert res.original_unit == "uM"
        assert res.normalized_value == 200.0
        assert res.normalized_unit == "uM"
        assert res.conversion_factor == 1.0

    def test_millimolar_to_micromolar(self):
        # 1 mM = 1000 uM
        res = normalize_concentration(0.2, "mM", target_unit="uM")
        assert res.status == UnitConversionStatus.CONVERTED
        assert res.original_value == 0.2
        assert res.original_unit == "mM"
        assert res.normalized_value == pytest.approx(200.0)
        assert res.normalized_unit == "uM"
        assert res.conversion_factor == pytest.approx(1000.0)

    def test_micromolar_variations(self):
        # µM, uM, micromolar, umol/L
        for u in ["µM", "\u00b5M", "micromolar", "umol/L", "uM"]:
            res = normalize_concentration(250, u, target_unit="uM")
            assert res.normalized_value == pytest.approx(250.0)
            assert res.normalized_unit == "uM"

    def test_nanomolar_to_micromolar(self):
        # 500 nM = 0.5 uM
        res = normalize_concentration(500, "nM", target_unit="uM")
        assert res.status == UnitConversionStatus.CONVERTED
        assert res.normalized_value == pytest.approx(0.5)
        assert res.normalized_unit == "uM"

    def test_molar_to_micromolar(self):
        # 0.001 M = 1000 uM
        res = normalize_concentration(0.001, "M", target_unit="uM")
        assert res.status == UnitConversionStatus.CONVERTED
        assert res.normalized_value == pytest.approx(1000.0)

    def test_mass_concentration_conversion(self):
        # mg/L = ug/mL
        res = normalize_concentration(10, "mg/L", target_unit="ug/mL")
        assert res.status == UnitConversionStatus.CONVERTED
        assert res.normalized_value == pytest.approx(10.0)
        assert res.normalized_unit == "ug/mL"

    def test_prohibition_of_mass_to_molar_cross_conversion(self):
        # Converting 10 ug/mL to uM must be blocked as INCOMPATIBLE_DIMENSIONS
        res = normalize_concentration(10, "ug/mL", target_unit="uM")
        assert res.status == UnitConversionStatus.INCOMPATIBLE_DIMENSIONS
        assert res.normalized_unit == "ug/mL"
        assert "without molecular weight" in res.message

    def test_unknown_unit_preserves_original(self):
        res = normalize_concentration(50, "puffs", target_unit="uM")
        assert res.status == UnitConversionStatus.UNKNOWN_UNIT
        assert res.original_value == 50
        assert res.original_unit == "puffs"
        assert res.normalized_value == 50.0
        assert res.normalized_unit is None

    def test_missing_unit_preserves_original(self):
        res = normalize_concentration(100, None, target_unit="uM")
        assert res.status == UnitConversionStatus.UNKNOWN_UNIT
        assert res.original_value == 100
        assert res.normalized_value == 100.0


class TestExposureDurationNormalization:
    def test_identical_hour_unit(self):
        res = normalize_duration(48, "h", target_unit="h")
        assert res.status == UnitConversionStatus.IDENTICAL
        assert res.normalized_value == 48.0
        assert res.normalized_unit == "h"

    def test_days_to_hours(self):
        # 2 days = 48 h
        res = normalize_duration(2, "days", target_unit="h")
        assert res.status == UnitConversionStatus.CONVERTED
        assert res.normalized_value == pytest.approx(48.0)
        assert res.normalized_unit == "h"
        assert res.conversion_factor == pytest.approx(24.0)

    def test_weeks_to_hours(self):
        # 1 week = 168 h
        res = normalize_duration(1, "week", target_unit="h")
        assert res.status == UnitConversionStatus.CONVERTED
        assert res.normalized_value == pytest.approx(168.0)
        assert res.normalized_unit == "h"

    def test_minutes_to_hours(self):
        # 90 minutes = 1.5 h
        res = normalize_duration(90, "min", target_unit="h")
        assert res.status == UnitConversionStatus.CONVERTED
        assert res.normalized_value == pytest.approx(1.5)
        assert res.normalized_unit == "h"

    def test_unknown_duration_unit(self):
        res = normalize_duration(3, "cycles", target_unit="h")
        assert res.status == UnitConversionStatus.UNKNOWN_UNIT
        assert res.original_value == 3
        assert res.original_unit == "cycles"
        assert res.normalized_value == 3.0
        assert res.normalized_unit is None
