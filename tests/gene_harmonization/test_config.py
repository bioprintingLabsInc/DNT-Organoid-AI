"""Unit tests for gene harmonization configuration loading and validation."""

import tempfile
import unittest
from pathlib import Path
import yaml

from src.gene_harmonization.config import (
    DEFAULT_CONFIG_PATH,
    config_checksum,
    load_config,
)


class ConfigTests(unittest.TestCase):
    def test_default_config_loads_cleanly(self):
        cfg = load_config()
        self.assertEqual(cfg["config_version"], "0.1.0")
        self.assertEqual(cfg["stage_name"], "gene_harmonization")
        self.assertIn("reference", cfg)
        self.assertIn("supported_identifier_types", cfg)
        self.assertIn("structural_rules", cfg)

    def test_config_checksum_is_deterministic(self):
        cfg1 = load_config()
        cfg2 = load_config()
        cs1 = config_checksum(cfg1)
        cs2 = config_checksum(cfg2)
        self.assertEqual(cs1, cs2)
        self.assertEqual(len(cs1), 64)

    def test_missing_required_section_raises(self):
        cfg = load_config()
        del cfg["reference"]
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "bad_config.yaml"
            p.write_text(yaml.safe_dump(cfg), encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                load_config(p)
            self.assertIn("Missing required configuration sections", str(ctx.exception))

    def test_invalid_stage_name_raises(self):
        cfg = load_config()
        cfg["stage_name"] = "wrong_stage"
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "bad_config.yaml"
            p.write_text(yaml.safe_dump(cfg), encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                load_config(p)
            self.assertIn("Invalid stage_name", str(ctx.exception))

    def test_unknown_section_raises(self):
        cfg = load_config()
        cfg["unexpected_section"] = {}
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "bad_config.yaml"
            p.write_text(yaml.safe_dump(cfg), encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                load_config(p)
            self.assertIn("Unknown configuration sections", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
