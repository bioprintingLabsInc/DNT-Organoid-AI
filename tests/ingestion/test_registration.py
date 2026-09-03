"""Asset registration tests."""

import unittest
from copy import deepcopy

from metadata.fixtures import canonical_bulk_study
from src.ingestion import register_assets


class RegistrationTests(unittest.TestCase):
    def test_registration_preserves_identity_reference_and_modality(self):
        document = canonical_bulk_study()
        asset = register_assets(document)[0]
        self.assertEqual("counts_asset", asset.asset_id)
        self.assertEqual("counts.tsv", asset.source_reference)
        self.assertEqual("raw_counts", asset.declared_format)
        self.assertEqual(("bulk_rna_seq",), asset.modalities)

    def test_registration_does_not_mutate_metadata(self):
        document = canonical_bulk_study()
        before = deepcopy(document)
        register_assets(document)
        self.assertEqual(before, document)


if __name__ == "__main__":
    unittest.main()
