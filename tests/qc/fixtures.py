from pathlib import Path
from tests.metadata.fixtures import canonical_bulk_study,reported,with_provenance

def write_matrix(path:Path,text=None):
    path.write_text(text or "gene\tsample_c1\tsample_c2\tsample_c3\tsample_t1\tsample_t2\tsample_t3\nG1\t1\t2\t3\t4\t5\t6\nG2\t0\t1\t2\t3\t4\t5\nG3\t2\t2\t2\t2\t2\t2\n",encoding="utf-8"); return path

def metadata(path:Path):
    d=canonical_bulk_study(); a=d["input_data_assets"][0]; a["source_file_identifier_or_reference"]=reported(str(path)); return with_provenance(d)
