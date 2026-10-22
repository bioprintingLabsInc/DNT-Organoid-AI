"""Approved 12-run pilot manifest and metadata validation for PRJNA802072."""

from dataclasses import dataclass

BIOLOGICAL_REPLICATE_UNIT = "UNVERIFIED"
TARGET_BIOPROJECT = "PRJNA802072"
TARGET_SRA_STUDY = "SRP357550"
TARGET_DEVELOPMENTAL_DAY = 28
TARGET_ORGANISM = "Homo sapiens"

@dataclass(frozen=True, slots=True)
class PilotSampleRecord:
    sample_title: str
    run_accession: str
    experiment_accession: str
    biosample_accession: str
    condition_group: str  # "control" or "treatment"
    vpa_dose_um: float
    developmental_day: int
    paired_reads_expected: int
    total_bases_expected: int
    fastq_1_url: str
    fastq_2_url: str
    fastq_1_md5: str
    fastq_2_md5: str
    biological_replicate_unit: str = BIOLOGICAL_REPLICATE_UNIT


# Authoritative live ENA/SRA records for Day 28 VPA 250 uM vs 0 uM vehicle control
APPROVED_PILOT_RUNS: tuple[PilotSampleRecord, ...] = (
    # Controls (0 uM, N=6)
    PilotSampleRecord(
        sample_title="VPA_0_1_1",
        run_accession="SRR17834920",
        experiment_accession="SRX13996121",
        biosample_accession="SAMN25510021",
        condition_group="control",
        vpa_dose_um=0.0,
        developmental_day=28,
        paired_reads_expected=21899521,
        total_bases_expected=6569856300,
        fastq_1_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/020/SRR17834920/SRR17834920_1.fastq.gz",
        fastq_2_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/020/SRR17834920/SRR17834920_2.fastq.gz",
        fastq_1_md5="d0a584f6e89d2418cbe624f4f0e7a243",
        fastq_2_md5="95cfbb92e05bfe9efcb2f2bc63f15513",
    ),
    PilotSampleRecord(
        sample_title="VPA_0_2_1",
        run_accession="SRR17834908",
        experiment_accession="SRX13996133",
        biosample_accession="SAMN25510023",
        condition_group="control",
        vpa_dose_um=0.0,
        developmental_day=28,
        paired_reads_expected=23152105,
        total_bases_expected=6945631500,
        fastq_1_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/008/SRR17834908/SRR17834908_1.fastq.gz",
        fastq_2_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/008/SRR17834908/SRR17834908_2.fastq.gz",
        fastq_1_md5="824923b1fb97fc1ceda41079dbbbb718",
        fastq_2_md5="d7bd53ab567bf090cc83194692ccaee7",
    ),
    PilotSampleRecord(
        sample_title="VPA_0_3_1",
        run_accession="SRR17834897",
        experiment_accession="SRX13996144",
        biosample_accession="SAMN25510025",
        condition_group="control",
        vpa_dose_um=0.0,
        developmental_day=28,
        paired_reads_expected=22237894,
        total_bases_expected=6671368200,
        fastq_1_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/097/SRR17834897/SRR17834897_1.fastq.gz",
        fastq_2_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/097/SRR17834897/SRR17834897_2.fastq.gz",
        fastq_1_md5="491209dcb939f7402af897badf5bb048",
        fastq_2_md5="57240ed81485018318e2349d108e3d73",
    ),
    PilotSampleRecord(
        sample_title="VPA_0_4_1",
        run_accession="SRR17834886",
        experiment_accession="SRX13996155",
        biosample_accession="SAMN25510027",
        condition_group="control",
        vpa_dose_um=0.0,
        developmental_day=28,
        paired_reads_expected=22600697,
        total_bases_expected=6780209100,
        fastq_1_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/086/SRR17834886/SRR17834886_1.fastq.gz",
        fastq_2_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/086/SRR17834886/SRR17834886_2.fastq.gz",
        fastq_1_md5="fe11fc754f72c19802c76ade0705b8bb",
        fastq_2_md5="2ca2665fb98b0d25a09f995596ff0cc5",
    ),
    PilotSampleRecord(
        sample_title="VPA_0_5_1",
        run_accession="SRR17834875",
        experiment_accession="SRX13996166",
        biosample_accession="SAMN25510029",
        condition_group="control",
        vpa_dose_um=0.0,
        developmental_day=28,
        paired_reads_expected=22052958,
        total_bases_expected=6615887400,
        fastq_1_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/075/SRR17834875/SRR17834875_1.fastq.gz",
        fastq_2_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/075/SRR17834875/SRR17834875_2.fastq.gz",
        fastq_1_md5="75fe7bc4a03bad5d38b8af2bba622872",
        fastq_2_md5="2584255212e8ba3ff54566bf53fea47e",
    ),
    PilotSampleRecord(
        sample_title="VPA_0_6_1",
        run_accession="SRR17834919",
        experiment_accession="SRX13996122",
        biosample_accession="SAMN25510031",
        condition_group="control",
        vpa_dose_um=0.0,
        developmental_day=28,
        paired_reads_expected=23085831,
        total_bases_expected=6925749300,
        fastq_1_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/019/SRR17834919/SRR17834919_1.fastq.gz",
        fastq_2_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/019/SRR17834919/SRR17834919_2.fastq.gz",
        fastq_1_md5="9a0790588e9ccacfacec62bc0645c250",
        fastq_2_md5="c24299c2fe123577bd8476f47c0cfe1a",
    ),
    # Treated (250 uM, N=6)
    PilotSampleRecord(
        sample_title="VPA_250_1_1",
        run_accession="SRR17834905",
        experiment_accession="SRX13996136",
        biosample_accession="SAMN25510045",
        condition_group="treatment",
        vpa_dose_um=250.0,
        developmental_day=28,
        paired_reads_expected=22634745,
        total_bases_expected=6790423500,
        fastq_1_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/005/SRR17834905/SRR17834905_1.fastq.gz",
        fastq_2_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/005/SRR17834905/SRR17834905_2.fastq.gz",
        fastq_1_md5="24b0074eac52303827cf315395419daa",
        fastq_2_md5="d385acfe7134fd7bd8f09e6a779b37a5",
    ),
    PilotSampleRecord(
        sample_title="VPA_250_2_1",
        run_accession="SRR17834904",
        experiment_accession="SRX13996137",
        biosample_accession="SAMN25510047",
        condition_group="treatment",
        vpa_dose_um=250.0,
        developmental_day=28,
        paired_reads_expected=20783074,
        total_bases_expected=6234922200,
        fastq_1_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/004/SRR17834904/SRR17834904_1.fastq.gz",
        fastq_2_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/004/SRR17834904/SRR17834904_2.fastq.gz",
        fastq_1_md5="e9b9378c782e2167233b7c48048b3487",
        fastq_2_md5="a4c87f381ad8f78fde0dc9b0bb88ebba",
    ),
    PilotSampleRecord(
        sample_title="VPA_250_3_1",
        run_accession="SRR17834903",
        experiment_accession="SRX13996138",
        biosample_accession="SAMN25510049",
        condition_group="treatment",
        vpa_dose_um=250.0,
        developmental_day=28,
        paired_reads_expected=21073965,
        total_bases_expected=6322189500,
        fastq_1_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/003/SRR17834903/SRR17834903_1.fastq.gz",
        fastq_2_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/003/SRR17834903/SRR17834903_2.fastq.gz",
        fastq_1_md5="b996ee69ae7803f2fbd725dbce8cede0",
        fastq_2_md5="bbc6810d91b6e6d553e419c5e51a70a8",
    ),
    PilotSampleRecord(
        sample_title="VPA_250_4_1",
        run_accession="SRR17834902",
        experiment_accession="SRX13996139",
        biosample_accession="SAMN25510051",
        condition_group="treatment",
        vpa_dose_um=250.0,
        developmental_day=28,
        paired_reads_expected=21110376,
        total_bases_expected=6333112800,
        fastq_1_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/002/SRR17834902/SRR17834902_1.fastq.gz",
        fastq_2_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/002/SRR17834902/SRR17834902_2.fastq.gz",
        fastq_1_md5="fb9e2d43e507d643786ff344cc8a3f22",
        fastq_2_md5="a45630402124ef0cb8699eb2a9c0e4ec",
    ),
    PilotSampleRecord(
        sample_title="VPA_250_5_1",
        run_accession="SRR17834901",
        experiment_accession="SRX13996140",
        biosample_accession="SAMN25510053",
        condition_group="treatment",
        vpa_dose_um=250.0,
        developmental_day=28,
        paired_reads_expected=20605954,
        total_bases_expected=6181786200,
        fastq_1_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/001/SRR17834901/SRR17834901_1.fastq.gz",
        fastq_2_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/001/SRR17834901/SRR17834901_2.fastq.gz",
        fastq_1_md5="dd3830d6a8a23413d925e85ec89a45c7",
        fastq_2_md5="7d769c4ee443ad0f2e62cc1611c2c7d6",
    ),
    PilotSampleRecord(
        sample_title="VPA_250_6_1",
        run_accession="SRR17834900",
        experiment_accession="SRX13996141",
        biosample_accession="SAMN25510055",
        condition_group="treatment",
        vpa_dose_um=250.0,
        developmental_day=28,
        paired_reads_expected=22878898,
        total_bases_expected=6863669400,
        fastq_1_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/000/SRR17834900/SRR17834900_1.fastq.gz",
        fastq_2_url="ftp.sra.ebi.ac.uk/vol1/fastq/SRR178/000/SRR17834900/SRR17834900_2.fastq.gz",
        fastq_1_md5="f5d4ac1e4c932b90ee0bf6d5bb125a28",
        fastq_2_md5="641ab146209ce8c112bb026526ce9885",
    ),
)

APPROVED_RUN_ACCESSIONS: tuple[str, ...] = tuple(r.run_accession for r in APPROVED_PILOT_RUNS)
RUN_TO_SAMPLE_MAP: dict[str, PilotSampleRecord] = {r.run_accession: r for r in APPROVED_PILOT_RUNS}
