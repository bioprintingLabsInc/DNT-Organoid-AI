# Metadata Data Dictionary

## Contract conventions

`config/metadata_schema.yaml` is the authoritative, versioned metadata contract. It standardizes identity, relationships, explicit missingness, and provenance for public and in-house human brain-organoid datasets. It does not assert that RNA-seq modalities or source formats are biologically equivalent and does not define processing, QC, normalization, features, or DNT classification logic.

Requirement terms mean:

- **Required:** the field must be present. For a `value_record`, its scientific value may still be unavailable if `value_status` explicitly records why.
- **Conditionally required:** required when the stated condition applies; otherwise use `not_applicable` where the field is represented as a `value_record`.
- **Optional:** may be omitted, but reported source information should be preserved when available and appropriate.

### Shared structures

| Structure/field | Owner | Purpose | Type | Requirement | Permitted values / missingness | Example |
|---|---|---|---|---|---|---|
| `schema_version` | Contract | Versions the contract independently of data records. | String | Required | Semantic version | `0.2.0` |
| `contract_name` | Contract | Stable contract identifier. | String | Required | `dnt_organoid_experiment_metadata` | Same as permitted value |
| `value_record.original_value` | Any source-derived assertion | Preserves the exact value supplied by the source before normalization. | String, number, boolean, or null | Optional when not reported | Null is allowed only with an explicit non-`reported` status; empty strings and tokens such as `NA` remain literal original values until curated. | `10 uM` |
| `value_record.normalized_value` | Any source-derived assertion | Stores a future approved canonical representation without overwriting the source value. | String, number, boolean, or null | Optional | Null until normalization occurs or when no normalized value applies. | `10` |
| `value_record.value_status` | Any source-derived assertion | Distinguishes actual reporting from distinct absence states. | String | Required | `reported`, `missing_not_reported`, `unknown`, `not_applicable` | `reported` |
| `value_record.assertion_provenance_ids` | Any source-derived assertion | Links an assertion to one or more provenance records. | Array of strings | Required | May be empty only while provenance is unresolved and validation reflects that condition. | `[prov_gsm_1_dose]` |
| `metadata_validation_status` | Every metadata entity except provenance | Records the validation outcome without deleting invalid data. | String | Required | `not_validated`, `valid`, `valid_with_warnings`, `invalid` | `not_validated` |
| `accession.namespace` | Accession | Names the accession system. | String | Required | Open vocabulary; examples include `GSE`, `GSM`, `SRA`, `BioProject`, `ArrayExpress`, `internal`. | `GSE` |
| `accession.identifier` | Accession | Preserves the identifier within its namespace. | String | Required | Non-empty source identifier | `GSE12345` |
| `accession.url` | Accession | Provides a resolvable source link when available. | String or null | Optional | Null means no URL was recorded, not that the accession is unknown. | `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12345` |
| `evidence_scope_assertion.dimension` | DNT evidence scope assertion | Names the exposure-context dimension constrained by independent evidence. | String | Required for each assertion | `species`, `biological_source_type`, `organoid_type`, `brain_region_or_model_identity`, `developmental_age_or_stage_at_exposure`, `concentration_or_dose`, `concentration_or_dose_unit`, `exposure_duration`, `exposure_duration_unit`, `other` | `species` |
| `evidence_scope_assertion.value` | DNT evidence scope assertion | Preserves the source-supported scope value and its provenance without creating a new label. | `value_record` | Required for each assertion | Explicit four-state missingness | `Homo sapiens` |

Identifiers are non-empty strings, unique within their entity collection, and stable across revisions. Foreign keys must resolve to the named entity collection. The schema places no biological limits on ages, stages, doses, durations, replicate counts, or sequencing depth.

## Study

| Field | Purpose | Type | Requirement | Permitted values / missingness | Example |
|---|---|---|---|---|---|
| `study_id` | Stable study identity and study-level split grouping key. | String identifier | Required | No missing value | `study_GSE12345` |
| `study_title` | Human-readable source title. | `value_record` | Optional | Explicit four-state missingness | `Public cortical organoid exposure study` |
| `source_type` | Distinguishes public datasets from future in-house work. | String | Required | `public`, `in_house` | `public` |
| `public_accessions` | Collects study-level public identifiers, including equivalent repositories. | Array of `accession` | Conditionally required for public studies | Empty only when validation records an unresolved public accession | `[{namespace: GSE, identifier: GSE12345}]` |
| `metadata_validation_status` | Study metadata validation outcome. | String | Required | Shared validation vocabulary | `valid_with_warnings` |

## Experiment

| Field | Purpose | Type | Requirement | Permitted values / missingness | Example |
|---|---|---|---|---|---|
| `experiment_id` | Stable experiment identity. | String identifier | Required | No missing value | `exp_001` |
| `study_id` | Links the experiment to one study and supports study-level splitting. | Foreign-key string | Required | Existing `studies.study_id` | `study_GSE12345` |
| `internal_experiment_id` | Preserves a separate source/internal operational identifier when one exists; `experiment_id` remains the platform-controlled stable identifier. | `value_record` | Optional | Explicit four-state missingness; public studies need not invent one | `DNT-EXP-0001` |
| `experiment_title` | Human-readable experiment name from the source or internal record. | `value_record` | Optional | Explicit four-state missingness | `Compound A 10 uM exposure` |
| `experimental_batch_id` | Preserves a reported batch/grouping identity for future confounding assessment. | `value_record` | Optional | Explicit four-state missingness | `batch_03` |
| `metadata_validation_status` | Experiment metadata validation outcome. | String | Required | Shared validation vocabulary | `valid` |

## Condition

A condition is a scientifically defined experimental group shared by one or more samples or biological replicates. Distinct concentrations, durations, stages, or planned exposure combinations can be represented as distinct conditions without imposing which distinctions are scientifically meaningful.

| Field | Purpose | Type | Requirement | Permitted values / missingness | Example |
|---|---|---|---|---|---|
| `condition_id` | Stable identity for a planned experimental group and condition-level split/comparison key. | String identifier | Required | No missing value | `condition_vpa_10um_24h` |
| `experiment_id` | Links the condition to its experiment. | Foreign-key string | Required | Existing `experiments.experiment_id` | `exp_001` |
| `condition_label` | Preserves a human-readable source or curated group label. | `value_record` | Required | Explicit four-state missingness | `VPA 10 µM, 24 h` |
| `treatment_control_status` | Declares the planned group as treatment or control. | `value_record` | Required | `treatment`, `control` | `treatment` |
| `control_type` | Distinguishes control-group designs without treating them as equivalent. | `value_record` | Conditionally required for control conditions | `vehicle`, `untreated`, `positive_reference`, `negative_reference`, `other` | `vehicle` |
| `condition_notes` | Preserves group-level qualifications, factorial design notes, or source descriptions. | `value_record` | Optional | Explicit four-state missingness | `Shared vehicle control across three doses` |
| `metadata_validation_status` | Condition metadata validation outcome. | String | Required | Shared validation vocabulary | `valid` |

## Biological source / cell line

| Field | Purpose | Type | Requirement | Permitted values / missingness | Example |
|---|---|---|---|---|---|
| `biological_source_id` | Stable identity for a cell line, donor-derived source, or other biological source. | String identifier | Required | No missing value | `source_line_01` |
| `species` | Preserves reported species without imposing a taxonomy normalization yet. | `value_record` | Required | Explicit four-state missingness | `Homo sapiens` |
| `biological_sex` | Preserves biological sex when reported. | `value_record` | Optional | Normalized vocabulary: `female`, `male`, `mixed`, `intersex`, `other`; absence is represented only by `value_status`. Original wording remains preserved. | `female` |
| `source_material_type` | Identifies the biological source class. | `value_record` | Required | Normalized vocabulary: `iPSC`, `hESC`, `other`; details may remain in the original value. | `iPSC` |
| `cell_line_or_source_identifier` | Links samples derived from the same cell line or source. | `value_record` | Required | Explicit four-state missingness | `WTC11` |
| `donor_identifier` | Preserves a donor grouping key when available, appropriate, and permitted. | `value_record` | Optional | Explicit missingness; privacy-safe coded identifiers only | `donor_07` |
| `metadata_validation_status` | Biological-source metadata validation outcome. | String | Required | Shared validation vocabulary | `valid` |

## Organoid context

| Field | Purpose | Type | Requirement | Permitted values / missingness | Example |
|---|---|---|---|---|---|
| `organoid_context_id` | Stable identity for a reusable organoid context. | String identifier | Required | No missing value | `orgctx_cortical_v1` |
| `organoid_type` | Preserves the reported organoid type. | `value_record` | Required | Open vocabulary with explicit missingness; no types are declared equivalent. | `cerebral organoid` |
| `brain_region_or_model_identity` | Preserves regional or model identity separately from organoid type. | `value_record` | Required | Open vocabulary with explicit missingness | `dorsal forebrain` |
| `differentiation_protocol_identifier` | Identifies a reported or internal differentiation protocol. | `value_record` | Optional when available | Explicit four-state missingness | `protocol_DO_01` |
| `differentiation_protocol_version` | Distinguishes protocol revisions. | `value_record` | Optional when available | Explicit four-state missingness | `2.1` |
| `metadata_validation_status` | Organoid-context metadata validation outcome. | String | Required | Shared validation vocabulary | `valid` |

## Sample

| Field | Purpose | Type | Requirement | Permitted values / missingness | Example |
|---|---|---|---|---|---|
| `sample_id` | Stable identity for the biological/collected sample. | String identifier | Required | No missing value | `sample_001` |
| `experiment_id` | Assigns the sample to its experiment. | Foreign-key string | Required | Existing `experiments.experiment_id` | `exp_001` |
| `condition_id` | Assigns the sample to its scientifically defined planned condition. | Foreign-key string | Required | Existing `conditions.condition_id`; the condition must belong to the same experiment | `condition_vpa_10um_24h` |
| `public_accessions` | Preserves sample-level identifiers such as GSM or SRA accessions. | Array of `accession` | Optional | Empty when none are reported | `[{namespace: GSM, identifier: GSM123456}]` |
| `biological_source_id` | Links the sample to its cell line/source and donor grouping context. | Foreign-key string | Required | Existing `biological_sources.biological_source_id` | `source_line_01` |
| `organoid_context_id` | Links the sample to its organoid identity and protocol context. | Foreign-key string | Required | Existing `organoid_contexts.organoid_context_id` | `orgctx_cortical_v1` |
| `biological_replicate_id` | Preserves biological replicate grouping without prescribing replicate counts. | `value_record` | Required | Explicit four-state missingness | `bio_rep_2` |
| `technical_replicate_id` | Preserves technical replicate grouping where applicable. | `value_record` | Conditionally required when applicable | Use `not_applicable` when the sample has no technical-replicate concept. | `tech_rep_1` |
| `developmental_age_or_stage_at_collection` | Retains collection age/stage for developmental comparisons. | `value_record` | Required | Source wording and units preserved; no accepted range is imposed. | `day 60` |
| `collection_time_relative_to_exposure` | Records collection timing relative to exposure. | `value_record` | Conditionally required for exposed samples | Explicit four-state missingness | `24 h after exposure start` |
| `metadata_validation_status` | Sample metadata validation outcome. | String | Required | Shared validation vocabulary | `valid` |

## Exposure

Planned group-level exposures use `exposure_scope: condition_planned` and omit `sample_id`. Factorial or multi-agent conditions use multiple planned exposure records for the same condition. A sample-specific departure is an additional `sample_deviation` record; it does not overwrite or replace the condition's planned exposure.

| Field | Purpose | Type | Requirement | Permitted values / missingness | Example |
|---|---|---|---|---|---|
| `exposure_id` | Stable identity for one treatment or control exposure record. | String identifier | Required | No missing value | `exposure_001` |
| `condition_id` | Links the planned exposure or sample deviation to its owning condition. | Foreign-key string | Required | Existing `conditions.condition_id` | `condition_vpa_10um_24h` |
| `exposure_scope` | Distinguishes the planned condition definition from an explicit sample deviation. | String | Required | `condition_planned`, `sample_deviation` | `condition_planned` |
| `sample_id` | Identifies the deviating sample without changing its planned condition. | Foreign-key string or null | Conditionally required for `sample_deviation` | Existing `samples.sample_id` belonging to `condition_id`; null for planned exposure | `sample_007` |
| `deviation_from_exposure_id` | Identifies the planned exposure from which the sample deviated when a specific record applies. | Foreign-key string or null | Conditionally required for a sample deviation when applicable | Existing condition-planned `exposures.exposure_id` from the same condition | `exposure_001` |
| `deviation_reason` | Preserves what changed or why the sample departed from plan. | `value_record` | Conditionally required for `sample_deviation` | Explicit four-state missingness with provenance | `Exposure ended two hours early` |
| `agent_name` | Preserves the compound or other agent name for this planned exposure or deviation. | `value_record` | Conditionally required when an agent applies | Explicit four-state missingness | `valproic acid` |
| `agent_identifier` | Provides a stable compound/agent identity for cross-study grouping and future compound-level splitting. | `value_record` | Optional when available | Any documented identifier namespace/value; original retained | `CAS:99-66-1` |
| `vehicle` | Preserves the reported vehicle where applicable. | `value_record` | Conditionally required when applicable | Explicit four-state missingness | `DMSO` |
| `concentration_or_dose` | Preserves the reported numeric or textual dose independently of its unit. | `value_record` | Conditionally required when an agent applies | No accepted range is imposed | `10` |
| `concentration_or_dose_unit` | Preserves the dose unit needed to interpret dose. | `value_record` | Conditionally required when dose is reported | Open vocabulary pending approved unit normalization | `µM` |
| `exposure_start_time_or_stage` | Records when exposure began. | `value_record` | Conditionally required if reported | Explicit four-state missingness | `day 45` |
| `developmental_age_or_stage_at_exposure` | Retains developmental context at exposure. | `value_record` | Required | Source wording preserved; no accepted stages imposed | `day 45 organoid` |
| `exposure_duration` | Preserves reported exposure duration separately from its unit. | `value_record` | Conditionally required when exposure duration applies | No accepted range is imposed | `24` |
| `exposure_duration_unit` | Provides the unit for exposure duration. | `value_record` | Conditionally required when duration is reported | Open vocabulary pending approved unit normalization | `hour` |
| `washout_or_recovery_duration` | Preserves a post-exposure washout/recovery interval when applicable. | `value_record` | Conditionally required when applicable | `not_applicable` is distinct from unreported | `48` |
| `washout_or_recovery_duration_unit` | Provides the washout/recovery duration unit. | `value_record` | Conditionally required when that duration is reported | Open vocabulary pending approved unit normalization | `hour` |
| `metadata_validation_status` | Exposure metadata validation outcome. | String | Required | Shared validation vocabulary | `valid_with_warnings` |

## Treatment/control relationship

| Field | Purpose | Type | Requirement | Permitted values / missingness | Example |
|---|---|---|---|---|---|
| `relationship_id` | Stable identity for an explicit condition-level matching assertion. | String identifier | Required | No missing value | `match_001` |
| `treatment_condition_ids` | Identifies one or more treatment conditions in the comparison. | Array of foreign-key strings | Required | Existing treatment `conditions.condition_id`; non-empty | `[condition_vpa_10um_24h]` |
| `matched_control_condition_ids` | Identifies scientifically appropriate matched control conditions; one control may be shared by multiple comparisons. | Array of foreign-key strings | Required | Existing control `conditions.condition_id`; non-empty | `[condition_vehicle_24h]` |
| `matching_basis` | Preserves why the conditions are considered matched; validation code must not invent the match. | `value_record` | Required | Open text/source assertion with provenance | `same line, batch, stage, duration, and vehicle` |
| `relationship_notes` | Holds additional limitations or context. | `value_record` | Optional | Explicit four-state missingness | `shared vehicle-control pool` |
| `metadata_validation_status` | Relationship validation outcome. | String | Required | Shared validation vocabulary | `not_validated` |

## Sample comparison exception

These records document exceptions to a condition-level comparison. They do not silently change condition membership or redefine the primary comparison.

| Field | Purpose | Type | Requirement | Permitted values / missingness | Example |
|---|---|---|---|---|---|
| `exception_id` | Stable identity for a sample-level comparison exception. | String identifier | Required | No missing value | `comparison_exception_001` |
| `relationship_id` | Identifies the condition-level comparison being qualified. | Foreign-key string | Required | Existing `treatment_control_relationships.relationship_id` | `match_001` |
| `affected_sample_ids` | Identifies samples that cannot follow the primary condition-level comparison as written. | Array of foreign-key strings | Required | Existing `samples.sample_id`; non-empty | `[sample_t_03]` |
| `replacement_control_sample_ids` | Records specifically approved replacement control samples when applicable. | Array of foreign-key strings | Optional when applicable | Existing `samples.sample_id`; omission does not imply a replacement | `[sample_c_02]` |
| `exception_reason` | Preserves the explicit scientific or operational reason for the exception. | `value_record` | Required | Explicit four-state missingness with provenance | `Only same-batch control is sample_c_02` |
| `metadata_validation_status` | Exception validation outcome. | String | Required | Shared validation vocabulary | `not_validated` |

## Sequencing assay

| Field | Purpose | Type | Requirement | Permitted values / missingness | Example |
|---|---|---|---|---|---|
| `assay_id` | Stable identity for a sequencing assay. | String identifier | Required | No missing value | `assay_001` |
| `sample_ids` | Links one assay to one or more represented samples. | Array of foreign-key strings | Required | Existing `samples.sample_id`; non-empty | `[sample_001]` |
| `rna_seq_modality` | Routes data without implying equivalence across modalities. | `value_record` | Required | `bulk_rna_seq`, `scRNA_seq`, `snRNA_seq` | `scRNA_seq` |
| `library_strategy` | Preserves reported library strategy or preparation description. | `value_record` | Optional when available | Open vocabulary | `10x 3' gene expression` |
| `sequencing_platform` | Preserves instrument/platform when reported. | `value_record` | Optional when available | Open vocabulary | `Illumina NovaSeq 6000` |
| `strandedness` | Preserves strandedness where meaningful. | `value_record` | Conditionally required when applicable | `unstranded`, `forward`, `reverse`; otherwise explicit missingness | `reverse` |
| `read_layout` | Preserves raw-read layout where meaningful. | `value_record` | Conditionally required when applicable | `single_end`, `paired_end`; processed matrices may be `not_applicable` | `paired_end` |
| `metadata_validation_status` | Assay metadata validation outcome. | String | Required | Shared validation vocabulary | `valid` |

## Input data asset

| Field | Purpose | Type | Requirement | Permitted values / missingness | Example |
|---|---|---|---|---|---|
| `asset_id` | Stable identity for a source file, repository object, or matrix. | String identifier | Required | No missing value | `asset_001` |
| `assay_ids` | Links the asset to all represented sequencing assays without assuming sample granularity. | Array of foreign-key strings | Required | Existing `sequencing_assays.assay_id`; non-empty | `[assay_001]` |
| `input_data_format` | Enables format routing while preserving heterogeneous inputs. | `value_record` | Required | `FASTQ`, `SRA`, `raw_counts`, `UMI_counts`, `H5AD`, `RDS`, `processed_expression_matrix`, `other` | `H5AD` |
| `gene_identifier_type` | Identifies matrix feature identifiers when known. | `value_record` | Conditionally required for expression matrices | Open vocabulary with original value preserved | `Ensembl Gene ID` |
| `source_file_identifier_or_reference` | Locates or names the source asset without requiring a local path. | `value_record` | Required | Path, filename, accession, object key, or equivalent reference | `GSE12345_counts.tsv.gz` |
| `source_url_or_accession_reference` | Provides a remote URL or repository reference when applicable. | `value_record` | Optional when applicable | Explicit four-state missingness | `SRA:SRR123456` |
| `checksum_algorithm` | Identifies the checksum algorithm. | `value_record` | Conditionally required when a checksum is available | Open vocabulary | `SHA-256` |
| `checksum_value` | Supports immutability and integrity verification. | `value_record` | Conditionally required when a checksum is available | Must correspond to `checksum_algorithm` | `3a7b...` |
| `availability` | Records whether the asset can currently be accessed. | String | Required | `available`, `restricted`, `unavailable` | `available` |
| `metadata_validation_status` | Asset metadata validation outcome. | String | Required | Shared validation vocabulary | `valid` |

## Provenance

| Field | Purpose | Type | Requirement | Permitted values / missingness | Example |
|---|---|---|---|---|---|
| `provenance_id` | Stable identity referenced by `value_record` assertions. | String identifier | Required | No missing value | `prov_gsm_1_dose` |
| `entity_type` | Names the owning entity collection. | String | Required | An entity name defined by the schema | `exposures` |
| `entity_id` | Identifies the specific owning record. | String | Required | Must resolve within `entity_type` | `exposure_001` |
| `field_path` | Locates the supported field within the entity. | String | Required | Dot/path notation | `concentration_or_dose` |
| `source_kind` | Describes the provenance source category. | String | Required | `public_repository`, `publication`, `source_file`, `in_house_record`, `curator_assertion`, `other` | `public_repository` |
| `source_reference` | Identifies the exact source record, file, document, or locator. | String | Required | Non-empty source reference | `GSM123456 characteristics: dose` |
| `original_metadata_key` | Preserves the source field/key name when applicable. | String or null | Optional | Null when the source has no keyed structure | `dose_ch1` |
| `recorded_by` | Identifies the curator or automated process that recorded the assertion. | String or null | Optional | Organization-approved user/process identifier | `curator_02` |
| `recorded_at` | Records when provenance was captured. | ISO 8601 datetime string or null | Optional | ISO 8601 when present | `2026-09-02T15:30:00Z` |
| `notes` | Preserves provenance qualifications. | String or null | Optional | Free text | `Unit supplied in publication supplement` |

## Independent DNT reference evidence

These records describe external reference evidence about an agent. They are conceptually and operationally separate from experimental gene-expression results and must never be generated by interpreting those results.

| Field | Purpose | Type | Requirement | Permitted values / missingness | Example |
|---|---|---|---|---|---|
| `evidence_id` | Stable identity for a versioned evidence assertion. | String identifier | Required | No missing value | `dnt_evidence_001` |
| `agent_name` | Names the referenced compound/agent. | `value_record` | Required | Explicit four-state missingness | `valproic acid` |
| `agent_identifier` | Links evidence to experimental agents and supports compound grouping. | `value_record` | Optional when available | Documented identifier namespace/value | `CAS:99-66-1` |
| `reference_state` | Represents the reference evidence state without forcing a binary label. | String | Required | `positive`, `negative_reference`, `unresolved_uncertain` | `positive` |
| `evidence_source` | Identifies the independent authority, dataset, review, or curated source. | String | Required | Non-empty source identity | `curated reference set` |
| `evidence_version` | Versions the evidence source or assertion. | String | Required | Non-empty version identifier | `2026.1` |
| `evidence_reference` | Gives a citation, URL, accession, or internal controlled reference. | String or null | Optional | Null only when unavailable and validation reflects the gap | `doi:10.xxxx/example` |
| `scope_type` | States whether the cited evidence is presented as general, context-limited, or insufficiently specified. | String | Required | `agent_general`, `context_limited`, `unspecified` | `context_limited` |
| `scope_assertions` | Structurally records only source-supported applicability constraints; it does not create context-specific labels. | Array of `evidence_scope_assertion` | Required | Empty for `agent_general` or `unspecified`; populated from the source for `context_limited` | `[{dimension: species, value: ...}]` |
| `scope_notes` | Records limitations, population/model scope, or unresolved interpretation not captured structurally. | String or null | Optional | Free text; does not override `reference_state` or invent a new label | `Evidence scope is limited to the cited context.` |
| `metadata_validation_status` | Evidence-record validation outcome. | String | Required | Shared validation vocabulary | `not_validated` |

## Relationship model

A `study` owns one or more `experiments`; each experiment owns scientifically defined `conditions`, and every sample references one condition as well as its experiment. A sample separately references its biological source and organoid context, preventing repeated text from becoming accidental identity. Planned exposure records belong to conditions, so multiple biological replicates share one explicit group definition. Multiple exposure records support multi-agent or factorial conditions. Sample-specific deviations are additional exposure records scoped to the affected sample and never silently overwrite the planned condition.

`treatment_control_relationships` compare treatment conditions with matched control conditions, allowing shared controls and multi-condition comparisons. `sample_comparison_exceptions` document any scientifically approved sample-level departure. Sequencing assays link to represented samples, and input assets link to represented assays. This supports sample-level FASTQ/SRA files, multi-sample count matrices, UMI matrices, H5AD/RDS objects, and processed expression matrices without pretending those formats or modalities are interchangeable. Provenance records can support any important field assertion. Independent DNT reference evidence links to compounds by documented identity but remains separate from expression results; scope assertions only preserve limitations stated by the independent evidence source.

Referential validation must confirm that a sample and its condition belong to the same experiment; a sample-deviation exposure refers to a sample in its condition; any referenced planned exposure belongs to that same condition; and treatment/control comparison IDs have the declared condition roles. Deviations and sample-level comparison exceptions supplement their planned records and never overwrite them.

Separate conditions can represent multiple concentrations, durations, or developmental stages for the same agent. Multiple planned exposure records on one condition represent reported factorial or multi-agent designs. A single control condition can participate in multiple relationships, and every condition can contain any number of biological-replicate samples without imposing a required count.

## Representational checks

The contract can represent the required hypothetical cases without dataset-specific fields:

1. Public bulk treatment and vehicle-control count matrices: treatment and control conditions with multiple sample replicates, bulk assays, `raw_counts` assets, and a condition-level matched-control relationship.
2. Public bulk FASTQ/SRA only: bulk assays linked to `FASTQ` or `SRA` assets, with public accessions and no matrix requirement.
3. Public scRNA-seq with UMI/H5AD: `scRNA_seq` assays linked to `UMI_counts` and/or `H5AD` assets.
4. Processed-expression-only public data: a `processed_expression_matrix` asset marked with availability and explicit missingness for unavailable raw assets or raw-read attributes.
5. In-house bulk with a positive reference compound: `in_house` study, treatment and control conditions, bulk assay, condition-level relationship, and a separate `positive` DNT evidence record.
6. In-house unresolved compound: the same experimental structure with an independent `unresolved_uncertain` evidence record.

These checks demonstrate structural coverage only; they do not approve data quality, scientific comparability, processing choices, or label interpretation.
