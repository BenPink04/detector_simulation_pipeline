# Detector Simulation Pipeline

Automated Geant4-based pipeline for kaon reconstruction analysis in configurable detector geometries.

## Quick Start

### Basic Usage
```bash
cd JOBSCRIPTS_TESTS
./detector_simulation_master_FIXED.sh T1 T2 T3 T4 P1 P2 F1 F2 E1
```

**Example**:
```bash
./detector_simulation_master_FIXED.sh 240 250 570 580 215 230 260 270 600
```

### Parameters
- **T1-T4**: Tracker positions (cm from z=0)
- **P1-P2**: Pizza detector positions
- **F1-F2**: FRI wall positions  
- **E1**: TOF wall position
- **Minimum**: All positions >210 cm (dipole extends 0-200 cm)

### Disabling Detectors
Set position to `0` to disable any detector:
```bash
./detector_simulation_master_FIXED.sh 240 250 680 690 0 0 260 270 700  # No pizza detectors
```

## Pipeline Workflow

```
Build (17 min) → Simulate (4 hrs, 100 parallel) → Parse → Combine → Analyze
```

1. **Build**: Compile Geant4 with detector geometry (100 unique executables)
2. **Simulate**: Generate kaon events (100 jobs, ~630k events total)
3. **Parse**: Extract physics data (momentum, vertices)
4. **Combine**: Merge results with `hadd`
5. **Analyze**: Generate acceptance and resolution histograms

## Requirements

- **Cluster**: SLURM-based HPC
- **Software**: Geant4, ROOT 6.30.06+, CMake 3.26.5+
- **Resources**: ~5GB disk space per configuration
- **Runtime**: 4-6 hours total (wall time)

## Output Files

```
SIMULATION_RESULTS/{CONFIG}/
├── {CONFIG}_combined_acceptance.root  # Reconstruction efficiency data
└── {CONFIG}_combined_vectors.root     # Momentum resolution data
```

## Analysis

### Single Configuration
```bash
cd VIKING_FOLDER/HISTOGRAM_MAKING
module load ROOT/6.30.06-foss-2023a

root -l 'KLong_plot_acceptance.C()'        # Efficiency vs momentum
root -l 'KLong_plot_resolution_histbar.C()' # Resolution analysis
```

### Multiple Configuration Comparison
```bash
root -l 'KLong_plot_acceptance_comparison.C()'  # Overlay multiple configs
root -l 'KLong_plot_resolution_comparison.C()'
```

## Monitoring

```bash
squeue -u $USER              # Check job status
./monitor_jobs.sh            # Pipeline monitoring tool
./monitor_jobs.sh watch      # Continuous updates
```

## Troubleshooting

### Empty Histograms
**Cause**: Detector positions too close to dipole magnet  
**Fix**: Ensure all positions >210 cm, use `detector_simulation_master_FIXED.sh`

### Job Failures
```bash
# Verify modules
module avail Geant4 ROOT

# Check disk space
df -h $HOME

# Review logs
tail -f temp_jobs/*.log
```

### Position Verification
```bash
# Confirm detector positions in generated geometry
grep "trkrPosZ\|pizzaPosZ\|friPosZ\|tofPosZ" \
  VIKING_FOLDER/SIMULATION_RUNNING/T1-*_T2-*/src/JLabKDetectorConstruction.cc
```

## Key Features

- **Parallel Processing**: 100x speedup vs sequential execution
- **Automated Dependencies**: SLURM job chains ensure proper execution order
- **Configurable Geometry**: Easy detector position adjustment
- **Comprehensive Analysis**: Acceptance, resolution, and comparison plots
- **Production Ready**: Validated on 100+ simulation campaigns

## Documentation

- **Full Manual**: [`JOBSCRIPTS_TESTS/README_detector_pipeline.md`](JOBSCRIPTS_TESTS/README_detector_pipeline.md)
- **Development History**: [`CHANGELOG.md`](CHANGELOG.md)
- **Status Reports**: [`PIPELINE_STATUS_REPORT.md`](PIPELINE_STATUS_REPORT.md)

## Project Structure

```
JOBSCRIPTS_TESTS/
├── detector_simulation_master_FIXED.sh  # Main pipeline (use this)
├── overnight_pipeline.sh                # Full automation with histograms
└── monitor_jobs.sh                      # Job monitoring tool

VIKING_FOLDER/
├── SIMULATION_RUNNING/BASE_SCENARIO_DIR/  # Geant4 source template
├── SIMULATION_RESULTS/                    # Generated data
├── DATA_PARSING/                          # ROOT analysis macros
└── HISTOGRAM_MAKING/                      # Plotting scripts
```

## Performance

- **Throughput**: 630,000 kaon events in ~4 hours
- **Efficiency**: 15-25% reconstruction rate (geometry dependent)
- **Resolution**: σ/p ~ 5-15% across momentum range
- **Scalability**: Tested up to 100 concurrent SLURM jobs

---

**Author**: bp969 | **Cluster**: Viking2 HPC | **Updated**: February 2026