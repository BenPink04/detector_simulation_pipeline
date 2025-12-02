# Detector Simulation Pipeline

A comprehensive Geant4-based detector simulation pipeline for kaon reconstruction analysis. This pipeline automates the complete process from detector geometry configuration through data analysis and histogram generation.

## 🚨 **IMPORTANT: Critical Fixes Applied (December 2024)**

**If you're getting empty histograms or zero reconstructable events, you MUST use the fixed version:**

- ✅ **Use**: `detector_simulation_master_FIXED.sh` 
- ❌ **Avoid**: `detector_simulation_master.sh` (has detector positioning bugs)

**Key fix**: Corrected hardcoded detector positions that were too close to dipole magnet, causing particles to never reach detectors.

## Overview

This project provides an automated pipeline for:
- **Detector Configuration**: Parameterized detector geometry setup
- **Parallel Simulation**: Efficient Geant4 simulation with job dependency management  
- **Data Processing**: ROOT-based analysis of simulation results
- **Visualization**: Comparison plots for different detector configurations

## Key Features

- **Configurable Detector Geometry**: Easy parameter adjustment for tracker, pizza detector, FRI wall, and TOF positions
- **Scalable Parallel Processing**: 100 simulations run in parallel with optimal resource utilization
- **Automated Job Dependencies**: SLURM-based pipeline with proper build → run → parse → combine workflow
- **Comprehensive Analysis Tools**: ROOT macros for acceptance analysis and momentum resolution studies
- **Multi-Configuration Comparison**: Overlay plots to compare different detector setups

## Quick Start

### ⚠️ **CRITICAL: Use Fixed Version**

**Use the FIXED master script** to avoid detector positioning issues:
```bash
./JOBSCRIPTS_TESTS/detector_simulation_master_FIXED.sh T1 T2 T3 T4 P1 P2 F1 F2 E1
```

Example:
```bash
./JOBSCRIPTS_TESTS/detector_simulation_master_FIXED.sh 250 260 330 340 215 230 290 300 400
```

This creates configuration `T1-250_T2-260_T3-330_T4-340_P1-215_P2-230_F1-290_F2-300_E1-400` and runs the complete pipeline.

### Detector Position Requirements
- **All positions in cm from dipole origin (z=0)**
- **Minimum safe distance: >210 cm** (dipole extends 0-200 cm)
- **Recommended ranges**:
  - Trackers (T1-T4): 240-600 cm
  - Pizza detectors (P1-P2): 215-300 cm  
  - FRI walls (F1-F2): 260-400 cm
  - TOF wall (E1): 400-800 cm

### Verification Before Running

**Quick check to ensure you have the fixed version:**
```bash
# 1. Check for fixed script
ls -la JOBSCRIPTS_TESTS/detector_simulation_master_FIXED.sh

# 2. Verify template has correct positions  
grep "240.*cm.*250.*cm.*570.*cm.*580.*cm" VIKING_FOLDER/SIMULATION_RUNNING/SCENARIO_5_SIM/src/JLabKDetectorConstruction.cc

# 3. Test with known working data
root -l -b -e '.L test_vectors.C+' -e 'KLong_save_vectors("VIKING_FOLDER/DATA_PARSING/T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600_DATA_READ/out_seed_Scen5_1.root")' -e '.q'
```

**Expected result**: Should show `Successfully reconstructed: 45` events (not 0).

### Pipeline Stages

1. **Build Phase** (17 minutes): Sequential compilation with unique executables
2. **Run Phase** (4 hours): 100 parallel Geant4 simulations  
3. **Parse Phase** (varies): ROOT analysis of simulation outputs
4. **Combine Phase** (minutes): Merge results into final datasets

## Project Structure

```
├── JOBSCRIPTS_TESTS/
│   ├── detector_simulation_master_FIXED.sh  # ✅ FIXED pipeline script (USE THIS)
│   ├── detector_simulation_master.sh        # ❌ Original (has positioning bugs)
│   ├── README_detector_pipeline.md          # Detailed documentation
│   └── test_*.job                           # Testing and verification scripts
├── VIKING_FOLDER/
│   ├── SIMULATION_RUNNING/              # Geant4 simulation code and scenarios
│   │   ├── SCENARIO_5_SIM/             # ✅ FIXED template (corrected positions)
│   │   └── T1-*_T2-*_*/                # Generated configuration directories
│   ├── SIMULATION_RESULTS/              # Generated data outputs  
│   ├── DATA_PARSING/                    # ROOT analysis macros
│   ├── HISTOGRAM_MAKING/                # Plotting and visualization tools
│   └── SCENARIOS/                       # Detector configuration files
├── test_vectors.C                       # ✅ FIXED analysis script (added includes)
└── README.md                           # This file
```

## Requirements

- **Geant4**: Physics simulation framework
- **ROOT**: Data analysis and visualization (v6.30.06+)  
- **SLURM**: Job scheduling system
- **CMake**: Build system (v3.26.5+)
- **GCC**: C++ compiler with C++17 support

## Documentation

For complete documentation including:
- Detailed parameter descriptions
- Resource requirements and optimization
- Troubleshooting guides  
- Analysis workflow examples
- Performance benchmarking

See: [`JOBSCRIPTS_TESTS/README_detector_pipeline.md`](JOBSCRIPTS_TESTS/README_detector_pipeline.md)

## Results Analysis

The pipeline generates two key output files per configuration:
- `{CONFIG}_combined_acceptance.root`: Kaon reconstruction acceptance data
- `{CONFIG}_combined_vectors.root`: Momentum reconstruction analysis

### Plotting Results
```bash
cd VIKING_FOLDER/HISTOGRAM_MAKING
module load ROOT/6.30.06-foss-2023a

# Single configuration analysis
root -l 'KLong_plot_resolution_histbar.C()'
root -l 'KLong_plot_acceptance.C()'

# Multi-configuration comparison  
root -l 'KLong_plot_resolution_comparison.C()'
root -l 'KLong_plot_acceptance_comparison.C()'
```

## Performance

- **Sequential baseline**: ~400 hours for 100 simulations
- **Optimized pipeline**: ~4.5 hours total (100x speedup)
- **Resource efficiency**: Parallel execution with dependency management
- **Scalability**: Tested up to 100 concurrent jobs

## Critical Fixes (December 2024)

### ⚠️ **Detector Positioning Issue - RESOLVED**

**Problem**: Original pipeline was generating empty histograms due to incorrect detector positions.

**Root Cause**: 
- SCENARIO_5_SIM template had hardcoded detector positions too close to dipole magnet (50, 60, 130, 140 cm)
- `sed` patterns in master script failed to update positions correctly
- Particles couldn't reach detectors, resulting in 0 reconstructable events

**Solution**: 
- **Fixed template**: Updated SCENARIO_5_SIM with correct baseline positions (240, 250, 570, 580 cm)
- **Fixed script patterns**: Corrected `sed` regex patterns using `|` delimiters  
- **Verification**: Added position verification output to confirm changes applied

**Impact**: 
- **Before fix**: 0 reconstructable events, empty histograms
- **After fix**: ~45 reconstructable events per 1M simulated events
- **Full pipeline**: Expected ~450 total events from 10M simulations

### ✅ **Verification Results**

Tested configuration T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600:
```
Total triple-pion events: 6,380
Successfully reconstructed: 45 events
Failed due to insufficient hits: 6,165 (reduced from 6,323)
Failed due to timing: 166 events
```

**Key detector positions verified**:
- Trackers: 240, 250, 570, 580 cm ✓
- Pizza detectors: 215, 230 cm ✓  
- FRI walls: 260, 270 cm ✓
- TOF wall: 600 cm ✓

## Development Status

This pipeline has been successfully tested and validated with:
- Multiple detector configurations
- Full 100-simulation production runs  
- Resource optimization on HPC clusters
- Comprehensive error handling and recovery
- **Critical detector positioning fixes (Dec 2024)**

## Troubleshooting

### Empty Histograms / Zero Reconstructable Events

**Symptoms**: Pipeline completes successfully but generates empty histograms or reports 0 reconstructable events.

**Diagnosis**:
```bash
# Check if you're using the fixed script
ls JOBSCRIPTS_TESTS/detector_simulation_master_FIXED.sh

# Verify detector positions in active configuration
grep -E "(trkr|pizza|fri|tof)PosZ.*=" VIKING_FOLDER/SIMULATION_RUNNING/*/src/JLabKDetectorConstruction.cc

# Test analysis on sample data
root -l -b -e '.L test_vectors.C+' -e 'KLong_save_vectors("your_file.root")' -e '.q'
```

**Solution**: 
1. Use `detector_simulation_master_FIXED.sh` instead of the original script
2. Ensure detector positions are >210 cm from dipole origin
3. Verify SCENARIO_5_SIM template has correct baseline positions

### Compilation Errors in Analysis Scripts

**Symptoms**: ROOT analysis fails with "unknown type name" errors.

**Solution**: Use the fixed `test_vectors.C` which includes required headers:
```cpp
#include "TRandom3.h"
#include <set>
```

### Job Submission Issues

**Symptoms**: SLURM jobs fail during build phase.

**Solution**:
```bash
# Check module availability
module avail Geant4 ROOT

# Verify SLURM account access  
sacctmgr show assoc user=$USER

# Test individual job submission
sbatch JOBSCRIPTS_TESTS/temp_jobs/build_1_*.job
```

## Contributing

The pipeline is designed for detector optimization studies. Key areas for extension:
- Additional detector types and configurations
- Enhanced analysis macros and visualization
- Integration with machine learning optimization
- Extended parameter sweep capabilities

## License

This project is developed for academic research purposes in particle physics detector optimization.

---

**Author**: bp969  
**Environment**: Viking2 HPC Cluster  
**Framework**: Geant4 + ROOT + SLURM  
**Last Updated**: December 2024 (Critical detector positioning fixes applied)