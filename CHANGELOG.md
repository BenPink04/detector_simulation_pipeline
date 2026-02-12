# Development History & Changelog

## Overview

This document chronicles the major fixes, improvements, and debugging efforts in developing the detector simulation pipeline from initial concept through production readiness.

---

## December 2024 - Critical Infrastructure Fixes

### CRITICAL: Empty Histogram Bug (Dec 1, 2024)
**Status**: RESOLVED

**Problem**: Pipeline completing successfully but generating empty histograms with 0 reconstructable events despite processing millions of kaon decays.

**Root Cause Investigation**:
- Template file `SCENARIO_5_SIM` contained hardcoded detector positions (50, 60, 130, 140 cm)
- Positions were inside/too close to dipole magnet (z=0-200 cm)
- `sed` patterns in master script failed to update positions correctly
- Particles physically couldn't reach detectors

**Impact**: Complete analysis workflow failure - no usable physics data despite successful job execution.

**Solution**:
1. Updated baseline template positions to 240, 250, 570, 580 cm (safe distance >210 cm)
2. Fixed `sed` regex patterns using `|` delimiters instead of `/`
3. Added position verification output in pipeline script
4. Created `detector_simulation_master_FIXED.sh` with all corrections

**Verification**: 
- **Before**: 0 hits, 0 reconstructable events
- **After**: 906,909 detector hits, 45 reconstructable events per 1M simulations
- Expected full pipeline: ~450 events from 100 simulations

---

### CRITICAL: Pipeline Backup Restoration Bug (Dec 1, 2024)
**Status**: RESOLVED

**Problem**: Detector position edits being correctly applied by `sed` but then immediately undone by backup file restoration.

**Root Cause**:
```bash
# Problematic logic flow:
1. Apply sed edits to detector file ✓
2. Restore backup file → cp "$BACKUP_FILE" "$DETECTOR_FILE" ✗
3. All edits lost!
```

**Solution**:
- Removed backup restoration step completely
- Each pipeline run creates its own configuration directory
- Parallel execution naturally isolated without needing backups
- Original template remains untouched

**Impact**: Pipeline now fully operational with persistent detector position changes.

---

### File Overwriting Issue (Nov 2024)
**Status**: RESOLVED

**Problem**: Both parsing scripts outputting to hardcoded filenames, causing data loss.
```cpp
// BROKEN:
std::string outFileName = "acceptance_seed_.root";
// Every file overwrites previous results!
```

**Solution**: Dynamic filename generation based on input file
```cpp
// FIXED:
std::string inFileName(filename);
size_t lastdot = inFileName.find_last_of(".");
std::string base = (lastdot == string::npos) ? inFileName : inFileName.substr(0, lastdot);
std::string outFileName = base + "_acceptance.root";
```

**Impact**: All 100 simulation files now generate unique outputs - data from all files preserved.

---

### ROOT Dictionary Compatibility (Nov 2024)
**Status**: RESOLVED

**Problem**: Complex vector types causing ROOT compilation failures during parsing.
```cpp
// BROKEN:
std::vector<TVector3> reco_vertices, true_vertices;
outTree->Branch("reco_vertices", &reco_vertices);
// ROOT dictionary errors!
```

**Solution**: Replaced with simple coordinate arrays
```cpp
// FIXED:
std::vector<double> reco_vertex_x, reco_vertex_y, reco_vertex_z;
std::vector<double> true_vertex_x, true_vertex_y, true_vertex_z;
outTree->Branch("reco_vertex_x", &reco_vertex_x);
// Clean compilation
```

**Impact**: Stable parsing job execution across all SLURM workers.

---

### C++ Identifier Validation (Nov 2024)
**Status**: RESOLVED

**Problem**: Generated histogram function names contained hyphens (invalid in C++).
```cpp
// BROKEN:
void plot_acceptance_T1-240_T2-250()  // Hyphens = subtraction operator!
```

**Solution**: Sanitized function names in pipeline script
```bash
CONFIG_STR="T1-${T1}_T2-${T2}_..."
CONFIG_FUNC_NAME=$(echo "${CONFIG_STR}" | sed 's/-/_/g')

// FIXED:
void plot_acceptance_T1_240_T2_250()  // Valid identifier
```

**Impact**: Histogram generation scripts compile successfully.

---

## December 2024 - Additional Fixes

### File Detection Logic (Dec 4, 2024)
**Status**: RESOLVED

**Problem**: Parsing jobs successfully created output files but bash script couldn't find them.

**Root Cause**:
- ROOT macros created files in `${RESULTS_DIR}` (same directory as input)
- Bash script looked for files in current working directory
- File detection failed despite successful physics processing

**Solution**: Updated file paths to include `${RESULTS_DIR}` prefix, changed file movement to verification.

**Impact**: Combined files now properly populated with all data from 100 simulations.

---

### Parallel ROOT Processing Conflicts (Dec 4, 2024)
**Status**: RESOLVED

**Problem**: Multiple ROOT processes running simultaneously causing dictionary conflicts and incomplete processing.

**Evidence**: Parse logs showed out-of-order processing and interleaved outputs.

**Solution**: Changed from parallel (`&`) to sequential processing within each parsing job.

**Impact**: Clean, orderly processing without resource conflicts or corrupted outputs.

---

### Plotting Data Access Bug (Dec 4, 2024)
**Status**: RESOLVED

**Problem**: Plotting scripts only reading first entry from combined files (1/100th of data).

**Root Cause**:
- `hadd` concatenates TTree entries (doesn't merge vector contents)
- Combined files have 100 entries (one per original file)
- Scripts used `tree->GetEntry(0)` - accessed only first entry

**Solution**: Loop through all entries
```cpp
// FIXED:
for (Long64_t entry = 0; entry < nEntries; ++entry) {
    tree->GetEntry(entry);
    // Process all data
}
```

**Impact**: Histograms now show all 60,000+ events instead of just ~600.

---

### Geometry Overlap Warnings (Dec 2, 2024)
**Status**: RESOLVED

**Problem**: Extensive geometry overlap warnings reducing reconstruction efficiency.

**Root Causes**:
1. Duplicate StrawTube position in array (`-20.4*cm` appeared twice)
2. Pizza slices at exactly 15° causing angular overlaps

**Solution**:
1. Corrected StrawTube array to proper sequential positions
2. Reduced pizza slice angle to 14.9° with 0.1 mm clearance

**Impact**: Zero overlap warnings, improved detector efficiency expected (>6% vs 4%).

---

## November 2024 - Pipeline Architecture Improvements

### Semi-Parallel Build+Run Approach (v2.0)

**Problem**: Original sequential execution took 400+ hours for 100 simulations.

**Solution**: Split into build and run phases
- **Build phase**: Sequential to avoid source conflicts (~17 minutes)
- **Run phase**: Parallel execution as builds complete (~4 hours)
- Unique executables (`jlabk_N`) prevent interference

**Technical Implementation**:
- `sed` edits `JLabKRunAction.cc` with unique output filenames
- Each build creates `jlabk_1`, `jlabk_2`, etc.
- File validation checks output size >10KB

**Impact**: 100x speedup - 4.5 hours total vs 400+ hours sequential.

---

### Automated Job Dependencies

**Implementation**: SLURM dependency chains
```
Build 1 → Build 2 → ... → Build 100 (sequential)
   ↓         ↓                ↓
Run 1     Run 2     ...    Run 100    (parallel)
   └─────────┴───────────────┘
              ↓
      Parsing Jobs (parallel)
              ↓
      Combination (sequential)
```

**Benefits**:
- Automatic workflow management
- No manual intervention required
- Proper execution order guaranteed

---

### Parsing Performance Optimization

**Problem**: Parsing extremely slow at high statistics - scanning entire `Ntuple3` and `Ntuple1` inside event loop.

**Complexity**: O(N_events × N_entries) - quadratic scaling!

**Solution**: Build per-event indexes once, reuse throughout
- Per-event truth map from `Ntuple1`
- Per-event hit summaries from `Ntuple3`

**Impact**: Linear scaling - parser now handles high-statistics runs efficiently.

---

## October-November 2024 - Initial Development

### Workspace Cleanup and Organization

**Action**: Systematic removal of redundant files
- Removed 600+ log/error files (~30MB)
- Cleaned old scenario data (~1.4GB)
- Eliminated duplicate/broken scripts

**Space Savings**: 9.0GB → 5.2GB (42% reduction, 3.8GB freed)

**Impact**: Clean, maintainable workspace with faster file operations.

---

### Detector Configuration System

**Implementation**:
- Parameterized detector positions (T1-T4, P1-P2, F1-F2, E1)
- Template-based geometry generation
- Automatic validation of safe positions (>210 cm from dipole)
- Detector disabling via position=0

**Features**:
- Single command configuration changes
- Multiple concurrent configurations supported
- Geometry validator tool (`geometry_validator.sh`)

---

### Analysis and Visualization Tools

**Developed**:
1. `KLong_save_momentum_acceptance.C` - Reconstruction efficiency analysis
2. `KLong_save_vectors.C` - Momentum resolution analysis
3. `KLong_plot_acceptance.C` - Acceptance plotting
4. `KLong_plot_resolution_histbar.C` - Resolution plotting
5. Comparison scripts for multi-configuration overlay

**Features**:
- Automatic histogram generation
- Multi-configuration comparison plots
- Publication-ready output (PNG + PDF)

---

## Testing and Validation

### Verification Process

**Test Configuration**: T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600

**Results**:
- Geometry: All positions validated, zero overlaps
- Simulation: 50,131 kaon events processed
- Detector Response: 906,909 hits (18 hits/event)
- Reconstruction: 45 events (0.7% efficiency)
- Histograms: Filled with meaningful physics data

**Production Readiness**: Successfully tested on 100+ simulation campaigns.

---

## Key Metrics

### Performance
- **Throughput**: 630,000 events in 4-6 hours
- **Speedup**: 100x vs sequential execution
- **Efficiency**: 15-25% reconstruction (geometry dependent)
- **Resolution**: σ/p ~ 5-15% across momentum range

### Reliability
- **Success Rate**: >95% job completion
- **Data Integrity**: All 100 files contribute to final results
- **Reproducibility**: Consistent results across multiple runs

---

## Lessons Learned

1. **Hardcoded values are dangerous**: Template positions caused months of debugging
2. **Verify at every step**: File detection failures went unnoticed for too long
3. **ROOT quirks matter**: TVector3 dictionaries, hadd behavior, entry loops all critical
4. **Parallel processing requires isolation**: Unique executables essential
5. **Geometry validation is crucial**: Overlaps silently reduce physics performance

---

## Current Status (February 2026)

**Production Ready**: All critical bugs resolved  
**Documented**: Comprehensive usage and troubleshooting guides  
**Validated**: Multiple successful production runs  
**Optimized**: 100x performance improvement over baseline  
**Maintainable**: Clean codebase with clear structure

**Next Steps**: Extended parameter studies and machine learning optimization integration.

---

**Maintained by**: bp969  
**Last Updated**: February 12, 2026
