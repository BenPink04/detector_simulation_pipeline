# Detector Simulation Pipeline Documentation

## Overview
This pipeline automates the complete detector simulation and data analysis process with configurable detector positions.

## Usage

### Basic Usage
```bash
./detector_simulation_master.sh T1 T2 T3 T4 P1 P2 F1 F2 E1
```

You can prevent cleanup of intermediate files by passing `--no-cleanup` before the positions:

```bash
./detector_simulation_master.sh --no-cleanup T1 T2 T3 T4 P1 P2 F1 F2 E1
```

Where:
- **T1-T4**: Tracker positions (cm from dipole end)
- **P1-P2**: Pizza detector positions (cm from dipole end)  
- **F1-F2**: FRI wall positions (cm from dipole end)
- **E1**: TOF wall (End detector) position (cm from dipole end)

### Example
```bash
./detector_simulation_master.sh 50 60 130 140 15 30 90 100 200
```
This creates configuration: `T1-50_T2-60_T3-130_T4-140_P1-15_P2-30_F1-90_F2-100_E1-200`

## What the Pipeline Does

### Step 1: Detector Configuration
- Backs up the original `JLabKDetectorConstruction.cc` file
- Edits detector positions based on input parameters
- Builds the Geant4 simulation with new geometry

### Step 2: Build Jobs (100 builds, sequential)
- Creates 100 build jobscripts that run sequentially to avoid source conflicts
- Each build job:
  - Uses `sed` to edit `JLabKRunAction.cc` with unique output filename
  - Rebuilds the simulation with `make -j4`
  - Creates a unique executable: `jlabk_1`, `jlabk_2`, etc.
- Sequential execution prevents race conditions during source editing and compilation
- Total time: ~17 minutes for all 100 builds

### Step 3: Run Jobs (100 simulations, parallel)
- Creates 100 run jobscripts that execute in parallel
- Each run job:
  - Starts immediately after its corresponding build job completes
  - Uses its unique executable (`jlabk_N`)
  - Runs simulation with unique random seed
  - Validates output file size (>10KB) to detect failures
  - Output files named: `{CONFIG}_{1-100}.root`
- All 100 simulations can run simultaneously (if nodes available)
- Total time: ~4 hours for all runs (same as single run!)

### Step 4: Data Parsing Jobs (100 files)
- Creates 100 parsing jobscripts, each processing 1 file
- Runs `KLong_save_momentum_acceptance.C` and `KLong_save_vectors.C`
- Creates acceptance and vector output files for each simulation
- Waits for all run jobs to complete before starting

### Step 5: File Combination
- Waits for all parsing to complete
- Combines all 100 acceptance files into one ROOT file
- Combines all 100 vector files into one ROOT file
- Uses ROOT's `hadd` command for efficient merging

## File Structure

```
VIKING_FOLDER/
├── SIMULATION_RESULTS/
│   └── {CONFIG_STR}/
│       ├── {CONFIG}_{1-100}.root              # Raw simulation files
│       ├── {CONFIG}_{1-100}_acceptance.root   # Acceptance analysis files
│       ├── {CONFIG}_{1-100}_vectors.root      # Vector analysis files
│       ├── {CONFIG}_combined_acceptance.root  # Final combined acceptance
│       └── {CONFIG}_combined_vectors.root     # Final combined vectors
├── DATA_PARSING/
│   └── {CONFIG_STR}_DATA_READ/
│       ├── KLong_save_momentum_acceptance.C
│       └── KLong_save_vectors.C
└── JOBSCRIPTS_TESTS/
    ├── temp_jobs/                             # Generated jobscripts
    └── submit_{CONFIG}.sh                     # Job submission script
```

## Job Dependencies

The pipeline uses SLURM job dependencies to ensure proper execution order:

```
Build 1 → Build 2 → ... → Build 100 (sequential chain)
   ↓         ↓                ↓
Run 1     Run 2     ...    Run 100    (parallel, each depends on its own build)
   └─────────┴───────────────┘
              ↓
  Parsing Jobs (1-100, parallel, wait for all runs)
              ↓
      Combination Job (sequential)
```

### Key Features:
- **Build jobs** run sequentially to avoid source code conflicts during sed editing and compilation
- **Run jobs** start in parallel as soon as their individual build completes
- **Total wall time**: ~17 minutes (builds) + ~4 hours (runs) instead of 400+ hours sequential

## Monitoring

### Check job status
```bash
squeue -u $USER
```

### Check specific job details
```bash
scontrol show job <job_id>
```

### View job logs
```bash
# In the directory where jobs were submitted
ls *.log *.err
```

## Resource Requirements

| Job Type | Time | CPUs | Memory | Count | Notes |
|----------|------|------|--------|-------|-------|
| Build | 10 min | 4 | 2GB | 100 | Sequential chain |
| Run | 4 hours | 4 | 4GB | 100 | Parallel execution |
| Parsing | 12 hours | 4 | 4GB | 100 | 1 file per job |
| Combination | 1 hour | 1 | 4GB | 1 | Sequential |

### Geometry Constraints
- **Dipole magnet**: extends from z=0 to z=200cm
- **Minimum detector position**: z > 210cm to avoid geometry overlaps
- The master script validates all positions before proceeding

## Troubleshooting

### If build jobs fail:
1. Check that Geant4 module loads correctly
2. Verify source file `JLabKRunAction.cc` exists and is writable
3. Check for sufficient disk space
4. Review build job logs: `Build_*-*.log`

### If run jobs fail:
1. Verify the corresponding build job completed successfully
2. Check that unique executable (`jlabk_N`) was created
3. Ensure `batch.mac` file exists in scenario directory
4. Review run job logs: `Run_*-*.log`
5. Check for output file size validation errors (files <10KB indicate early failure)

### If parsing jobs fail:
1. Ensure ROOT is available with the specified version
2. Check that simulation output files exist in results directory
3. Verify parsing scripts are in the parsing directory
4. Review parsing job logs: `Parse_*-*.log`

### If parsing jobs are very slow (high statistics):
**Cause**: Older parsing macros scanned `Ntuple3` and `Ntuple1` *inside* the per‑event loop, which scales roughly as $O(N_{events} \times N_{entries})$ and becomes extremely slow at high statistics.

**Fix**: The macros now build per‑event indexes once and reuse them:
- Per‑event truth map from `Ntuple1`
- Per‑event hit summaries from `Ntuple3`

This reduces the runtime to roughly linear in the number of entries.

**Important**: Parsing jobs now copy the latest macros from [VIKING_FOLDER/DATA_PARSING](VIKING_FOLDER/DATA_PARSING) into the configuration’s parsing directory at job start, so updated performance fixes are always used.

### If combination fails:
1. Check that all acceptance and vector files exist
2. Ensure ROOT's `hadd` command is available
3. Verify sufficient disk space for combined files
4. Review combination job log: `Combine_*-*.log`

### Common Issues:
- **Empty output files (174 bytes)**: Simulation crashed early, check run logs for segfaults
- **Build jobs stuck**: Check if previous build is still running
- **Run jobs pending**: Wait for corresponding build job to complete
- **Source conflicts**: Ensure only one build runs at a time (handled automatically by dependencies)

## Manual Intervention

If you need to run parts manually:

### Cancel all running jobs:
```bash
scancel -u $USER
```

### Submit only build jobs:
```bash
# Submit first build
BUILD_ID=$(sbatch temp_jobs/build_1_{CONFIG}.job | awk '{print $4}')

# Submit remaining builds in chain
for run in {2..100}; do
    BUILD_ID=$(sbatch --dependency=afterok:$BUILD_ID temp_jobs/build_${run}_{CONFIG}.job | awk '{print $4}')
done
```

### Submit run jobs after builds complete:
```bash
# Each run depends on its corresponding build
for run in {1..100}; do
    # Get the build job ID (if still in queue/running)
    BUILD_ID=$(squeue -u $USER -n Build_{CONFIG}_${run} -h -o "%i")
    
    if [ -n "$BUILD_ID" ]; then
        sbatch --dependency=afterok:$BUILD_ID temp_jobs/run_${run}_{CONFIG}.job
    else
        # Build already completed, submit without dependency
        sbatch temp_jobs/run_${run}_{CONFIG}.job
    fi
done
```

### Run parsing after simulations complete:
```bash
# Submit parsing jobs manually with appropriate dependencies
for file_num in {1..100}; do
  sbatch temp_jobs/parse_file_${file_num}_{CONFIG}.job
done
```

### Combine files manually:
```bash
module load ROOT/6.30.06-foss-2023a
cd /users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS/{CONFIG}/

# Combine acceptance files
hadd -f {CONFIG}_combined_acceptance.root {CONFIG}_*_acceptance.root

# Combine vector files  
hadd -f {CONFIG}_combined_vectors.root {CONFIG}_*_vectors.root
```

### Clean up unique executables if needed:
```bash
cd /users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/SCENARIO_5_SIM
rm -f jlabk_*
```

### Keep or delete intermediate files
By default the pipeline will clean up intermediate files and logs after successful completion. To keep files for debugging, pass `--no-cleanup` before the positions:

```bash
./detector_simulation_master.sh --no-cleanup T1 T2 T3 T4 P1 P2 F1 F2 E1
```

## Configuration Examples

### Standard production run:
```bash
./detector_simulation_master.sh 250 260 330 340 215 230 290 300 400
```

### Close detectors (high resolution):
```bash
./detector_simulation_master.sh 230 240 310 320 215 225 270 280 350
```

### Spread out detectors (wider acceptance):
```bash
./detector_simulation_master.sh 300 350 400 450 220 250 350 400 500
```

### Test run with debugging (keeps all intermediate files):
```bash
./detector_simulation_master.sh --no-cleanup 250 260 330 340 215 230 290 300 400
```

## Pipeline Improvements

### v2.0 - Semi-Parallel Build+Run Approach
- **Problem**: Previous version ran simulations sequentially (400+ hours total)
- **Solution**: Split into build and run phases:
  - Builds run sequentially to avoid source conflicts (~17 minutes total)
  - Runs execute in parallel as soon as their build completes (~4 hours total)
  - Uses unique executables (`jlabk_N`) to prevent interference
- **Result**: 100x speedup for total pipeline execution

### Key Technical Details:
- **sed editing**: Each build uses `sed` to modify `JLabKRunAction.cc` with unique output filename
- **Unique executables**: Prevents parallel runs from interfering with each other
- **File validation**: Checks output file size >10KB to detect early failures
- **Automatic cleanup**: Removes unique executables after successful runs

## Analyzing Results

Once the pipeline completes, you'll have two combined output files for each configuration:
- `{CONFIG}_combined_acceptance.root` - Contains kaon reconstruction acceptance data
- `{CONFIG}_combined_vectors.root` - Contains kaon momentum reconstruction data

### Plotting Individual Configuration Results

#### Resolution Histogram (original):
```bash
cd /users/bp969/scratch/VIKING_FOLDER/HISTOGRAM_MAKING
module load ROOT/6.30.06-foss-2023a

# Plot mean kaon momentum resolution
root -l 'KLong_plot_resolution_histbar.C()'
```

#### Acceptance Histogram (original):
```bash
# Plot kaon reconstruction acceptance
root -l 'KLong_plot_acceptance.C()'
```

**Note**: You'll need to edit these macros to point to your combined files instead of the hardcoded individual seed files.

### Comparing Multiple Configurations

To overlay and compare results from different detector configurations:

#### Resolution Comparison:
```bash
cd /users/bp969/scratch/VIKING_FOLDER/HISTOGRAM_MAKING
module load ROOT/6.30.06-foss-2023a

 Method 1: Edit the configuration list in the macro
 Edit lines 90-96 in KLong_plot_resolution_comparison.C to add your files:
   {"T1-250_T2-260_T3-330_T4-340_P1-215_P2-230_F1-290_F2-300_E1-400_combined_vectors.root", 
    "Config 1: T1=250, T2=260, T3=330, T4=340", kBlue},
   {"T1-300_T2-310_T3-380_T4-390_P1-215_P2-230_F1-290_F2-300_E1-400_combined_vectors.root", 
    "Config 2: T1=300, T2=310, T3=380, T4=390", kRed},

root -l 'KLong_plot_resolution_comparison.C()'

# Method 2: Pass files as arguments in ROOT
root -l
root [0] .L KLong_plot_resolution_comparison.C
root [1] vector<string> files = {"config1_combined_vectors.root", "config2_combined_vectors.root"};
root [2] vector<string> labels = {"Config 1", "Config 2"};
root [3] vector<int> colors = {kBlue, kRed};
root [4] KLong_plot_resolution_comparison(files, labels, colors);
```

Output: `kaon_resolution_comparison.png` and `kaon_resolution_comparison.pdf`

#### Acceptance Comparison:
```bash
 Method 1: Edit the configuration list in the macro
 Edit lines 87-93 in KLong_plot_acceptance_comparison.C to add your files:
   {"T1-250_T2-260_T3-330_T4-340_P1-215_P2-230_F1-290_F2-300_E1-400_combined_acceptance.root", 
    "Config 1: T1=250, T2=260, T3=330, T4=340", kBlue},
   {"T1-300_T2-310_T3-380_T4-390_P1-215_P2-230_F1-290_F2-300_E1-400_combined_acceptance.root", 
    "Config 2: T1=300, T2=310, T3=380, T4=390", kRed},

root -l 'KLong_plot_acceptance_comparison.C()'

# Method 2: Pass files as arguments in ROOT
root -l
root [0] .L KLong_plot_acceptance_comparison.C
root [1] vector<string> files = {"config1_combined_acceptance.root", "config2_combined_acceptance.root"};
root [2] vector<string> labels = {"Config 1", "Config 2"};
root [3] vector<int> colors = {kBlue, kRed};
root [4] KLong_plot_acceptance_comparison(files, labels, colors);
```

Output: `kaon_acceptance_comparison.png` and `kaon_acceptance_comparison.pdf`

### Comparison Macro Features:
- **Overlays multiple configurations** with different colors for easy comparison
- **Automatic legend** generation with custom labels
- **Smart y-axis scaling** to fit all histograms
- **Processes all entries** in combined ROOT files (not just first entry)
- **Saves as PNG and PDF** for presentations and papers
- **Default color palette**: kBlue, kRed, kGreen+2, kMagenta+2, kOrange+1, kCyan+2, kViolet, kPink+1

### Typical Workflow:

1. **Run multiple configurations**:
   ```bash
   ./detector_simulation_master.sh 250 260 330 340 215 230 290 300 400
   ./detector_simulation_master.sh 300 310 380 390 215 230 290 300 400
   ./detector_simulation_master.sh 270 280 350 360 215 230 290 300 400
   ```

2. **Wait for all pipelines to complete** (check with `squeue -u $USER`)

3. **Compare results**:
   ```bash
   cd /users/bp969/scratch/VIKING_FOLDER/HISTOGRAM_MAKING
   module load ROOT/6.30.06-foss-2023a
   
   # Create comparison plots for both resolution and acceptance
   root -l 'KLong_plot_resolution_comparison.C()'
   root -l 'KLong_plot_acceptance_comparison.C()'
   ```

4. **Analyze plots** to determine optimal detector configuration based on:
   - **Resolution**: Lower is better (shows better momentum reconstruction)
   - **Acceptance**: Higher is better (shows more kaons are reconstructed)