# Detector Simulation Pipeline

A comprehensive Geant4-based detector simulation pipeline for kaon reconstruction analysis. This pipeline automates the complete process from detector geometry configuration through data analysis and histogram generation.

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

### Basic Usage
```bash
./JOBSCRIPTS_TESTS/detector_simulation_master.sh T1 T2 T3 T4 P1 P2 F1 F2 E1
```

Example:
```bash
./JOBSCRIPTS_TESTS/detector_simulation_master.sh 250 260 330 340 215 230 290 300 400
```

This creates configuration `T1-250_T2-260_T3-330_T4-340_P1-215_P2-230_F1-290_F2-300_E1-400` and runs the complete pipeline.

### Pipeline Stages

1. **Build Phase** (17 minutes): Sequential compilation with unique executables
2. **Run Phase** (4 hours): 100 parallel Geant4 simulations  
3. **Parse Phase** (varies): ROOT analysis of simulation outputs
4. **Combine Phase** (minutes): Merge results into final datasets

## Project Structure

```
├── JOBSCRIPTS_TESTS/
│   ├── detector_simulation_master.sh    # Main pipeline script
│   ├── README_detector_pipeline.md      # Detailed documentation
│   └── test_*.job                       # Testing and verification scripts
├── VIKING_FOLDER/
│   ├── SIMULATION_RUNNING/              # Geant4 simulation code and scenarios
│   ├── SIMULATION_RESULTS/              # Generated data outputs
│   ├── DATA_PARSING/                    # ROOT analysis macros
│   ├── HISTOGRAM_MAKING/                # Plotting and visualization tools
│   └── SCENARIOS/                       # Detector configuration files
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

## Development Status

This pipeline has been successfully tested and validated with:
- Multiple detector configurations
- Full 100-simulation production runs
- Resource optimization on HPC clusters
- Comprehensive error handling and recovery

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
**Last Updated**: November 2024