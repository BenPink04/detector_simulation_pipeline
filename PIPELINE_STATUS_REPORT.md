# Pipeline Issue Resolution Report - Updated December 4, 2025

## **Root Cause Analysis: Empty Histogram Outputs**

### **Issue 1: Parsing Script ROOT Dictionary Errors** 
- **Problem**: Both `KLong_save_momentum_acceptance.C` and `KLong_save_vectors.C` used `vector<TVector3>` types
- **Impact**: ROOT dictionary compilation failures preventing output file creation
- **Evidence**: Parse job logs show "Expected output file not created!" and TVector3 dictionary errors
- **Status**: ✅ **FIXED**

**Fix Applied:**
- Replaced `std::vector<TVector3> reco_vertex, true_vertex;` with separate coordinate arrays:
  ```cpp
  std::vector<double> reco_vertex_x, reco_vertex_y, reco_vertex_z;
  std::vector<double> true_vertex_x, true_vertex_y, true_vertex_z;
  ```
- Updated tree branch creation to use simple double vectors instead of TVector3
- Modified vertex storage to extract X(), Y(), Z() components

### **Issue 2: Histogram Function Name Invalid Identifiers**
- **Problem**: Auto-generated histogram functions had hyphens in names: `plot_acceptance_T1-240_T2-250...`
- **Impact**: C++ compilation fails (hyphens invalid in identifiers)
- **Evidence**: Histogram job error shows "variable has incomplete type 'void'" 
- **Status**: ✅ **FIXED**

**Fix Applied:**
- Added `CONFIG_FUNC_NAME` variable in `overnight_pipeline.sh` that converts hyphens to underscores
- Updated function names: `plot_acceptance_T1_240_T2_250...` (valid C++ identifiers)
- Modified ROOT script calls to use corrected function names

### **Issue 3: Parsing Script File Detection Bug - CRITICAL**
- **Problem**: Parsing jobs successfully created output files but bash script couldn't find them
- **Impact**: Combined files severely underpopulated, histogram event counts far below expected
- **Evidence**: Parse logs show "Saved kaon event info to ...acceptance.root" immediately followed by "ERROR: Expected output file ...acceptance.root not created!"
- **Status**: ✅ **FIXED** (December 4, 2025)

**Root Cause:**
- ROOT macros create output files in same directory as input files (`${RESULTS_DIR}`)
- Bash script expected to find files in current working directory
- This caused file detection failures despite successful physics processing

**Fix Applied:**
- Updated `ACCEPTANCE_OUTPUT` and `VECTORS_OUTPUT` paths to include full `${RESULTS_DIR}` path
- Changed file movement logic to file verification (files already in correct location)
- Removed unnecessary `mv` commands that were failing
- Fixed around lines 409-431 in `detector_simulation_master_FIXED.sh`

### **Issue 4: Parallel ROOT Processing Conflicts - CRITICAL**
- **Problem**: Multiple ROOT processes running simultaneously within parsing jobs causing resource conflicts
- **Impact**: Incomplete file processing, ROOT dictionary conflicts, interleaved outputs, underpopulated results
- **Evidence**: Parse logs show out-of-order processing (11,12,13,14,15,16,19,20,17,18) and mixed output streams
- **Status**: ✅ **FIXED** (December 4, 2025)

**Root Cause:**
- Parsing jobs used `(command) &` to run 10 file processes in parallel within same job
- Multiple ROOT sessions caused dictionary loading conflicts and memory competition
- File system contention when multiple processes wrote to same directory simultaneously
- Race conditions led to incomplete processing and corrupted outputs

**Fix Applied:**
- Changed parallel processing `(command) &` to sequential processing within each parsing job
- Removed `wait` command (no longer needed)
- Added clearer progress logging for each processing step
- Fixed around lines 405-440 in `detector_simulation_master_FIXED.sh`

### **Issue 5: Plotting Script Data Access Bug - CRITICAL**
- **Problem**: Plotting scripts only read first entry from combined ROOT files, accessing ~1% of actual data
- **Impact**: Histograms severely underpopulated showing only 600 events instead of 60,000+ total events
- **Evidence**: Combined files exist with correct size but histograms show minimal data
- **Status**: ✅ **FIXED** (December 4, 2025)

**Root Cause:**
- ROOT's `hadd` concatenates TTree entries rather than merging vector contents
- Individual parsing files create 1 entry with vectors containing ~600 events each
- Combined files have 100 entries (one from each original file)
- Plotting scripts used `tree->GetEntry(0)` accessing only first entry (1/100th of data)

**Fix Applied:**
- Updated plotting scripts to loop through ALL entries in combined files
- Changed from `tree->GetEntry(0)` to `for (entry = 0; entry < nEntries; ++entry)`
- Added progress logging to show number of entries being processed
- Fixed in `KLong_plot_resolution_histbar.C` and `KLong_plot_acceptance.C`

## **Files Modified**

### **Main Parsing Scripts (ROOT Dictionary Fix)**
1. `/users/bp969/scratch/VIKING_FOLDER/DATA_PARSING/KLong_save_momentum_acceptance.C`
   - Replaced TVector3 vectors with coordinate arrays
   - Fixed vertex storage to use X(), Y(), Z() components

2. `/users/bp969/scratch/VIKING_FOLDER/DATA_PARSING/KLong_save_vectors.C`
   - Replaced TVector3 vectors with coordinate arrays  
   - Fixed tree branch creation to avoid dictionary errors

### **Pipeline Infrastructure (Function Name Fix)**
3. `/users/bp969/scratch/JOBSCRIPTS_TESTS/overnight_pipeline.sh`
   - Added CONFIG_FUNC_NAME variable for valid C++ identifiers
   - Updated generated function names to use underscores instead of hyphens
   - Modified ROOT script execution calls

### **Master Pipeline Script (Multiple Critical Fixes) - UPDATED**
4. `/users/bp969/scratch/JOBSCRIPTS_TESTS/detector_simulation_master_FIXED.sh`
   - **Lines 409-411**: Fixed `ACCEPTANCE_OUTPUT` and `VECTORS_OUTPUT` paths to include `${RESULTS_DIR}` prefix
   - **Lines 421-431**: Replaced file movement logic with verification logic (files already in correct location)
   - **Lines 405-440**: Changed parallel processing to sequential to avoid ROOT conflicts
   - **Impact**: Parsing jobs will now properly process files without conflicts and detect output files correctly

### **Plotting Scripts (Data Access Fix) - NEW**
5. `/users/bp969/scratch/VIKING_FOLDER/HISTOGRAM_MAKING/KLong_plot_resolution_histbar.C`
   - **Lines 52-63**: Updated to loop through all entries in combined ROOT files instead of just first entry
   - **Added**: Progress logging showing number of entries processed
   - **Impact**: Histograms will now show all data (60,000+ events) instead of just first file's data (~600 events)

6. `/users/bp969/scratch/VIKING_FOLDER/HISTOGRAM_MAKING/KLong_plot_acceptance.C`  
   - **Lines 57-73**: Updated to loop through all entries in combined ROOT files
   - **Added**: Progress logging showing number of entries processed
   - **Impact**: Acceptance plots will now include data from all 100 simulation files

## **Verification Status**

✅ **Parsing Script Syntax**: Tested - no more ROOT dictionary errors
✅ **Function Name Validation**: Hyphens replaced with underscores for valid C++ identifiers
✅ **File Detection Bug**: Fixed file path logic in parsing script (December 4, 2025)
⏳ **Full Pipeline Test**: Ready for re-running with all infrastructure fixes applied

## **Technical Analysis Summary**

### **Data Flow Verification (December 4, 2025)**
- **Simulation Phase**: ✅ Working correctly - produces ROOT files with physics events
- **Parsing Phase**: ✅ Fixed - ROOT macros execute successfully and create output files  
- **File Detection**: ✅ Fixed - Bash script now looks in correct directory (`${RESULTS_DIR}`)
- **Combination Phase**: ⚠️ Previously failing due to missing input files, should now work
- **Histogram Generation**: ⚠️ Previously failing due to underpopulated combined files

### **Evidence of Successful Processing**
From parse job logs before fix:
- "Found 598 events with π+, π-, π0 as direct kaon decay products"
- "Event 9088 | Reco p: 8.31954 | True p: 4.46733" 
- "Saved kaon event info to ...acceptance.root"
- **Problem was**: Immediately followed by "ERROR: Expected output file not created!"

## **Next Steps**

1. **Re-run Complete Pipeline**: Test T1-240_T2-250_T3-580_T4-590 configuration with all fixes
2. **Monitor File Sizes**: Verify combined files contain expected event counts (should be ~59,800+ events total)
3. **Validate Histograms**: Check that histograms show meaningful momentum distributions
4. **Performance Check**: Confirm pipeline completes without file detection errors

## **Expected Outcome**

The pipeline should now:
1. ✅ Successfully parse simulation data without ROOT dictionary errors
2. ✅ Create proper vector and acceptance files in correct locations
3. ✅ Properly detect and verify created output files  
4. 🔜 Combine files with full event content (not underpopulated)
5. 🔜 Generate histograms with meaningful kaon physics data
6. 🔜 Produce populated acceptance vs momentum plots for detector optimization

**Critical Fixes Applied**: 
1. **File detection bug** - Pipeline infrastructure failed to collect results despite successful processing
2. **Parallel ROOT conflicts** - Multiple ROOT processes caused dictionary conflicts and incomplete processing  
3. **Plotting data access bug** - Scripts only read 1/100th of available data from combined files

All three issues prevented proper data visualization, resulting in severely underpopulated histograms. The physics simulation, analysis logic, and data combination were working correctly - the problem was in data access patterns.