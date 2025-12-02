# Pipeline Issue Resolution Report

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

## **Verification Status**

✅ **Parsing Script Syntax**: Tested - no more ROOT dictionary errors
✅ **Function Name Validation**: Hyphens replaced with underscores for valid C++ identifiers  
⏳ **Full Pipeline Test**: Ready for re-running with fixed infrastructure

## **Next Steps**

1. **Test Fixed Pipeline**: Re-run a complete pipeline with one detector configuration
2. **Verify Histogram Generation**: Ensure histograms are created successfully without compilation errors
3. **Check Data Flow**: Confirm that parsing creates output files and combination works correctly

## **Expected Outcome**

The pipeline should now:
1. Successfully parse simulation data without ROOT dictionary errors
2. Create proper vector and acceptance files 
3. Combine files correctly
4. Generate histograms without C++ compilation failures
5. Produce meaningful physics analysis results

The empty histogram issue was caused by infrastructure problems, not physics analysis issues. With these fixes, the pipeline should produce populated histograms showing kaon reconstruction acceptance vs momentum.