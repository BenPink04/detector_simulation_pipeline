// ============================================================================
// KLong_visualize_config.C — multi-seed event display for K_L -> pi+pi-pi0
//
// PURPOSE:
//   Merges all raw seed files in a configuration directory into a single ROOT
//   file using TFileMerger (equivalent to hadd), then runs the full
//   KLong_visualize_event reconstruction on that merged file.
//
//   The merged file is saved as <config_dir>/<stem>_raw_merged.root and is
//   REUSED on subsequent calls (delete it manually to force a re-merge).
//
// INPUT:  path to a configuration directory, e.g.:
//           .../SIMULATION_RESULTS/T1-240_..._E1-700/
//         Must contain raw *.root files (not *_vectors*, *_acceptance*,
//         *_combined*, *_raw_merged*).
// OUTPUT: <stem>_all_seeds_delta_p.txt  — all events sorted by |delta_p|
//         KLong_event<N>_view.png        — 4-panel figure for event N (optional)
//         <stem>_raw_merged.root         — cached merged file (in config_dir)
//
// USAGE (from VISUALISER/ directory):
//   root -l -q 'KLong_visualize_config.C("/path/to/config_dir")'
//   root -l -q 'KLong_visualize_config.C("/path/to/config_dir", 46480)'
//   root -l -q 'KLong_visualize_config.C("/path/to/config_dir", 46480, "out.txt")'
// ============================================================================

#include "KLong_visualize_event.C"

void KLong_visualize_config(
    const char* config_dir   = ".",
    int         target_event = -1,
    const char* output_txt   = "")
{
    // --- Collect all raw seed files from the directory ---
    TSystemDirectory tdir("d", config_dir);
    TList* files = tdir.GetListOfFiles();
    if (!files) { std::cout << "Cannot list directory: " << config_dir << "\n"; return; }

    std::vector<std::string> seed_files;
    TIter next(files);
    TObject* obj;
    while ((obj = next())) {
        std::string name(obj->GetName());
        if (name.size() < 5 || name.substr(name.size()-5) != ".root") continue;
        if (name.find("_vectors")    != std::string::npos) continue;
        if (name.find("_acceptance") != std::string::npos) continue;
        if (name.find("_combined")   != std::string::npos) continue;
        if (name.find("_raw_merged") != std::string::npos) continue;
        seed_files.push_back(std::string(config_dir) + "/" + name);
    }
    std::sort(seed_files.begin(), seed_files.end());

    if (seed_files.empty()) {
        std::cout << "No raw seed files found in: " << config_dir << "\n";
        std::cout << "(Looking for *.root not matching *_vectors*, *_acceptance*, "
                     "*_combined*, *_raw_merged*)\n";
        return;
    }
    std::cout << "Found " << seed_files.size() << " seed files in " << config_dir << "\n";

    // Directory base name (for filenames and geometry parsing)
    auto dir_basename = [](const std::string& p) -> std::string {
        size_t sl = p.find_last_of("/\\");
        return (sl == std::string::npos) ? p : p.substr(sl+1);
    };
    std::string stem = dir_basename(std::string(config_dir));

    // Merged file lives alongside the seed files
    std::string merged_path = std::string(config_dir) + "/" + stem + "_raw_merged.root";

    // --- Merge seed files with TFileMerger (skip if already done) ---
    bool need_merge = true;
    {
        TFile* probe = TFile::Open(merged_path.c_str(), "READ");
        if (probe && !probe->IsZombie()) { need_merge = false; probe->Close(); delete probe; }
    }
    if (!need_merge) {
        std::cout << "Reusing cached merged file: " << merged_path << "\n";
        std::cout << "  (delete it to force a re-merge)\n";
    } else {
        std::cout << "Merging " << seed_files.size() << " seed files -> "
                  << merged_path << "\n";
        TFileMerger merger(/*isLocal=*/kFALSE);
        merger.OutputFile(merged_path.c_str(), "RECREATE");
        for (const auto& sf : seed_files) merger.AddFile(sf.c_str());
        if (!merger.Merge()) {
            std::cout << "ERROR: TFileMerger failed.\n";
            return;
        }
        std::cout << "Merge complete.\n";
    }

    // Output text file: default to <stem>_all_seeds_delta_p.txt
    std::string txt_path = (std::string(output_txt).empty())
                           ? (stem + "_all_seeds_delta_p.txt") : std::string(output_txt);

    // --- Hand off to the standard single-file visualiser ---
    std::cout << "Running reconstruction on merged file...\n";
    KLong_visualize_event(merged_path.c_str(), target_event, txt_path.c_str());
}
