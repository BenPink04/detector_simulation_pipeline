#!/bin/bash

echo "Applying comprehensive fix to all parsing macros..."

# Use the backup files we created earlier as the source
find VIKING_FOLDER/DATA_PARSING -name "*.C.backup" | while read backup_file; do
    original_file="${backup_file%.backup}"
    echo "Restoring and fixing $original_file from backup"
    
    # Copy from backup and apply the correct fix
    cp "$backup_file" "$original_file"
    
    # Now apply the correct verbosity fix - comment out entire problematic sections
    sed -i '/std::cout << "Selected Event/,/std::endl;/s/^/\/\/ /' "$original_file"
    
    # Also comment out any standalone cout lines with Reco p:
    sed -i 's/std::cout << " | Reco p:/\/\/ std::cout << " | Reco p:/' "$original_file"
done

# Also fix the main parsing files (not in config directories)
for main_file in VIKING_FOLDER/DATA_PARSING/KLong_save_momentum_acceptance.C VIKING_FOLDER/DATA_PARSING/KLong_save_vectors.C; do
    if [ -f "${main_file}.backup" ]; then
        echo "Restoring and fixing $main_file from backup"
        cp "${main_file}.backup" "$main_file"
        sed -i '/std::cout << "Selected Event/,/std::endl;/s/^/\/\/ /' "$main_file"
        sed -i 's/std::cout << " | Reco p:/\/\/ std::cout << " | Reco p:/' "$main_file"
    fi
done

echo "Comprehensive fix applied!"
echo "Verifying syntax..."