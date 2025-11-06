#!/usr/bin/env bash

# Script to fix the verbosity issue in parsing macros
echo "Fixing parsing macro verbosity issues..."

# Create backup copies first
cp VIKING_FOLDER/DATA_PARSING/KLong_save_momentum_acceptance.C VIKING_FOLDER/DATA_PARSING/KLong_save_momentum_acceptance.C.backup
cp VIKING_FOLDER/DATA_PARSING/KLong_save_vectors.C VIKING_FOLDER/DATA_PARSING/KLong_save_vectors.C.backup

echo "Created backups..."

# Comment out the verbose per-event output in acceptance macro
sed -i 's/std::cout << "Selected Event #"/\/\/ std::cout << "Selected Event #"/' VIKING_FOLDER/DATA_PARSING/KLong_save_momentum_acceptance.C
sed -i 's/std::cout << " | Reco p: " << kaon_p << std::endl;/\/\/ std::cout << " | Reco p: " << kaon_p << std::endl;/' VIKING_FOLDER/DATA_PARSING/KLong_save_momentum_acceptance.C
sed -i 's/std::cout << " | Reco p: not reconstructable" << std::endl;/\/\/ std::cout << " | Reco p: not reconstructable" << std::endl;/' VIKING_FOLDER/DATA_PARSING/KLong_save_momentum_acceptance.C

# Do the same for vectors macro if it has similar issues
if grep -q "std::cout.*Event" VIKING_FOLDER/DATA_PARSING/KLong_save_vectors.C; then
    sed -i 's/std::cout.*Event.*/\/\/ &/' VIKING_FOLDER/DATA_PARSING/KLong_save_vectors.C
fi

echo "Fixed verbosity issues in parsing macros"
echo "Summary output will still be printed, but per-event output is now commented out"