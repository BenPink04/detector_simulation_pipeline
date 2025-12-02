#!/bin/bash

# Comprehensive Geometry Overlap Fix for All Pipeline Configurations
# This script provides multiple strategies to fix geometry overlaps across all configurations

set -e

show_usage() {
    echo "Usage: $0 [OPTION]"
    echo ""
    echo "Options:"
    echo "  --fix-templates     Fix all template files (recommended first step)"
    echo "  --fix-running ALL   Fix all existing running configurations"
    echo "  --fix-running NAME  Fix specific running configuration"
    echo "  --list-configs      List all available configurations"
    echo "  --check-overlaps    Check for overlap issues in templates and configs"
    echo "  --help             Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 --fix-templates"
    echo "  $0 --fix-running T1-240_T2-250_T3-280_T4-290_P1-215_P2-230_F1-260_F2-270_E1-600"
    echo "  $0 --fix-running ALL"
    echo "  $0 --check-overlaps"
}

# Function to apply geometry fixes to a file
apply_geometry_fixes() {
    local detector_file="$1"
    local name="$2"
    
    if [ ! -f "$detector_file" ]; then
        echo "  ❌ File not found: $detector_file"
        return 1
    fi
    
    echo "  🔧 Processing: $name"
    
    # Create backup
    local backup_file="${detector_file}.backup_$(date +%Y%m%d_%H%M%S)"
    cp "$detector_file" "$backup_file"
    echo "    📁 Backup: $(basename "$backup_file")"
    
    # Count initial issues
    local initial_duplicates=$(grep -o "\-20\.4\*cm" "$detector_file" | wc -l)
    
    # Fix 1: Remove duplicate -20.4*cm StrawTube position
    sed -i 's/-20\.4\*cm, -20\.4\*cm/-20.4*cm/g' "$detector_file"
    
    # Fix 2: Adjust pizza slice angles to 14.9° (from 15.0°)
    sed -i 's/15\.0\*deg/14.9*deg/g' "$detector_file"
    sed -i 's/15\*deg/14.9*deg/g' "$detector_file"
    
    # Verify fixes
    local final_duplicates=$(grep -o "\-20\.4\*cm" "$detector_file" | wc -l)
    
    if [ $initial_duplicates -gt 1 ] && [ $final_duplicates -eq 1 ]; then
        echo "    ✅ Fixed StrawTube duplicates ($initial_duplicates → $final_duplicates)"
    elif [ $initial_duplicates -gt 1 ]; then
        echo "    ⚠️  Partial fix: StrawTube duplicates ($initial_duplicates → $final_duplicates)"
    else
        echo "    ✅ StrawTube positions already clean"
    fi
    
    if ! grep -q "15\.0\*deg\|15\*deg" "$detector_file"; then
        echo "    ✅ Pizza slice angles optimized to 14.9°"
    else
        echo "    ⚠️  Some pizza slice angles may still need adjustment"
    fi
    
    return 0
}

# Function to check for overlap issues
check_overlaps() {
    local detector_file="$1"
    local name="$2"
    
    if [ ! -f "$detector_file" ]; then
        echo "  ❌ Not found: $name"
        return 1
    fi
    
    local duplicates=$(grep -o "\-20\.4\*cm" "$detector_file" | wc -l)
    local has_15deg=$(grep -c "15\.0\*deg" "$detector_file" 2>/dev/null || echo 0)
    
    if [ $duplicates -gt 1 ] || [ $has_15deg -gt 0 ]; then
        echo "  ⚠️  $name: Issues detected"
        [ $duplicates -gt 1 ] && echo "     - StrawTube duplicates: $duplicates occurrences of -20.4*cm"
        [ $has_15deg -gt 0 ] && echo "     - Pizza slice angles: $has_15deg occurrences of 15.0°"
        return 1
    else
        echo "  ✅ $name: Clean geometry"
        return 0
    fi
}

# Main execution based on command line argument
case "${1:-}" in
    --fix-templates)
        echo "=== Fixing All Template Files ==="
        echo "This will apply geometry fixes to all base template files"
        echo ""
        
        # Define template directories
        TEMPLATE_DIRS=(
            "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/SCENARIO_3_SIM"
            "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/SCENARIO_4_SIM" 
            "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/SCENARIO_5_SIM"
        )
        
        fixed_count=0
        total_count=0
        
        for template_dir in "${TEMPLATE_DIRS[@]}"; do
            detector_file="$template_dir/src/JLabKDetectorConstruction.cc"
            template_name=$(basename "$template_dir")
            
            total_count=$((total_count + 1))
            
            if apply_geometry_fixes "$detector_file" "$template_name"; then
                fixed_count=$((fixed_count + 1))
            fi
            echo ""
        done
        
        echo "📊 Summary: Fixed $fixed_count/$total_count templates"
        echo "✅ Future pipeline runs will use clean geometry automatically"
        ;;
        
    --fix-running)
        if [ -z "$2" ]; then
            echo "Error: Please specify configuration name or 'ALL'"
            echo "Usage: $0 --fix-running [CONFIG_NAME|ALL]"
            exit 1
        fi
        
        if [ "$2" = "ALL" ]; then
            echo "=== Fixing All Running Configurations ==="
            echo ""
            
            # Find all configuration directories
            CONFIG_BASE="/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING"
            CONFIGS=($(find "$CONFIG_BASE" -maxdepth 1 -name "T1-*" -type d | sort))
            
            if [ ${#CONFIGS[@]} -eq 0 ]; then
                echo "❌ No running configurations found in $CONFIG_BASE"
                exit 1
            fi
            
            echo "Found ${#CONFIGS[@]} configurations to fix..."
            echo ""
            
            fixed_count=0
            
            for config_dir in "${CONFIGS[@]}"; do
                config_name=$(basename "$config_dir")
                detector_file="$config_dir/src/JLabKDetectorConstruction.cc"
                
                if apply_geometry_fixes "$detector_file" "$config_name"; then
                    fixed_count=$((fixed_count + 1))
                fi
                echo ""
            done
            
            echo "📊 Summary: Fixed $fixed_count/${#CONFIGS[@]} configurations"
            
        else
            # Fix specific configuration
            CONFIG_NAME="$2"
            CONFIG_DIR="/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/$CONFIG_NAME"
            DETECTOR_FILE="$CONFIG_DIR/src/JLabKDetectorConstruction.cc"
            
            echo "=== Fixing Specific Configuration ==="
            echo "Configuration: $CONFIG_NAME"
            echo ""
            
            apply_geometry_fixes "$DETECTOR_FILE" "$CONFIG_NAME"
        fi
        ;;
        
    --list-configs)
        echo "=== Available Configurations ==="
        CONFIG_BASE="/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING"
        
        echo "Templates:"
        find "$CONFIG_BASE" -maxdepth 1 -name "SCENARIO_*" -type d | sort | while read dir; do
            echo "  📁 $(basename "$dir")"
        done
        
        echo ""
        echo "Running Configurations:"
        find "$CONFIG_BASE" -maxdepth 1 -name "T1-*" -type d | sort | while read dir; do
            config_name=$(basename "$dir")
            if [ -f "$dir/src/JLabKDetectorConstruction.cc" ]; then
                echo "  🔧 $config_name"
            else
                echo "  ❓ $config_name (no geometry file)"
            fi
        done
        ;;
        
    --check-overlaps)
        echo "=== Checking for Geometry Overlap Issues ==="
        echo ""
        
        # Check templates
        echo "📁 Template Files:"
        TEMPLATE_DIRS=(
            "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/SCENARIO_3_SIM"
            "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/SCENARIO_4_SIM"
            "/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING/SCENARIO_5_SIM"
        )
        
        template_issues=0
        for template_dir in "${TEMPLATE_DIRS[@]}"; do
            detector_file="$template_dir/src/JLabKDetectorConstruction.cc"
            template_name=$(basename "$template_dir")
            
            if ! check_overlaps "$detector_file" "$template_name"; then
                template_issues=$((template_issues + 1))
            fi
        done
        
        echo ""
        echo "🔧 Running Configurations:"
        
        # Check running configs
        CONFIG_BASE="/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RUNNING"
        CONFIGS=($(find "$CONFIG_BASE" -maxdepth 1 -name "T1-*" -type d | sort))
        
        config_issues=0
        for config_dir in "${CONFIGS[@]}"; do
            config_name=$(basename "$config_dir")
            detector_file="$config_dir/src/JLabKDetectorConstruction.cc"
            
            if ! check_overlaps "$detector_file" "$config_name"; then
                config_issues=$((config_issues + 1))
            fi
        done
        
        echo ""
        echo "📊 Summary:"
        echo "  Templates with issues: $template_issues/${#TEMPLATE_DIRS[@]}"
        echo "  Configs with issues: $config_issues/${#CONFIGS[@]}"
        
        if [ $template_issues -gt 0 ] || [ $config_issues -gt 0 ]; then
            echo ""
            echo "💡 Recommended actions:"
            [ $template_issues -gt 0 ] && echo "  1. Run: $0 --fix-templates"
            [ $config_issues -gt 0 ] && echo "  2. Run: $0 --fix-running ALL"
        else
            echo "  ✅ All files have clean geometry!"
        fi
        ;;
        
    --help|"")
        show_usage
        ;;
        
    *)
        echo "Error: Unknown option '$1'"
        echo ""
        show_usage
        exit 1
        ;;
esac