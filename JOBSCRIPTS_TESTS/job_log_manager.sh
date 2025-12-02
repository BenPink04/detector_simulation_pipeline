#!/bin/bash

# Job Log Management Script
# Provides utilities for monitoring and managing organized job logs

CONFIG=""
COMMAND=""

# Function to display usage
show_usage() {
    echo "=== Job Log Management Utility ==="
    echo "Usage: $0 <command> [configuration]"
    echo ""
    echo "Commands:"
    echo "  list [CONFIG]           - List all log directories or specific configuration"
    echo "  monitor CONFIG          - Monitor active jobs for configuration"  
    echo "  summary CONFIG          - Show job completion summary"
    echo "  errors CONFIG           - Show recent errors for configuration"
    echo "  clean [CONFIG]          - Clean old log files (older than 7 days)"
    echo "  tail CONFIG TYPE        - Tail latest log of type (build|simulation|parsing|combination|histograms)"
    echo ""
    echo "Examples:"
    echo "  $0 list"
    echo "  $0 monitor T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600"
    echo "  $0 summary T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600"
    echo "  $0 tail T1-240_T2-250_T3-570_T4-580_P1-215_P2-230_F1-260_F2-270_E1-600 simulation"
}

# Function to list log directories
list_logs() {
    local config=$1
    echo "=== Job Log Directories ==="
    
    if [ -z "$config" ]; then
        echo "Available configurations:"
        ls -la /users/bp969/scratch/JOB_LOGS/ 2>/dev/null || echo "No log directories found"
    else
        local log_dir="/users/bp969/scratch/JOB_LOGS/${config}"
        if [ -d "$log_dir" ]; then
            echo "Log structure for configuration: $config"
            echo "📁 $log_dir"
            for subdir in build simulation parsing combination histograms; do
                if [ -d "$log_dir/$subdir" ]; then
                    local count=$(ls -1 "$log_dir/$subdir"/*.log 2>/dev/null | wc -l)
                    echo "  📂 $subdir/ ($count log files)"
                fi
            done
        else
            echo "❌ No logs found for configuration: $config"
        fi
    fi
}

# Function to monitor active jobs
monitor_jobs() {
    local config=$1
    if [ -z "$config" ]; then
        echo "❌ Configuration required for monitoring"
        return 1
    fi
    
    echo "=== Active Jobs for Configuration: $config ==="
    echo "🔍 Checking SLURM queue..."
    
    # Check for jobs matching this configuration
    squeue -u bp969 --format="%.8i %.12j %.8T %.10M %.6D %R" | grep "$config" || echo "No active jobs found for $config"
    
    echo ""
    echo "📊 Job Status Summary:"
    local total_jobs=$(squeue -u bp969 | grep -c "$config" || echo "0")
    local running_jobs=$(squeue -u bp969 | grep "$config" | grep -c "RUNNING" || echo "0")
    local pending_jobs=$(squeue -u bp969 | grep "$config" | grep -c "PENDING" || echo "0")
    
    echo "  Total jobs: $total_jobs"
    echo "  Running: $running_jobs"
    echo "  Pending: $pending_jobs"
}

# Function to show job completion summary
show_summary() {
    local config=$1
    if [ -z "$config" ]; then
        echo "❌ Configuration required for summary"
        return 1
    fi
    
    local log_dir="/users/bp969/scratch/JOB_LOGS/${config}"
    
    echo "=== Job Completion Summary for: $config ==="
    
    if [ ! -d "$log_dir" ]; then
        echo "❌ No log directory found for configuration: $config"
        return 1
    fi
    
    for job_type in build simulation parsing combination histograms; do
        if [ -d "$log_dir/$job_type" ]; then
            local log_count=$(ls -1 "$log_dir/$job_type"/*.log 2>/dev/null | wc -l)
            local err_count=$(ls -1 "$log_dir/$job_type"/*.err 2>/dev/null | wc -l)
            local error_files=$(find "$log_dir/$job_type" -name "*.err" -size +0c 2>/dev/null | wc -l)
            
            printf "📊 %-12s: %3d logs, %3d errors" "$job_type" "$log_count" "$error_files"
            
            if [ "$error_files" -gt 0 ]; then
                echo " ❌"
            else
                echo " ✅"
            fi
        fi
    done
}

# Function to show recent errors
show_errors() {
    local config=$1
    if [ -z "$config" ]; then
        echo "❌ Configuration required for error checking"
        return 1
    fi
    
    local log_dir="/users/bp969/scratch/JOB_LOGS/${config}"
    
    echo "=== Recent Errors for Configuration: $config ==="
    
    if [ ! -d "$log_dir" ]; then
        echo "❌ No log directory found for configuration: $config"
        return 1
    fi
    
    # Find error files with content
    find "$log_dir" -name "*.err" -size +0c -exec ls -la {} \; 2>/dev/null | head -10
    
    echo ""
    echo "🔍 Sample error content:"
    find "$log_dir" -name "*.err" -size +0c -exec head -10 {} \; 2>/dev/null | head -20
}

# Function to tail latest log
tail_log() {
    local config=$1
    local job_type=$2
    
    if [ -z "$config" ] || [ -z "$job_type" ]; then
        echo "❌ Both configuration and job type required"
        echo "Job types: build, simulation, parsing, combination, histograms"
        return 1
    fi
    
    local log_dir="/users/bp969/scratch/JOB_LOGS/${config}/${job_type}"
    
    if [ ! -d "$log_dir" ]; then
        echo "❌ No log directory found: $log_dir"
        return 1
    fi
    
    # Find most recent log file
    local latest_log=$(ls -t "$log_dir"/*.log 2>/dev/null | head -1)
    
    if [ -z "$latest_log" ]; then
        echo "❌ No log files found in $log_dir"
        return 1
    fi
    
    echo "=== Tailing latest $job_type log ==="
    echo "📄 File: $latest_log"
    echo ""
    tail -f "$latest_log"
}

# Function to clean old logs
clean_logs() {
    local config=$1
    local days=7
    
    echo "=== Cleaning Old Log Files ==="
    
    if [ -z "$config" ]; then
        echo "🧹 Cleaning all configurations (files older than $days days)..."
        find /users/bp969/scratch/JOB_LOGS -name "*.log" -o -name "*.err" -mtime +$days -exec rm -f {} \;
    else
        local log_dir="/users/bp969/scratch/JOB_LOGS/${config}"
        if [ -d "$log_dir" ]; then
            echo "🧹 Cleaning configuration: $config (files older than $days days)..."
            find "$log_dir" -name "*.log" -o -name "*.err" -mtime +$days -exec rm -f {} \;
        else
            echo "❌ No log directory found for configuration: $config"
            return 1
        fi
    fi
    
    echo "✅ Cleanup complete"
}

# Parse command line arguments
if [ $# -eq 0 ]; then
    show_usage
    exit 1
fi

COMMAND=$1
CONFIG=$2
JOB_TYPE=$3

# Execute requested command
case $COMMAND in
    "list")
        list_logs "$CONFIG"
        ;;
    "monitor")
        monitor_jobs "$CONFIG"
        ;;
    "summary")
        show_summary "$CONFIG"
        ;;
    "errors")
        show_errors "$CONFIG"
        ;;
    "tail")
        tail_log "$CONFIG" "$JOB_TYPE"
        ;;
    "clean")
        clean_logs "$CONFIG"
        ;;
    *)
        echo "❌ Unknown command: $COMMAND"
        show_usage
        exit 1
        ;;
esac