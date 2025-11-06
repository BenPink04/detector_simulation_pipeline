#!/bin/bash

# Job Monitoring Script for Detector Pipeline

echo "=== Detector Pipeline Job Monitor ==="
echo "Monitoring jobs for user: $USER"
echo

# Function to check job status
check_job_status() {
    local job_id=$1
    if [ -n "$job_id" ]; then
        echo "Checking job $job_id:"
        scontrol show job "$job_id" 2>/dev/null | grep -E "JobId|JobName|JobState|RunTime|ExitCode" || echo "Job $job_id not found (may have completed)"
        echo
    fi
}

# Function to show recent jobs
show_recent_jobs() {
    echo "Recent jobs for $USER:"
    squeue -u "$USER" --format="%.10i %.20j %.8T %.10M %.6D %.20S %.20e" || echo "No active jobs found"
    echo
}

# Function to show job logs
show_logs() {
    local pattern=$1
    echo "Recent job logs matching '$pattern':"
    find /users/bp969/scratch/JOBSCRIPTS_TESTS -name "*${pattern}*.log" -o -name "*${pattern}*.err" -mtime -1 | sort -t- -k2 -n | tail -10
    echo
}

# Function to tail active job logs
tail_logs() {
    local job_name=$1
    local log_file=$(find /users/bp969/scratch/JOBSCRIPTS_TESTS -name "*${job_name}*.log" -mtime -1 | sort | tail -1)
    local err_file=$(find /users/bp969/scratch/JOBSCRIPTS_TESTS -name "*${job_name}*.err" -mtime -1 | sort | tail -1)
    
    if [ -n "$log_file" ] || [ -n "$err_file" ]; then
        echo "=== Latest log files for '$job_name' ==="
        if [ -n "$log_file" ]; then
            echo "Log file: $log_file"
            echo "--- Last 20 lines ---"
            tail -20 "$log_file"
            echo
        fi
        if [ -n "$err_file" ] && [ -s "$err_file" ]; then
            echo "Error file: $err_file"
            echo "--- Last 10 lines ---"
            tail -10 "$err_file"
            echo
        fi
    else
        echo "No recent log files found for '$job_name'"
    fi
}

# Function to check test simulation status
check_test_status() {
    echo "=== Test Simulation Status ==="
    
    # Look for test job logs
    TEST_LOGS=$(find /users/bp969/scratch/JOBSCRIPTS_TESTS -name "TEST_Sim_*" -name "*.log" -o -name "*.err" | sort -t- -k3 -n | tail -5)
    
    if [ -n "$TEST_LOGS" ]; then
        echo "Found test simulation logs:"
        echo "$TEST_LOGS"
        echo
        
        # Check latest test
        LATEST_LOG=$(echo "$TEST_LOGS" | tail -1)
        echo "Latest test log: $LATEST_LOG"
        
        if [[ "$LATEST_LOG" == *.err ]]; then
            if [ -s "$LATEST_LOG" ]; then
                echo "❌ Test failed - error file has content:"
                tail -10 "$LATEST_LOG"
            else
                echo "✅ Test error file is empty (good sign)"
            fi
        fi
        
        # Look for corresponding .log file
        LOG_FILE="${LATEST_LOG%.err}.log"
        if [ -f "$LOG_FILE" ]; then
            echo "Checking output log:"
            if grep -q "SUCCESS: Test simulation completed" "$LOG_FILE"; then
                echo "✅ Test simulation completed successfully!"
            elif grep -q "ERROR" "$LOG_FILE"; then
                echo "❌ Test simulation had errors:"
                grep "ERROR" "$LOG_FILE"
            else
                echo "⏳ Test simulation may still be running or incomplete"
                echo "Last few lines:"
                tail -5 "$LOG_FILE"
            fi
        fi
        echo
        
        # Check for output files
        echo "Checking for output files..."
        RESULTS_DIR="/users/bp969/scratch/VIKING_FOLDER/SIMULATION_RESULTS"
        TEST_RESULTS=$(find "$RESULTS_DIR" -name "TEST_*" -name "*.root" 2>/dev/null | head -5)
        if [ -n "$TEST_RESULTS" ]; then
            echo "✅ Found test output files:"
            echo "$TEST_RESULTS"
        else
            echo "⏳ No test output files found yet"
        fi
        
    else
        echo "No test simulation logs found"
    fi
}

# Main menu
case "${1:-status}" in
    "status"|"")
        show_recent_jobs
        check_test_status
        ;;
    "job")
        if [ -n "$2" ]; then
            check_job_status "$2"
        else
            echo "Usage: $0 job <job_id>"
        fi
        ;;
    "logs")
        if [ -n "$2" ]; then
            tail_logs "$2"
        else
            show_logs "TEST"
        fi
        ;;
    "test")
        check_test_status
        ;;
    "watch")
        watch -n 30 "$0 status"
        ;;
    "clean-logs")
        # Remove job logs older than a given number of days (default: 7)
        DAYS=${2:-7}
        echo "Cleaning job logs older than ${DAYS} days in temp_jobs..."
        find /users/bp969/scratch/JOBSCRIPTS_TESTS -name "*.log" -o -name "*.err" -mtime +${DAYS} -delete || true
        echo "Cleanup complete."
        ;;
    *)
        echo "Usage: $0 [command] [options]"
        echo "Commands:"
        echo "  status          - Show current job status and test results (default)"
        echo "  job <job_id>    - Show details for specific job"
        echo "  logs [pattern]  - Show recent log files"
        echo "  test            - Check test simulation status"
        echo "  watch           - Monitor status continuously (refresh every 30s)"
        echo
        echo "Examples:"
        echo "  $0                    # Show current status"
        echo "  $0 job 25529883      # Check specific job"
        echo "  $0 logs TEST_Sim     # Show test simulation logs"
        echo "  $0 watch             # Monitor continuously"
        ;;
esac