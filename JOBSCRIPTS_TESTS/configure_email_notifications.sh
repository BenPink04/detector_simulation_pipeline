#!/bin/bash

# Enhanced Email Notification Options for Detector Pipeline
# This script allows you to customize email notifications

# Email notification options:
# NONE     - No email notifications
# END      - Email when job ends (success only)
# FAIL     - Email when job fails
# END_FAIL - Email when job ends or fails (recommended)
# ALL      - Email for all job state changes (BEGIN, END, FAIL, REQUEUE, etc.)

# Set your email notification preference here:
MAIL_TYPE="END,FAIL"  # Change this to your preference
MAIL_USER="bp969@york.ac.uk"  # Change this to your email address

echo "=== Email Notification Configuration ==="
echo "Mail Type: $MAIL_TYPE"
echo "Mail User: $MAIL_USER"
echo

# Available options:
echo "Available email notification options:"
echo "  NONE     - No emails"
echo "  END      - Email only when job completes successfully"
echo "  FAIL     - Email only when job fails"
echo "  END,FAIL - Email when job ends (success or failure) - RECOMMENDED"
echo "  ALL      - Email for all job state changes (can be many emails)"
echo

# Function to update email settings in a file
update_email_settings() {
    local file="$1"
    local job_type="$2"
    
    if [ ! -f "$file" ]; then
        echo "Warning: File $file not found"
        return 1
    fi
    
    # Check if file already has email settings
    if grep -q "#SBATCH --mail-type" "$file"; then
        # Update existing settings
        sed -i "s/#SBATCH --mail-type=.*/#SBATCH --mail-type=$MAIL_TYPE/" "$file"
        sed -i "s/#SBATCH --mail-user=.*/#SBATCH --mail-user=$MAIL_USER/" "$file"
        echo "Updated email settings in $file ($job_type)"
    else
        echo "Warning: No email settings found in $file"
    fi
}

# Check if scripts exist and update them
MAIN_SCRIPT="/users/bp969/scratch/JOBSCRIPTS_TESTS/detector_simulation_master.sh"
TEST_SCRIPT="/users/bp969/scratch/JOBSCRIPTS_TESTS/test_detector_pipeline.sh"

if [ -f "$MAIN_SCRIPT" ]; then
    echo "Updating email settings in main pipeline script..."
    # The main script generates job scripts dynamically, so we need to update the template
    sed -i "s/#SBATCH --mail-type=.*/#SBATCH --mail-type=$MAIL_TYPE/" "$MAIN_SCRIPT"
    sed -i "s/#SBATCH --mail-user=.*/#SBATCH --mail-user=$MAIL_USER/" "$MAIN_SCRIPT"
    echo "✓ Main script updated"
else
    echo "✗ Main script not found: $MAIN_SCRIPT"
fi

if [ -f "$TEST_SCRIPT" ]; then
    echo "Updating email settings in test script..."
    sed -i "s/#SBATCH --mail-type=.*/#SBATCH --mail-type=$MAIL_TYPE/" "$TEST_SCRIPT"
    sed -i "s/#SBATCH --mail-user=.*/#SBATCH --mail-user=$MAIL_USER/" "$TEST_SCRIPT"
    echo "✓ Test script updated"
else
    echo "✗ Test script not found: $TEST_SCRIPT"
fi

# Update any existing generated job scripts
echo
echo "Checking for existing generated job scripts..."
JOB_DIR="/users/bp969/scratch/JOBSCRIPTS_TESTS/temp_jobs"
if [ -d "$JOB_DIR" ]; then
    JOB_COUNT=$(find "$JOB_DIR" -name "*.job" | wc -l)
    if [ $JOB_COUNT -gt 0 ]; then
        echo "Found $JOB_COUNT existing job scripts. Updating email settings..."
        find "$JOB_DIR" -name "*.job" -exec sed -i "s/#SBATCH --mail-type=.*/#SBATCH --mail-type=$MAIL_TYPE/" {} \;
        find "$JOB_DIR" -name "*.job" -exec sed -i "s/#SBATCH --mail-user=.*/#SBATCH --mail-user=$MAIL_USER/" {} \;
        echo "✓ Updated $JOB_COUNT job scripts"
    else
        echo "No existing job scripts found"
    fi
else
    echo "No temp_jobs directory found"
fi

echo
echo "=== Email Configuration Complete ==="
echo
echo "Your jobs will now send emails with the following settings:"
echo "  When: $MAIL_TYPE"
echo "  To: $MAIL_USER"
echo
echo "Email content will include:"
echo "  - Job ID and name"
echo "  - Start and end times"
echo "  - Exit status"
echo "  - Resource usage summary"
echo
echo "To change these settings in the future:"
echo "1. Edit this script to change MAIL_TYPE and MAIL_USER variables"
echo "2. Run this script again to apply changes"
echo
echo "Recommended setting: END,FAIL (notifies on completion and failures only)"