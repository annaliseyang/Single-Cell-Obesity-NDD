#!/bin/bash

# Define the folder to move files to
backup_dir="/home/anna_y/data/log"
completed_key="Job Completed"
cancelled_key="CANCELLED"

# Loop through all text files in the current log directory
for file in log/*.out; do
    if grep -q "$completed_key" "$file"; then
        # Move the file to the backup folder
        mv "$file" "$backup_dir"
        echo "Moved $file to $backup_dir"
    elif grep -q "$cancelled_key" "$file"; then
        rm $file
        echo "Deleted $file"
    fi
done
