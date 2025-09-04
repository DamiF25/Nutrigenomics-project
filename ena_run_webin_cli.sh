#!/bin/bash

# === CONFIGURATION ===
JAVA_HOME="/mnt/scratch2/users/40319381/java/jdk-17.0.8+7-jre"          # where you installed Java
WEBIN_CLI_JAR="/mnt/scratch2/users/40319381/webin-cli/webin-cli-9.0.1.jar"   # path to Webin-CLI jar
INPUT_DIR="/mnt/scratch2/users/40319381/New_Kraken/Read_Files/Read_Files"      # folder containing FASTQs and manifests
WEBIN_USER="Webin-xxxxx"                        # your Webin username
WEBIN_PASS='xxxxxxxxxxxx'                    # your Webin password

# === ACTION: validate or submit ===
ACTION="-validate"
# To submit after successful validation, change ACTION to "-submit"

cd "$INPUT_DIR" || { echo "ERROR: Input directory not found: $INPUT_DIR"; exit 1; }

for manifest in *_manifest.txt; do
    echo "=== Processing manifest: $manifest ==="
    "$JAVA_HOME/bin/java" -jar "$WEBIN_CLI_JAR" \
        -context reads \
        -manifest "$manifest" \
        -inputDir . \
        $ACTION \
        -userName "$WEBIN_USER" \
        -password "$WEBIN_PASS"

    if [ $? -ne 0 ]; then
        echo "ERROR: Webin-CLI failed for $manifest"
    else
        echo "SUCCESS: $manifest processed"
    fi
    echo
done
