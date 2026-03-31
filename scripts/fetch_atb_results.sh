#!/usr/bin/env bash
# fetch_atb_results.sh — Fetch LexicMap results from EC2 and terminate instance
#
# Usage:
#   bash scripts/fetch_atb_results.sh <INSTANCE_ID> <PUBLIC_IP>
#
# Or reads from DB/atb_lexicmap_results/instance_info.txt if no args given.

set -euo pipefail

RESULTS_DIR="DB/atb_lexicmap_results"
REGION="eu-west-2"
KEY_FILE="$HOME/.ssh/atb_search_key"

if [ $# -eq 2 ]; then
    INSTANCE_ID="$1"
    PUBLIC_IP="$2"
else
    source "$RESULTS_DIR/instance_info.txt"
fi

echo "Fetching results from $PUBLIC_IP (instance $INSTANCE_ID)..."

# Fetch results and log
scp -i "$KEY_FILE" -o StrictHostKeyChecking=no \
    "ec2-user@${PUBLIC_IP}:results.tsv.gz" \
    "$RESULTS_DIR/atb_lexicmap_raw.tsv.gz"

scp -i "$KEY_FILE" -o StrictHostKeyChecking=no \
    "ec2-user@${PUBLIC_IP}:/home/ec2-user/lexicmap_run.log" \
    "$RESULTS_DIR/lexicmap_run.log"

echo "Results saved to $RESULTS_DIR/atb_lexicmap_raw.tsv.gz"

# Terminate instance
echo "Terminating instance $INSTANCE_ID..."
aws ec2 terminate-instances --region "$REGION" --instance-ids "$INSTANCE_ID" > /dev/null
echo "Instance terminated."

echo ""
echo "Next: python scripts/analyse_atb_results.py"
