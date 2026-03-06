#!/usr/bin/env bash
# =============================================================================
# run_atb_lexicmap_ec2.sh
#
# Launch a spot EC2 instance in eu-west-2, mount the AllTheBacteria LexicMap
# index from S3, search KL300/303/306/307 query sequences, and return results.
#
# Prerequisites (local):
#   aws configure         (valid credentials + region eu-west-2)
#   ssh-keygen            (creates ~/.ssh/atb_search_key if needed)
#
# Cost: ~$0.30–0.80 total for a c7g.8xlarge spot instance (~20–30 min)
#
# Usage:
#   bash scripts/run_atb_lexicmap_ec2.sh
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
REGION="eu-west-2"           # Same region as s3://allthebacteria-lexicmap
INSTANCE_TYPE="c7g.8xlarge"   # ARM64, 32 vCPU, 64 GiB — fits default spot vCPU limit
KEY_NAME="atb-search-key"
KEY_FILE="$HOME/.ssh/atb_search_key"
SECURITY_GROUP_NAME="atb-lexicmap-sg"
QUERY_FASTA="DB/kl_failing_query.fasta"
RESULTS_DIR="DB/atb_lexicmap_results"
REMOTE_RESULTS="/home/ec2-user/results.tsv.gz"

# Amazon Linux 2023 ARM64 in eu-west-2 (graviton2)
# Verify with: aws ec2 describe-images --region eu-west-2 \
#   --filters "Name=name,Values=al2023-ami-2023*-kernel-6*-arm64" \
#   "Name=owner-alias,Values=amazon" --query 'sort_by(Images,&CreationDate)[-1].ImageId'
AMI_ID="ami-0a3c640e8f9c74b61"   # Amazon Linux 2023 ARM64, eu-west-2 (update if needed)

mkdir -p "$RESULTS_DIR"

echo "============================================================"
echo " ATB LexicMap Search — KL300/303/306/307"
echo " Region: $REGION  |  Instance: $INSTANCE_TYPE (spot)"
echo "============================================================"

# ---------------------------------------------------------------------------
# Step 1: SSH key
# ---------------------------------------------------------------------------
if [ ! -f "$KEY_FILE" ]; then
    echo "[1] Generating SSH key: $KEY_FILE"
    ssh-keygen -t ed25519 -f "$KEY_FILE" -N "" -C "atb-search"
else
    echo "[1] SSH key exists: $KEY_FILE"
fi
PUB_KEY=$(cat "${KEY_FILE}.pub")

# Import key pair to EC2 (idempotent)
aws ec2 import-key-pair \
    --region "$REGION" \
    --key-name "$KEY_NAME" \
    --public-key-material "$(base64 < "${KEY_FILE}.pub")" \
    2>/dev/null || true
echo "    Key pair '$KEY_NAME' ready in $REGION"

# ---------------------------------------------------------------------------
# Step 2: Security group (SSH only)
# ---------------------------------------------------------------------------
echo "[2] Setting up security group..."
SG_ID=$(aws ec2 describe-security-groups \
    --region "$REGION" \
    --filters "Name=group-name,Values=$SECURITY_GROUP_NAME" \
    --query "SecurityGroups[0].GroupId" \
    --output text 2>/dev/null || echo "None")

if [ "$SG_ID" = "None" ] || [ -z "$SG_ID" ]; then
    SG_ID=$(aws ec2 create-security-group \
        --region "$REGION" \
        --group-name "$SECURITY_GROUP_NAME" \
        --description "Temporary SG for ATB LexicMap search" \
        --query "GroupId" --output text)
    aws ec2 authorize-security-group-ingress \
        --region "$REGION" \
        --group-id "$SG_ID" \
        --protocol tcp --port 22 --cidr 0.0.0.0/0
    echo "    Created security group: $SG_ID"
else
    echo "    Security group exists: $SG_ID"
fi

# ---------------------------------------------------------------------------
# Step 3: Verify AMI (offer to auto-detect if needed)
# ---------------------------------------------------------------------------
echo "[3] Verifying AMI $AMI_ID in $REGION..."
AMI_CHECK=$(aws ec2 describe-images \
    --region "$REGION" \
    --image-ids "$AMI_ID" \
    --query "Images[0].Name" --output text 2>/dev/null || echo "NOT_FOUND")
if [ "$AMI_CHECK" = "NOT_FOUND" ] || [ -z "$AMI_CHECK" ]; then
    echo "    AMI $AMI_ID not found — looking up latest Amazon Linux 2023 ARM64..."
    AMI_ID=$(aws ec2 describe-images \
        --region "$REGION" \
        --filters \
            "Name=name,Values=al2023-ami-2023*-kernel-6*-arm64" \
            "Name=owner-alias,Values=amazon" \
            "Name=state,Values=available" \
        --query "sort_by(Images, &CreationDate)[-1].ImageId" \
        --output text)
    echo "    Using AMI: $AMI_ID"
else
    echo "    AMI found: $AMI_CHECK"
fi

# ---------------------------------------------------------------------------
# Step 4: User-data script (runs on boot inside EC2)
# ---------------------------------------------------------------------------
# Embed query FASTA inline so we don't need a separate upload step
QUERY_B64=$(base64 < "$QUERY_FASTA")

USER_DATA=$(cat <<USERDATA_EOF
#!/bin/bash
set -euxo pipefail
exec > /var/log/lexicmap_run.log 2>&1

echo "=== Starting ATB LexicMap search ==="
cd /home/ec2-user

# Install dependencies
dnf install -y wget tar fuse fuse-libs amazon-cloudwatch-agent 2>/dev/null || true

# Install mountpoint-s3
wget -q "https://s3.amazonaws.com/mountpoint-s3-release/latest/arm64/mount-s3.rpm"
dnf install -y ./mount-s3.rpm || yum install -y ./mount-s3.rpm

# Install LexicMap v0.8.1 (ARM64)
wget -q https://github.com/shenwei356/LexicMap/releases/download/v0.8.1/lexicmap_linux_arm64.tar.gz
tar xzf lexicmap_linux_arm64.tar.gz
mv lexicmap /usr/local/bin/lexicmap
lexicmap version

# Mount ATB index from S3 (read-only, no credentials needed)
mkdir -p /mnt/atb_index
UNSTABLE_MOUNTPOINT_MAX_PREFETCH_WINDOW_SIZE=65536 \
    mount-s3 --read-only --prefix 202408/ allthebacteria-lexicmap \
    /mnt/atb_index --no-sign-request
echo "=== Index mounted ==="
ls /mnt/atb_index/

# Write query FASTA
echo "${QUERY_B64}" | base64 -d > /home/ec2-user/query.fasta
echo "=== Query sequences ==="
grep "^>" /home/ec2-user/query.fasta

# Run LexicMap search
# --min-qcov-per-genome 70: genome must cover >=70% of query locus
# --align-min-match-pident 90: >=90% nucleotide identity
# --top-n-genomes 0: return all matches (no cap)
echo "=== Running LexicMap search ==="
lexicmap search \
    -d /mnt/atb_index \
    /home/ec2-user/query.fasta \
    -o /home/ec2-user/results.tsv.gz \
    --align-min-match-pident 90 \
    --min-qcov-per-genome 70 \
    --top-n-genomes 0 \
    -j 32 \
    --debug
echo "=== Search complete ==="

# Signal completion
touch /home/ec2-user/DONE
USERDATA_EOF
)

# ---------------------------------------------------------------------------
# Step 5: Launch spot instance
# ---------------------------------------------------------------------------
echo "[4] Launching spot instance ($INSTANCE_TYPE)..."

LAUNCH_SPEC=$(cat <<LAUNCH_EOF
{
  "ImageId": "$AMI_ID",
  "InstanceType": "$INSTANCE_TYPE",
  "KeyName": "$KEY_NAME",
  "SecurityGroupIds": ["$SG_ID"],
  "BlockDeviceMappings": [
    {
      "DeviceName": "/dev/xvda",
      "Ebs": {
        "VolumeSize": 30,
        "VolumeType": "gp3",
        "DeleteOnTermination": true
      }
    }
  ],
  "UserData": "$(echo "$USER_DATA" | base64 | tr -d '\n')"
}
LAUNCH_EOF
)

INSTANCE_ID=$(aws ec2 run-instances \
    --region "$REGION" \
    --instance-market-options '{"MarketType":"spot","SpotOptions":{"SpotInstanceType":"one-time"}}' \
    --cli-input-json "$LAUNCH_SPEC" \
    --query "Instances[0].InstanceId" \
    --output text)

echo "    Instance launched: $INSTANCE_ID"

# ---------------------------------------------------------------------------
# Step 6: Wait for instance to be running and get IP
# ---------------------------------------------------------------------------
echo "[5] Waiting for instance to reach running state..."
aws ec2 wait instance-running --region "$REGION" --instance-ids "$INSTANCE_ID"

PUBLIC_IP=$(aws ec2 describe-instances \
    --region "$REGION" \
    --instance-ids "$INSTANCE_ID" \
    --query "Reservations[0].Instances[0].PublicIpAddress" \
    --output text)
echo "    Public IP: $PUBLIC_IP"

echo ""
echo "============================================================"
echo " Instance is running. The LexicMap search is running in the"
echo " background via user-data. Monitor with:"
echo ""
echo "   ssh -i $KEY_FILE -o StrictHostKeyChecking=no ec2-user@$PUBLIC_IP"
echo "   tail -f /var/log/lexicmap_run.log"
echo ""
echo " When search completes (/home/ec2-user/DONE exists), run:"
echo ""
echo "   bash scripts/fetch_atb_results.sh $INSTANCE_ID $PUBLIC_IP"
echo ""
echo " Or wait automatically with: --wait flag"
echo ""
echo " INSTANCE ID (save this): $INSTANCE_ID"
echo "============================================================"

# Save instance details for fetch script
cat > "$RESULTS_DIR/instance_info.txt" <<INFO_EOF
INSTANCE_ID=$INSTANCE_ID
PUBLIC_IP=$PUBLIC_IP
REGION=$REGION
KEY_FILE=$KEY_FILE
LAUNCHED=$(date -u)
INFO_EOF

echo ""
echo "[6] Waiting for search to complete (polls every 60s)..."
echo "    (Ctrl+C to stop waiting and monitor manually)"
echo ""

MAX_WAIT=7200   # 2 hours max
ELAPSED=0
while [ $ELAPSED -lt $MAX_WAIT ]; do
    sleep 60
    ELAPSED=$((ELAPSED + 60))

    # Check if DONE file exists
    DONE=$(ssh -i "$KEY_FILE" \
        -o StrictHostKeyChecking=no \
        -o ConnectTimeout=10 \
        -o BatchMode=yes \
        "ec2-user@$PUBLIC_IP" \
        "test -f /home/ec2-user/DONE && echo yes || echo no" 2>/dev/null || echo "connecting")

    if [ "$DONE" = "yes" ]; then
        echo "    [${ELAPSED}s] Search complete!"
        break
    else
        echo "    [${ELAPSED}s] Still running... ($DONE)"
    fi
done

# ---------------------------------------------------------------------------
# Step 7: Fetch results
# ---------------------------------------------------------------------------
echo "[7] Fetching results..."
scp -i "$KEY_FILE" \
    -o StrictHostKeyChecking=no \
    "ec2-user@$PUBLIC_IP:results.tsv.gz" \
    "$RESULTS_DIR/atb_lexicmap_raw.tsv.gz"

scp -i "$KEY_FILE" \
    -o StrictHostKeyChecking=no \
    "ec2-user@$PUBLIC_IP:/var/log/lexicmap_run.log" \
    "$RESULTS_DIR/lexicmap_run.log"

echo "    Results saved to $RESULTS_DIR/atb_lexicmap_raw.tsv.gz"

# ---------------------------------------------------------------------------
# Step 8: Terminate instance
# ---------------------------------------------------------------------------
echo "[8] Terminating instance $INSTANCE_ID..."
aws ec2 terminate-instances --region "$REGION" --instance-ids "$INSTANCE_ID"
echo "    Instance termination initiated."

echo ""
echo "============================================================"
echo " Done! Results in: $RESULTS_DIR/atb_lexicmap_raw.tsv.gz"
echo " Next: python scripts/analyse_atb_results.py"
echo "============================================================"
