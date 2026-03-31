#!/usr/bin/env bash
# =============================================================================
# find_subset_reps_ec2.sh
#
# Launch a spot EC2 instance to search for better ATB representatives for
# high-priority subset loci using the AllTheBacteria LexicMap index.
#
# Targets:
#   KL713 — superset with 100-N cross-contig artefact; find clean same-contig rep
#   KL601 — subset, 3'-truncated assembly; find longer rep
#   KL742 — subset, 3'-truncated assembly; find longer rep
#
# Strategy:
#   For KL601 / KL742: search with current (truncated) representative; select
#   hits where the subject region EXTENDS beyond the query end — these are
#   assemblies that contain a more complete version of the locus.
#
#   For KL713: search and flag hits where the full match is on a SINGLE CONTIG
#   (no 100-N spacer in the matched region) — these are clean representatives.
#
# Prerequisites (local):
#   aws configure  (valid credentials, region eu-west-2)
#   Query FASTAs produced by scripts/extract_subset_queries.py:
#     DB/subset_queries/KL601_query.fasta
#     DB/subset_queries/KL713_query.fasta
#     DB/subset_queries/KL742_query.fasta
#
# Cost: ~$0.30–0.80 total (c7g.8xlarge spot, ~20–30 min)
#
# Usage:
#   python scripts/extract_subset_queries.py   # generates query FASTAs
#   bash scripts/find_subset_reps_ec2.sh
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
REGION="eu-west-2"
INSTANCE_TYPE="c7g.8xlarge"    # ARM64, 32 vCPU, 64 GiB
KEY_NAME="atb-search-key"
KEY_FILE="$HOME/.ssh/atb_search_key"
SECURITY_GROUP_NAME="atb-lexicmap-sg"
RESULTS_DIR="DB/atb_lexicmap_results/subset_fix"
AMI_ID="ami-0a3c640e8f9c74b61"   # Amazon Linux 2023 ARM64, eu-west-2

QUERY_DIR="DB/subset_queries"

# Verify query FASTAs exist
for kl in KL601 KL713 KL742; do
    qf="${QUERY_DIR}/${kl}_query.fasta"
    if [ ! -f "$qf" ]; then
        echo "ERROR: Query FASTA not found: $qf"
        echo "Run: python scripts/extract_subset_queries.py"
        exit 1
    fi
done

mkdir -p "$RESULTS_DIR"

echo "============================================================"
echo " ATB LexicMap Search — subset locus representative fix"
echo " Targets: KL601, KL713, KL742"
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

aws ec2 import-key-pair \
    --region "$REGION" \
    --key-name "$KEY_NAME" \
    --public-key-material "$(base64 < "${KEY_FILE}.pub")" \
    2>/dev/null || true
echo "    Key pair '$KEY_NAME' ready in $REGION"

# ---------------------------------------------------------------------------
# Step 2: Security group
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
# Step 3: Verify AMI
# ---------------------------------------------------------------------------
echo "[3] Verifying AMI $AMI_ID..."
AMI_CHECK=$(aws ec2 describe-images \
    --region "$REGION" \
    --image-ids "$AMI_ID" \
    --query "Images[0].Name" --output text 2>/dev/null || echo "NOT_FOUND")
if [ "$AMI_CHECK" = "NOT_FOUND" ] || [ -z "$AMI_CHECK" ]; then
    echo "    AMI not found — looking up latest Amazon Linux 2023 ARM64..."
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
# Step 4: Combine query FASTAs locally (will be SCP'd after instance is up)
# ---------------------------------------------------------------------------
COMBINED_QUERY=$(mktemp /tmp/subset_query_XXXXXX.fasta)
cat "$QUERY_DIR/KL601_query.fasta" \
    "$QUERY_DIR/KL713_query.fasta" \
    "$QUERY_DIR/KL742_query.fasta" > "$COMBINED_QUERY"
echo "[4] Combined query FASTA: $COMBINED_QUERY ($(wc -c < "$COMBINED_QUERY") bytes)"

# ---------------------------------------------------------------------------
# Step 5: User-data script — install deps + mount index only (no query embed)
#         Query FASTA is SCP'd after instance is running to avoid 25KB limit.
# ---------------------------------------------------------------------------
USER_DATA=$(cat <<'USERDATA_EOF'
#!/bin/bash
set -euxo pipefail
exec > /var/log/lexicmap_run.log 2>&1

echo "=== Setting up ATB LexicMap environment — subset locus fix ==="
cd /home/ec2-user

# Install dependencies
dnf install -y wget tar fuse fuse-libs 2>/dev/null || true

# Install mountpoint-s3
wget -q "https://s3.amazonaws.com/mountpoint-s3-release/latest/arm64/mount-s3.rpm"
dnf install -y ./mount-s3.rpm || yum install -y ./mount-s3.rpm

# Install LexicMap v0.8.1 (ARM64)
wget -q https://github.com/shenwei356/LexicMap/releases/download/v0.8.1/lexicmap_linux_arm64.tar.gz
tar xzf lexicmap_linux_arm64.tar.gz
mv lexicmap /usr/local/bin/lexicmap
lexicmap version

# Mount ATB index
mkdir -p /mnt/atb_index
UNSTABLE_MOUNTPOINT_MAX_PREFETCH_WINDOW_SIZE=65536 \
    mount-s3 --read-only --prefix 202408/ allthebacteria-lexicmap \
    /mnt/atb_index --no-sign-request
echo "=== Index mounted — signalling READY ==="
touch /home/ec2-user/READY
USERDATA_EOF
)

# ---------------------------------------------------------------------------
# Step 6: Launch spot instance
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
      "Ebs": { "VolumeSize": 30, "VolumeType": "gp3", "DeleteOnTermination": true }
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
# Step 7: Wait, fetch, terminate
# ---------------------------------------------------------------------------
echo "[5] Waiting for instance to reach running state..."
aws ec2 wait instance-running --region "$REGION" --instance-ids "$INSTANCE_ID"

PUBLIC_IP=$(aws ec2 describe-instances \
    --region "$REGION" \
    --instance-ids "$INSTANCE_ID" \
    --query "Reservations[0].Instances[0].PublicIpAddress" \
    --output text)
echo "    Public IP: $PUBLIC_IP"

cat > "$RESULTS_DIR/instance_info.txt" <<INFO_EOF
INSTANCE_ID=$INSTANCE_ID
PUBLIC_IP=$PUBLIC_IP
REGION=$REGION
KEY_FILE=$KEY_FILE
LAUNCHED=$(date -u)
INFO_EOF

echo "[6] Waiting for instance setup to complete (polls every 30s)..."
MAX_SETUP=1800
ELAPSED=0
while [ $ELAPSED -lt $MAX_SETUP ]; do
    sleep 30
    ELAPSED=$((ELAPSED + 30))
    READY=$(ssh -i "$KEY_FILE" \
        -o StrictHostKeyChecking=no \
        -o ConnectTimeout=10 \
        -o BatchMode=yes \
        "ec2-user@$PUBLIC_IP" \
        "test -f /home/ec2-user/READY && echo yes || echo no" 2>/dev/null || echo "connecting")
    if [ "$READY" = "yes" ]; then
        echo "    [${ELAPSED}s] Instance ready!"
        break
    else
        echo "    [${ELAPSED}s] Still setting up..."
    fi
done

echo "[7] Uploading query FASTA..."
scp -i "$KEY_FILE" -o StrictHostKeyChecking=no \
    "$COMBINED_QUERY" "ec2-user@$PUBLIC_IP:/home/ec2-user/query.fasta"
rm -f "$COMBINED_QUERY"

echo "    Triggering LexicMap search..."
ssh -i "$KEY_FILE" -o StrictHostKeyChecking=no -o BatchMode=yes \
    "ec2-user@$PUBLIC_IP" bash <<'REMOTE_EOF'
grep "^>" /home/ec2-user/query.fasta
nohup bash -c '
    set -euxo pipefail
    exec >> /var/log/lexicmap_run.log 2>&1
    echo "=== Running LexicMap search ==="
    lexicmap search \
        -d /mnt/atb_index \
        /home/ec2-user/query.fasta \
        -o /home/ec2-user/results.tsv.gz \
        --align-min-match-pident 90 \
        --min-qcov-per-genome 70 \
        --top-n-genomes 0 \
        -j 32
    echo "=== Search complete ==="
    zcat /home/ec2-user/results.tsv.gz | wc -l
    touch /home/ec2-user/DONE
' &>/dev/null &
echo "    LexicMap launched (PID $!)"
REMOTE_EOF

echo "[8] Waiting for search to complete (polls every 60s)..."
MAX_WAIT=7200
ELAPSED=0
while [ $ELAPSED -lt $MAX_WAIT ]; do
    sleep 60
    ELAPSED=$((ELAPSED + 60))
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
        echo "    [${ELAPSED}s] Still running..."
    fi
done

echo "[9] Fetching results..."
scp -i "$KEY_FILE" -o StrictHostKeyChecking=no \
    "ec2-user@$PUBLIC_IP:results.tsv.gz" \
    "$RESULTS_DIR/atb_subset_raw.tsv.gz"
scp -i "$KEY_FILE" -o StrictHostKeyChecking=no \
    "ec2-user@$PUBLIC_IP:/var/log/lexicmap_run.log" \
    "$RESULTS_DIR/lexicmap_run.log"

echo "    Results: $RESULTS_DIR/atb_subset_raw.tsv.gz"

echo "[10] Terminating instance $INSTANCE_ID..."
aws ec2 terminate-instances --region "$REGION" --instance-ids "$INSTANCE_ID"

echo ""
echo "============================================================"
echo " Done! Next: python scripts/analyse_subset_hits.py"
echo "============================================================"
