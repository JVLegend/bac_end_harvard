#!/usr/bin/env bash
# Pack outputs into a single tarball + write MANIFEST.json describing them.
set -euo pipefail
cd "$(dirname "$0")/.."

TS=$(date -u +%Y%m%dT%H%M%SZ)
ARCHIVE="outputs/smartsepsis_5090_${TS}.tar.gz"
MANIFEST="outputs/MANIFEST.json"

# Build manifest
python3 - <<PY
import os, json, hashlib, glob
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda: f.read(1<<16), b''): h.update(b)
    return h.hexdigest()[:16]
files=[]
for path in sorted(glob.glob("outputs/**/*", recursive=True)):
    if os.path.isfile(path):
        files.append({"path":path,"bytes":os.path.getsize(path),"sha256_16":sha(path)})
json.dump({"generated_at_utc":"${TS}","files":files,"job_count":len(files)},
          open("$MANIFEST","w"), indent=2)
print(f"MANIFEST has {len(files)} files")
PY

# Tarball (excluding any FAILED.txt? keep them — honest)
tar -czf "$ARCHIVE" -C outputs .
echo ""
echo "Packed: $ARCHIVE"
du -h "$ARCHIVE" | awk '{print "Size:", $1}'
echo ""
echo "Send this tarball back to the user."
echo "Suggested transports: scp, rsync, GitHub release attachment, Google Drive."
