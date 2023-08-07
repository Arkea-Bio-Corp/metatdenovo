#!/bin/bash --login

# Temporarily disable strict mode and activate conda:
set +euo pipefail
conda activate kreation

# Re-enable strict mode:
set -euo pipefail

# exec the final command:
exec python2.7 /kreation/KREATION/KREATION.py