#!/bin/bash
set -e

for root_setup in /opt/root/bin/thisroot.sh /usr/local/bin/thisroot.sh /opt/conda/bin/thisroot.sh; do
    if [ -f "$root_setup" ]; then
        # ROOT needs its runtime environment sourced before Python can import ROOT.
        source "$root_setup"
        break
    fi
done

exec "$@"
