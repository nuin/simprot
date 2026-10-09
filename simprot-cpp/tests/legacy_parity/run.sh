#!/usr/bin/env bash
# Run the C++20 port and the original SIMPROT 1.04 (legacy/) with the same
# seeds over a matrix of trees and options, and check the outputs agree.
#
#   tests/legacy_parity/run.sh path/to/simprot [seeds-per-config]   (default 20)
#
# Legacy has no seed option and seeds from getpid(), so a scratch copy of
# legacy/simprot.cpp is built with SetSeed() reading $SIMPROT_SEED, and with
# the port's model fixes (rebuilt PAM and JTT, source-frequency divisor; see
# below). Needs a
# C++ compiler, popt (brew install popt / apt install libpopt-dev) and python3.
# Flags differ between the versions, so each config lists both spellings.
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
repo="$(cd "$here/../../.." && pwd)"
cpp="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
seeds="${2:-20}"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# Build legacy with the seed override
cp "$repo"/legacy/{simprot.cpp,random.c,random.h,eigen.h} "$tmp/"
needle='SetSeed(getpid());'
[ "$(grep -c "$needle" "$tmp/simprot.cpp")" = 1 ] || { echo "expected one '$needle' in legacy/simprot.cpp" >&2; exit 1; }
python3 - "$tmp/simprot.cpp" "$needle" <<'EOF'
import sys
path, needle = sys.argv[1], sys.argv[2]
s = open(path, errors='surrogateescape').read()
s = s.replace(needle, 'SetSeed(getenv("SIMPROT_SEED") ? atoi(getenv("SIMPROT_SEED")) : getpid());')
open(path, 'w', errors='surrogateescape').write(s)
EOF
# The C++ port ships rebuilt PAM and JTT data (tools/make_eigen.py) and divides
# substitution probabilities by the source residue's frequency; legacy/ keeps
# the 1.04 data and divisor. Give the scratch legacy copy the same model data
# and divisor so every run still compares the two code paths.
python3 - "$tmp/eigen.h" "$repo/simprot-cpp/include/simprot/evolution/matrix_data.hpp" <<'EOF'
import re, sys
eigen_h, data_hpp = sys.argv[1], sys.argv[2]
d = open(data_hpp).read()
def numbers(name):
    i = d.index(name + ' =')
    body = d[d.index('{', i):d.index(';', i)]
    return re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', body)
s = open(eigen_h, errors='surrogateescape').read()
def replace(s, name, body):
    i = s.index(name)
    return s[:s.index('{', i)] + body + s[s.index(';', i):]
for model in ('pam', 'jtt'):
    lam = numbers(model + '_eigenvalues')
    vec = numbers(model + '_eigenvectors')
    assert len(lam) == 20 and len(vec) == 400, (model, len(lam), len(vec))
    s = replace(s, model + 'eigmat[20]', '{' + ', '.join(lam) + '}')
    s = replace(s, model + 'probmat[20][20]',
                '{\n' + ',\n'.join('{' + ', '.join(vec[k*20:(k+1)*20]) + '}' for k in range(20)) + '\n}')
open(eigen_h, 'w', errors='surrogateescape').write(s)
EOF
grep -q 'sum += p/freqaa\[j\];' "$tmp/simprot.cpp" || { echo "expected 'sum += p/freqaa[j];' in legacy" >&2; exit 1; }
sed -i.bak 's#sum += p/freqaa\[j\];#sum += p/freqaa[i];#' "$tmp/simprot.cpp"
popt="$(brew --prefix popt 2>/dev/null || echo /usr)"
c++ -O2 -w -o "$tmp/legacy" "$tmp/simprot.cpp" "$tmp/random.c" -I"$popt/include" -L"$popt/lib" -lpopt -lm

# Trees: two small ones (one with nested zero-length branches) and two shipped ones
echo '(((A:0.1,B:0.2):0.15,(C:0.3,D:0.05):0.2):0.1,((E:0.2,F:0.1):0.25,(G:0.15,H:0.3):0.1):0.2);' > "$tmp/t8.nwk"
echo '(((A:0.2,B:0.2):0,(C:0.2,D:0.2):0):0.3,E:0.2);' > "$tmp/nz.nwk"
cp "$repo/data/bigree40.txt" "$tmp/b40.nwk"
cp "$repo/data/bigree80.txt" "$tmp/b80.nwk"
printf '0.5\n0.2\n0.1\n0.1\n0.05\n0.05\n' > "$tmp/custom.txt"

# tree | legacy flags | C++ flags
configs=(
  "t8.nwk|-r 100|-r 100"
  "nz.nwk|-r 100|-r 100"
  "b40.nwk|-r 100|-r 100"
  "b80.nwk|-r 200|-r 200"
  "b40.nwk|-r 300 -g 0.05|-r 300 -g 0.05"
  "t8.nwk|-r 300 -g 0.1|-r 300 -g 0.1"
  "t8.nwk|-r 20 -g 0.2|-r 20 -g 0.2"
  "t8.nwk|-r 100 -x 0.5|-r 100 -x 0.5"
  "t8.nwk|-r 100 -x 2.5|-r 100 -x 2.5"
  "t8.nwk|-r 100 -x -1|-r 100 -x -1"
  "t8.nwk|-r 100 -g 0|-r 100 -g 0"
  "t8.nwk|-r 150 -p 0|-r 150 -p 0"
  "t8.nwk|-r 150 -p 1|-r 150 -p 1"
  "b40.nwk|-r 500 -p 1 -x 0.3 -g 0.06|-r 500 -p 1 -x 0.3 -g 0.06"
  "t8.nwk|-r 200 -y 1 -g 0.08|-r 200 -b 1 -g 0.08"
  "t8.nwk|-r 200 -y 1 -k -3 -g 0.08|-r 200 -b 1 -k -3 -g 0.08"
  "t8.nwk|-r 200 -u $tmp/custom.txt -g 0.08|-r 200 -u $tmp/custom.txt -g 0.08"
  "t8.nwk|-r 150 -b 2.5|-r 150 -t 2.5"
  "b40.nwk|-r 150 -e 2|-r 150 -v 2"
  "b40.nwk|-r 150 -m 0.2|-r 150 -e 0.2"
  "b40.nwk|-r 150 -e 2 -m 0.15|-r 150 -v 2 -e 0.15"
)

fail=0
for spec in "${configs[@]}"; do
  IFS='|' read -r tree la ca <<<"$spec"
  ok=0; bad=""; crash=""
  for s in $(seq 1 "$seeds"); do
    rm -rf "$tmp/L" "$tmp/C"; mkdir "$tmp/L" "$tmp/C"
    # shellcheck disable=SC2086
    # A legacy run that dies (it occasionally gets "Killed: 9") can leave a
    # partial seq.fa, so judge it by its exit status
    if ! (cd "$tmp/L" && SIMPROT_SEED=$s "$tmp/legacy" -f "$tmp/$tree" -a aln.fa -s seq.fa -o indel.log $la >/dev/null 2>&1) 2>/dev/null \
        || [ ! -s "$tmp/L/seq.fa" ]; then
      crash="$crash $s"; continue
    fi
    # shellcheck disable=SC2086
    (cd "$tmp/C" && "$cpp" -f "$tmp/$tree" -a aln.fa -s seq.fa -o indel.log -S "$s" $ca >/dev/null 2>&1)
    if python3 -I "$here/compare_legacy.py" "$tmp/L" "$tmp/C" >/dev/null; then
      ok=$((ok + 1))
    else
      bad="$bad $s"
    fi
  done
  ran=$((seeds - $(echo $crash | wc -w)))
  echo "$tree [$ca]: $ok/$ran match${bad:+; MISMATCH seeds:$bad}${crash:+; legacy produced no output for seeds:$crash}"
  [ -z "$bad" ] || fail=1
done

[ "$fail" = 0 ] && echo "LEGACY PARITY OK" || echo "LEGACY PARITY FAILED"
exit "$fail"
