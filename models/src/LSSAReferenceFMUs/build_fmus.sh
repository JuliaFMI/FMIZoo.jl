#!/usr/bin/env bash
# Builds the Linux binaries (x86_64-linux) of the LS-SA NBouncingBalls FMUs for all ball counts and adds them to
#   models/bin/LSSAReferenceFMUs/<version>/3.0/NBouncingBalls_<N>.fmu
# Binaries of other platforms already in the FMU (x86_64-windows from build_fmus.ps1) are kept. If an existing FMU
# has a different modelDescription.xml, its binaries were built from different sources: the script stops.
#
# usage (from any directory):
#   ./build_fmus.sh [-v VERSION] [-n "1 3 10 30 100 300 1000"] [-b BUILD_DIR] [-o OUT_DIR] [-- CMAKE_ARGS...]
#
# Needs cmake (>= 3.15) and gcc. All FMUs of one version must be built with the same compiler / flags per platform
# (they are used for benchmarks). Variants go to a separate build and output directory, e.g.
#   ./build_fmus.sh -b build-linux-nocache -o /tmp/nocache -- -DLSSA_PLAN_CACHE=OFF
set -euo pipefail

src="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
version="0.2.0"
counts="1 3 10 30 100 300 1000"
build_dir="build-linux"
out=""
while getopts "v:n:b:o:" opt; do
    case "$opt" in
        v) version="$OPTARG" ;;
        n) counts="$OPTARG" ;;
        b) build_dir="$OPTARG" ;;
        o) out="$OPTARG" ;;
        *) echo "usage: $0 [-v VERSION] [-n COUNTS] [-b BUILD_DIR] [-o OUT_DIR] [-- CMAKE_ARGS...]" >&2; exit 2 ;;
    esac
done
shift $((OPTIND - 1))
cmake_args=("$@")
out="${out:-$src/../../bin/LSSAReferenceFMUs/$version/3.0}"
mkdir -p "$out"
out="$(cd "$out" && pwd)"

platform="x86_64-linux"
[[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]] || { echo "builds $platform binaries: run on x86_64 Linux" >&2; exit 1; }

# adds the content of the staging directory $1 (modelDescription.xml, extra/, binaries/<platform>/) to the FMU $2
merge_into() {
    local staging="$1" target="$2" tmp
    tmp="$(mktemp -d)"
    if [[ -f "$target" ]]; then
        (cd "$tmp" && cmake -E tar xf "$target")
        if ! cmp -s <(tr -d '\r' < "$tmp/modelDescription.xml") <(tr -d '\r' < "$staging/modelDescription.xml"); then
            echo "error: $target has a different modelDescription.xml (built from different sources);" >&2
            echo "       rebuild all platforms or remove the FMU first" >&2
            rm -rf "$tmp"
            exit 1
        fi
        rm -rf "$tmp/binaries/$platform"
    fi
    cp -r "$staging"/. "$tmp"/
    (cd "$tmp" && cmake -E tar cf "$target.tmp" --format=zip -- $(ls -A))
    mv -f "$target.tmp" "$target"
    rm -rf "$tmp"
}

for n in $counts; do
    build="$src/$build_dir/n$n"
    echo "=== N = $n ==="
    cmake -S "$src" -B "$build" -DCMAKE_BUILD_TYPE=Release "-DNUM_BALLS=$n" "${cmake_args[@]}"
    cmake --build "$build"
    target="$out/NBouncingBalls_$n.fmu"
    merge_into "$build/fmu_staging" "$target"
    echo "-> $target ($(cmake -E tar tf "$target" | grep -c '^binaries/.*/NBouncingBalls\.' ) platform binaries)"
done
