#!/usr/bin/env sh
# =============================================================================
#  SolarSim - demo script: runs every command once and collects the outputs
#  in ./out (images, animations, CSV, SVG). Takes about a minute.
#
#    ./scripts/demo.sh [path/to/solarsim]
# =============================================================================
set -e
cd "$(dirname "$0")/.."
BIN="${1:-./build/solarsim}"
[ -x "$BIN" ] || BIN=./solarsim
if [ ! -x "$BIN" ]; then
  echo "solarsim binary not found - run ./build.sh first" >&2
  exit 1
fi
mkdir -p out
EPOCH=2026-10-21

echo "== selftest";   "$BIN" selftest
echo "== info";       "$BIN" info --epoch $EPOCH --moons --minor > out/info.txt
echo "== ephemeris";  "$BIN" ephem --date "$EPOCH 03:00" --lat 51.4779 --lon 0.0 --csv out/ephem.csv > out/ephem.txt
echo "== events";     "$BIN" events --from $EPOCH --to 2028-01-01 --csv out/events.csv > out/events.txt
echo "== validate";   "$BIN" validate > out/validate.txt || true
echo "== compare";    "$BIN" compare --years 20 --plot out/energy.svg --csv out/energy.csv > out/compare.txt
echo "== precession"; "$BIN" precession > out/precession.txt
echo "== lagrange";   "$BIN" lagrange --pair Sun,Jupiter --test-years 150 --epoch $EPOCH > out/lagrange.txt
echo "== mission";    "$BIN" mission --epoch $EPOCH --porkchop out/porkchop.png --plot out/transfer.svg --csv out/transfer.csv > out/mission.txt
echo "== export";     "$BIN" export --epoch $EPOCH --duration 365 --every 5 --out out/states.csv --csv out/diagnostics.csv
echo "== renders"
"$BIN" render --epoch $EPOCH --out out/overview.png --width 1600 --height 900 --grid --asteroids 4000 --trojans 800 --warmup 400
"$BIN" render --epoch $EPOCH --out out/saturn.png --focus Saturn --moons --width 1600 --height 900
"$BIN" render --epoch $EPOCH --out out/jupiter.png --focus Jupiter --moons --width 1600 --height 900 --pitch 12
"$BIN" render --epoch $EPOCH --out out/earth_moon.png --focus Earth --width 1600 --height 900
"$BIN" render --epoch $EPOCH --out out/orbits.svg --duration 4000
"$BIN" render --epoch $EPOCH --out out/solar_system.gif --frames 180 --days-per-frame 4 --width 720 --height 405 \
       --asteroids 3000 --trojans 600 --warmup 300 --orbit-camera 0.4
"$BIN" render --epoch $EPOCH --out out/inner.gif --view inner --frames 150 --days-per-frame 2 --width 640 --height 400 --minor
echo
echo "Outputs written to ./out"
