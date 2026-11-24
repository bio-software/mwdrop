#!/bin/bash
set -e

# MWDROP needs the patched MobyWat fork (exits non-zero on errors, accepts
# MOBYWAT_RMS_ATOM, supports DNA/RNA targets), not the original 1.1 release
# from mobywat.com. MOBYWAT_REF selects the branch, tag or commit to build
# (default: main), e.g. MOBYWAT_REF=7495847 ./install.sh
MOBYWAT_REPO="https://github.com/bio-software/MobyWat.git"
MOBYWAT_REF="${MOBYWAT_REF:-main}"
TMPDIR="$(mktemp -d /tmp/mobywat.XXXXXX)"
SRC_DIR="$TMPDIR/mobywat"

cleanup() {
    rm -rf "$TMPDIR"
}

trap cleanup EXIT INT TERM

if ! command -v git >/dev/null 2>&1; then
    echo "Error: git is not installed."
    exit 1
fi

echo "Cloning MobyWat from $MOBYWAT_REPO ..."
git clone --quiet "$MOBYWAT_REPO" "$SRC_DIR"

echo "Checking out $MOBYWAT_REF ..."
git -C "$SRC_DIR" checkout --quiet "$MOBYWAT_REF"

# Compiling MobyWat
echo "Compiling MobyWat"
make -C "$SRC_DIR/src"

# Installing MobyWat (make install copies into ~/bin, which must already be a directory)
echo "Installing MobyWat"
mkdir -p "$HOME/bin"
make -C "$SRC_DIR/src" install

echo "Installation completed."
