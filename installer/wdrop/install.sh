#!/bin/bash
set -e

# wdrop lives in its own repository. WDROP_REF selects the branch, tag or
# commit to build (default: master), e.g. WDROP_REF=0.6 ./install.sh
WDROP_REPO="https://github.com/bio-software/wdrop.git"
WDROP_REF="${WDROP_REF:-master}"
TMPDIR="$(mktemp -d /tmp/wdrop.XXXXXX)"
SRC_DIR="$TMPDIR/wdrop"

cleanup() {
    rm -rf "$TMPDIR"
}

trap cleanup EXIT INT TERM

if ! command -v git >/dev/null 2>&1; then
    echo "Error: git is not installed."
    exit 1
fi

echo "Cloning wdrop from $WDROP_REPO ..."
git clone --quiet "$WDROP_REPO" "$SRC_DIR"

echo "Checking out $WDROP_REF ..."
git -C "$SRC_DIR" checkout --quiet "$WDROP_REF"

# Compiling wdrop
echo "Compiling wdrop"
make -C "$SRC_DIR/src"

# Installing wdrop (make install copies into ~/bin, which must already be a directory)
echo "Installing wdrop"
mkdir -p "$HOME/bin"
make -C "$SRC_DIR/src" install

echo "Installation completed."
