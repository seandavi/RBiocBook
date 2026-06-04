#!/usr/bin/env bash
# Runs once when the devcontainer is created. Dev-only tooling installed here
# (not baked into the public image, so CI/Orchestra stay lean).
set -euo pipefail

# Project-local titlepage extension (needed for the PDF format; _extensions/ is
# gitignored, so it isn't in the image).
quarto install extension --no-prompt nmfs-opensci/quarto_titlepages

# Claude Code CLI for vibe coding inside the container (Node comes from the image).
npm install -g @anthropic-ai/claude-code

# --- antigravity-cli (`agy`) — disabled pending review ---------------------
# Google's Go-based terminal agent. Installs via a remote script piped to bash:
#   curl -fsSL https://antigravity.google/cli/install.sh | bash
# Uncomment to enable (and add /root/.local/bin to PATH via devcontainer.json's
# remoteEnv so the `agy` binary is found):
#   curl -fsSL https://antigravity.google/cli/install.sh | bash

echo "devcontainer ready: $(claude --version 2>/dev/null | head -1)"
