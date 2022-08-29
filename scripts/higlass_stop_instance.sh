#!/usr/bin/env bash

set -e
set -u

git_root="$(git rev-parse --show-toplevel)"
name="$(basename "$git_root")"
sudo -E higlass-manage stop "$name"
