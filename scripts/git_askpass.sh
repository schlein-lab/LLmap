#!/usr/bin/env bash
# Git askpass helper for the SCHLEIN_LAB_TOKEN.
#
# When git asks for a username or password to authenticate to github.com,
# GIT_ASKPASS invokes this script with the prompt text as $1.
# We respond with the username "x-access-token" and the token as password.
# The token itself is read from the SCHLEIN_LAB_TOKEN env var — never stored
# in .git/config, never logged.

case "$1" in
    *[Uu]sername*) echo "x-access-token" ;;
    *[Pp]assword*) echo "${SCHLEIN_LAB_TOKEN:-}" ;;
    *) echo "${SCHLEIN_LAB_TOKEN:-}" ;;
esac
