#!/bin/bash
# LLmap Release Script
#
# Usage:
#   ./scripts/release.sh [version]
#
# Examples:
#   ./scripts/release.sh 1.0.0      # Tag v1.0.0 release
#   ./scripts/release.sh            # Dry-run, shows what would be done

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

VERSION="${1:-}"
DRY_RUN=false

if [[ -z "$VERSION" ]]; then
    DRY_RUN=true
    VERSION="x.y.z"
    echo "=== LLmap Release Script (DRY RUN) ==="
    echo ""
    echo "No version specified. Showing what would be done."
    echo "To create a release, run: $0 <version>"
    echo ""
else
    echo "=== LLmap Release Script ==="
    echo ""
    echo "Version: $VERSION"
    echo ""
fi

# Validate version format
if [[ "$VERSION" != "x.y.z" ]] && ! [[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+(-[a-zA-Z0-9]+)?$ ]]; then
    echo "ERROR: Invalid version format: $VERSION"
    echo "Expected format: MAJOR.MINOR.PATCH or MAJOR.MINOR.PATCH-suffix"
    exit 1
fi

# Check for uncommitted changes
if ! git diff --quiet || ! git diff --staged --quiet; then
    echo "ERROR: Uncommitted changes detected. Please commit or stash before releasing."
    git status --short
    exit 1
fi

# Check we're on main branch
BRANCH=$(git rev-parse --abbrev-ref HEAD)
if [[ "$BRANCH" != "main" ]]; then
    echo "WARNING: Not on main branch (currently on: $BRANCH)"
    if [[ "$DRY_RUN" == "false" ]]; then
        read -p "Continue anyway? [y/N] " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            exit 1
        fi
    fi
fi

# Check build passes
echo "=== Checking build ==="
if [[ -d "build" ]]; then
    cd build
    if cmake --build . -j"$(nproc)" 2>&1 | tail -5; then
        echo "Build: PASS"
    else
        echo "ERROR: Build failed"
        exit 1
    fi
    cd ..
else
    echo "WARNING: No build directory found. Run cmake first."
    if [[ "$DRY_RUN" == "false" ]]; then
        exit 1
    fi
fi

# Check tests pass
echo ""
echo "=== Checking tests ==="
if [[ -d "build" ]]; then
    cd build
    # Run tests and capture full output
    TEST_OUTPUT=$(ctest --output-on-failure -j"$(nproc)" 2>&1)
    # Extract the summary line
    TEST_RESULT=$(echo "$TEST_OUTPUT" | grep "tests passed" || echo "")
    if [[ -n "$TEST_RESULT" ]]; then
        echo "$TEST_RESULT"
        if echo "$TEST_RESULT" | grep -q "100% tests passed"; then
            echo "Tests: PASS"
        else
            echo "ERROR: Some tests failed"
            exit 1
        fi
    else
        echo "WARNING: Could not parse test results"
        echo "$TEST_OUTPUT" | tail -10
    fi
    cd ..
else
    echo "WARNING: No build directory found. Skipping tests."
fi

# Check CHANGELOG has entry for this version
echo ""
echo "=== Checking CHANGELOG ==="
if [[ -f "CHANGELOG.md" ]]; then
    if grep -q "## \[$VERSION\]" CHANGELOG.md; then
        echo "CHANGELOG entry found for $VERSION: PASS"
    else
        echo "WARNING: No CHANGELOG entry found for [$VERSION]"
        if [[ "$DRY_RUN" == "false" ]]; then
            read -p "Continue anyway? [y/N] " -n 1 -r
            echo
            if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                exit 1
            fi
        fi
    fi
else
    echo "WARNING: CHANGELOG.md not found"
fi

# Check CMakeLists.txt version matches
echo ""
echo "=== Checking CMake version ==="
if [[ -f "CMakeLists.txt" ]]; then
    CMAKE_VERSION=$(grep -oP 'VERSION\s+\K[0-9]+\.[0-9]+\.[0-9]+' CMakeLists.txt | head -1 || echo "not found")
    echo "CMakeLists.txt VERSION: $CMAKE_VERSION"
    if [[ "$VERSION" != "x.y.z" ]] && [[ "$CMAKE_VERSION" != "$VERSION" ]]; then
        echo "WARNING: CMakeLists.txt ($CMAKE_VERSION) does not match release version ($VERSION)"
        if [[ "$DRY_RUN" == "false" ]]; then
            read -p "Update CMakeLists.txt? [y/N] " -n 1 -r
            echo
            if [[ $REPLY =~ ^[Yy]$ ]]; then
                sed -i "s/VERSION $CMAKE_VERSION/VERSION $VERSION/" CMakeLists.txt
                git add CMakeLists.txt
                git commit -m "release: bump version to $VERSION"
                echo "CMakeLists.txt updated and committed"
            fi
        fi
    fi
else
    echo "WARNING: CMakeLists.txt not found"
fi

echo ""
echo "=== Release Summary ==="
echo "Version:    $VERSION"
echo "Tag:        v$VERSION"
echo "Branch:     $BRANCH"
echo "Commit:     $(git rev-parse --short HEAD)"
echo ""

if [[ "$DRY_RUN" == "true" ]]; then
    echo "DRY RUN complete. To create the release:"
    echo "  $0 $VERSION"
    exit 0
fi

# Create tag
echo "=== Creating tag ==="
if git tag -l "v$VERSION" | grep -q "v$VERSION"; then
    echo "ERROR: Tag v$VERSION already exists"
    exit 1
fi

# Generate tag message from CHANGELOG
TAG_MESSAGE="LLmap $VERSION

$(sed -n "/## \[$VERSION\]/,/## \[/p" CHANGELOG.md | head -n -1 | tail -n +2)"

git tag -a "v$VERSION" -m "$TAG_MESSAGE"
echo "Tag created: v$VERSION"

# Push
echo ""
echo "=== Pushing to remote ==="
read -p "Push tag v$VERSION to origin? [y/N] " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    git push origin "v$VERSION"
    echo "Tag pushed successfully"
    echo ""
    echo "Release URL: https://github.com/schlein-lab/LLmap/releases/tag/v$VERSION"
else
    echo "Tag not pushed. To push manually:"
    echo "  git push origin v$VERSION"
fi

echo ""
echo "=== Release complete ==="
