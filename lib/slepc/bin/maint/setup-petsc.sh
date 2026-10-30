#!/bin/bash -e

if [ -z "${CI_PIPELINE_ID+x}" ]; then
  echo "Running locally outside GitLab CI..."
  GIT_CLEAN_FLAGS="${GIT_CLEAN_FLAGS:--f}"
fi

if [ -z "$SLEPC_DIR" ]; then
  echo "Error: SLEPC_DIR environment variable is not defined."
  exit 1
fi

if [ ! -d petsc ]; then
  for i in $(seq 5); do
    if git clone https://gitlab.com/petsc/petsc.git; then
      break
    fi
    echo "Failed to clone PETSc. Retrying in 10 seconds (attempt $i of 5)..."
    sleep 10
  done
fi

cd petsc
git clean ${GIT_CLEAN_FLAGS}

for i in $(seq 5); do
  if git fetch; then
    break
  fi
  sleep 10
done

# Determine the correct parent branch based on slepcversion.h
if [ $(awk '/define SLEPC_VERSION_RELEASE/{print $3}' "$SLEPC_DIR/include/slepcversion.h") == "1" ]; then
  PARENT_BRANCH="release"
else
  PARENT_BRANCH="main"
fi

echo "PARENT_BRANCH=${PARENT_BRANCH}"


if [ -z "${PETSC_BRANCH}" ]; then
  echo "PETSC_BRANCH is not defined. Using origin/${PARENT_BRANCH}..."
  git checkout --no-track -B runner/ci-branch origin/${PARENT_BRANCH}
else
  git checkout --no-track -B runner/ci-branch origin/${PETSC_BRANCH}

  if [ $(awk '/define PETSC_VERSION_RELEASE/{print $3}' "include/petscversion.h") == "1" ]; then
    PETSC_PARENT_BRANCH="release"
  else
    PETSC_PARENT_BRANCH="main"
  fi

  if [ "$PARENT_BRANCH" != "$PETSC_PARENT_BRANCH" ]; then
    echo "Error: Branch compatibility mismatch!"
    echo "SLEPc expects a '${PARENT_BRANCH}' branch, but PETSc '${PETSC_BRANCH}' is a '${PETSC_PARENT_BRANCH}' branch."
    exit 1
  fi

  MERGED_LOWER=$(echo "${PETSC_MERGED}" | tr '[:upper:]' '[:lower:]')

  if [ "$MERGED_LOWER" == "true" ] || [ "$MERGED_LOWER" == "1" ]; then
    echo "Testing PETSC_BRANCH=${PETSC_BRANCH} merged onto origin/${PARENT_BRANCH}..."
    git log -1 --format="HEAD is at %h (%as by %an): %s"
    git -c user.name="GitLab CI" -c user.email="gitlab-ci@slepc.org" merge --no-edit origin/${PARENT_BRANCH}
  else
    echo "Using PETSC_BRANCH=${PETSC_BRANCH} directly..."
  fi
fi

git log -1 --format="HEAD is at %h (%as by %an): %s"
