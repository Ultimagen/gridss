# GRIDSS DevContainer Configuration

This directory contains the VSCode DevContainer configuration for GRIDSS development.

## Overview

The devcontainer uses the **production GRIDSS Dockerfile directly** (no duplication!) and adds development tools via a post-create script. This provides a complete development environment with:

- **Java Development**: JDK 11 + Maven for building GRIDSS
- **R Development**: R 4.0+ with all required CRAN and Bioconductor packages
- **Python Development**: Python 3 with pysam, biopython, and other dependencies
- **Complete Toolchain**: All bioinformatics tools (samtools, bwa, bedtools, bcftools, etc.)
- **VSCode Extensions**: Java, R, Python language support with IntelliSense

## How It Works

1. **Builds the production `gridss` stage** from `../Dockerfile`
2. **Runs `post-create.sh`** after container creation to install dev tools (Maven, vim, nano, git)
3. **No Dockerfile duplication** - uses the production Dockerfile directly

## Getting Started

### Prerequisites

1. **VSCode** with the [Dev Containers extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)
2. **Docker** installed and running
3. **Host directories** must exist:
   - `/data/Runs/VariantCalling/work/260416` (test data)
   - `/data/Runs/genomes` (reference genomes)

### Opening in Container

1. Open this repository in VSCode
2. Press `F1` and select "Dev Containers: Reopen in Container"
3. Wait for the container to build (30-40 minutes first time, cached after)

## Configuration Details

### Image Size
- **Full image**: ~8-10GB (includes all production tools)
- **Build time**: 30-40 minutes first build, <1 minute on subsequent builds (Docker cache)

### Resource Allocation
- **CPU**: 4 cores (helps with R package compilation and Maven builds)
- **Memory**: 8GB (Maven and R can be memory-intensive)

Adjust these in `devcontainer.json` `runArgs` if needed.

### Mounted Directories

The devcontainer mounts:
- **Workspace**: `/workspace` (your GRIDSS repository)
- **Test data**: `/data/Runs/VariantCalling/work/260416`
- **Genomes**: `/data/Runs/genomes`

### Environment Variables

The container inherits all environment variables from the production Dockerfile:
- `GRIDSS_JAR`: Points to the GRIDSS JAR in `/opt/gridss/`
- `PATH`: Includes all bioinformatics tools
- `R_INSTALL_STAGED=false`: Required for R package installation
- `LC_ALL=C`: Locale setting

## Development Workflow

### 1. Building GRIDSS

```bash
# Clean build
mvn clean package -DskipTests

# With tests
mvn clean package

# Specific test
mvn test -Dtest=CallVariantsTest
```

### 2. Testing R Scripts

```bash
# Load R libraries
Rscript -e 'library(StructuralVariantAnnotation)'

# Run R script
Rscript scripts/libgridss.R
```

### 3. Testing Python Scripts

```bash
# Test imports
python3 -c 'import pysam'

# Run script
python3 scripts/align_long_homopolymers.py --help
```

### 4. Running GRIDSS

Since this is the full production image, you can run complete GRIDSS workflows:

```bash
/opt/gridss/gridss --help
```

## Customization

### Switching to Lightweight Image

If the full image is too large, you can target a different stage by modifying `devcontainer.json`:

```json
"build": {
    "target": "gridss_c_build_environment"
}
```

Then modify `post-create.sh` to install R and other missing dependencies. See the plan file for details.

### Adding VSCode Extensions

Edit `devcontainer.json` and add extension IDs to the `extensions` array:

```json
"extensions": [
    "vscjava.vscode-java-pack",
    "your.extension.id"
]
```

### Adjusting Resource Limits

Edit `runArgs` in `devcontainer.json`:

```json
"runArgs": [
    "--cpus=8",
    "--memory=16g"
]
```

## Synchronization with Production

**Important**: This devcontainer uses the production `Dockerfile` **directly** - no duplication! Changes to production dependencies automatically flow to the dev environment.

Development tools (Maven, vim, nano, git) are added via `post-create.sh` after the container is created, keeping the production Dockerfile unchanged.

## Troubleshooting

### Container Build Fails

- **Check Docker resources**: Ensure Docker has enough CPU/memory allocated
- **Network issues**: R package installation requires internet access
- **Build args**: Ensure `GRIDSS_VERSION` matches a valid version

### Git Submodule Issues

The `postCreateCommand` should initialize the htslib submodule automatically. If not:

```bash
git submodule update --init --recursive
```

### Permission Issues

The container runs as `root`. VSCode handles file permissions automatically. If you encounter issues:

```bash
git config --global --add safe.directory /workspace
```

### Maven Build Issues

First build downloads ~500MB of dependencies. Ensure:
- Internet connection is stable
- Enough disk space available

### R Package Loading Fails

The production image includes all required R packages. If a package is missing:

```R
# CRAN package
install.packages("package-name")

# Bioconductor package
BiocManager::install("package-name")
```

## Performance Tips

1. **Use Docker BuildKit**: Enable for faster builds
   ```bash
   export DOCKER_BUILDKIT=1
   ```

2. **Pre-cache Maven dependencies** (optional):
   Add to `postCreateCommand` in `devcontainer.json`:
   ```bash
   "postCreateCommand": "git config --global --add safe.directory /workspace && git submodule update --init --recursive && mvn dependency:go-offline"
   ```

3. **Exclude large directories**: The `.dockerignore` file already excludes `.git`, `target`, etc.

## Production Dockerfile Optimizations

See the plan file for potential optimizations to the production Dockerfile:
- BuildKit cache mounts for apt packages
- BuildKit cache mounts for R package compilation
- Multi-architecture support

These should be implemented in separate commits to avoid conflating changes.
