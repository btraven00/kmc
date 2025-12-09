# GitHub Actions Workflows

## Build and Upload Artifacts

This workflow automatically builds the `kmc` binary for multiple platforms and uploads them as artifacts.

### Trigger Events

The workflow runs on:
- **Push** to `main` or `master` branches
- **Pull requests** targeting `main` or `master` branches  
- **Manual trigger** via the Actions tab (workflow_dispatch)

### Build Matrix

Builds are performed on three platforms:

| Platform | Artifact Name | Binary Location |
|----------|---------------|-----------------|
| Ubuntu Linux | `kmc-linux` | `target/release/kmc` |
| macOS | `kmc-macos` | `target/release/kmc` |
| Windows | `kmc-windows` | `target/release/kmc.exe` |

### Workflow Steps

1. **Checkout code** - Clones the repository
2. **Setup Rust** - Installs the stable Rust toolchain
3. **Cache dependencies** - Caches Cargo registry, git dependencies, and build artifacts
4. **Run tests** - Executes `cargo test --verbose`
5. **Build release** - Compiles with `cargo build --release --verbose`
6. **Upload artifacts** - Uploads the release binary for each platform

### Downloading Artifacts

After a successful workflow run:

1. Go to the **Actions** tab in your repository
2. Click on the workflow run
3. Scroll to the **Artifacts** section at the bottom
4. Download the artifact for your platform:
   - `kmc-linux` for Linux
   - `kmc-macos` for macOS
   - `kmc-windows` for Windows

Artifacts are retained for **90 days**.

### Local Testing

To test the build locally before pushing:

```bash
# Run tests
cargo test --verbose

# Build release binary
cargo build --release --verbose

# The binary will be at:
# Linux/macOS: target/release/kmc
# Windows: target/release/kmc.exe
```

### Caching

The workflow uses caching to speed up builds:
- **Cargo registry** - Downloaded crate metadata
- **Cargo git** - Git dependencies
- **Target directory** - Compiled dependencies

This significantly reduces build times on subsequent runs.
