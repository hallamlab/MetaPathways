# Container validation and nested-build evidence

[Development](development.md) · [User container guide](containers.md)

For maintainers validating runtime changes. These observations describe particular tests, not a guarantee of portability to every host.

## Nested Pathway Tools containers

MP's release container includes an unlicensed Ubuntu dependency base for
`build_pt`. When present, MP uses it automatically instead of installing system
packages inside the nested build. Mamba installations without the bundled base download the digest-pinned public
dependency image automatically. The observations below predate that download route.

A local nested build passed installation, official patch startup, SIF creation,
BLAST validation and Pathway Tools readiness. Complete nested PGDB inference and
portability across HPC installations require separate validation. Use the Mamba
route when a host restricts nested user namespaces or container execution.

### What the nested checks established

- Java/Nextflow stalled on a FUSE filesystem request when executing from the
  outer SIF on the test host. Launching the outer image with `--unsquash`
  resolved this; extraction requires additional temporary space and startup time.
- Outer bind mounts were inherited by the inner build and targeted directories
  absent from its writable root. Clearing inherited bind variables inside the
  outer container let the build proceed. Explicit inner runtime mounts are still
  required; removing all mounts indiscriminately is not a general solution.
- The inner build fell back to a root-mapped namespace without full UID mappings
  or a fakeroot helper. Package installation failed when APT tried to switch to
  its `_apt` user. Preparing dependencies outside this restricted namespace
  removed that step; the licensed installer then completed successfully.
- The prepared-base build completed in approximately six minutes, including
  official patch startup, SIF packaging, BLAST database/search validation and
  the Pathway Tools readiness marker. This duration excludes building the base.
- SRI returned HTTP 503 twice during the experiment. The successful test reused
  107 unmodified official patch files from an existing private SIF, checking
  their hashes against its retained manifest. A fresh vendor download was
  therefore not validated in that successful attempt.

### Build behavior and provenance

The inner build subprocess clears inherited bind variables, while leaving the
parent environment unchanged. The built image metadata records its installer,
recipe, image and official-patch hashes, plus a package-manifest hash for the
bundled dependency base. Validation remains blocking: a failed build is never
registered. Transient vendor HTTP errors receive at most three attempts; missing
patches still stop the build.

Use `--unsquash` on the outer Apptainer invocation if Java stalls on its FUSE
mount. This is a host launch setting; code inside an already-mounted container
cannot remount that outer image. Keep configuration inside the mounted project
when building from a container, and pass the container-visible SIF path with
`--image` when running it. A `/work` registration is meaningful inside that mount,
not from an unrelated host shell.

Maintainers can repeat the licensed checks on a target host from a source checkout:

```bash
python scripts/check_nested_ptools.py \
  --mp-image /absolute/path/metapathways.sif \
  --pt-image /absolute/path/pathway-tools.sif \
  --installer /absolute/path/pathway-tools-29.5-linux-64-tier1-install \
  --output /absolute/path/new-nested-test
```

Activate the host Apptainer environment first so its helper tools are on `PATH`.
Omit `--installer` to check only execution of an existing licensed image. The
output directory must be new; it contains private logs and a machine-readable
result. Each check has a configurable timeout (`--timeout`, default 900 seconds).
The test does not upload licensed files. Its configuration directory is isolated
under the requested output, so it does not replace the normal MP image registry.
Use `--unsquash` when testing the extraction route described above; the default
continues to test the mounted SIF.
A passed startup check verifies Pathway Tools readiness and BLAST integration;
it does not replace a complete PGDB inference test.
