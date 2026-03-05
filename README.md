# NECCTON — BAMHBI (FABM submodules)

This repository packages the **BAMHBI biogeochemical model** as a set of **FABM-compatible Fortran submodules**, with documentation derived from the model “hats” (process-block specifications).

## What’s in here

- Fortran sources (FABM submodules): `fortran/`
- Example FABM configurations: `fortran/fabm.yaml`, `fortran/gotm.yaml`
- Process-block documentation: `docs/process/`
- Tracer, diagnostic, parameter registries: `fortran/Databases/`

## Build / compile (FABM)

This repo is intended to be compiled as part of a **FABM** build (e.g., coupled to NEMO-FABM or GOTM-FABM).

See:
- `docs/build/compiling.md` for a practical build checklist
- `fortran/CMakeLists.txt` if you integrate via CMake

## Documentation

- Local preview:
  ```bash
  pip install -r docs/requirements.txt
  mkdocs serve
  ```

## License

BSD-3-Clause — see [LICENSE](LICENSE).

## Citation

See [CITATION.cff](CITATION.cff).
