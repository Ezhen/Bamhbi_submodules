# Scope & conventions

## What this repo is (and is not)

This repository is a **FABM module collection**:
- It does **not** provide a standalone executable.
- It is designed to be compiled **inside a FABM host model** (e.g., NEMO-FABM, GOTM-FABM).

## Conventions (recommended)

These are conventions used throughout the “hats”. If your host model uses different unit conventions, document the translation here.

- **State units:** tracer units follow `fortran/Databases/Tracers.txt` (typically mmol m-3).
- **Tendencies:** FABM exchanges tendencies in the host model’s convention (usually per second).
- **Flux sign conventions:** document host-model conventions for air–sea and bottom fluxes when activating `do_surface` / `do_bottom`.

## Files of record

- Tracers registry: `fortran/Databases/Tracers.txt`
- Diagnostics registry: `fortran/Databases/Diagnostics.txt`
- Parameters registry: `fortran/Databases/Parameters.txt`
- Module registry: `fortran/Databases/Modules.txt`
