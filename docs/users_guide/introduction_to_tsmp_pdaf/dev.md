# Developer information

## List of TSMP-related changes to PDAF repository

1. the PDAF-library with slight changes in
   1. `make.arch/`: changed include files, in particular added the
include file `cmake.h` for TSMP2-PDAF
   2. `src/`: debug output of the observation ensemble for EnKF/LEnKF
2. `interface/` the interface routines for TSMP-PDAF
3. `.github/`: GitHub workflow files
   1. `docs.yml`: Deploying the user documentation.
4. `docs/`: Documentation sources.


## Updating the pre-patched PDAF repository

The pre-patched PDAF is updated with three types of changes. Each type
of change calls for diferent kind of 

1. Updating the interface and documentation.
   - Feel free to open a PR.
   - Try to follow the style of coding in the existing interface.

2. Version updates from PDAF
   - **Remark:** Only do this in collaboration with repository
     administrators.
   - **Explanation**: PDAF version updates are pulled from PDAF's
     branch `master` into the branch `tsmp-pdaf-patched`. For a PDAF
     version, f.e. `PDAF_V2.2.1` in the `master`, a corresponding
     version tag is defined, f.e. `PDAF_V2.2.1-tsmp` in
     `tsmp-pdaf-patched`.

3. TSMP-related changes in the PDAF library
   -  **Remark**: Avoid TSMP-related changes in the PDAF library. If
   necessary, collaborate with repository administrators and carefully
   document the changes in the list of changes above.
