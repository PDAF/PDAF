# Debugging TSMP2-PDAF

## Debug Tip 1: Build with debug version

The following command will build eCLM-PDAF in debug version
``` shell
./build_tsmp2.sh --eCLM --PDAF --build_type DEBUG
```

## Debug Tip 2: Building TSMP-PDAF after source code changes in component models

This is under construction in TSMP2

## Debug Tip 3: Debugging OASIS3-MCT

In `namcouple`, change the input `$NLOGPRT` to get more debug
output. More information:
<https://cerfacs.fr/oa4web/oasis3-mct_5.0/oasis3mct_UserGuide/node39.html>

## Debug Tip 4: Debugging CLM (3.5)

In `lnd.stdin`, change the input `wrtdia` to `.true.`.

```
If true, global average 2-m temperature written to standard out (ascii log file of the run) (see ex. 4). 
```

More information:
<https://www2.cgd.ucar.edu/tss/clm/distribution/clm3.0/UsersGuide/UsersGuide/node6.html>

## Debug Tip 5: Debugging PDAF

### Screen variables

Turn the screen variables on completely in order to get maximum output of PDAF and the PDAF-interface.
- [`screen`](cmd:screen)
- [`screen_wrapper`](enkfpf:da:screen_wrapper)

### PDAF Debugging Information

For TSMP-PDAF, debug output can be called via the preprocessor
variable [`PDAF_DEBUG`](precompiler:pdaf_debug).

Information about turning on PDAF debugging output by changing the
source code of TSMP-PDAF:

<https://pdaf.awi.de/trac/wiki/PDAF_debugging>


### PDAF's Timing and Memory output

`/intf_DA/pdaf/framework/finalize_pdaf.F90`: The input for the timing
information can be changed for adding more output.

See also <https://pdaf.awi.de/trac/wiki/PDAF_print_info>.

