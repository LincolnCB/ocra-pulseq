# ocra-pulseq
Pulseq interpreter for OCRA

# Usage:
`PSAssembler` is an object that initializes with timing and system specifications. After initializing, run `PSAssembler.assemble("[filepath]")` to return a list of `[tx_bytes, [gx_bytes, gy_bytes, gz_bytes], command_bytes`