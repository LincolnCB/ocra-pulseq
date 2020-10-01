# ocra-pulseq
Pulseq interpreter for OCRA -- upon testing, will be merged into the [OCRA repository](https://github.com/OpenMRI/ocra)

# Information:
Current OCRA software uses hand-written machine code to function -- this limits ease of use. [PulSeq](https://pulseq.github.io/) is an open-source library widely used in MR pulse design, available in MATLAB and Python. The library primarily functions by generating a `.seq` file, which scanners interpret. `pulseq_assembler.py` is an assembler that takes a `.seq` file and assembles OCRA machine code and data to the specifications of the file. 

# Usage
`PSAssembler` is an object that initializes with timing and system specifications.

After initializing, run `PSAssembler.assemble("[filepath]")` to return a list of:

`[tx_bytes, [gx_bytes, gy_bytes, gz_bytes], command_bytes]`
