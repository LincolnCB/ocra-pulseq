# ocra-pulseq
Pulseq interpreter for OCRA -- upon testing, will be merged into the [OCRA repository](https://github.com/OpenMRI/ocra)

# Information:
Current OCRA software uses hand-written machine code to function -- this limits ease of use. [PulSeq](https://pulseq.github.io/) is an open-source library widely used in MR pulse design, available in MATLAB and Python. The library primarily functions by generating a `.seq` file, which scanners interpret. `pulseq_assembler.py` is an assembler that takes a `.seq` file and assembles OCRA machine code and data to the specifications of the file. 

# Usage
`PSAssembler` is an object that initializes with timing and system specifications.

After initializing, run `PSAssembler.assemble("[filepath]")` to return a list of:

`[tx_bytes, [gx_bytes, gy_bytes, gz_bytes], command_bytes, readout_count]`

You can pass these to existing OCRA servers, such as [MaRCoS](https://github.com/vnegnev/marcos_extras/wiki/setting_marcos_up)

# Initialization Parameters

`rf_center` (int): RF center (local oscillator frequency) in Hz

`rf_amp_max` (int): RF amplitude max in Hz -- used to set the fractional amplitude from PulSeq's Hz amplitude

`grad_max` (int): Gradient max in Hz/m -- used to set the fractional amplitude from PulSeq's Hz/m amplitude

`clk_t` (float): Hardware clock period in us

`tx_t` (float): Transmit raster period in us -- should be a multiple of clk_t, will round otherwise

`grad_t` (float): Gradient raster period in us -- should be a multiple of clk_t, will round otherwise

`pulseq_t_match` (bool): Default True -- Set to False if PulSeq file transmit and gradient raster times do not match OCRA transmit and raster times.

`ps_tx_t` (float): PulSeq transmit raster period in us, if pulseq_t_match is False

`ps_grad_t` (float): PulSeq gradient raster period in us, if pulseq_t_match is False

`grad_pad` (int): Default 0 -- Padding zeros at the end of gradients to prevent maintained gradient levels. Currently may change timing.

`addresses_per_grad_sample` (int): Default 1 -- Memory offset step per gradient readout, to account for different DAC boards with byte interleaving

