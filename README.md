# ocra-pulseq
Pulseq interpreter for OCRA -- upon testing, will be merged into the [OCRA repository](https://github.com/OpenMRI/ocra)

# Information:
Current OCRA software uses hand-written machine code to function -- this limits ease of use. [PulSeq](https://pulseq.github.io/) is an open-source library widely used in MR pulse design, available in MATLAB and Python. The library primarily functions by generating a `.seq` file, which scanners interpret. `pulseq_assembler.py` is an assembler that takes a `.seq` file and assembles OCRA machine code and data to the specifications of the file. 

# Usage
`PSAssembler` is an object that initializes with timing and system specifications.

After initializing, run `PSAssembler.assemble("[filepath]")` to return a list of:

`[tx_data, [gx_data, gy_data, gz_data], command_bytes, output_dict]` with data in either `bytes` format or `numpy.ndarray`

You can pass these to existing OCRA servers, such as [MaRCoS](https://github.com/vnegnev/marcos_extras/wiki/setting_marcos_up)

# Initialization Parameters

`rf_center` (int): RF center (local oscillator frequency) in Hz.

`rf_amp_max` (int): Default 5e+3 -- System RF amplitude max in Hz.

`grad_max` (int): Default 1e+6 -- System gradient max in Hz/m.

`clk_t` (float): Default 7e-3 -- System clock period in us.

`tx_t` (float): Default 1.001 -- Transmit raster period in us.

`grad_t` (float): Default 10.003 -- Gradient raster period in us.

`pulseq_t_match` (bool): Default False -- If PulSeq file transmit and gradient raster times match OCRA transmit and raster times.

`ps_tx_t` (float): Default 1 -- PulSeq transmit raster period in us. Used only if pulseq_t_match is False.

`ps_grad_t` (float): Default 10 -- PulSeq gradient raster period in us. Used only if pulseq_t_match is False.

`tx_warmup` (float): Default 0 -- Padding delay at the beginning of RF to give TR warmup in us. PADDING WILL CHANGE TIMING.

`grad_pad` (int): Default 0 -- Padding zero samples at the end of gradients to prevent maintained gradient levels. PADDING WILL CHANGE TIMING.

`adc_pad` (int): Default 0 -- Padding samples in ADC to account for junk in system buffer. PADDING WILL CHANGE TIMING.

`addresses_per_grad_sample` (int): Default 1 -- Memory offset step per gradient readout, to account for different DAC boards.

