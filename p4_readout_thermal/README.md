# P4: readout noise and thermal design of the camera head (stub)

Not started. Scope: read-noise distribution versus CMS sampling and line time, dark current D(T)
and the cooled-head thermal design (see the P4 project prompt). Outputs: `read_noise`
(`GaussianRTSMixture` or `EmpiricalSigmaCDF`), `read_noise_median_e`, `read_noise_rms_e`,
`shutter.line_time_us` and `dark` (`DarkArrhenius`) entries of `camsim_common.SensorModel`.
