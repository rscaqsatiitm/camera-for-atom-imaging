# P2: thick BSI pixel optics and charge transport (stub)

Not started. Scope: QE(lambda, T), fringing, MTF and dark current versus silicon thickness, bias and
AR stack (see the P2 project prompt). Outputs go into `camsim_common.SensorModel` fields:
`qe` (QETable with temperature axis), `qe_780`, `qe_852`, `mtf` (table or `kernel_sigma_um`) and
`cosmic` (muon deposit versus thickness). Every value must be a `ProvValue`; P1 consumes the
resulting YAML presets unchanged.
