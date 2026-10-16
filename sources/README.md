# Source retrieval record

Sources retrieved and read on 2026-10-14 (session of the P1 validation stage). PDFs are not
committed (copyright); the SHA-256 identifies the exact file that was read. Re-fetch with the URL.

| Source | Version read | URL | SHA-256 of the PDF | Used for |
|---|---|---|---|---|
| Manetsch, Nomura, Bataille, Leung, Lv, Endres, "A tweezer array with 6100 highly coherent atomic qubits" | arXiv:2403.12021v4 (29 Jul 2025) | https://arxiv.org/abs/2403.12021v4 | abd0eef3e6be178d68037230a4e8ec8b43e88c93881d5771377cc825ea78bb70 | V1 (Methods "Single-atom imaging", "Imaging model and characterization"; Fig. 2, 3; Ext. Data Fig. 4, 6) |
| Picken, Legaie, Pritchard, "Single Atom Imaging with an sCMOS camera" | arXiv:1709.00946v2 (5 Oct 2017) | https://arxiv.org/abs/1709.00946v2 | db88cbadf67323d06e61439d744707afef4e51a535f083cde90e8be0e3b89d2a | V2 (Sec. II-V, Eq. 1-5, Tables I-II) |
| Basden, Haniff, Mackay, "Photon counting strategies with low light level CCDs", MNRAS 345, 985 (2003) | arXiv:astro-ph/0307305v1 (16 Jul 2003) | https://arxiv.org/abs/astro-ph/0307305v1 | 34324539b46d7d1bddcb03190782318c95b5a658dc52b3f588ab31c369b12aab | EMCCD register model (Eq. 1, 5, App. A) |
| Steck, "Cesium D Line Data" | rev. 2.3.4 (8 Aug 2025) | https://steck.us/alkalidata/cesiumnumbers.pdf | f3f6fefb10447fd26ec2076a7c2049a5043b45063a66ffbfec4d04367f2023a4 | Cs D2 constants (Tables 2, 3, 7) |
| Steck, "Rubidium 87 D Line Data" | rev. 2.3.4 (8 Aug 2025) | https://steck.us/alkalidata/rubidium87numbers.pdf | 0d7ecd708e44ad950db280cbf980b16090660646f5e530e544ae78ccd6385c3b | Rb87 D2 constants (Tables 2, 3, 7) |
| "Achieving higher photoabsorption than group III-V semiconductors in silicon using photon-trapping surface structures" | arXiv:2306.05170v1 | https://arxiv.org/abs/2306.05170v1 | a8cef95e3f7c74be4eee65fde3cbc1a52e89252ffe5e87d423b23d5ceb19a83f | retrieved for P3; not yet read |
| "Technology advancement of the CCD201-20 EMCCD for the WFIRST coronagraph instrument" | arXiv:1601.01761v1 | https://arxiv.org/abs/1601.01761v1 | 52d4e1e3301ec03bcee0d378c970119550d9121ca95f9b902fc6aae780876ce5 | retrieved for P2; not yet read |
| PDG Review "Cosmic Rays" (2024) | rpp2024-rev-cosmic-rays | https://pdg.lbl.gov/2024/reviews/rpp2024-rev-cosmic-rays.pdf | c8f0620d58d3d61a7b0eae5d2606ce65bbe581a9000c5435299d88ca9ea0125e | searched for the sea-level muon flux; not found in this review, value remains a placeholder |

Journal reference of Basden et al. verified via Crossref on 2026-10-14: MNRAS 345, 985-991 (2003),
doi:10.1046/j.1365-8711.2003.07020.x. The content read is the arXiv preprint version.

Published versions found via Crossref on 2026-10-14 (metadata only; the text read is the arXiv version
listed above, and numbers were not checked against the journal version):
Manetsch et al., Nature 647, 60-67 (2025), doi:10.1038/s41586-025-09641-4;
Picken et al., Appl. Phys. Lett. 111 (2017), doi:10.1063/1.5003304.

## Figure digitisation (DATA)

`data/manetsch_ed4b_inset_digitised.csv` and `data/manetsch_ed4d_rate_hist_digitised.csv` were
digitised from 400 dpi renders of page 14 of the v4 PDF (pdftoppm), by locating the top edge of
the filled histogram in each pixel column and calibrating against the printed axis ticks.
Moments of Ext. Data Fig. 4b (main linear panel) used in V1 were obtained the same way. The
estimated resolution is 0.2 photon on the x axis and about 1 px (0.25% of the peak) on linear
y axes; 0.04 decades on the log inset.

## Production stage (2026-10-15): camera documents

USER-supplied PDFs (uploaded by the PI; read in this session):

| Document | SHA-256 (prefix) | Used for |
|---|---|---|
| Teledyne Photometrics, Prime BSI Express datasheet, revision date 2024-08-20 | 1c867b6f09affb24 | modes table (read noise, line time, full well, bit depth, dark current), QE curve (digitised) |
| Hamamatsu ORCA-Quest 2 C15550-22UP, Cat. No. SCAS0166E03 (Jun. 2026) | 8497f943da734aad | spec table (read noise median/rms per scan, full well, PRNU, dark current, frame rates), QE curve p.6 (digitised) |
| Andor iXon Ultra EMCCD specifications (888 and 897) | d0b88cd4c37144ed | read noise per speed, EM gain, spurious background, dark current, well depths, crop-mode frame rates, QE chart (digitised) |
| Andor Marana sCMOS specifications, MARANAFAMILY SS 0926 R1 | 7863aba3351af711 | read noise (median), well depth, dark current, frame rates, QE chart (digitised) |

Vendor documents retrieved on 2026-10-15 by retrieval sub-agents of this session and spot-checked by
direct reading of the quoted lines (files kept in the session scratchpad, not committed):

| Document | URL | SHA-256 (prefix) | Used for |
|---|---|---|---|
| Hamamatsu C15550-22UP instruction manual B9200301-02 Ver.1.2 (Mar. 2025) | hamamatsu.com/.../manual/C15550-22UP_IM_En.pdf | 5417560f644925ff | conversion factor 0.107 e-/count, 1H line times, PNR output range, dark offset |
| Hamamatsu ORCA-Fusion BT instruction manual B0850301-05 Ver.1.5 (Oct. 2025) | hamamatsu.com/.../manual/C15440-20UP,-20UP01_IM_En.pdf | ffa23b21a177e454 | read noise per scan (p.80, checked), conversion factor, dark current, 1H |
| Hamamatsu ORCA-Fusion BT brochure SCAS0143E03 (Jun. 2026) | hamamatsu.com/.../SCAS0143E_C15440-20UP.pdf | e481b889d2851dba | QE curve (digitised, +-2%) |
| Hamamatsu ORCA-Quest technical note SCAS0154E06 (Apr. 2025) | hamamatsu.com/.../SCAS0154E_C15550-20UP_tec.pdf | 91784acb9abfc31a | original Quest QE curve: 39.6% at 852 nm (cf. 44% stated by Manetsch et al.) |
| Andor "Zyla for Physical Sciences sCMOS Specifications" LZYLASPECTSS 1020 R1 | oxinst.com/learning/uploads/Andor_Zyla_for_Physical_Sciences_Specifications.pdf | a791a99950be137d | read noise median [rms] rolling and global shutter (p.6, checked), dark current, QE chart (digitised) |
| FLIR EMVA 1288 data, BFS-U3-17S7 (11/15/2022) | softwareservices.flir.com/BFS-U3-17S7/latest/EMVA/EMVA.html | (HTML) | read noise 4.69 e-, gain 0.31 e-/ADU, saturation 19557 e- (checked); QE curve PNG (digitised +-1%) |
| FLIR BFS-U3-17S7 specification page | softwareservices.flir.com/BFS-U3-17S7/latest/Model/spec.html | (HTML) | global shutter, 9 um pixels, ADC (checked) |
| QHYCCD QHY268 product page (modified 2026-06-16) and QE article (2021-05-24) | qhyccd.com/astronomical-camera-qhy268/ | (HTML) | 3.76 um, rolling shutter, dark current, read-noise chart (1.58 e- HG-2CMS gain 56, checked); absolute QE curve (digitised; +-3% rel. scale) |
| Gpixel GSENSE2020BSI flyer GP-PR250801 V1.1 and QE curve | gpixel.com | db3bc5559493e0d6 | cross-check: 48.1% at 852 nm vs 48.5% from the Prime BSI curve |

QE at 780 and 852 nm is not printed in any of these documents; every such value is DATA,
digitised from a typical-performance curve, with the method and uncertainty in the preset pointer.
