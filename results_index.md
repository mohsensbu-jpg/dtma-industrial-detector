# Results (text output of the author's MATLAB Online runs)

| File | Script | Note |
|---|---|---|
| `01_first_test_output.txt` | `PQ_first_test.m` | First exploratory run. Exp B = synthetic grid with noise whose level changes: the local noise rule against a global threshold, wavelet and cycle comparison. Exp A used the meter data before the timing fix (superseded by 03–05). |
| `02_timing_check_output.txt` | `PQ_timing_check.m` | Sample-clock jitter and recording gaps of the low-cost meter |
| `03_four_loads_output.txt` | `PQ_four_loads.m` | Four loads, gaps masked; false alarms on clean voltage, calibration |
| `04_lowcost_equal_FA_output.txt` | `PQ_lowcost_equal_FA.m` | Detection rate at equal false-alarm budget, meter noise only |
| `05_online_integer_output.txt` | `PQ_online_integer.m` | Online integer version: detection rate, delay, operations, memory, overflow check |
| `06_calibration_fix_output.txt` | `PQ_calibration_fix.m` | Counting vs tail-fit calibration, 16 calibration/test pairs |
| `07_adaptive_threshold_output.txt` | `PQ_adaptive_threshold.m` | Noise trackers and threshold strategies; **`static/q99` = final detector** (specifications in the main README) |

Script names inside the output headers are the working names used during development (`PQ_TMA_step…`).
