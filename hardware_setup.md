# Hardware set-up: what the low-cost data taught

## 1. The meter used for the tests

The recordings come from the **BidyutPi** open-source meter: a Raspberry Pi with an external ADC, voltage and current sensing, nominally 10 kHz, and 60 s per load.

Data: Mendeley Data, doi:[10.17632/5jp4r7nnpw.2](https://doi.org/10.17632/5jp4r7nnpw.2).

What the data showed ([`results/02_timing_check_output.txt`](02_timing_check_output.txt)):

| Problem | Evidence | Effect on detection |
|---|---|---|
| Sample-clock jitter | sampling interval 85–146 µs (nominal 100 µs) | looks like noise: apparent noise fell from 0.087 to 0.037 p.u. after resampling on the real time stamps |
| Recording gaps | 278 gaps, the largest 17.2 ms | the strongest "events" in the raw file were these gaps |
| Heavy-tailed noise | default threshold gave 11–14 false alarms/s on a clean voltage | the threshold factor has to be calibrated (k ≈ 7.5–11.5 instead of 3.9) |

The first two problems come from **operating-system timing**: the Pi reads the ADC in software. They are not grid phenomena.

## 2. Recommended set-up for an industrial module

### Signal chain

```
sensor → isolation / conditioning → anti-alias filter → ADC (timer + DMA) → MCU → outputs
```

- **Sampling.** The ADC conversion is started by a hardware timer, and samples are moved to memory by DMA. This removes jitter and gaps.
- **ADC.** 12 bit is enough for the tests here: 1 p.u. = 2048 counts, and the largest internal value is 2.7·10⁸, which fits in 32-bit integers. Use 16 bit for small events.
- **Anti-alias filter.** An RC or second-order filter at about 0.4 × sampling rate.
- **Isolation.** Isolated amplifier or PT/CT for mains; IEPE supply for accelerometers; 4–20 mA / 0–10 V inputs for process sensors.
- **MCU.** Any 32-bit microcontroller with a timer-triggered ADC and DMA is enough. The algorithm needs about 20 integer operations per sample, which is about 1 MOPS at 50 kHz.

### Outputs

- Relay or transistor alarm output and a status LED
- Modbus RTU (RS-485) or Modbus TCP: event count, last event time, peak size, noise level
- An event log with time stamps

### Set-up procedure (what the customer does)

1. Choose a profile: power quality, bearing, pulse counting or pressure. The profile sets W and Tmin.
2. Set the false-alarm budget, for example 0.1 /s.
3. Press **Calibrate** while the machine runs normally. The tests here used 15 s; about 60 s is recommended in a product.
4. Run.

## 3. Before selling (to check with a test lab)

For an industrial product, electromagnetic compatibility (EMC) and electrical safety tests are normally needed. Which standards apply depends on the product and the market, so check this with an accredited lab or a standards consultant.
