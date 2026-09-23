# 03 - W2O6: chemical-path selection

```bash
bash examples/paper/03_w2o6_chemistry/run.sh
```

The coarse electronic reduction uses O `2p` and W `5d` channels. At the 3.4 A graph cutoff used in the paper there are W-O and O-O edges but no W-W edges. The code reports edge-resolved couplings and the raw W(5d)-O(2p) level difference.

The raw level contrast is a screening scale, not a DFT-quality trap energy. `reproduction/paper/atomic_cases.py` compares zero site contrast with the raw EHT screening contrast in first-passage transport.
