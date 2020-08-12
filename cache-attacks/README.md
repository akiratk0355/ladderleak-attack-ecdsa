## Scripts

### `build.sh`
- Download OpenSSL `1.0.2u` and build with the default configuration.

### `ecdsa_sign.sh`
- Generage ECDSA signatures with a given curve parameter, save output log under the specified directory.
- Usage `ecdsa_sign.sh <#signatures> <curve_parameter> <result_directory>`

### `fr-trace-gf2m.sh` and `fr-trace-gfp.sh`
- Wrapper scripts for `FR-trace`.

### `analysis-gf2m.py` and `analysis-gfp.sh`
- Detect the second most significant bit of nonces given traces and log files stored in the result directory, calculate the accurary.
- Usage `python analysis-*.py <result_directory>`

### `traces-gf2m.py` and `traces-gfp.py`
- Utility scripts to visualize traces.
- Usage `python traces-gf*.py <result_directory>`

### Combined attack script
- Run `./run-gf2m.sh` or `./run-gfp.sh` to generate key pairs (of sect163r1 and secp192r1 respectively), run Flush+Reload, and analyze traces, and calculate the success rate.

## Experimental results
### Binary curves
- Run `python analysis-gf2m.py results-gf2m-1000` to see the experimental results for 1000 sect163r1 signatures.

### Prime curves
- Run `python analysis-gfp.py results-gfp-1000-r1` or `-k1` to see the experimental results for 1000 secp192r1 or secp192k1 signatures, respectively. The same applies to the respective folders with 10000 signatures in both curves.
 
