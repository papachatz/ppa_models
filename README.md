# File Description

MATLAB scripts are provided that evaluate the delay yield error for various adder topologies.
The following scripts are provided:

[cdf_ks_spice.m](https://github.com/papachatz/ppa_models/blob/main/matlab/cdf_ks_spice.m): Delay yield estimation for 16-bit Kogge-Stone adder <br />
[cdf_ks_32_spice.m](https://github.com/papachatz/ppa_models/blob/main/matlab/cdf_ks_32_spice.m): Delay yield estimation for 32-bit Kogge-Stone adder <br />
[cdf_sl_spice.m](https://github.com/papachatz/ppa_models/blob/main/matlab/cdf_sl_spice.m): Delay yield estimation for 16-bit Sklansky adder <br />
[cdf_kn1112_spice.m](https://github.com/papachatz/ppa_models/blob/main/matlab/cdf_kn1112_spice.m): Delay yield estimation for 16-bit Knowles (1,1,1,2) adder <br />
[cdf_kn1222_spice.m](https://github.com/papachatz/ppa_models/blob/main/matlab/cdf_kn1222_spice.m): Delay yield estimation for 16-bit Knowles (1,2,2,2) adder <br />
[cdf_bk_spice.m](https://github.com/papachatz/ppa_models/blob/main/matlab/cdf_bk_spice.m): Delay yield estimation for 16-bit Brenk-Kung adder <br />
[cdf_hc_spice.m](https://github.com/papachatz/ppa_models/blob/main/matlab/cdf_hc_spice.m): Delay yield estimation for 16-bit Han-Carlson adder <br />
[cdf_lf_spice.m](https://github.com/papachatz/ppa_models/blob/main/matlab/cdf_lf_spice.m): Delay yield estimation for 16-bit Ladner-Fischer adder <br />
[cdf_ks_udm.m](https://github.com/papachatz/ppa_models/blob/main/matlab/unit_delay/ks_adder.m): Delay yield estimation based on unit delay model for Kogge-Stone adder 

The delay models are constructed based on Spice mean maximum delay and standard deviation 
for each cell.

The following files comprise experimental Monte-Carlo-based data, where the first column lists 
maximu-delay data and the second list the dissipated total power, for parallel-prefix adders:

[KOGGE_STONE_VTHINTRA_hist.txt](https://github.com/papachatz/ppa_models/blob/main/matlab/KOGGE_STONE_VTHINTRA_hist.txt) <br />
[KOGGE_STONE_32_VTHINTRA_hist.txt](https://github.com/papachatz/ppa_models/blob/main/matlab/KOGGE_STONE_32_VTHINTRA_hist.txt) <br />
[SKLANSKY_VTHINTRA_hist.txt](https://github.com/papachatz/ppa_models/blob/main/matlab/SKLANSKY_VTHINTRA_hist.txt) <br />
[KNOWLES1112_VTHINTRA_hist.txt](https://github.com/papachatz/ppa_models/blob/main/matlab/KNOWLES1112_VTHINTRA_hist.txt) <br />
[KNOWLES1222_VTHINTRA_hist.txt](https://github.com/papachatz/ppa_models/blob/main/matlab/KNOWLES1222_VTHINTRA_hist.txt) <br />
[BRENT_KUNG_VTHINTRA_hist.txt](https://github.com/papachatz/ppa_models/blob/main/matlab/BRENT_KUNG_VTHINTRA_hist.txt) <br />
[HAN_CARLSON_VTHINTRA_hist.txt](https://github.com/papachatz/ppa_models/blob/main/matlab/HAN_CARLSON_VTHINTRA_hist.txt) <br />
[LADNER_FISCHER_VTHINTRA_hist.txt](https://github.com/papachatz/ppa_models/blob/main/matlab/LADNER_FISCHER_VTHINTRA_hist.txt) <br />

Python scripts for the proposed algorithmic construction of transformation matrix A for parallel-prefix nodes are provided: <br />
[KS16.py](https://github.com/papachatz/ppa_models/blob/main/python/KS16.py) <br />
[LF16.py](https://github.com/papachatz/ppa_models/blob/main/python/LF16.py) <br />
[HC16.py](https://github.com/papachatz/ppa_models/blob/main/python/HC16.py) <br />
[BK16.py](https://github.com/papachatz/ppa_models/blob/main/python/BK16.py) <br />


# Contributors

Kleanthis Papachatzopoulos, and Vassilis Paliouras

# License 

The software is licensed under the [MIT License](https://github.com/papachatz/ppa_models/blob/main/LICENSE). 


