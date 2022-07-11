# ppa_models

The following MATLAB scripts are provided:

cdf_ks_spice.m: Delay yield estimation for 16-bit Kogge-Stone adder 
cdf_ks_32_spice.m: Delay yield estimation for 32-bit Kogge-Stone adder 
cdf_sl_spice.m: Delay yield estimation for 16-bit Sklansky adder 
cdf_kn1112_spice.m: Delay yield estimation for 16-bit Knowles (1,1,1,2) adder 
cdf_kn1222_spice.m: Delay yield estimation for 16-bit Knowles (1,2,2,2) adder 
cdf_bk_spice.m: Delay yield estimation for 16-bit Brenk-Kung adder 
cdf_hc_spice.m: Delay yield estimation for 16-bit Han-Carlson adder 
cdf_lf_spice.m: Delay yield estimation for 16-bit Ladner-Fischer adder 

Models are constructed based on Spice mean maximum delay and standard deviation for 
each cell.

The following files comprise experimental Monte-Carlo based maximum-delay data for 
parallel-prefix adders:

KOGGE_STONE_VTHINTRA_hist.txt
KOGGE_STONE_32_VTHINTRA_hist.txt 
SKLANSKY_VTHINTRA_hist.txt 
KNOWLES1112_VTHINTRA_hist.txt 
KNOWLES1222_VTHINTRA_hist.txt 
BRENT_KUNG_VTHINTRA_hist.txt 
HAN_CARLSON_VTHINTRA_hist.txt 
LADNER_FISCHER_VTHINTRA_hist.txt 
