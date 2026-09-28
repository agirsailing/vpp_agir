function [Fn, a] = delft_coefficients()
%DELFT_COEFFICIENTS Rows: Fn; columns: a0,...,a7 (no extra 10^3 factor).
% J.A. Keuning & M. Katgert (2008), A bare hull resistance prediction method
% derived from the results of the Delft Systematic Yacht Hull Series
% extended to higher speeds, INNOVSAIL, pp. 13-21. Table 2, PDF page 6.
% https://research.tudelft.nl/en/publications/a-bare-hull-resistance-prediction-method-derived-from-the-results/
% Visually checked against the original typeset table at:
% https://www.scribd.com/document/346520042/Keuning-2008-pdf
Fn = (15:5:75)' / 100;
a = [ ...
    -0.0005  0.0023 -0.0086 -0.0015  0.0061  0.0010  0.0001  0.0052;
    -0.0003  0.0059 -0.0064  0.0070  0.0014  0.0013  0.0005 -0.0020;
    -0.0002 -0.0156  0.0031 -0.0021 -0.0070  0.0148  0.0010 -0.0043;
    -0.0009  0.0016  0.0337 -0.0285 -0.0367  0.0218  0.0015 -0.0172;
    -0.0026 -0.0567  0.0446 -0.1091 -0.0707  0.0914  0.0021 -0.0078;
    -0.0064 -0.4034 -0.1250  0.0273 -0.1341  0.3578  0.0045  0.1115;
    -0.0218 -0.5261 -0.2945  0.2485 -0.2428  0.6293  0.0081  0.2086;
    -0.0388 -0.5986 -0.3038  0.6033 -0.0430  0.8332  0.0106  0.1336;
    -0.0347 -0.4764 -0.2361  0.8726  0.4219  0.8990  0.0096 -0.2272;
    -0.0361  0.0037 -0.2960  0.9661  0.6123  0.7534  0.0100 -0.3352;
     0.0008  0.3728 -0.3667  1.3957  1.0343  0.3230  0.0072 -0.4632;
     0.0108 -0.1238 -0.2026  1.1282  1.1836  0.4973  0.0038 -0.4477;
     0.1023  0.7726  0.5040  1.7867  2.1934 -1.5479 -0.0115 -0.0977];
end
