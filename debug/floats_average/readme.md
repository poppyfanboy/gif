Results of computing an average of one million random floats in the interval between -500000 and 1500000.

```
average_naive
================
Mean absolute error: 5.6223437499999998
================
437.6 micros per call

average_naive_in_f64
================
Mean absolute error: 0.0092187499999999995
================
430.1 micros per call

average_multi_x4_in_f64
================
Mean absolute error: 0.0092187499999999995
================
109.3 micros per call

average_multi_x16_avx_in_f64
================
Mean absolute error: 0.0092187499999999995
================
62.2 micros per call

average_kahan
================
Mean absolute error: 0.00015625
================
1.7 millis per call

average_block_kahan
================
Mean absolute error: 0
================
126.3 micros per call

average_block_kahan_avx
================
Mean absolute error: 0.00296875
================
54.5 micros per call

average_pairwise
================
Mean absolute error: 0.0096874999999999999
================
2.3 millis per call

average_block_pairwise
================
Mean absolute error: 0.010625000000000001
================
158.7 micros per call

average_python_fsum
================
Mean absolute error: 0
================
11.1 millis per call
```

Results of computing an average of the colors from an actual image. The concrete values can be found in `input.txt`. The file contains a list of values for the color component "a" (from the "Lab" colorspace), so there might be negative values or values with an absolute value greater than 1.

```
average_naive
================
expected value: 6.98441601
actual value: 6.98442507
================
Mean absolute error: 9.059906005859375e-06
================
3.5 micros per call

average_naive_in_f64
================
expected value: 6.98441601
actual value: 6.98441601
================
Mean absolute error: 0
================
3.5 micros per call

average_multi_x4_in_f64
================
expected value: 6.98441601
actual value: 6.98441601
================
Mean absolute error: 0
================
873.2 nanos per call

average_multi_x16_avx_in_f64
================
expected value: 6.98441601
actual value: 6.98441601
================
Mean absolute error: 0
================
312.4 nanos per call

average_kahan
================
expected value: 6.98441601
actual value: 6.98441601
================
Mean absolute error: 0
================
14.0 micros per call

average_block_kahan
================
expected value: 6.98441601
actual value: 6.98441601
================
Mean absolute error: 0
================
1.0 micros per call

average_block_kahan_avx
================
expected value: 6.98441601
actual value: 6.98441648
================
Mean absolute error: 4.76837158203125e-07
================
420.6 nanos per call

average_pairwise
================
expected value: 6.98441601
actual value: 6.98441601
================
Mean absolute error: 0
================
18.4 micros per call

average_block_pairwise
================
expected value: 6.98441601
actual value: 6.98441601
================
Mean absolute error: 0
================
1.3 micros per call

average_python_fsum
================
expected value: 6.98441601
actual value: 6.98441601
================
Mean absolute error: 0
================
56.2 micros per call
```

So, when considering a real-world example, even the naive algorithm which uses single-precision floats might be enough to compute an average of a bunch of values.
