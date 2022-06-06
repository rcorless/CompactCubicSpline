# CompactCubicSpline

Code supporting the paper at https://arxiv.org/abs/1805.07659

file compactslopesA.m generates the matrix given the diff(tau) = h, ie the widths of each interval.  Assumed positive.

file Randocondsdist.m generates 10,000 matrices with random h's of dimension 10946 (a Fibonacci number, which helps in doing tests of increasing dimension, but not really used for anything here) and shows a histogram of the condition numbers of the matrices.  We see that usually it is less than the dimension, and sometimes as large as the dimension squared.
