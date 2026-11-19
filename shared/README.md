# Shared code policy

This folder is reserved for routines whose mathematical definition and
numerical behavior are identical across datasets.  No numerical routines have
been moved here in the initial release because the two simulations use
different covariance references, spatial multipliers, and reporting paths.

Prematurely merging similarly named functions can hide scientifically
important differences.  A routine should move here only after the dataset
implementations have matching interfaces and regression tests.

