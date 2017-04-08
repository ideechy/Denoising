# Denoising with subspace learning

Bowen Liu, Renfei Gong, Haoyu Chen

## Description

Denoising is a real world problem that embed low-dimensional structure into a high-dimensional manifold. Unlike ordinary subspace learning or information reduction problems, which have more attributes (rank r) than the dimension of target subspace (d), denoising have r << d. This requires new treatment of the data, and hence an interesting problem from a computational perspective. To handle subspace learning with partial information, we will mainly exploit Partially Observed PCA and Matrix Bandit Exponentiated Gradient algorithm by Gonen et al. (2016), and Parallel Estimation and Tracking by REcursive Least Squares algorithm by Chi et al. (2013).
Our main contribution is the application of these algorithms to the denoising problem. We will also compare of the denoising result and computational price across the algorithms. As our project goes along, we may give some revision to these general algorithms such that they can be better applied to the specific denoising problem.

## Weekly goal

- [x] 3/15 - 3/21: Understand the theory part of the two paper.
- [x] 3/22 - 3/28: Realize the algorithms.
- [x] 3/29 - 4/4 : Debug and data test.
- [ ] 4/5  - 4/11: prepare the milestone report.
- [ ] 4/12 - 4/18: Try to improve existed algorithms.
- [ ] 4/19 - 4/25: Make poster and prepare for presentation.
- [ ] 4/26 - 5/1 : Finish final report.

## Miscellaneous

Clone it to your local folder:
```
git clone https://github.ncsu.edu/hchen36/Denoising.git
```
