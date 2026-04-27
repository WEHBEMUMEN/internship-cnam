# DEIM Error Analysis & Feedback

Based on investigating the `dynamicslab/databook_matlab` repository (specifically the Chapter 12 companion code `CH12_SEC06_1_DEIM.m`), here is the mathematical feedback on where the massive DEIM errors originate and how to fix them.

## 1. Architectural Verification (We are doing it right!)
The Matlab script clearly shows that DEIM requires **two separate SVDs**:
```matlab
[U,S,W]=svd(X,0);          % SVD of displacements (Galerkin Basis)
...
[XI,S_NL,W]=svd(NL,0);     % SVD of the NON-LINEAR term (DEIM Basis)
```
In our Phase 3.4a JavaScript implementation, we are successfully following this architecture. We compute the SVD of the displacements in `romEngine.computePOD()` and a completely separate SVD of the internal forces in `deimEngine.train()`. The error is **not** coming from a flawed architecture.

## 2. The Root Cause of Massive Errors: Ill-Conditioning
The massive error spikes (divergences and 15%+ errors) in DEIM are almost always caused by an **ill-conditioned interpolation matrix**.
In DEIM, the non-linear force is reconstructed using:
$F_{int} \approx \Phi_F (P^T \Phi_F)^{-1} F_{sampled}$

The standard Greedy Algorithm selects points by finding the maximum residual. However, occasionally, the greedy algorithm selects a point that is "mathematically close" to an already selected point in the feature space. When this happens, the matrix $(P^T \Phi_F)$ becomes nearly singular (its determinant approaches zero). 
When you invert a nearly singular matrix, small numerical noises explode to infinity, causing the "massive error" or `Divergence detected` logs you are seeing.

## 3. The Solution: Q-DEIM (QR Factorization)
At the very bottom of the `CH12_SEC06_1_DEIM.m` script, Brunton explicitly introduces an alternative to the Greedy loop:
```matlab
%% QR DEIM
[Q,R,pivot]=qr(NL.');
P_qr=pivot(:,1:3);
```
**Why does he include this?**
Because of the exact problem you are facing! Researchers discovered that doing a **QR Factorization with Column Pivoting** on the transposed force basis perfectly selects the optimal interpolation indices. 

Unlike the Greedy loop, **Q-DEIM** mathematically guarantees that the selected points are as orthogonal as possible, preventing the interpolation matrix from ever becoming ill-conditioned. 

### Recommended Next Steps
If you want to completely eliminate these massive error spikes, our next task should be to upgrade the point selection logic in `deim-engine.js` from the **Standard Greedy Algorithm** to the **Q-DEIM (QR with Column Pivoting)** algorithm!
