import numpy as np
from scipy.stats import pearsonr

def vstrans(nk):
        return np.var(np.arange(nk)) / nk**4

# stratified correlation coefficient used for HiC prediction evaluation
def scc(mat1: np.ndarray, mat2: np.ndarray, max_bins: int = None):
        assert mat1.shape == mat2.shape
        if max_bins is None:
                max_bins = int(mat1.shape[0] - 5)
        else:
                max_bins = min(max_bins, mat1.shape[0] - 5)
        corr_diag = np.zeros(max_bins)
        weight_diag = np.zeros_like(corr_diag)
        for d in range(max_bins):
                d1 = mat1.diagonal(d)
                d2 = mat2.diagonal(d)
                mask = ~np.isnan(d1) & ~np.isnan(d2)
                d1 = d1[mask]
                d2 = d2[mask]
                corr_diag[d] = pearsonr(d1, d2)[0] if np.var(d1) > 0 and np.var(d2) > 0 else np.nan
                r2k = vstrans(len(d1))
                weight_diag[d] = len(d1) * r2k
        corr_diag, weight_diag = corr_diag[1:], weight_diag[1:]
        mask = ~np.isnan(corr_diag)
        corr_diag, weight_diag = corr_diag[mask], weight_diag[mask]
        # Normalize weights
        weight_diag /= sum(weight_diag)
        # Weighted sum of coefficients to get SCCs
        scc = np.nansum(corr_diag * weight_diag)
        return scc