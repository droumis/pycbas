# Algorithm

This describes the CBAS procedure as implemented in pycbas, matching the original Igor implementation by David Kastner.

## Overview

1. Enumerate all behavioral subsequences up to length L for each subject
2. Build a count matrix (subjects x sequences)
3. Compute a one-sided test statistic per sequence per direction
4. Generate a bootstrap null distribution
5. Apply Romano-Wolf step-down to get adjusted p-values
6. Iterate k-FWER until the false discovery proportion converges

## Step 1. Sequence enumeration

For each subject, extract the choice stream (filtered by contingency block if applicable) and count all subsequences of lengths 1 through L starting within the first `criterion` positions (inclusive).

A sliding window at each position `i` from 0 to `min(criterion, stream_length - seq_len)` extracts the tuple `stream[i : i + seq_len]` and increments its count. The inclusive upper bound on start position matches Igor's `seqTrialWv[p][0] <= critWv` check.

## Step 2. Count matrix

Union all sequences observed across subjects. Sort by total frequency descending. Build the matrix C of shape (N, S) where `C[n, s]` is the count of sequence s for subject n.

## Step 3. Test statistics

For each sequence s, compute two one-sided studentized statistics:

```
t_pos = (mean_grp0 - mean_grp1) / sqrt(sem0^2 + sem1^2)    if mean_grp0 > mean_grp1
t_neg = (mean_grp1 - mean_grp0) / sqrt(sem0^2 + sem1^2)    if mean_grp1 > mean_grp0
```

where `sem = std(ddof=1) / sqrt(n_group)`. This doubles the hypothesis count to 2S, with one entry per pair set to NaN (the unobserved direction). This is Type III error handling from the paper.

For correlative mode, the test statistic is the studentized Pearson correlation using the robust tau estimator from DiCiccio & Romano (2017).

## Step 4. Bootstrap null

Generate M resamples by drawing subjects with replacement from the pooled sample, ignoring group labels. For each resample m and sequence s:

1. Compute the bootstrap statistic the same way as the observed (studentized difference)
2. Store the magnitude `|t*|` and which direction it went (positive=0, negative=1)

With `centering=False` (default, matches Igor), the raw bootstrap delta is used. With `centering=True`, the observed delta is subtracted before dividing by sigma.

For correlative mode, the covariate is permuted rather than resampling subjects.

The result is a null matrix of shape (M, S) containing magnitudes, and a direction matrix of shape (M, S) containing int8 direction indicators.

## Step 5. Romano-Wolf step-down

This procedure computes adjusted p-values that control the family-wise error rate while gaining power from the step-down structure.

Sort the observed test statistics descending. Starting with all sequences active:

1. For each bootstrap row, compute the k-th largest value among active columns
2. The p-value for the current (largest remaining) statistic is the fraction of rows where the k-th largest exceeds it: `p = (count + 1) / (M + 1)`
3. Enforce monotonicity: if `p < p_previous`, set `p = p_previous`
4. Remove the current sequence from the active set
5. Repeat for the next largest statistic

**Direction-conditional removal.** When removing a sequence from a bootstrap row, only actually remove it if the bootstrap went the same direction as the observed statistic. If the bootstrap went the opposite direction, that magnitude stays in the active set for that row. This matches Igor's behavior and makes the procedure slightly more conservative.

**Optimization.** After each removal step, only recompute the k-th largest for rows where (a) the direction matched and the sequence was removed, AND (b) the removed value was at least as large as the previous k-th largest for that row. Rows where the removed value couldn't have been in the top-k are skipped. This gives a 10-30x speedup depending on the dataset.

## Step 6. k-FWER iteration

The standard step-down (k=1) uses the row maximum, controlling FWER. To control FDP instead, iterate k:

1. Start with k=1
2. Run the step-down, count rejections R (sequences with p < alpha)
3. Check convergence: if `R < k/gamma - 1`, stop
4. Otherwise set `k = ceil((R + 1) * gamma)` and repeat

After convergence, run a final full step-down at the converged k to get the reported g-values. The converged k represents the number of false positives you're willing to tolerate among the rejections.

## Significance

A sequence is significant if its g-value (adjusted p-value from the converged step-down) is below alpha (default 0.5, which controls median FDP at level gamma=0.05).

## Chunked pipeline

The chunked variant generates bootstrap statistics directly into the sorted null submatrix in row-chunks, rather than computing the full M x S null matrix first and then extracting valid columns. This avoids the intermediate allocation and reduces peak memory by about 40%. The statistical output is identical.
