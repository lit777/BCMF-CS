# Joint Confounder Selection for Causal Mediation Analysis in High Dimensions

Code and data for the NeurIPS 2026 paper by Chanmin Kim and Hyunwoo Kim (Department of Statistics, Sungkyunkwan University). The proposed method is BCMF-CS: Bayesian causal mediation forests with joint confounder selection. Lab homepage: https://bci.creatorlink.net/

## Simulations Scenarios 1-7

To run each simulation scenario, run test_XXXX.R for BCMFCS, BCMF, DCMF, LSEM in each scenario folder:

```scenario
R CMD BATCH test_BCMFCS.R
R CMD BATCH test_BCMF.R
R CMD BATCH test_DCMF.R
R CMD BATCH test_LSEM.R
```

## Results

To summarize the outcomes, run this command:

```summary
R CMD BATCH summary.R
```


## Sensitivity Analysis

For the sensitivity analysis, run this command for each case:

```summary
R CMD BATCH test_sim_50.R
R CMD BATCH test_sim_100.R
```
Do the same step for other additional scenarios

## Real Application

To run the proposed model for the real data, run this command in the Application folder:

```summary
R CMD BATCH main_application.R
```
