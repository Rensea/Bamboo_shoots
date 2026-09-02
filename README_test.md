# Best Models for Review

## Files
- best_model_<trait>.rds : full model bundle (model object + preprocessing + feature selection + train/test split IDs)
- coefficients_<trait>.csv : ridge regression coefficients

## Reproduce test-set R2 in R
```r
bundle <- readRDS('best_model_soluble_sugar.rds')
# bundle$test_ids contains the sample IDs used as held-out test set
# Preprocess data with apply_preprocess_pair() + build_feature_matrix() from main script
# then predict with bundle$model
```

## Traits: soluble_sugar, protein, moisture, lignin, cellulose
Generated: 2026-09-02 23:35:38.497521
