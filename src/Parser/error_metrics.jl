using Configurations

@option struct ErrorMetricsParams
    MAE::Bool = false
    eigenvalue_error::Bool = false
    el_entropy_error::Bool = false
end
