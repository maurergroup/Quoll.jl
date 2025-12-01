### OPERATOR CONVERSION INTERFACE ###

# Could be specialised for particular operator conversions if the conversion can be done in a more efficient way
# e.g. if metadata is always shared between operators
function convert_operators(in_operators, out_operator_type::Type{<:AbstractOperator};
    radii = nothing, hermitian = nothing, float = nothing)
    return convert_operator.(in_operators, out_operator_type; radii = radii, hermitian = hermitian, float = float)
end

# Requires defining appropriate constructor for T2.
# This function is an alias for the constructor, but also
# selects appropriate kwargs based on T1 at runtime
function convert_operator(in_operator::T₁, ::Type{T₂};
    radii = nothing, hermitian = nothing, float = nothing) where {T₁<:AbstractOperator, T₂<:AbstractOperator}

    isnothing(hermitian) && (hermitian = get_sparsity(in_operator).hermitian)
    isnothing(float) && (float = get_float(in_operator))

    return T₂(in_operator, radii = radii, hermitian = hermitian, float = float)
end

### FOURIER TRANSFORM INTERFACE ###

function fourier_transform(in_operator::T₁, kpoints, ::Type{T₂}) where {T₁<:AbstractOperator, T₂<:AbstractOperator}
    if !isa(kpoints, AbstractVector{<:AbstractVector})
        kpoints = [kpoints]
    end
    phases = eachcol(precompute_phases(kpoints, get_sparsity(in_operator).images))
    return fourier_transform(in_operator, kpoints, phases, T₂)
end

### SPARSITY CONVERSION INTERFACE ###

function convert_sparsity(in_metadata::AbstractOperatorMetadata, radii::SpeciesAnyDict, T::Type{<:AbstractSparsity}; hermitian = false)
    @debug "Computing sparsity manually from radii"
    return T(get_atoms(in_metadata), radii, hermitian = hermitian)
end

function convert_sparsity(in_metadata::AbstractOperatorMetadata, ::Nothing, out_sparsity_type::Type{<:AbstractSparsity}; hermitian = false)
    @debug "Converting sparsity"
    return convert_sparsity(in_metadata, out_sparsity_type, hermitian = hermitian)
end

# Works out of the box for in_sparsity_type == out_sparsity_type. Also would work for sparsity conversions
# which define convert_sparsity(::AbstractSparsity, ::Type{<:AbstractSparsity}; hermitian).
# However for some sparsity conversions in_sparsity is not enough to convert to out_sparsity.
# In that case this method needs to be specialised for the OperatorMetadata in question.
function convert_sparsity(in_metadata::AbstractOperatorMetadata, out_sparsity_type::Type{<:AbstractSparsity}; hermitian = false)
    in_sparsity = get_sparsity(in_metadata)
    return convert_sparsity(in_sparsity, out_sparsity_type, hermitian = hermitian)
end
