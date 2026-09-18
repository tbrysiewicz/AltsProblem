include("CertifyLower.jl")

"""
    certify_alt_lower(; kwargs...)

Compatibility wrapper matching the paper's naming. This file is intentionally
named `certify_alt_lower.jl` so the lower-bound certification can be referenced
by that exact name.
"""
certify_alt_lower(; kwargs...) = certify_lower(; kwargs...)

if Base.abspath(PROGRAM_FILE) == @__FILE__
    certify_alt_lower()
end
