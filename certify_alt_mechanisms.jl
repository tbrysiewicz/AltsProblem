include("CertifyLower.jl")

"""
    certify_alt_mechanisms(; kwargs...)

Compatibility wrapper matching the paper's naming. This file is intentionally
named `certify_alt_mechanisms.jl` so the certification of alternative
mechanisms can be referenced by that exact name.
"""
certify_alt_mechanisms(; kwargs...) = certify_lower(; kwargs...)

if Base.abspath(PROGRAM_FILE) == @__FILE__
    certify_alt_mechanisms()
end
