module SelectionTutorial

using CairoMakie
using Colors
using Cosmology
using DataFrames
using DimensionalData
using Distributions
using GaussianKDEs
using HDF5
using PolyLog
using ProgressLogging
using Unitful
using UnitfulAstro

export lvk_cosmology, categorical_palette, load_pe_from_dir, load_selection, chi_eff_marginal
export li_nocosmo_prior_logwt_m1qzchie, kdecontour, kdecontour!, traceplot, pop_logwt_m1qz
export traceplot

# From https://dcc.ligo.org/LIGO-P2000318/public/
const h_lvk = 0.679
const Ω_M_lvk = 0.3065
const lvk_cosmology = cosmology(h=h_lvk, OmegaM=Ω_M_lvk, w0=-1, wa=0)

"""
    categorical_palette(n; l=65, c=90)

Generate a palette of `n` colors cycling through the LCHuv color space with
lightness `l` and chroma `c`.  Roughly equivalent to the `HUSL` palette of
Seaborn.
"""
function categorical_palette(n; l=65, c=90)
    [LCHuv(l, c, h) for h in range(0, stop=360, length=n+1)[1:n]]
end

"""
    load_pe_from_dir(dir)

Load parameter estimation samples from all the `_nocosmo.h5` files in the given
directory.

Returns a DataFrame with all the samples, with an additional column `gwname`
giving the full (GWYYMMDD_NNNNNNN) name of the event to which the samples
belong.
"""
function load_pe_from_dir(dir)
    dfs = []
    @progress name="Loading Directory" for file in readdir(dir)
        if occursin(r"GW[0-9]+.*_nocosmo.h5", file)
            gwname = match(r".*(GW[0-9]+_[0-9]+).*", file).captures[1]
            h5open(joinpath(dir, file), "r") do f
                k = keys(f)
                if "C01:Mixed" in keys(f)
                    samps = read(f["C01:Mixed/posterior_samples"])

                    d = DataFrame(samps)
                    d[!, :gwname] .= gwname

                    push!(dfs, d)
                else
                    @info "Could not read $file"
                end
            end
        end
    end
    vcat(dfs...; cols=:intersect)
end

"""
    load_selection(file)

Given a path to a file containing mock injections, load it into a data frame.

Returns `(df, T_in_years, N)` where `df` is the data frame, `T_in_years` is the
analysis time in years over which the injections have been generated, and `N` is
the total number of injections generated.

The data frame will include generated columns `q` and `chi_eff` and
`sampling_pdf_qchieff` which gives the (marginal) sampling PDF for parameters
`mass1`, `q`, `chi_eff`, `z`.
"""
function load_selection(file)
    h5open(file, "r") do f
        df = DataFrame(read(f["injections"]))
        T_yr = attrs(f)["analysis_time_s"] / (3600 * 24 * 365.25)
        N = attrs(f)["total_generated"]

        df.mass_1 = @. df.mass1_source * (1 + df.redshift)
        df.q = @. df.mass2_source ./ df.mass1_source
        df.chi_eff = @. (df.spin1z + df.q * df.spin2z) / (1 + df.q)
        df.luminosity_distance = @. ustrip(u"Gpc", luminosity_dist((lvk_cosmology, ), df.redshift))

        a1 = @. sqrt(df.spin1x^2 + df.spin1y^2 + df.spin1z^2)
        a2 = @. sqrt(df.spin2x^2 + df.spin2y^2 + df.spin2z^2)

        spin_sampling_pdf = @. 1 / (16 * π^2 * a1^2 * a2^2)

        df.sampling_pdf_qchieff = @. df.sampling_pdf / spin_sampling_pdf * df.mass1_source * chi_eff_marginal(df.chi_eff, df.q)
        df.sampling_pdf_q = @. df.sampling_pdf * df.mass1_source

        dc = @. ustrip(u"Gpc", comoving_transverse_dist((lvk_cosmology, ), df.redshift))
        dh_z = @. ustrip(u"Gpc", 2.99792e8*u"m"/u"s" / Cosmology.H((lvk_cosmology, ), df.redshift))

        df.sampling_pdf_m1dqdlchieff = @. df.sampling_pdf_qchieff / (1 + df.redshift) / (dc + (1+df.redshift)*dh_z)

        (df, T_yr, N)
    end
end

"""
    chi_eff_marginal(chi_eff, q; amax=1)

Returns the marginal prior density on `chi_eff` conditional on `q` assuming an
isotropic, independent spin prior with flat priors on the component spins for
`a_i < amax`.

Taken from [Callister (2021)](https://arxiv.org/abs/2104.09508).
"""
function chi_eff_marginal(chi_eff, q; amax=1)
    abs_chi_eff = abs(chi_eff)
    l1 = amax*((1-q)/(1+q))
    l2 = q*amax/(1+q)
    l3 = amax/(1+q)

    if chi_eff == 0
        (1+q)/(2*amax)*(2 - log(q))
    elseif abs_chi_eff < l1 && abs_chi_eff < l2
        (1+q)/(4*q*amax^2)*( 
            q*amax*(4 + 2*log(amax) - log(q^2*amax^2 - (1+q)^2*abs_chi_eff^2)) 
            - 2*(1+q)*abs_chi_eff*atanh((1+q)*abs_chi_eff/(q*amax)) 
            + (1+q)*abs_chi_eff*(reli2(-q*amax/((1+q)*abs_chi_eff)) - reli2(q*amax/((1+q)*abs_chi_eff)))
        )
    elseif abs_chi_eff < l1 && abs_chi_eff > l2
        (1+q)/(4*q*amax^2)*(
            4*q*amax + 2*q*amax*log(amax)
            - 2*(1+q)*abs_chi_eff*atanh(q*amax/((1+q)*abs_chi_eff))
            - q*amax*log((1+q)^2*abs_chi_eff^2 - q^2*amax^2)
            + (1+q)*abs_chi_eff*(reli2(-q*amax/((1+q)*abs_chi_eff)) - reli2(q*amax/((1+q)*abs_chi_eff)))
        )
    elseif abs_chi_eff > l1 && abs_chi_eff < l2
        (1+q)/(4*q*amax^2)*(
            2*(1+q)*(amax - abs_chi_eff) - (1+q)*abs_chi_eff*log(amax)^2
            + (amax + (1+q)*abs_chi_eff*log((1+q)*abs_chi_eff))*log(q*amax/(amax - (1+q)*abs_chi_eff))
            - (1+q)*abs_chi_eff*log(amax)*(2 + log(q) - log(amax - (1+q)*abs_chi_eff))
            + q*amax*log(amax / (q*amax - (1+q)*abs_chi_eff))
            + (1+q)*abs_chi_eff*log((amax - (1+q)*abs_chi_eff)*(q*amax - (1+q)*abs_chi_eff)/q)
            + (1+q)*abs_chi_eff*(reli2(1-amax/((1+q)*abs_chi_eff)) - reli2(q*amax/((1+q)*abs_chi_eff)))
        )
    elseif abs_chi_eff > l1 && abs_chi_eff > l2 && abs_chi_eff < l3
        (1+q)/(4*q*amax^2)*(
            - abs_chi_eff*log(amax)^2 + 2*(1+q)*(amax - abs_chi_eff)
            + q*amax*log(amax/((1+q)*abs_chi_eff - q*amax)) + amax*log(q*amax/(amax - (1+q)*abs_chi_eff))
            - abs_chi_eff*log(amax)*(2*(1+q) - log((1+q)*abs_chi_eff) - q*log((1+q)*abs_chi_eff/amax))
            + (1+q)*abs_chi_eff*log((-q*amax + (1+q)*abs_chi_eff)*(amax - (1+q)*abs_chi_eff)/q)
            + (1+q)*abs_chi_eff*log(amax/((1+q)*abs_chi_eff))*log((amax - (1+q)*abs_chi_eff)/q)
            + (1+q)*abs_chi_eff*(reli2(1-amax/((1+q)*abs_chi_eff)) - reli2(q*amax/((1+q)*abs_chi_eff)))
        )
    elseif abs_chi_eff > l3 && abs_chi_eff < amax
        (1+q)/(4*q*amax^2)*(
            2*(1+q)*(amax - abs_chi_eff) - (1+q)*abs_chi_eff*log(amax)^2
            + log(amax)*(amax - 2*(1+q)*abs_chi_eff - (1+q)*abs_chi_eff*log(q/((1+q)*abs_chi_eff - amax)))
            - amax*log(((1+q)*abs_chi_eff - amax)/q)
            + (1+q)*abs_chi_eff*log(((1+q)*abs_chi_eff - amax)*((1+q)*abs_chi_eff - q*amax)/q)
            + (1+q)*abs_chi_eff*log((1+q)*abs_chi_eff)*log(q*amax / ((1+q)*abs_chi_eff - amax))
            - q*amax*log(((1+q)*abs_chi_eff - q*amax)/amax)
            + (1+q)*abs_chi_eff*(reli2(1 - amax/((1+q)*abs_chi_eff)) - reli2(q*amax/((1+q)*abs_chi_eff)))
        )
    else abs_chi_eff > amax
        zero(chi_eff)
    end
end

"""
    li_nocosmo_prior_logwt_m1qzchie(df)

Returns the LALInference "nocosmo" prior weight (i.e. un-normalized prior
density) over `m1`, `q`, `chi_eff`, `z` for each row in the given data frame.
"""
function li_nocosmo_prior_logwt_m1qzchie(df)
    z = df[!, :redshift]
    m1 = df[!, :mass_1_source]
    q = df[!, :mass_ratio]
    chi_eff = df[!, :chi_eff]

    log_opz = log1p.(z)

    dl = @. ustrip(u"Gpc", luminosity_dist((lvk_cosmology, ), z))
    dc = @. ustrip(u"Gpc", comoving_transverse_dist((lvk_cosmology, ), z))
    dh_z = @. ustrip(u"Gpc", 2.99792e8*u"m"/u"s" / Cosmology.H((lvk_cosmology, ), z))

    # According to [Callister (2021)](https://arxiv.org/abs/2104.09508), the
    # prior on m1_source, m2_source, z is (1+z)^2*dl^2*(dc + (1+z)*c/H(z))
    # the final log(m1) comes from d(m2_source)/d(q) = m1_source
    m1_m2_z_logwt = @. 2*log_opz + 2*log(dl) + log(dc + (1+z)*dh_z)
    m1_q_z_logwt = @. m1_m2_z_logwt + log(m1)

    # Also from [Callister (2021)](https://arxiv.org/abs/2104.09508)
    chi_eff_logwt = @. log(chi_eff_marginal(chi_eff, q))

    m1_q_z_logwt .+ chi_eff_logwt
end

"""
    pop_logwt_m1qz(df)

Returns the log of the population weight over m1, q, z for each row of the data
frame.

Our fiducial population is a power law with exponent -2.35 in m1, flat in q, and
follows a Madau-Dickinson merger rate over redshift.
"""
function pop_logwt_m1qz(df)
    alpha_m = -2.35

    z = df[!, :redshift]
    dVdz = @. ustrip(u"Gpc^3", comoving_volume_element((lvk_cosmology,), z))
    
    @. alpha_m * log(df[!, :mass_1_source]) + 2.7*log1p(z) - log1p(((1+z)/2.9)^5.6) + log(dVdz) - log1p(z)
end

@recipe(KDEContour, x, y) do scene
    Theme(
        levels = [0.1, 0.5],
        xexpansionfactor = 0.1,
        yexpansionfactor = 0.1,
        xgridsize=128,
        ygridsize=129
    )
end

function Makie.plot!(kdecontour::KDEContour)
    xexpansionfactor = pop!(kdecontour.attributes, :xexpansionfactor)
    yexpansionfactor = pop!(kdecontour.attributes, :yexpansionfactor)
    xgridsize = pop!(kdecontour.attributes, :xgridsize)
    ygridsize = pop!(kdecontour.attributes, :ygridsize)
    levels = pop!(kdecontour.attributes, :levels)

    x = kdecontour.x
    y = kdecontour.y

    xgrid = lift(x, xexpansionfactor, xgridsize) do x, xexpansionfactor, xgridsize
        dx = maximum(x)-minimum(x)
        fx = xexpansionfactor/2
        range(minimum(x)-fx*dx, maximum(x)+fx*dx, length=xgridsize)
    end

    ygrid = lift(y, yexpansionfactor, ygridsize) do y, yexpansionfactor, ygridsize
        dy = maximum(y)-minimum(y)
        fy = yexpansionfactor/2
        range(minimum(y)-fy*dy, maximum(y)+fy*dy, length=ygridsize)
    end

    k = lift(x, y) do x, y
        pts = vcat(x', y')
        KDE(pts)
    end

    p_levels = lift(x, y, k, levels) do x, y, k, levels
        pts = vcat(x', y')
        p_kde_pts = [pdf(k, pts[:,i]) for i in axes(pts, 2)]
        [quantile(p_kde_pts, l) for l in levels]
    end

    z = lift(x, y, k, xgrid, ygrid) do x, y, k, xgrid, ygrid
        [pdf(k, [x, y]) for x in xgrid, y in ygrid]
    end

    contour!(kdecontour, xgrid, ygrid, z; levels=p_levels, kdecontour.attributes...)
end

function traceplot(trace, params)
    f = Figure()

    # Nevt = length(dims(trace.posterior, :event))

    for (i, p) in enumerate(params)
        a_trace = Axis(f[i,1]; ylabel=string(p), xlabel="Iter")
        a_dens = Axis(f[i,2]; ylabel="Density", xlabel=string(p))

        for (ic, _) in enumerate(dims(trace.posterior, :chain))
            lines!(a_trace, trace.posterior[p][chain=ic])
            density!(a_dens, vec(trace.posterior[p][chain=ic]))
        end

        # if p == :Neff_sel
        #     hlines!(a_trace, 4*Nevt, color=:black)
        #     vlines!(a_dens, 4*Nevt, color=:black)
        # end
    end
    f
end

end # module SelectionTutorial
