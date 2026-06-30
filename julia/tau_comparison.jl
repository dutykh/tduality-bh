# Authorship:
#   Davide Batic  (Mathematics Department, Khalifa University of Science and
#                  Technology, Abu Dhabi, UAE)
#   Denys Dutykh  (Mathematics Department, Khalifa University of Science and
#                  Technology, Abu Dhabi, UAE)
#
# tau_comparison.jl
# Comparison between the exact proper time along the PG radial congruence
# and the endpoint expression obtained from the JN integrated prescription.
# Output: tau_comparison.pdf
#
# It only requires Plots.jl (no dependency on LaTeXStrings).

using Plots

# Parameters used in the figure. In the dimensionless ratio r/ell0,
# changing M only rescales the vertical axis by ell0^(3/2)/sqrt(2M).
M = 5.0
ell0 = 1.0
ri = 4.0 * ell0

# Expansion scalar for the selected PG congruence.
theta(r) = -3.0 * sqrt(2.0*M) * (r^2 + 2.0*ell0^2) / (2.0 * (r^2 + ell0^2)^(7.0/4.0))

# Function appearing in d(1/Theta)/d tau = F(r(tau)).
Ffun(r) = r^2 * (3.0*r^2 + 10.0*ell0^2) / (3.0 * (r^2 + 2.0*ell0^2)^2)

# Antiderivative for int (r^2 + ell0^2)^(3/4)/r dr after the change
# u = (1 + r^2/ell0^2)^(1/4).
u_of_r(r) = (1.0 + (r/ell0)^2)^(1.0/4.0)
Ifun(u) = (2.0/3.0)*u^3 + atan(u) + 0.5*log((u - 1.0)/(u + 1.0))

# Exact PG proper time from ri to r.
function tau_exact(r)
    ui = u_of_r(ri)
    u = u_of_r(r)
    return ell0^(3.0/2.0) / sqrt(2.0*M) * (Ifun(ui) - Ifun(u))
end

# Endpoint expression corresponding to replacing int F(r(tau)) d tau by tau F(r).
theta_i = theta(ri)
tau_JN(r) = (1.0/theta(r) - 1.0/theta_i) / Ffun(r)

# Logarithmic sampling resolves the small-r region where the two asymptotics differ.
rmin = 1.0e-3 * ell0
rmax = 0.999 * ri
rs = exp.(range(log(rmin), log(rmax), length=900))

p = plot(
    rs ./ ell0,
    tau_exact.(rs),
    label = "tau_exact",
    lw = 2.2,
    color = :black,
    xscale = :log10,
    yscale = :log10,
    xlabel = "r/ell0",
    ylabel = "",
    framestyle = :box,
    legend = :topright,
    size = (520, 360),
    dpi = 300,
    minorgrid = false,
)

plot!(
    p,
    rs ./ ell0,
    tau_JN.(rs),
    label = "tau_JN",
    lw = 2.2,
    ls = :dash,
    color = :black,
)

savefig(p, "tau_comparison.pdf")
