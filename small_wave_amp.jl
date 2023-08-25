using Printf
using LaTeXStrings
using Plots

const g = 9.81  # gravity acceleration m/s²
const ρ = 1.013e3 # kg/m³

## wave conditions
A = 0.2 # amplitude
L = 5.0 # wavelength: deep
#L = 20.0 # wavelength: shallow
#L = 60.0 # wavelength: very shallow (long wave)
k = 2π/L # wavenumber
h = 3.0 # depth (constant)
kh = k*h
ω = sqrt(g*k*tanh(k*h)) # angular frequency
T = 2π/ω # wave period
c = ω/k # celerity

## horizontal
Δx = 0.25
x = collect(Float64, -10.0:Δx:10.0)
nx = length(x)
## vertical
Δy = 0.2
y = collect(Float64, -h:Δy:0.0)
ny = length(y)
## for plotting
nxmid = ceil(Int64,0.50nx)
nymid = ceil(Int64,0.50ny)
ny3qt = ceil(Int64,0.75ny)


## time specification
nt = 21
t = LinRange(0,T,nt)
Δt = T/(nt-1)

## formulae
wx(x₀,y₀,t₀) = -A*cosh(k*(h+y₀))/sinh(k*h) * sin(k*x₀-ω*t₀)
wy(x₀,y₀,t₀) =  A*sinh(k*(h+y₀))/sinh(k*h) * cos(k*x₀-ω*t₀)
dp(x₀,y₀,t₀) = ρ*g*A*cosh(k*(h+y₀))/cosh(k*h) * sin(k*x₀-ω*t₀)
u(x₀,y₀,t₀) = A*ω*cosh(k*(h+y₀))/sinh(k*h) * cos(k*x₀-ω*t₀)
v(x₀,y₀,t₀) = A*ω*sinh(k*(h+y₀))/sinh(k*h) * sin(k*x₀-ω*t₀)

## particle location
# preallocate
X = zeros(nx,ny,nt)
Y = zeros(nx,ny,nt)
ΔP = zeros(nx,ny,nt)
U = zeros(nx,ny,nt)
for i = 1:nx
    for j = 1:ny
        for it = 1:nt
            X[i,j,it] = wx(x[i],y[j],t[it]) + x[i]
            Y[i,j,it] = wy(x[i],y[j],t[it]) + y[j]
            ΔP[i,j,it] = dp(x[i],y[j],t[it])
            U[i,j,it] = u(x[i],y[j],t[it])
        end
    end
end

umax = maximum(U[nxmid,end,:])

## for annotations
str1 = @sprintf("\$ h/L= %0.2f ~~~~ kh= %0.3f ~~~~ \\tanh~kh=%0.3f \$", h/L, kh, tanh(kh))
str2 = @sprintf("\$ T= %0.2f ~\\mathrm{s} ~~~~ c=%0.2f ~\\mathrm{m/s} ~~~~ u_{\\max} = %0.2f ~\\mathrm{m/s}\$", T, c, umax)


## functions
function snapshot(X,Y,it; kwargs...)
    #plt = scatter(X[:,:,it], Y[:,:,it], legend=false, zcolor=ΔP[:,:,it],
    plt = scatter(X[:,:,it], Y[:,:,it]; axis_ratio=:equal, legend=false, kwargs...)
    return plt
end
location_particle!(plt,x,y; kwargs...) = scatter!(plt, x, y; kwargs...)
trajectory_particle!(plt,xv,yv; kwargs...) = plot!(plt, xv, yv; kwargs...)


## surface
anim = @animate for it = 1:nt-1
    plt = plot(X[:,end,it], Y[:,end,it], ls=:solid, lw=4, axis_ratio=2, xlims=extrema(x), ylims=(-2A,2A), label=false, guide="", ticks=([]), margin=1*Plots.PlotMeasures.mm)
end
gifname = "surface.gif"
if isfile(gifname); rm(gifname); end
gif(anim, gifname, fps=5)


## whole domain
# preallocate
plts = Vector{Plots.Plot{Plots.GRBackend}}(undef, nt-1)
# plot
anim = @animate for it = 1:nt-1
    plts[it] = snapshot(X,Y,it; color=:cyan, xlims=extrema(x), ylims=(-h,2A), xlabel=L"x ~~\mathrm{(m)}", ylabel=L"z ~~\mathrm{(m)}")
    plts[it] = trajectory_particle!(plts[it], X[nxmid,  end,:], Y[nxmid,  end,:]; lc=:black)
    plts[it] = trajectory_particle!(plts[it], X[nxmid,ny3qt,:], Y[nxmid,ny3qt,:]; lc=:black)
    plts[it] = trajectory_particle!(plts[it], X[nxmid,nymid,:], Y[nxmid,nymid,:]; lc=:black)
    plts[it] = location_particle!(plts[it], [X[nxmid,  end,it]], [Y[nxmid,  end,it]]; color=:magenta)
    plts[it] = location_particle!(plts[it], [X[nxmid,ny3qt,it]], [Y[nxmid,ny3qt,it]]; color=:magenta)
    plts[it] = location_particle!(plts[it], [X[nxmid,nymid,it]], [Y[nxmid,nymid,it]]; color=:magenta)
end

## save
gifname = "deep_nozoom.gif"
if isfile(gifname); rm(gifname); end
gif(anim, gifname, fps=5)


