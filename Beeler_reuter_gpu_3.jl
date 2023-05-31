using CUDA
using DifferentialEquations, Sundials
using LinearAlgebra, SparseArrays
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using Debugger

device!(1)

const v0 = -84.624
const v1 = 10.0
const C_K1 = 1.0f0
const C_x1 = 1.0f0
const C_Na = 1.0f0
const C_s = 1.0f0
const D_Ca = 0.0f0
const D_Na = 0.0f0
const g_s = 0.09f0
const g_Na = 4.0f0
const g_NaC = 0.005f0
const ENa = 50.0f0 + D_Na
const γ = 0.5f0
const C_m = 1.0f0

mutable struct BeelerReuterGpu <: Function
    t::Float64                  # the last timestep time to calculate Δt
    diff_coef::Float64          # the diffusion-coefficient (coupling strength)

    d_C::CuArray{Float32, 2}    # intracellular calcium concentration
    d_M::CuArray{Float32, 2}    # sodium current activation gate (m)
    d_H::CuArray{Float32, 2}    # sodium current inactivation gate (h)
    d_J::CuArray{Float32, 2}    # sodium current slow inactivaiton gate (j)
    d_D::CuArray{Float32, 2}    # calcium current activaiton gate (d)
    d_F::CuArray{Float32, 2}    # calcium current inactivation gate (f)
    d_XI::CuArray{Float32, 2}   # inward-rectifying potassium current (iK1)

    d_u::CuArray{Float64, 2}    # place-holder for u in the device memory
    d_du::CuArray{Float64, 2}   # place-holder for d_u in the device memory
    
    Δv::Array{Float64, 2}       # place-holder for voltage gradient

    # crank nicolson rhs
    a_rhs::Float64
    c_rhs::Float64

    # block tridiagonal matrix
    bM::SparseMatrixCSC{Float64, Int64}

    function BeelerReuterGpu(u0, diff_coef, dt)
        self = new()

        ny, nx = size(u0)

        @assert (nx % 16 == 0) && (ny % 16 == 0)
        self.t = 0.0
        self.diff_coef = diff_coef

        self.d_C = CuArray(fill(0.0001f0, (ny,nx)))
        self.d_M = CuArray(fill(0.01f0, (ny,nx)))
        self.d_H = CuArray(fill(0.988f0, (ny,nx)))
        self.d_J = CuArray(fill(0.975f0, (ny,nx)))
        self.d_D = CuArray(fill(0.003f0, (ny,nx)))
        self.d_F = CuArray(fill(0.994f0, (ny,nx)))
        self.d_XI = CuArray(fill(0.0001f0, (ny,nx)))

        self.d_u = CuArray(u0)
        self.d_du = CuArray(zeros(ny,nx))

        self.Δv = zeros(ny,nx)

        # 5-point crank nicolson stencil
        dh = 1/(nx-1)

        a = 1+( diff_coef *2*dt/dh^2)
        c = -diff_coef *dt/(2*dh^2)

        self.a_rhs = 1-( diff_coef *2*dt/dh^2)
        self.c_rhs = diff_coef *dt/(2*dh^2)

        A1 = vec(fill(a, (1,nx-1)))
        A2 = vec(fill(c, (1,nx)))
        C = vec(fill(c, (1,nx)))
        diag = Tridiagonal(A1,A2,A1)
        upper_lower = Diagonal(C)

        # put together the larger matrix for the block tridiagonal
        μ = [sparse(zeros( nx, nx )) for i in 1:nx, j in 1:nx];
        μ[diagind(μ)] .= [diag];
        μ[diagind(μ, -1)] .= [upper_lower];
        μ[diagind(μ, 1)] .= [upper_lower'];
        bM = reduce(vcat, [reduce(hcat, μ[i, :]) for i in 1:nx]);

        # Boundary conditions, neumann?
        bM[2,1] = 2.0
        bM[end-1,end] = 2.0
        bM[1,2] = 2.0
        bM[end,end-1] = 2.0

        self.bM = bM

        return self
    end
end

function crank_nicolson(Δv, u, a, c, bM)
    n1, n2 = size(u)

    # RHS
    # internal nodes
    for j = 2:n2-1
        for i = 2:n1-1
            @inbounds  Δv[i,j] = a * u[i,j] + c * ( u[i+1,j] + u[i-1,j] + u[i,j+1] + u[i,j-1] )
        end
    end

    # left/right edges
    for i = 2:n1-1
        @inbounds Δv[i,1] = a * u[i,1] + c * ( u[i+1,1] + u[i-1,1] + u[i,2] )
        @inbounds Δv[i,n2] = a * u[i,n2] + c * ( u[i+1,n2] + u[i-1,n2] + u[i,n2-1] )
    end

    # top/bottom edges
    for j = 2:n2-1
        @inbounds Δv[1,j] = a * u[1,j] + c * ( u[1,j+1] + u[1,j-1] + u[2,j] )
        @inbounds Δv[n1,j] = a * u[n1,j] + c * ( u[n1,j+1] + u[n1,j-1] + u[n1-1,j] )
    end

    # corners
    @inbounds Δv[1,1] = a * u[1,1] + c * ( (u[2,1] + u[1,2]) )
    @inbounds Δv[n1,1] = a * u[n1,1] + c * ( u[n1-1,1] + u[n1,2] ) 
    @inbounds Δv[1,n2] = a * u[1,n2] + c * ( u[2,n2] + u[1,n2-1] ) 
    @inbounds Δv[n1,n2] = a * u[n1,n2] + c * ( u[n1-1,n2] + u[n1,n2-1] ) 
    
    # LHS
    return reshape(vec( Δv ) \ bM, (n1,n2)) #.- u
end


function rush_larsen_gpu(g, α, β, Δt)
    inf = α/(α+β)
    τ = 1.0/(α+β)
    return clamp(g + (g - inf) * CUDA.expm1(-Δt/τ), 0f0, 1f0)
end

function update_M_gpu(g, v, Δt)
    # the condition is needed here to prevent NaN when v == 47.0
    α = CUDA.isapprox(v, -47.0f0, atol=1e-9) ? 10.0f0 : -(v+47.0f0) / (CUDA.exp(-0.1f0*(v+47.0f0)) - 1.0f0)
    β = (40.0f0 * CUDA.exp(-0.056f0*(v+72.0f0)))
    return rush_larsen_gpu(g, α, β, Δt)
end

function update_H_gpu(g, v, Δt)
    α = 0.126f0 * CUDA.exp(-0.25f0*(v+77.0f0))
    β = 1.7f0 / (CUDA.exp(-0.082f0*(v+22.5f0)) + 1.0f0)
    return rush_larsen_gpu(g, α, β, Δt)
end

function update_J_gpu(g, v, Δt)
    α = (0.55f0 * CUDA.exp(-0.25f0*(v+78.0f0))) / (CUDA.exp(-0.2f0*(v+78.0f0)) + 1.0f0)
    β = 0.3f0 / (CUDA.exp(-0.1f0*(v+32.0f0)) + 1.0f0)
    return rush_larsen_gpu(g, α, β, Δt)
end

function update_D_gpu(g, v, Δt)
    α = γ * (0.095f0 * CUDA.exp(-0.01f0*(v-5.0f0))) / (CUDA.exp(-0.072f0*(v-5.0f0)) + 1.0f0)
    β = γ * (0.07f0 * CUDA.exp(-0.017f0*(v+44.0f0))) / (CUDA.exp(0.05f0*(v+44.0f0)) + 1.0f0)
    return rush_larsen_gpu(g, α, β, Δt)
end

function update_F_gpu(g, v, Δt)
    α = γ * (0.012f0 * CUDA.exp(-0.008f0*(v+28.0f0))) / (CUDA.exp(0.15f0*(v+28.0f0)) + 1.0f0)
    β = γ * (0.0065f0 * CUDA.exp(-0.02f0*(v+30.0f0))) / (CUDA.exp(-0.2f0*(v+30.0f0)) + 1.0f0)
    return rush_larsen_gpu(g, α, β, Δt)
end

function update_XI_gpu(g, v, Δt)
    α = (0.0005f0 * CUDA.exp(0.083f0*(v+50.0f0))) / (CUDA.exp(0.057f0*(v+50.0f0)) + 1.0f0)
    β = (0.0013f0 * CUDA.exp(-0.06f0*(v+20.0f0))) / (CUDA.exp(-0.04f0*(v+20.0f0)) + 1.0f0)
    return rush_larsen_gpu(g, α, β, Δt)
end

function update_C_gpu(c, d, f, v, Δt)
    ECa = D_Ca - 82.3f0 - 13.0278f0 * CUDA.log(c)
    kCa = C_s * g_s * d * f
    iCa = kCa * (v - ECa)
    inf = 1.0f-7 * (0.07f0 - c)
    τ = 1f0 / 0.07f0
    return c + (c - inf) * CUDA.expm1(-Δt/τ)
end

# iK1 is the inward-rectifying potassium current
function calc_iK1(v)
    ea = CUDA.exp(0.04f0*(v+85f0))
    eb = CUDA.exp(0.08f0*(v+53f0))
    ec = CUDA.exp(0.04f0*(v+53f0))
    ed = CUDA.exp(-0.04f0*(v+23f0))
    return 0.35f0 * (4f0*(ea-1f0)/(eb + ec)
            + 0.2f0 * (CUDA.isapprox(v, -23f0, atol=1e-9) ? 25f0 : (v+23f0) / (1f0-ed)))
end

# ix1 is the time-independent background potassium current
function calc_ix1(v, xi)
    ea = CUDA.exp(0.04f0*(v+77f0))
    eb = CUDA.exp(0.04f0*(v+35f0))
    return xi * 0.8f0 * (ea-1f0) / eb
end

# iNa is the sodium current (similar to the classic Hodgkin-Huxley model)
function calc_iNa(v, m, h, j)
    return C_Na * (g_Na * m^3 * h * j + g_NaC) * (v - ENa)
end

# iCa is the calcium current
function calc_iCa(v, d, f, c)
    ECa = D_Ca - 82.3f0 - 13.0278f0 * CUDA.log(c)    # ECa is the calcium reversal potential
    return C_s * g_s * d * f * (v - ECa)
end

function update_gates_gpu(u, XI, M, H, J, D, F, C, Δt)
    i = (blockIdx().x-UInt32(1)) * blockDim().x + threadIdx().x
    j = (blockIdx().y-UInt32(1)) * blockDim().y + threadIdx().y

    v = Float32(u[i,j])

    let Δt = Float32(Δt)
        XI[i,j] = update_XI_gpu(XI[i,j], v, Δt)
        M[i,j] = update_M_gpu(M[i,j], v, Δt)
        H[i,j] = update_H_gpu(H[i,j], v, Δt)
        J[i,j] = update_J_gpu(J[i,j], v, Δt)
        D[i,j] = update_D_gpu(D[i,j], v, Δt)
        F[i,j] = update_F_gpu(F[i,j], v, Δt)

        C[i,j] = update_C_gpu(C[i,j], D[i,j], F[i,j], v, Δt)
    end
    nothing
end

function update_du_gpu(du, u, XI, M, H, J, D, F, C)#, Istim)
    i = (blockIdx().x-UInt32(1)) * blockDim().x + threadIdx().x
    j = (blockIdx().y-UInt32(1)) * blockDim().y + threadIdx().y

    v = Float32(u[i,j])

    # calculating individual currents
    iK1 = calc_iK1(v)
    ix1 = calc_ix1(v, XI[i,j])
    iNa = calc_iNa(v, M[i,j], H[i,j], J[i,j])
    iCa = calc_iCa(v, D[i,j], F[i,j], C[i,j])

    # total current
    I_sum = iK1 + ix1 + iNa + iCa

    # the reaction part of the reaction-diffusion equation
    du[i,j] = -I_sum / C_m #+ Istim[i,j]
    nothing
end

function (f::BeelerReuterGpu)(du, u, p, t)
    Δt = t - f.t

    copyto!(f.d_u, u)
    ny, nx = size(u)

    if Δt != 0 || t == 0
        @cuda blocks=(12,12) threads=(16,16) update_gates_gpu(
            f.d_u, f.d_XI, f.d_M, f.d_H, f.d_J, f.d_D, f.d_F, f.d_C, Δt)
        f.t = t
    end

    # calculate the reaction portion
    @cuda blocks=(12,12) threads=(16,16) update_du_gpu(
        f.d_du, f.d_u, f.d_XI, f.d_M, f.d_H, f.d_J, f.d_D, f.d_F, f.d_C)#, f.d_Istim)

    copyto!(du, f.d_du)

    crank_nicolson(f.Δv, u, f.a_rhs, f.c_rhs, f.bM) #, f.d_rhs, f.a_lhs, f.c_lhs, f.d_lhs, f.d_C, Δt, f.Δx, f.Δy, nx, ny)

    # f.diff_coef
    du .+=  f.Δv .+ p  
    nothing
end

### S2

const tstop1 = [1150. ]#, 300.]
const tstop2 = [1155. ]#, 301.]


# starts
function condition(u,t,integrator)
    t in tstop1
end

# stops
function condition2(u,t,integrator)
    t in tstop2
end

function affect!(integrator)
    # integrator.u[end:72,end:100] .+= -53
    # println("affect!")
    # println(integrator.u[72:end,100:end])
    
    #integrator.p = zeros(N,N);
    integrator.p[1:100,100:end] .+= 43.0
    
    # integrator.p[72:end,100:end] .= integrator.p[72:end,100:end]
    #integrator.p[72:end,100:end] .+= -
end

function affect2!(integrator)
    integrator.p .= 0
end

save_positions = (true,true)
cb = DiscreteCallback(condition, affect!, save_positions=save_positions)
save_positions = (true,true)
cb2 = DiscreteCallback(condition2, affect2!, save_positions=save_positions)
cbs = CallbackSet(cb ,cb2)

# solve

const N = 16*12;
u0 = fill(v0, (N, N));
u0[1:2,1:N] .= v1;   # a small square in the middle of the domain
p = zeros(N,N);     # the external stimulus, 0 initially

deriv_gpu = BeelerReuterGpu(u0, 0.1, 0.001);
prob = ODEProblem(deriv_gpu, u0, (0.0, 4000.0), p);

i = init(prob, CVODE_BDF(linear_solver=:GMRES), adaptive=false, dt=0.001)

# @time
# sol = solve(prob, CVODE_BDF(linear_solver=:GMRES), adaptive=false, dt=0.001, maxiters=Int(1e15), # callback=cbs, tstops=[tstop1[1],tstop2[1]],
#     saveat=1, progress = true, progress_steps = 100); 
# @ diff_coef=0.1
# results: LARGE spiral wave, no fast breakup
# Time: 0:01:39


# using Plots

# ENV["GKSwstype"]="nul"

# n = size(sol.u, 1)
# anim = @animate for t ∈ collect(range(1,n,step=10)) # 1:n#
#     heatmap(sol.u[t], c = cgrad(:rainbow, scale=(-100,30)), clim=(-100,20));
# end

# mp4(anim, "./heatmap.mp4", fps=15)

# using UnicodePlots
# # heatmap(sol.u[end], colormap_lim=(-180,90))

# function move_up(s::AbstractString)
#     move_up_n_lines(n) = "\u1b[$(n)F"
#     # actually string_height - 1, but we're assuming cursor is on the last line
#     string_height = length(collect(eachmatch(r"\n", s)))
#     print(move_up_n_lines(string_height))
#     nothing
# end

# function animate(frames; frame_delay = 0)
#     print("\u001B[?25l") # hide cursor
#     for frame in frames[1:end-1]
#         print(frame)
#         sleep(frame_delay)
#         move_up(string(frame))
#     end
#     print(frames[end])
#     print("\u001B[?25h") # visible cursor
#     nothing
# end

# # frames = [heatmap(sol.u[i], colorbar_lim=(-180,90), ylabel=string(i))  for i in range(1,length(sol.u))]
# using ThreadsX
# frames = ThreadsX.collect([heatmap(sol.u[i], colorbar_lim=(-180,90), ylabel=string(i))  for i in range(1,length(sol.u))]);

# animate(frames; frame_delay = 0)

# #frames2 = [densityplot(collect(range(1,192;step=1)),vec(sol.u[i][:,[80]]), ylim=(-100,80)) for i in range(1,length(sol.u))]
# frames2 = [lineplot(collect(range(1,192;step=1)),vec(sol.u[i][:,[80]]), canvas=BrailleCanvas, ylim=(-100,80)) for i in range(1,length(sol.u))]

# animate(frames2; frame_delay = 0.001)

# global sol = []
# global u = fill(v0, (N, N)
# global du = zeros(N,N)
# for i in range(1,100;step=0.001)
#     push!(sol, deriv_gpu(du, u, p, i))



# using DiffEqOperators
# hx = 1.0 / (1 .+ N)
# hy = 1.0 / (1 .+ N)
# Qx, Qy = Dirichlet0BC(Float64, size(u0));
# Dxx = CenteredDifference{1}(2,2,hx,N);
# Dyy = CenteredDifference{2}(2,2,hx,N);
# A = (Dxx + Dyy)*compose(Qx, Qy);
# @btime A*u0 # Error! can't figure it out, SLOW
# @btime laplacian(deriv_gpu.Δv, u0) # Fast!