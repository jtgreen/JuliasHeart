# using CUDA
# using DifferentialEquations, Sundials
# using Logging: global_logger
# using TerminalLoggers: TerminalLogger
# global_logger(TerminalLogger())
# using Debugger
using LinearAlgebra, SparseArrays

#device!(1)

#####################################################################################################################
### Simple matrix inversion, native solution; small matrix
# const v0 = 0.0f0
# const D = 0.01f0
# const v1 = 100.0f0

# const h = 10#16*12;
# const dh = 1/(h-1)
# u0 = fill(v0, (h, h));
# u0[3:6,3:6] .= v1;   # a small square in the middle of the domain

# dt = 0.001
# a = 1+(D*2*dt/dh^2)
# c = -dt/(D*2*dh^2)

# A1 = vec(fill(a, (1,h-1)))
# A2 = vec(fill(c, (1,h)))
# C = vec(fill(c, (1,h)))
# diag = Tridiagonal(A1,A2,A1)
# upper_lower = Diagonal(C)

# # sparse(
# # put together the larger matrix
# μ = [zeros(Float32, h, h) for i in 1:h, j in 1:h];

# μ[diagind(μ)] .= [diag];

# μ[diagind(μ, -1)] .= [upper_lower];

# μ[diagind(μ, 1)] .= [upper_lower'];

# bM = reduce(vcat, [reduce(hcat, μ[i, :]) for i in 1:h]);
# bM[1,2] = bM[1,2]*2
# bM[end,end-1] = bM[end,end-1]*2

# du = []
# push!(du, reshape(vec(u0) \ bM, (h,h)))
# for i in range(1,100;step=dt)
#     push!(du, reshape(vec(du[end]) \ bM, (h,h)))
# end

#############################################################################################################


#####################################################################################################################
### Simple matrix inversion, native solution; large matrix, SPARSE
const v0 = 0.0f0
const D = 0.01f0 # diffusion coefficient
const v1 = 100.0f0

const h = 202#16*12;
const dh = 1/(h-1)
u0 = fill(v0, (h, h));
u0[95:105,95:105] .= v1;   # a small square in the middle of the domain

dt = 0.001

a = 1+(D*2*dt/dh^2)
c = -D*dt/(2*dh^2)
a_ = 1-(D*2*dt/dh^2)
c_ = D*dt/(2*dh^2)

A1 = vec(fill(a, (1,h-1)))
A2 = vec(fill(c, (1,h)))
C = vec(fill(c, (1,h)))
diag = Tridiagonal(A1,A2,A1)
upper_lower = Diagonal(C)

# sparse(
# put together the larger matrix
μ = [sparse(zeros( h, h)) for i in 1:h, j in 1:h];

μ[diagind(μ)] .= [diag];

μ[diagind(μ, -1)] .= [upper_lower];

μ[diagind(μ, 1)] .= [upper_lower'];

bM = reduce(vcat, [reduce(hcat, μ[i, :]) for i in 1:h]);
# bM[1,2] = bM[1,2]*2
# bM[end,end-1] = bM[end,end-1]*2
bM[2,1] = 0#2.0
bM[end-1,end] = 0#2.0
bM[1,2] = 0#2.0
bM[end,end-1] = 0#2.0

#b
#bM = kron(My, sparse(I, h, h)) + kron(sparse(I, h, h), Mx)
# function rhs(u)
#     u = reshape(u, (h,h))
#     u[1,:] .= 0
#     u[end,:] .= 0
#     u[:,1] .= 0
#     u[:,end] .= 0
#     u = vec(u)
#     u[1] = 0
#     u[end] = 0
#     u[h] = 0
#     u[end-h+1] = 0
#     return u
# end

function rhs(u, a, c)
    n1, n2 = size(u)
    du = zeros(Float64, n1, n2)

    # internal nodes
    for j = 2:n2-1
        for i = 2:n1-1
            @inbounds  du[i,j] = a * u[i,j] + c * ( u[i+1,j] + u[i-1,j] + u[i,j+1] + u[i,j-1] )
        end
    end

    # left/right edges
    for i = 2:n1-1
        @inbounds du[i,1] = a * u[i,1] + c * ( u[i+1,1] + u[i-1,1] + u[i,2] )
        @inbounds du[i,n2] = a * u[i,n2] + c * ( u[i+1,n2] + u[i-1,n2] + u[i,n2-1] )
    end

    # top/bottom edges
    for j = 2:n2-1
        @inbounds du[1,j] = a * u[1,j] + c * ( u[1,j+1] + u[1,j-1] + u[2,j] )
        @inbounds du[n1,j] = a * u[n1,j] + c * ( u[n1,j+1] + u[n1,j-1] + u[n1-1,j] )
    end

    # corners
    @inbounds du[1,1] = a * u[1,1] + c * ( (u[2,1] + u[1,2]) )
    @inbounds du[n1,1] = a * u[n1,1] + c * ( u[n1-1,1] + u[n1,2] ) 
    @inbounds du[1,n2] = a * u[1,n2] + c * ( u[2,n2] + u[1,n2-1] ) 
    @inbounds du[n1,n2] = a * u[n1,n2] + c * ( u[n1-1,n2] + u[n1,n2-1] ) 

    return du
end

du = []
push!(du, u0)
push!(du, reshape(vec(rhs(u0, a_, c_)) \ bM, (h,h)))
for i in range(1,100;step=dt)

    push!(du, reshape(vec(rhs(du[end], a_, c_)) \ bM, (h,h)))
end

#############################################################################################################

######
# High res plot
# using Plots

# ENV["GKSwstype"]="nul"

# n = size(du, 1)
# anim = @animate for t ∈ collect(range(2,n,step=2)) # 1:n#
#     heatmap(du[t], c = cgrad(:rainbow, scale=(0,100)), clim=(0,100));
# end

# mp4(anim, "./heatmap.mp4", fps=15)

#########################################################################################################################################
## Low res plot

using UnicodePlots
# heatmap(sol.u[end], colormap_lim=(-180,90))

function move_up(s::AbstractString)
    move_up_n_lines(n) = "\u1b[$(n)F"
    # actually string_height - 1, but we're assuming cursor is on the last line
    string_height = length(collect(eachmatch(r"\n", s)))
    print(move_up_n_lines(string_height))
    nothing
end

function animate(frames; frame_delay = 0)
    print("\u001B[?25l") # hide cursor
    for frame in frames[1:end-1]
        print(frame)
        sleep(frame_delay)
        move_up(string(frame))
    end
    print(frames[end])
    print("\u001B[?25h") # visible cursor
    nothing
end

# frames = [heatmap(sol.u[i], colorbar_lim=(-180,90), ylabel=string(i))  for i in range(1,length(sol.u))]
using ThreadsX
#frames = ThreadsX.collect([heatmap(sol.u[i], colorbar_lim=(-180,90), ylabel=string(i))  for i in range(1,length(sol.u))]);
frames = ThreadsX.collect([heatmap(du[i], colorbar_lim=(0,100), ylabel=string(i))  for i in range(1,length(du); step=80)]);

animate(frames; frame_delay = 0)

#########################################################################################################################################
