# SHAPE FROM SHADING DEMO

# Based on: https://gist.github.com/Balaje/409276c9501e29dc9a79e9e3e41bfd4a
# Papers on this topic: 
# https://doi.org/10.1137/0729053
# https://doi.org/10.1109/ICCV.2003.1238433

# Call required packages
using LinearAlgebra
using GLMakie
using Images, Colors

# Visualization parameters
fontsize=24

# Build Discretization
N = 200;
x = LinRange(-1,1,N+1)
y = LinRange(-1,1,N+1)
Δx = x[2] - x[1]
Δy = y[2] - y[1]

# Load image and convert to Grayscale
imgRaw = load("/home/kevin/DATA/Julia/JuliaAdventures/moz.png");
imgResampled = imresize(imgRaw, 201, 201)
Imatrix = Gray.(imgResampled);

# Compute functions for Perspective model
I(i,j) = Imatrix[i,j]
Q(x,y) = √(f^2/(f^2 + x^2 + y^2))
f = 2.32

# Rouy and Tourin/Godunov Scheme (Return Hamiltonian)
function H(i,j,U)
  #Iᵢⱼ = I(x[i], y[j])
  Iᵢⱼ = I(i,j)
  Iᵢⱼ = (Iᵢⱼ ≈ 0) ? 0.01 : Iᵢⱼ
  Qᵢⱼ = Q(x[i], y[j])

  D₊ˣ = (U[i+1,j] - U[i,j])/Δx
  D₋ˣ = (U[i,j] - U[i-1,j])/Δx
  D₊ʸ = (U[i,j+1] - U[i,j])/Δy
  D₋ʸ = (U[i,j] - U[i,j-1])/Δy
  Dˣ = max(0, -D₊ˣ, D₋ˣ)
  Dʸ = max(0, -D₊ʸ, D₋ʸ)

  ∇u = (Dˣ^2 + Dʸ^2)
  ∇ux = [Dˣ,Dʸ]⋅[x[i],y[j]]

  H = -exp(-2*U[i,j]) + (Iᵢⱼ*f^2/Qᵢⱼ)*√(f^2*∇u + (∇ux)^2 + Qᵢⱼ^2)
end

# Build Arrays
U = zeros(N+1,N+1)
U₁ = zeros(N+1,N+1)
U[1,:] .= 100
U[N+1,:] .= 100
U[:,1] .= 100
U[:,N+1] .= 100
U₁ .= U

# Choose Δt
Δt = 0.1*Δx
τ = 1e-5

# Solve. Should not take more than a minute ...
while 1 > 0
  for i = 2:N
    for j = 2:N
      U₁[i,j] = U[i,j] - Δt*H(i,j,U)
    end
  end
  δ = abs.(U₁[2:N,2:N] - U[2:N,2:N])
  ϵ = norm(δ,Inf)
  if(ϵ < τ)
    @show ϵ, norm(U, Inf)
    break
  end
  U .= U₁
end

# Visualization
fig = Figure(fontsize=fontsize, resolution=(1000,1000))
ax1 = GLMakie.Axis(fig[1, 1],title = "Intensity Image", xlabel = "j", ylabel = "i", aspect=DataAspect());
GLMakie.image!(ax1,rotr90(imgRaw))

ax2 = GLMakie.Axis3(fig[1, 2],title = "Reconstructed Surface", xlabel = "x", ylabel = "y", zlabel = "z");
GLMakie.surface!(ax2, x[2:end-1], y[2:end-1], -U[2:end-1,2:end-1])

fig
