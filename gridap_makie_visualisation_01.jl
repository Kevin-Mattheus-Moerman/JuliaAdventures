# This code is based on the gridap hyperelasticity demo. Here I expanded it to 3D and added Makie based model visualisation. 

# Note this code currently requires: ] add Makie@0.15.2 GLMakie@0.4.6

using Gridap 
using FileIO
using LineSearches: BackTracking
using GLMakie, GeometryBasics

# Material parameters
const λ = 100.0
const μ = 1.0

# Deformation Gradient
F(∇u) = one(∇u) + ∇u'

J(F) = sqrt(det(C(F)))

#Green strain

#E(F) = 0.5*( F'*F - one(F) )

dE(∇du,∇u) = 0.5*( ∇du⋅F(∇u) + (∇du⋅F(∇u))' )

# Right Cauchy-green deformation tensor

C(F) = (F')⋅F

# Constitutive law (Neo hookean)

function S(∇u)
  Cinv = inv(C(F(∇u)))
  μ*(one(∇u)-Cinv) + λ*log(J(F(∇u)))*Cinv
end

function dS(∇du,∇u)
  Cinv = inv(C(F(∇u)))
  _dE = dE(∇du,∇u)
  λ*(Cinv⊙_dE)*Cinv + 2*(μ-λ*log(J(F(∇u))))*Cinv⋅_dE⋅(Cinv')
end

# Cauchy stress tensor

σ(∇u) = (1.0/J(F(∇u)))*F(∇u)⋅S(∇u)⋅(F(∇u))'

# Model
domain = (0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition)

# Define new boundaries
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"diri_0",[1,3,5,7,13,15,17,19,25])
add_tag_from_tags!(labels,"diri_1",[2,4,6,8,14,16,18,20,26])

# Setup integration
degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

# Weak form

res(u,v) = ∫( (dE∘(∇(v),∇(u))) ⊙ (S∘∇(u)) )*dΩ

jac_mat(u,du,v) =  ∫( (dE∘(∇(v),∇(u))) ⊙ (dS∘(∇(du),∇(u))) )*dΩ

jac_geo(u,du,v) = ∫( ∇(v) ⊙ ( (S∘∇(u))⋅∇(du) ) )*dΩ

jac(u,du,v) = jac_mat(u,v,du) + jac_geo(u,v,du)

# Construct the FEspace
reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},1)
V = TestFESpace(model,reffe,conformity=:H1,dirichlet_tags = ["diri_0", "diri_1"])

# Setup non-linear solver
nls = NLSolver(
  show_trace=true,
  method=:newton,
  linesearch=BackTracking())

solver = FESolver(nls)

function run(x0,disp_x,step,nsteps,cache)

  g0 = VectorValue(0.0,0.0,0.0)
  g1 = VectorValue(disp_x,0.0,0.0)
  U = TrialFESpace(V,[g0,g1])

  #FE problem
  op = FEOperator(res,jac,U,V)

  println("\n+++ Solving for disp_x $disp_x in step $step of $nsteps +++\n")

  uh = FEFunction(U,x0)

  uh, cache = solve!(uh,solver,op,cache)

  #writevtk(Ω,"results_$(lpad(step,3,'0'))",cellfields=["uh"=>uh,"sigma"=>σ∘∇(uh)])

  return get_free_dof_values(uh), cache

end

function runs()

 disp_max = 0.75
 disp_inc = 0.02
 nsteps = ceil(Int,abs(disp_max)/disp_inc)

 x0 = zeros(Float64,num_free_dofs(V))

 cache = nothing
 for step in 1:nsteps
   disp_x = step * disp_max / nsteps
   x0, cache = run(x0,disp_x,step,nsteps,cache)
 end

end

#Do the work!
runs()


function getMakieMesh(model)
  #TO DO: Implement other element types, hex->quads only shown here

  #Get gridap element and node descriptions
  E=model.grid.cell_node_ids[:] #Elements
  V=model.grid.node_coords[:] #Nodes/Vertices
  
  #Get faces and convert to QuadFace type
  F=Vector{QuadFace{Int64}}(undef,size(E,1)*6)
  for q=1:1:size(E,1)    
    F[q]=convert(QuadFace{Int64},E[q][[1,2,4,3],1])             #top
    F[q+size(E,1)*1]=convert(QuadFace{Int64},E[q][[5,6,8,7],1]) #bottom
    F[q+size(E,1)*2]=convert(QuadFace{Int64},E[q][[1,2,6,5],1]) #side 1
    F[q+size(E,1)*3]=convert(QuadFace{Int64},E[q][[3,4,8,7],1]) #side 2
    F[q+size(E,1)*4]=convert(QuadFace{Int64},E[q][[2,4,8,6],1]) #front
    F[q+size(E,1)*5]=convert(QuadFace{Int64},E[q][[2,4,8,6],1]) #back
  end

  #Create face type labels
  faceTypeLabel=[ones(Int64,size(E,1)); 
                 ones(Int64,size(E,1))*2;
                 ones(Int64,size(E,1))*3;  
                 ones(Int64,size(E,1))*4;  
                 ones(Int64,size(E,1))*5;  
                 ones(Int64,size(E,1))*6;]

  #Convert gridap coordinate type to Makie Point3f type
  P=Vector{GeometryBasics.Point{3, Float32}}(undef,size(V,1))
  for q=1:1:size(V,1)
    P[q]=convert(Point3f,convert(Tuple{Float64, Float64, Float64},V[q]))
  end

  #Gather face and point set as GeometryBasics mesh
  M=GeometryBasics.Mesh(pointSet, faceSet)

  return M, faceTypeLabel
end

#Create makie compatible face and point set
M,faceTypeLabel=getMakieMesh(model)

#Visualize mesh
fig = Figure()
ax=Axis3(fig[1, 1], aspect = :data, xlabel = "x label", ylabel = "y label", zlabel = "z label", title = "Title")
poly!(M, strokewidth=2,shading=false,color=RGBAf0(0.5, 0.5, 1, 0.5), transparency=true, overdraw=false)
fig
