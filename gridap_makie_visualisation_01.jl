# This code is based on the gridap hyperelasticity demo: https://gridap.github.io/Tutorials/dev/pages/t005_hyperelasticity/
# Here I expanded it to 3D and added Makie based model visualisation. 

# Note this code currently requires: ] add Makie@0.15.2 GLMakie@0.4.6

using Gridap 
using Gridap.Visualization
using Gridap.ReferenceFEs
using Gridap.Geometry
using FileIO
using LineSearches: BackTracking
using GLMakie, GeometryBasics

# Material parameters
const λ = 100.0
const μ = 1.0

# BC parameters
disp_max = 0.8
disp_inc = 0.1

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
nElem=8
domain = (0,1,0,1,0,1)
partition = (nElem,nElem,nElem)
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

  return get_free_dof_values(uh), uh, cache

end

function runs(disp_max,disp_inc)

 nsteps = ceil(Int,abs(disp_max)/disp_inc)

 x0 = zeros(Float64,num_free_dofs(V))
 nodalDisplacements = Vector{Vector{VectorValue{3, Float64}}}(undef,nsteps)

 cache = nothing
 for step in 1:nsteps
  disp_x = step * disp_max / nsteps
  x0, uh, cache = run(x0,disp_x,step,nsteps,cache)

  vd = visualization_data(Ω,"",cellfields=["u"=>uh])
  nodalDisplacements[step] = vd[1].nodaldata["u"]
 end
 
 return nodalDisplacements
end


function nodesToPointset(V)
  #Convert gridap coordinate type to Makie Point3f type
  P=Vector{GeometryBasics.Point{3, Float32}}(undef,size(V,1))
  for q=1:1:size(V,1)
    P[q]=convert(GeometryBasics.Point,convert(Tuple{Float64, Float64, Float64},V[q]))
  end
  return P
end

function convertToFacePointSet(Ω)
  #TO DO: Implement other element types, hex->quads only shown here

  #Get gridap element and node descriptions
  # E=Ω.cell_node_ids[:] #Elements
  # V=Ω.node_coords[:] #Nodes/Vertices
  vd=visualization_data(Ω,"");
  grid = vd[1].grid
  E = get_cell_node_ids(grid)
  V = get_node_coordinates(grid)

  #Get faces and convert to QuadFace type
  F=Vector{QuadFace{Int64}}(undef,size(E,1)*6)
  for q=1:1:size(E,1)    
    F[q]=convert(QuadFace{Int64},E[q][[1,2,4,3],1])             #top
    F[q+size(E,1)*1]=convert(QuadFace{Int64},E[q][[5,6,8,7],1]) #bottom
    F[q+size(E,1)*2]=convert(QuadFace{Int64},E[q][[1,2,6,5],1]) #side 1
    F[q+size(E,1)*3]=convert(QuadFace{Int64},E[q][[4,3,7,8],1]) #side 2
    F[q+size(E,1)*4]=convert(QuadFace{Int64},E[q][[2,4,8,6],1]) #front
    F[q+size(E,1)*5]=convert(QuadFace{Int64},E[q][[3,1,5,7],1]) #back
  end

  #Create face type labels
  faceTypeLabel=[ones(Int64,size(E,1))*1; 
                 ones(Int64,size(E,1))*2;
                 ones(Int64,size(E,1))*3;  
                 ones(Int64,size(E,1))*4;  
                 ones(Int64,size(E,1))*5;  
                 ones(Int64,size(E,1))*6;]

  P=nodesToPointset(V)

  return P,F, faceTypeLabel
end


#Do the work!
nodalDisplacements = runs(disp_max,disp_inc)

#Create makie compatible face and point set
pointSet,faceSet,faceTypeLabel=convertToFacePointSet(Ω)


function getCoordStep(Ω,nodalDisplacement)
  vd = visualization_data(Ω,"")
  grid = vd[1].grid
  V = get_node_coordinates(grid)
  V2= V + nodalDisplacement
  pointSet2=nodesToPointset(V2)
  return pointSet2
end

function getMagnitude(U)
  M=zeros(size(U,1))
  for q=1:1:size(U,1)
    M[q]=sqrt(U[q][1]^2 + U[q][2]^2 + U[q][3]^2)
  end
  return M
end

pointSet=getCoordStep(Ω,nodalDisplacements[1])

#Gather face and point set as GeometryBasics mesh
M=GeometryBasics.Mesh(pointSet,faceSet)

nodalColor=getMagnitude(nodalDisplacements[1])

#Visualize mesh
fig = Figure()

sl_step = Slider(fig[2, 1], range = 1:1:size(nodalDisplacements,1), startvalue = size(nodalDisplacements,1))

nodalColor = lift(sl_step.value) do stepIndex
  getMagnitude(nodalDisplacements[stepIndex])
end

M = lift(sl_step.value) do stepIndex
  GeometryBasics.Mesh(getCoordStep(Ω,nodalDisplacements[stepIndex]),faceSet)
end

titleString = lift(sl_step.value) do stepIndex
  "Step: "*string(stepIndex)
end

ax=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = titleString)
poly!(M, strokewidth=2,shading=true,color=nodalColor, transparency=false, overdraw=false)
fig

# Colorbar(fig[1, 2], limits = (0, 10), colormap = :viridis,vertical = true)