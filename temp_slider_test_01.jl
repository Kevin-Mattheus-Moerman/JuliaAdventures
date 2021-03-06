
using GLMakie

xs = range(0, 10, length = 25)
ys = range(0, 15, length = 25)

fig = Figure()

sl = Slider(fig[2, 1], range = 0:0.1:3, startvalue = 1)

zs = lift(sl.value) do s
    convert(Float64,s).*[cos(x) * sin(y) for x in xs, y in ys]
end

titleString = lift(sl.value) do s
    "Step: "*string(s)
end

ax = Axis(fig[1, 1], title = titleString)

hm=heatmap!(xs, ys, zs,colormap = Reverse(:Spectral))
Colorbar(fig[1, 2],hm,label = "Values")
fig
