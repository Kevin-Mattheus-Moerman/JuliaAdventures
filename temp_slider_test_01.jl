
using GLMakie, GeometryBasics
fig = Figure()

ax = Axis(fig[1, 1])

sl_x = Slider(fig[2, 1], range = 0:0.01:10, startvalue = 3)
sl_y = Slider(fig[1, 2], range = 0:0.01:10, horizontal = false, startvalue = 6)

point = lift(sl_x.value, sl_y.value) do x, y
    convert(GeometryBasics.Point,(x, y))
end

scatter!(point, color = :red, markersize = 20)

limits!(ax, 0, 10, 0, 10)

fig
