using GLMakie
using GLMakie.GLFW
using GLMakie: to_native

#imitialize plot
rx = 0:0.5:600
ry = rand(length(rx)) .+ rx .^ 2 .+ 3

fig = Figure(resolution = (1440,900))

#add Axis

ax1 = Axis(fig[1, 1],
        #borders
        aspect = 1,
        #title
        title = "Determine Zero Point",
        titlegap = 20, titlesize = 20,
        #x-axis 
        xautolimitmargin =(0,0), xgridwidth = 2, xticklabelsize = 16,
        xticks = LinearTicks(10), xticksize = 18,
        #y-axis
        yautolimitmargin = (0,0), ygridwidth = 2, yticklabelpad = 14,
        yticklabelsize = 16, yticks = LinearTicks(10), yticksize = 18
        )

#create Sliders
lsgrid = labelslidergrid!(fig, 
    ["Center Area"],
    Ref(LinRange(1:0.1:length(rx)));
    formats = [x-> "$(round(Int,x))"],
    labelkw = Dict([(:textsize, 20)]),
    sliderkw = Dict([(:linewidth, 18)]),
    valuekw = Dict([(:textsize, 20)])
)

#set starting ost for slope

set_close_to!(lsgrid.sliders[1], 1.0)

#Layout sliders:

sl_sublayout = GridLayout(height = 150)
fig[2,1] = sl_sublayout
fig[2,1] = lsgrid.layout

#plot:
scatter1 = scatter!(ax1, rx, ry, color =:red, makersize = 15)

#make listener axes limits:

limits = lift(lsgrid.sliders[1].value) do x
    xlimit = xlims!(ax1, rx[round(Int,x)]-25, rx[round(Int,x)]+25)
    ylimit = ylims!(ax1, ry[round(Int,x)]-25, ry[round(Int,x)]+25)
end

fig[1,2] = buttongrid = GridLayout(tellwidth = false)

buttonlabels = ["Continue Calculation"]

buttons = buttongrid[1, 1] = [
    Button(fig,
     label = l, height = 30, width = 200, textsize = 16) 
     for l in buttonlabels
    ]
tb = Textbox(fig[2, 2], placeholder = "Select  zero point")
bt_sublayout = GridLayout(width = 200, height = 150)
fig[1, 2] = bt_sublayout

fig[1,1] = GridLayout(width = 500, height = 500)
glfw_window = to_native(display(fig))
on(buttons[1].clicks) do click
    GLFW.SetWindowShouldClose(glfw_window, true)
    vlines!(ax1, rx[returnindex], linewidth = 2)
    hlines!(ax1, ry[returnindex], linewidth = 2)
    xlims!(ax1,minimum(rx),maximum(rx))
    ylims!(ax1,minimum(ry),maximum(ry))
    schedule(ret)
    #get zero index finally here!!
    #safe close figure
    #return index continue calculations .. DataFrame..
end
returnindex=1
spoint = select_point(ax1.scene)
on(spoint) do z
    x, y = z
    inXax = findall(p -> p < round(x,digits=2)+0.2 &&  p > round(x,digits=2)-0.2, rx)
    inYax = findall(p -> p < round(y,digits=2)+0.2 &&  p > round(y,digits=2)-0.2, ry)
    if !(isempty(inXax) || isempty(inYax))
        index = indexin(inXax, inYax)
        if !isempty(index)
            println("Selection is $(rx[inXax[index[1]]]),$(ry[inXax[index[1]]])")
            tb.displayed_string="Selection is $(round(rx[inXax[index[1]]],digits=2)),$(round(ry[inXax[index[1]]],digits=2))"
            global returnindex=inXax[index[1]]
        else 
            tb.displayed_string="Please try to select a valid data point!"
        end
    else 
        tb.displayed_string="Please try to select a valid data point!"
    end
end
#t = @task begin; ; end
#ret = @task begin; println("done"); end
#schedule(t); wait(t)
