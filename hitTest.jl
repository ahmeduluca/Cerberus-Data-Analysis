using GLMakie
using GLMakie.GLFW
using GLMakie: to_native

#
#
global tempind = Observable(0)
global returnindex = Observable(0)

function zeroPt(dataF, expNumber)
    drawGraph(dataF[:, "Z-Disp Total (nm)"], dataF[:, "Z-Load (mN)"], expNumber)
end
#listen = on(tempind) do
#    println(tempind[])
#end

function drawGraph(rx, ry, no)
    rx0=Observable(rx)
    ry0=Observable(ry)
    listry=on(ry0)do 
        delete!(ax1,scatter1)
        scatter1 = GLMakie.scatter!(ax1, rx0[], ry0[], color =:red, markersize = 15);
        println("Deleted and renewed")
    end
    fig = Figure(resolution = (1440,900));
    #scene = Scene(resolution = (1440,900))
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
            );
    #create Sliders
    lsgrid = labelslidergrid!(fig, 
        ["Center Area"],
        Ref(LinRange(1:0.1:length(rx)));
        formats = [x-> "$(round(Int,x))"],
        labelkw = Dict([(:fontsize, 20)]),
        sliderkw = Dict([(:linewidth, 18)]),
        valuekw = Dict([(:fontsize, 20)])
    );

    #Layout sliders:

    sl_sublayout = GLMakie.GridLayout(height = 150);
    fig[2,1] = sl_sublayout;
    fig[2,1] = lsgrid.layout;

    #plot:
    scatter1 = GLMakie.scatter!(ax1, rx0[], ry0[], color =:red, markersize = 15);
    #make listener axes limits:
    #set starting ost for slope

    set_close_to!(lsgrid.sliders[1], experimentalSequence[no,"Load"]);

    limits = lift(lsgrid.sliders[1].value) do x
        x=round(Int,x)
        if(x<length(rx) && x>1 )
            xpm = abs(rx0[][x+1]-rx0[][x-1])
            ypm = abs(ry0[][x+1]-ry0[][x-1])
            if(xpm==0)
                xpm=1
            end
            if(ypm==0)
                ypm=0.01
            end
        else
            xpm = 1
            ypm = 0.01
        end
        xlimit = GLMakie.xlims!(ax1, rx0[][x]-xpm, rx0[][x]+xpm)
        ylimit = GLMakie.ylims!(ax1, ry0[][x]-ypm, ry0[][x]+ypm)
    end

    fig[1,2] = buttongrid = GLMakie.GridLayout(tellwidth = false);

    buttonlabels = ["Continue Calculation"];
    buttons = buttongrid[1, 1] = [
        Button(fig,
        label = l, height = 30, width = 200, fontsize = 16) 
        for l in buttonlabels
        ];
    tb = GLMakie.Textbox(fig[2, 2], placeholder = "Select  zero point");
    bt_sublayout = GLMakie.GridLayout(width = 200, height = 150);
    fig[1, 2] = bt_sublayout;

    fig[1,1] = GLMakie.GridLayout(width = 500, height = 500);
    glfw_window = to_native(display(fig))
    on(buttons[1].clicks) do click
        global returnindex[]= tempind[]
        global tempind[]=0
        GLMakie.vlines!(ax1, rx0[][returnindex[]], linewidth = 2)
        GLMakie.hlines!(ax1, ry0[][returnindex[]], linewidth = 2)
        GLMakie.xlims!(ax1,minimum(rx0[]),maximum(rx0[]))
        GLMakie.ylims!(ax1,minimum(ry0[]),maximum(ry0[]))
        GLFW.SetWindowShouldClose(glfw_window, true)
#        return returnindex[]
    # schedule(ret)
        #get zero index finally here!!
        #safe close figure
        #return index continue calculations .. DataFrame..
    end
    Makie.deactivate_interaction!(ax1, :rectanglezoom)
    spoint = select_points(ax1.scene); 
    on(spoint) do z
        x, y = z
        inXax = findall(p -> p < round(x,digits=2)+1 &&  p > round(x,digits=2)-1, rx0[])
        #println(inXax)
        inYax = findall(p -> p < round(y,digits=2)+0.2 &&  p > round(y,digits=2)-0.2, ry0[])
        #println(inYax)
        if !(isempty(inXax) || isempty(inYax))
            index = indexin(inXax, inYax)
            if !isempty(index)
                #println("Selection is $(rx[inXax[1]]),$(ry[inYax[index[1]]])")
                tb.displayed_string="Selection is $(round(rx0[][inXax[1]],digits=2)),$(round(ry0[][inYax[index[1]]],digits=2))"
                global tempind[] = inXax[1]
                global tempind[] = tempind[]
                for i in eachindex(rx)
                    rx0[][i] = rx[i]-rx[tempind[]]
                    ry0[][i] = ry[i]-ry[tempind[]]
                end
            else 
                tb.displayed_string="Please try to select a valid data point!"
         end
        else 
            tb.displayed_string="Please try to select a valid data point!"
        end
    end


end

function select_points(scene; blocking = false, priority=1, kwargs...)
    key = Mouse.left
    waspressed = Observable(false)
    point = Observable([Point2f(0,0)])
    point_ret = Observable(Point2f(0,0))

    on(events(scene).mousebutton, priority=priority) do event
        if event.button == key && is_mouseinside(scene)
            mp = mouseposition(scene)
            if event.action == Mouse.press
                waspressed[] = true
                point[][1] = mp
                point[] = point[]
                return Consume(blocking)
            end
        end
        if !(event.button == key && event.action == Mouse.press)
            if waspressed[] # User has selected the rectangle
                waspressed[] = false
                point_ret[] = copy(point[][1])
            end
            return Consume(blocking)
        end
        return Consume(false)
    end
    on(events(scene).mouseposition, priority=priority) do event
        if waspressed[]
            mp = mouseposition(scene)
            point[][1] = mp
            point[] = point[] # actually update observable
            return Consume(blocking)
        end
        return Consume(false)
    end

    return point_ret
end