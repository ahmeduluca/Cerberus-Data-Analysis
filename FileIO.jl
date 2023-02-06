using Gtk
using Glob
using DataFrames
using DelimitedFiles
using CSV
using GLMakie
using XLSX

# Choose whether all the txts or by selection of files -> hold filepath names in a container then read
# and process
global rawPaths=[]
global fin=0
function starter()
    #GTK WINDOW FOR SELECTION OF BATCH TYPE
    win=GtkWindow("Batch Process Type")
    hbox=GtkBox(:h)
    push!(win,hbox)
    allIn = GtkButton("All indents in a File")
    bySel = GtkButton("Selection of indents in a Folder")
    push!(hbox, allIn)
    push!(hbox, bySel)
    showall(win)
    function button_clicked_callback(widget)
    # open dialog depending on the selection -All txts in Main folder or by selection of txts..
        destroy(win)
        if (occursin("All",get_gtk_property(widget,:label,String)))
            ## get the file paths in container..
                dir = open_dialog("Select Main Folder", action=GtkFileChooserAction.SELECT_FOLDER)
                if isdir(dir)
                   #Get applicable text paths  as R_###-### etc..
                   #println(dir)
                  global rawPaths=glob([r"(?i)R\p{Lu}[[:xdigit:]]{3}\d{3}.TXT"],string(dir))
                end
        else
            # Test each selection whether it is in applicable format
                dirs = open_dialog("Pick some text files", GtkNullContainer(), (GtkFileFilter("*.TXT, *.CSV", name="TXT Files of Indents"),), select_multiple=true)
                global rawPaths=[i for i in dirs if occursin(r"(?i)R\p{Lu}[[:xdigit:]]{3}\d{3}.TXT",i)] # Maybe we can lift this control to give more functionality
        end
        sleep(1)
        global fin = 1
        #print(rawPaths)
    end
#CONNECTING BUTTONS TO FUNCTION (EVENT HANDLER)
    id0 = signal_connect(button_clicked_callback, allIn, "clicked")
    id1 = signal_connect(button_clicked_callback, bySel, "clicked")
if(!isinteractive())
    Gtk.waitforsignal(win, :destroy)
end
end
#BUTTON CLICKED EVENT CALLBACK

function writer(file, data)
    file=file*raw"\processed.xlsx"
    XLSX.openxlsx(file, mode="w") do xl
        XLSX.rename!(xl[1],basename(rawPaths[1]))
        for i in eachindex(data)
            if i>1
                XLSX.addsheet!(xl, basename(rawPaths[i]))
            end
            XLSX.writetable!(xl[i], collect(eachcol(data[i])), names(data[i]))
        end
    end
end



