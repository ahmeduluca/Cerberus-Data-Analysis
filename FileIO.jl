using Gtk
using Glob
# Choose whether all the txts or by selection of files -> hold filepath names in a container then read
# and process
rawPaths=[]
#GTK WINDOW FOR SELECTION OF BATCH TYPE
win=GtkWindow("Batch Process Type")
hbox=GtkBox(:h)
push!(win,hbox)
allIn = GtkButton("All indents in a File")
bySel = GtkButton("Selection of indents in a Folder")
push!(hbox, allIn)
push!(hbox, bySel)
showall(win)

#BUTTON CLICKED EVENT CALLBACK
function button_clicked_callback(widget)
# open dialog depending on the selection -All txts in Main folder or by selection of txts..
    if(occursin("All",get_gtk_property(widget,:label,String)))
        ## get the file paths in container..
        dir = open_dialog("Select Main Folder", action=GtkFileChooserAction.SELECT_FOLDER)
            if isdir(dir)
            #Get applicable text paths  as R_###-### etc..
            println(dir)
            rawPaths=glob([r"(?i)R\p{Lu}[[:xdigit:]]{3}\d{3}.TXT"],string(dir))
            end
    else
        # Test each selection whether it is in applicable format
        dirs = open_dialog("Pick some text files", GtkNullContainer(), (GtkFileFilter("*.TXT, *.CSV", name="TXT Files of Indents"),), select_multiple=true)
        rawPaths=[i for i in dirs if occursin(r"(?i)R\p{Lu}[[:xdigit:]]{3}\d{3}.TXT",i)] # Maybe we can lift this control to give more functionality
    end
    destroy(win)
    print(rawPaths)
end

#CONNECTING BUTTONS TO FUNCTION (EVENT HANDLER)
id0 = signal_connect(button_clicked_callback, allIn, "clicked")
id1 = signal_connect(button_clicked_callback, bySel, "clicked")

## Process with taken paths.. read file -> DataFrame/ or just arrays (16 inputs-max-) 
#-> calculations & new Calculated DataFrame (70 columns of calculation) 
#-> write to file as CSV /& txt -> draw graphs.. 
