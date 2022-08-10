using Gtk
# Choose whether all the txts or by selection of files -> hold filepath names in a container then read
# and process

win=GtkWindow("Batch Process Type")
hbox=GtkBox(:h)
push!(win,hbox)
allIn = GtkButton("All indents in a File")
bySel = GtkButton("Selection of indents in a Folder")
push!(hbox, allIn)
push!(hbox, bySel)
showall(win)

function button_clicked_callback(widget)
#await
# open dialog depending on the selection -All txts in Main folder or by selection of txts..
if(occursin("All",get_gtk_property(widget,:label,String)))
    ## get the file paths in container..
   dir = open_dialog("Select Main Folder", action=GtkFileChooserAction.SELECT_FOLDER)
        if isdir(dir)
            #Get applicable text paths  as R_###-### etc..
            println(dir)
        end
else
    # Test each selection whether it is in applicable format
    open_dialog("Pick some text files", GtkNullContainer(), ("*.txt, *.csv",), select_multiple=true)
end
    destroy(win)
end

id0 = signal_connect(button_clicked_callback, allIn, "clicked")
id1 = signal_connect(button_clicked_callback, bySel, "clicked")

