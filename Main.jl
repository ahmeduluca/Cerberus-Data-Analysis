## It is the main file
include("FileIO.jl")
include("CalculationFunctions.jl")
global rawPaths = []
global fin = 0
input=DataFrame()
function readInp(raws)
    ##READING FILES
    i=1
    global input=[DataFrame() for _ in 1:length(raws)]
        for ind in raws
            mat,head=readdlm(ind,',',Float64,header=true,)
            global input[i]=DataFrame(mat, vec(head))
            i+=1
                ##PUT DATA IN DATAFRAME
        end
    end
## Process with taken paths.. read file -> DataFrame/ or just arrays (16 inputs-max-) 
#-> calculations & new Calculated DataFrame (70 columns of calculation) 
#-> write to file as CSV /& txt -> draw graphs.. 
## data should be separated / indexed for segments -> Approach, Load, Hold at Peak Load, First Unload, Drift, Full Unload..
#starter()
#yield();
#while(!istaskdone(fileTask))
#    nothing
#end
starter()
while(fin==0)
    sleep(1)
    println("waiting to complete FileIO")
end
readInp(rawPaths)

