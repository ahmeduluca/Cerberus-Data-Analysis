## It is the main file
include("FileIO.jl")
include("CalculationFunctions.jl")
include("hitTest.jl")

global rawPaths = []
global fin = 0
global exps = 0
global reading=1
global processing=0
global zeroing=0
global experimentalSequence = DataFrame( Exp = Int[], Approach = Int[], Load = Int[], PeakHold = Int[], UnloadToDrift = Int[],
Drift = Int[], UnloadToZero = Int[])

function readInp(raws)
    ##READING FILES
    for i in eachindex(raws)
        lines = readlines(raws[i])
        seq=1
        push!(experimentalSequence, [i,2,1,1,1,1,1])
        for l in eachindex(lines)
            if isempty( filter(x -> !isspace(x), lines[l]) )
                if(seq==1)
                    experimentalSequence[i,"Load"] = l
                elseif(seq==2)
                    experimentalSequence[i,"PeakHold"] = l
                elseif(seq==3)
                    experimentalSequence[i,"UnloadToDrift"] = l
                elseif(seq==4)
                    experimentalSequence[i,"Drift"] = l
                elseif(seq==5)
                    experimentalSequence[i,"UnloadToZero"] = l
                end
                seq+=1
            end
        end
    end
    global input=[DataFrame() for _ in 1:length(raws)]
    global output=[DataFrame() for _ in 1:length(raws)]
    for ind in raws
        mat,head=readdlm(ind,',',Float64,header=true,)
        global input[exps+1]=DataFrame(mat, vec(head))
        global exps+=1
    end
    ##PUT DATA IN DATAFRAME
    sleep(1)
    global reading = 0
end

function process(rawDataFrames)
    axes3d=["Z","X","Y"]
    global processing=1
# 1. zeroing time
    for df in eachindex(rawDataFrames)
        t0 = rawDataFrames[df][2,"Time (s)"]
        output[df] = transform(rawDataFrames[df], "Time (s)" => ByRow(t -> t-t0) => "Time (s)")
        output[df][1,"Time (s)"]=0.
# 2. Scaling Displacement Measurement Range and Offsetting Displacements -Raw (Voltage) to Nanometer-
        for axis in axes3d
            output[df] = displacementVoltToNano(output[df], axis)
             #Amplitude data Voltage to Lenght conversion
            output[df] = displacementVoltToNano(output[df], axis, 1)
        end
# 3. Scaling Load Measurement Range and Offsetting to mN 
        for axis in axes3d
            output[df] = loadVoltToMili(output[df], axis, 0)
        end
        # Harmonic Force data Voltage to MicroNewtons
        output[df] = loadVoltToMili(output[df], "Z-Excite", 1)
        output[df] = loadVoltToMili(output[df], "X-Exc (V)", 1)        
        output[df] = loadVoltToMili(output[df], "Y-Exc (V)", 1)           
# 4. Zero Point Determination
        index, output[df] = zeroDetermination(output[df], df)
        while(zeroing==1)
            sleep(1)
        end
        
        # 5. Leaf Spring Correction on Vertical Load
        output[df] = leafspringCorrection(output[df], index)
        # 6. Thermal Drift Correction on Displacement Data
        output[df] = thermaldriftCorrection(output[df], experimentalSequence[df,:])
        
        ## Stiffness Calculations CSM
# (S2/Load) vs Depth and (1/S) vs (1/disp) frame compliance subtraction accuracy &/ calculation of frame compliances 
## this might be interactive too; change the number and see calibration etc..
## Machine Compliance Subtraction from Load
##Viscous damping and Loss/Storage Moduli - if available
        for axis in axes3d
            output[df] = CSM(axis, output[df])
        end
## *Contact Depth vs Load *on Sample --Finally--
## Tip Area Function -> Contact Area Calculation
## Calibration from the data can be an option via interactive curve fit application! by given moduli and load-disp data
        output[df] = contact(output[df])
## Modulus & Hardness CSM
## Shear Modulus - Interfacial Shear Strength - Mindlin calcs. - Static Friction
## Modulus & Hardness Static Method
        output[df] = mechanicalProperties(output[df])

## Interfacial Shear/Friction properties ; Poisson Ratio
        for ax in axes3d
            output[df] = harmonicPerStaticF(ax,output[df])
        end
        ## mindlin contact analysis functions..
#= 
        zeroedGraphs = Figure(resolution=(1440,900))
        axDisplacement = Axis(zeroedGraphs[1,1])
        axTime = Axis(zeroedGraphs[1,2])
        axForceAmp = Axis(zeroedGraphs[1,3])
        display(zeroedGraphs)
        scatter!(axDisplacement, output[df][:, "Time Zeroed (s)"], output[df][:, "Z-F/X (N/m)"])
        sleep(5)
        scatter!(axDisplacement, output[df][:, "Time Zeroed (s)"], output[df][:, "Z-F/XCosPhi-(kₛ-mω²) (N/m)"])
        sleep(5)
        scatter!(axDisplacement, output[df][:, "Time Zeroed (s)"], output[df][:, "Z-Stiffness (N/m)"])
        for axis in axes3d
        #    loadVsDisp = scatter!(axDisplacement,output[df][:, axis*"-Disp Zeroed (nm)"],output[df][:, axis*"-Load Zeroed (mN)"])
        #    ampVsDisp = scatter!(axDisplacement,output[df][:, axis*"-Disp Zeroed (nm)"],output[df][:, axis*"-Amplitude (nm)"])
            scatter!(axTime,output[df][:,"Contact Depth (nm)"],output[df][:, "Sample Elastic Modulus (MPa)"])
            scatter!(axTime,output[df][:,"Contact Depth (nm)"],output[df][:, "Hardness (GPa)"])
            scatter!(axForceAmp,output[df][:, "Time Zeroed (s)"],output[df][:, axis*"-Amplitude (nm)"])
            scatter!(axForceAmp,output[df][:, "Time Zeroed (s)"],output[df][:, axis*"-Force (μN)"])
        end =#
    end
    global processing=0
end

starter()
while(fin==0)
    sleep(1)
end
readInp(rawPaths)
while(reading==1)
    sleep(1)
end
process(input)
while(processing==1)
    sleep(1)
end
writer(dirname(rawPaths[1]),output)
