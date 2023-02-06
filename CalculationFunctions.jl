import Interpolations: interpolate, gradient, Gridded, Linear
import Statistics: mean
import CurveFit: linear_fit
include("hitTest.jl")
global returnindex = Observable(0)
global zeroing = 0
global epsilon = 0.75
global poissonSample = 0.2
global youngIndenter = 1e6
global poissonIndenter = 0.2
global oscFreq = 125
#=
Data Calculation Functions
#
#DataFrame/ or just arrays (16 inputs-max-) 
#-> calculations & new Calculated DataFrame (70 columns of calculation) 
=#


## ! Segment determination -via meta file or user defined..
## data has blank rows in between segments but does not indicates any name for segment

#=
To find range and voltage offset through the data ([1]st element)
=#
function dispOffset(dispArrayInput,range=3)
    return round(((dispArrayInput[1]-(range+1)*1e4)-5)/10)
end

function findRange(inputArray)
    return round((inputArray[1]/1e4)-1)
end

#= Range 3 is the range for Cerberus, but to make it more generic for any DCM system
  the function can take range argument too.
Raw displacement data firstly will be taken into dispRange function which subtracts
the offset of voltage measurement at that range from all data
=#
function dispRange(dispOffset,range=3)
    return 5 + dispOffset*10 + (range+1)*1e4
end

###### DC Displacement Calibrations #####

### Displacement Constant Term
z_dispCalDC0 = 9944.3
x_dispCalDC0 = 18506
y_dispCalDC0 = 11091

### Displacement Linear Term
z_dispCalDC1 = 2.0903
x_dispCalDC1 = 4.8506
y_dispCalDC1 = 1.6121

### Displacement Squared Term
z_dispCalDC2 = -0.3789
x_dispCalDC2 = -0.5316
y_dispCalDC2 = -0.26013

##### DC Displacement Calibration Values (calValue) for  Z, X and Y axes:
function dispCalDC(dispOff, Axis)
    if(Axis == "Z")
        return z_dispCalDC0 + z_dispCalDC1 * (dispOff - 128) + z_dispCalDC2 * (dispOff - 128)^2
    elseif(Axis == "X")
        return x_dispCalDC0 + x_dispCalDC1 * (dispOff - 128) + x_dispCalDC2 * (dispOff - 128)^2
    elseif(Axis == "Y")
        return y_dispCalDC0 + y_dispCalDC1 * (dispOff - 128) + y_dispCalDC2 * (dispOff - 128)^2
    end
end

###### AC Displacement Calibrations #####

### Displacement Constat Term
z_dispCalAC0 = 9935.9
x_dispCalAC0 = 18509
y_dispCalAC0 = 11088

### Displacement Linear Term
z_dispCalAC1 = 2.0639
x_dispCalAC1 = 4.8118
y_dispCalAC1 = 1.6513

### Displacement Squared Term
z_dispCalAC2 = -0.37801
x_dispCalAC2 = -0.53142
y_dispCalAC2 = -0.26013

##### AC Displacement Calibration Values (calValue) for  Z, X and Y axes:
function dispCalAC(dispOff, Axis)
    if(Axis == "Z")
        return z_dispCalAC0 + z_dispCalAC1 * (dispOff - 128) + z_dispCalAC2 * (dispOff - 128)^2
    elseif(Axis == "X")
        return x_dispCalAC0 + x_dispCalAC1 * (dispOff - 128) + x_dispCalAC2 * (dispOff - 128)^2
    elseif(Axis == "Y")
        return y_dispCalAC0 + y_dispCalAC1 * (dispOff - 128) + y_dispCalAC2 * (dispOff - 128)^2
    end
end

#Finally, create the voltage data to actual displacement by
#applying calibration values.

function totalDisp(dispVolt, Axis, dispOf, mode=0 ,range=3)
    if(mode == 0)
        calValue = dispCalDC(dispOf, Axis)
    elseif(mode == 1)
        calValue = dispCalAC(dispOf, Axis)
    end
    return (dispVolt * calValue) / 2^range
end

function displacementVoltToNano(rawData, Axis, mode=0)
    if(mode == 0)
        Range = findRange(rawData[:, Axis*"-Disp (V)"])
        dispOf = dispOffset(rawData[:, Axis*"-Disp (V)"], Range)
        d0 = dispRange(dispOf, Range)
        transform!(rawData, Axis*"-Disp (V)" => ByRow(d -> d - d0) => Axis*"-Disp R3 (V)")
        d_off = rawData[1,Axis*"-Disp R3 (V)"]
        transform!(rawData, Axis*"-Disp R3 (V)" => ByRow(d -> d-d_off) => Axis*"-Disp R3 Zeroed (V)")
        transform!(rawData, Axis*"-Disp R3 (V)" => ByRow(d -> totalDisp(d, Axis, dispOf, mode, Range)) => Axis*"-Disp Total (nm)")
    elseif(mode == 1)
        Range = findRange(rawData[:, Axis*"-Disp (V)"])
        dispOf = dispOffset(rawData[:, Axis*"-Disp (V)"], Range)
        transform!(rawData, Axis*"-Amp (V)" => ByRow(d -> totalDisp(d, Axis, dispOf, mode, Range)) => Axis*"-Amplitude (nm)")
    end
    return rawData
end

#=
Zero Point Determination; live GUI need here!
after determination return 3 displacement 3 load and time arrays with zeroed form
=#

function zeroDetermination(dataF, number)
    global zeroing=1
    zeroPt(dataF, number)
    while(returnindex[] == 0)
        sleep(1)
    end
    index = returnindex[]
    zeroingList = ["Z-Disp Total (nm)", "X-Disp Total (nm)", "Y-Disp Total (nm)", "Z-Load (mN)",
    "X-Load (mN)", "Y-Load (mN)", "Time (s)"]
    zeroedList = ["Z-Disp Zeroed (nm)", "X-Disp Zeroed (nm)", "Y-Disp Zeroed (nm)", "Z-Load Zeroed (mN)",
     "X-Load Zeroed (mN)", "Y-Load Zeroed (mN)", "Time Zeroed (s)"]
    for col in eachindex(zeroingList)
        d0 = dataF[index, zeroingList[col]]
        transform!(dataF, zeroingList[col] => ByRow(d -> d-d0) => zeroedList[col])
    end
    global returnindex[]=0
    global zeroing=0
    return index, dataF
end

#=
Thermal drift correction to displacement:
subtract from the "zeroed" displacement (depth into surface) 
slope of that displacement-time at the peak hold *times zeroed time[i]
=#
function thermaldriftCorrection(data, indexex)
    start=indexex["Drift"]
    stop=indexex["UnloadToZero"]
    cnst,slope = linear_fit(data[start:stop,"Time Zeroed (s)"],data[start:stop,"Z-Disp Zeroed (nm)"])
    thermalDriftCorr = convert(Array{Union{Float64,Missing}}, data[:,"Z-Disp Zeroed (nm)"])
    for i in eachindex(thermalDriftCorr)
        thermalDriftCorr[i] = data[i,"Z-Disp Zeroed (nm)"] - convert(Float64,slope) * data[i,"Time Zeroed (s)"] 
    end
    data[!,"Z-Disp Thermal Drift Corrected (nm)"] = thermalDriftCorr
    return data
end

##### Displacement Offset (dispOffset) Values for Z, X and Y axes: !These will be calculated for each test below are generic examples
#z_dispOff = 131
#x_dispOff = 129
#y_dispOff = 128
###############################

## DC Load Calibrations

### Load Constant Term
z_loadCalDC0 = 47624
x_loadCalDC0 = 43313
y_loadCalDC0 = 58032

### Load Linear Term
z_loadCalDC1 = 8.0558
x_loadCalDC1 = 21.588
y_loadCalDC1 = 3.5987

### Load Squared Term
z_loadCalDC2 = 0.003
x_loadCalDC2 = 0.0067645
y_loadCalDC2 = 0.00064854

##### DC Load Calibration Values (calValue) for  Z, X and Y axes:
function loadCalDC(dispOff, Axis)
    if(Axis == "Z")
        return z_loadCalDC0 + z_loadCalDC1 * (dispOff - 128) + z_loadCalDC2 * (dispOff - 128)^2
    elseif(Axis == "X")
        return x_loadCalDC0 + x_loadCalDC1 * (dispOff - 128) + x_loadCalDC2 * (dispOff - 128)^2
    elseif(Axis == "Y")
        return y_loadCalDC0 + y_loadCalDC1 * (dispOff - 128) + y_loadCalDC2 * (dispOff - 128)^2
    end
end

## AC Load Calibrations

### Load Constant Term
z_loadCalAC0 = 14023.5
x_loadCalAC0 = 12767
y_loadCalAC0 = 17098

### Load Linear Term
z_loadCalAC1 = 2.36
x_loadCalAC1 = 6.4062
y_loadCalAC1 = 1.1621

### Load Squared Term
z_loadCalAC2 = -0.0002
x_loadCalAC2 = 0.0019054
y_loadCalAC2 = -0.0013593

##### AC Load Calibration Values (calValue) for  Z, X and Y axes:
function loadCalAC(dispOff, Axis)
    if(startswith(Axis,"Z"))
        return ((z_loadCalAC0 + z_loadCalAC1 * (dispOff - 128) + z_loadCalAC2 * (dispOff - 128)^2)*z_modCal)
    elseif(startswith(Axis,"X"))
        return ((x_loadCalAC0 + x_loadCalAC1 * (dispOff - 128) + x_loadCalAC2 * (dispOff - 128)^2)*x_modCal)
    elseif(startswith(Axis,"Y"))
        return ((y_loadCalAC0 + y_loadCalAC1 * (dispOff - 128) + y_loadCalAC2 * (dispOff - 128)^2)*y_modCal)
    end
end

#=
Similar function for range of load
=#
function findRangeLoad(inputArray)
    return round((inputArray[1]/10)-1)
end

function loadAdcScaled(rawLoad,range=3)
    return rawLoad-(range+1)*10-5
end
#=
Convert Load data in voltage to mN
=#
function loadConvert(scaledLoad, Axis, dispOff,mode=0)
    if(mode == 0)
        loadCal = loadCalDC(dispOff, Axis)
        return -scaledLoad*loadCal/1000
    elseif(mode == 1)
        loadCal = loadCalAC(dispOff, Axis)
        return scaledLoad*loadCal
    end
end

function loadVoltToMili(rawData, Axis, mode)
    if(mode == 0) 
        Range = findRangeLoad(rawData[:, Axis*"-Load (V)"])
        d_off = dispOffset(rawData[:, Axis*"-Disp (V)"], Range)
        transform!(rawData, Axis*"-Load (V)" => ByRow(d -> loadAdcScaled(d,Range)) => Axis*"-Load (V)")
        transform!(rawData, Axis*"-Load (V)" => ByRow(d -> loadConvert(d, Axis, d_off, mode)) => Axis*"-Load (mN)")
    elseif(mode == 1)
        Range = findRange(rawData[:, Axis[1]*"-Disp (V)"])
        d_off = dispOffset(rawData[:, Axis[1]*"-Disp (V)"], Range)
        transform!(rawData, Axis => ByRow(d -> loadConvert(d, Axis, d_off, mode)) => Axis[1]*"-Force (μN)")
    end
    return rawData
end

    #=
subtract the load on the leaf spring by calculating it via
slope of load-disp at approach *times displacement[i] 
=#
function leafspringCorrection(data, surfaceIndex)
    fApproach = interpolate((data[2:surfaceIndex, "Z-Disp Zeroed (nm)"],),data[2:surfaceIndex, "Z-Load Zeroed (mN)"],Gridded(Linear()))
    slope =  mean(only.(gradient.(Ref(fApproach), data[2:surfaceIndex, "Z-Disp Zeroed (nm)"])))
    Z_Load_kLS_Corrected = zeros(length(data[:,"Z-Load Zeroed (mN)"]))
    for i in eachindex(Z_Load_kLS_Corrected)
        Z_Load_kLS_Corrected[i] = data[i, "Z-Load Zeroed (mN)"]-slope*data[i, "Z-Disp Zeroed (nm)"]
    end
    data[!,"Z-Load Leaf Spring Corrected (mN)"] = Z_Load_kLS_Corrected
    return  data
end

##### Modcal Factor 
z_modCal = 1.15e-8
x_modCal = 1.15e-8
y_modCal = 1.15e-8

### Spring Constant Constant Term
z_spring0 = 230
x_spring0 = 207
y_spring0 = 204

### Spring Constant Linear Term
z_spring1 = 0.0055042
x_spring1 = 0.
y_spring1 = 0.

### Spring Constant Squared Term
z_spring2 = 1.77e-5
x_spring2 = 0.
y_spring2 = 0.

### Spring Constant Cubic Term
z_spring3 = 1.77e-5
x_spring3 = 0.
y_spring3 = 0.

### DC Load Frame Stiffness
z_frameStiffDC = 2.2e5
x_frameStiffDC = 2.2e5
y_frameStiffDC = 2.2e5

### Mass of indenter = Fitted Masses

### Damping Coefficient Constant Term
z_damp0 = 0.019
x_damp0 = 0.0066396
y_damp0 = 0.017396

### Damping Coefficient Linear Term
z_damp1 = 2.61e-6
x_damp1 = 0.
y_damp1 = 2.61e-6

### Damping Coefficient Squared Term
z_damp2 = 7.94e-6
x_damp2 = 0.
y_damp2 = 7.94e-6

### Damping Coefficient Cubic Term
z_damp3 = -3.46e-9
x_damp3 = 0.
y_damp3 = -3.46e-9

### Damping Coefficient Quad Term
z_damp4 = 3.1e-9
x_damp4 = 0.
y_damp4 = 3.1e-9

### AC Load Frame Stiffness
z_frameStiffAC = 1.3e5
x_frameStiffAC = 2.9e4
y_frameStiffAC = 6.5e4


###### Individual Masses ###


### Mass From Uncoupled Dyncal
z_dynaMass = 1.04e-4
x_dynaMass = 9.65e-5
y_dynaMass = 9.95e-5

### Mass of Coupler
z_couplerMass = 2.1e-5
x_couplerMass = 2.1e-5
y_couplerMass = 2.1e-5

### Mass of Indenter (Ti Rod)
z_indMass = 1.47e-6
x_indMass = 1.47e-6
y_indMass = 1.47e-6

### Mass of Quartz Rod
z_qMass = 6.32e-6
x_qMass = 6.32e-6
y_qMass = 6.32e-6

### Total Mass
z_totMass = z_dynaMass + z_couplerMass + z_indMass + z_couplerMass
x_totMass = x_dynaMass + x_couplerMass + x_indMass + x_couplerMass
y_totMass = y_dynaMass + y_couplerMass + y_indMass + y_couplerMass

### Fitted Mass
z_fitMass = 1.65e-4
x_fitMass = 1.45e-4
y_fitMass = 1.4e-4

###### Area Function ######

area2 = 24.56
area1 = 272.89
area05 = 23376
area25 = -288280
area0125= 754610
area006125 = -492830
area0030625 = -1319.00673 ## ! 
area00153125 = -2428.99083 ## !
area00075625 =-2985.8581 ## !

#############################
#=
Phase Correction
=#
function phaseCorrection(rawPhase, co=1)
    return rawPhase*co- (-0.0844 * oscFreq - 0.0596)
end

#Cos Phase
cosPhase(correctedPhase) = cos(deg2rad(correctedPhase))

#Sin Phase
sinPhase(correctedPhase) = sin(deg2rad(correctedPhase))

#=
Harmonic Force / Harmonic Displacement (for getting stiffness from DCM)
=#
function forcePerAmp(forceUn, ampNm)
    return (forceUn * 1e-6) / (ampNm * 1e-9) # μN/nm
end

function fperAsubSpring(fperA, corrPhase, ax)
    mass=0
    spring=0
    if(ax=="Z")
        mass = z_fitMass
        spring = z_spring0
    elseif(ax=="X")
        mass = x_fitMass
        spring = x_spring0
    elseif(ax=="Y")
        mass = y_fitMass
        spring = y_spring0
    end
    return fperA*cosPhase(corrPhase)-spring+mass*(oscFreq*2pi)^2
end

function stiffnessAC(fperA, ax)
    kLoadFrame = Inf
    if(ax=="Z")
        kLoadFrame = z_frameStiffAC
    elseif(ax=="X")
        kLoadFrame = x_frameStiffAC
    elseif(ax=="Y")
        kLoadFrame = y_frameStiffAC
    end
    return 1/((1/(fperA)) - (1/kLoadFrame))
end

function StiffnesSqrPerLoad(lfCorrectedLoad, acStiffness)
    return  1000 * acStiffness^2 / lfCorrectedLoad
end

function dampingAC(fperA, correctedPhase, acStiffness, ax)
    dampingCoeff = 0
    if(ax=="Z")
        dampingCoeff = z_damp0
    elseif(ax=="X")
        dampingCoeff = z_damp0
    elseif(ax=="Y")
        dampingCoeff = z_damp0
    end
    return (fperA*sinPhase(correctedPhase)-2*pi*dampingCoeff*oscFreq)/acStiffness
end

function CSM(axis, data)
    correctPhase = zeros(length(data[:, axis*"-Pha (Dg)"])) 
    for i in eachindex(correctPhase)
        correctPhase[i] = phaseCorrection(data[i, axis*"-Pha (Dg)"], 1)
    end
    data[!, axis*"-Phase (°)"] = correctPhase

    FdivX = zeros(length(data[:, axis*"-Force (μN)"]))
    for i in eachindex(FdivX)
        FdivX[i] = forcePerAmp(data[i, axis*"-Force (μN)"], data[i, axis*"-Amplitude (nm)"])
    end
    data[!, axis*"-F/X (N/m)"] = FdivX

    FdivXSubS = zeros(length(FdivX))
    for i in eachindex(FdivXSubS)
        FdivXSubS[i] = fperAsubSpring(FdivX[i], correctPhase[i], axis)
    end
    data[!, axis*"-F/XCosPhi-(kₛ-mω²) (N/m)"] = FdivXSubS

    stiffness = zeros(length(FdivXSubS))
    compliance = zeros(length(stiffness))
    for i in eachindex(stiffness)
        stiffness[i] = stiffnessAC(FdivXSubS[i], axis)
        compliance[i] = 1/stiffness[i]
    end
    data[!, axis*"-Stiffness (N/m)"] = stiffness
    data[!, axis*"-Compliance (m/N)"] = compliance

    s2DivLoad = zeros(length(stiffness))
    for i in eachindex(stiffness)
        s2DivLoad[i] = StiffnesSqrPerLoad(data[i,"Z-Load Leaf Spring Corrected (mN)"], stiffness[i])
    end
    data[!, axis*"-Stiffness² / Load (N/m²)"] = s2DivLoad
    
    damp = zeros(length(stiffness))
    for i in eachindex(stiffness)
        damp[i] = dampingAC(FdivX[i], correctPhase[i], stiffness[i], axis)
    end
    data[!, axis*"-Damping"] = damp
    
    return  data
end

function displacementIntoSurf(tdCorrectedDisp, lfCorrectedLoad, lfStiffness)
    return tdCorrectedDisp - (lfCorrectedLoad * 0.001 / lfStiffness) * 1E9
end

function contactDepth(dispInSur, lfCorrectedLoad, acStiff, ϵ)
    return dispInSur - 1E9 * ϵ * (0.001*lfCorrectedLoad / acStiff)
end

function contactArea(contDepth, coQuad, coLin, co0_5, co0_25, co0_125, co_0_06125, co0_030625, co0_0153125, co0_00765625)
    if(contDepth>0)
    	return coQuad*contDepth^2 + coLin*contDepth + co0_5*contDepth^0.5 +
        co0_25*contDepth^0.25 + co0_125*contDepth^0.125 + co_0_06125*contDepth^0.06125 + co0_030625*contDepth^0.030625 +
        co0_0153125*contDepth^0.0153125 + co0_00765625*contDepth^0.00765625 
    else
        return 0.
    end
end

function contact(data)

    dispInSurf = zeros(length(data[:,"Z-Disp Zeroed (nm)"]))
    for i in eachindex(dispInSurf)
        dispInSurf[i] = displacementIntoSurf(data[i,"Z-Disp Thermal Drift Corrected (nm)"],
         data[i,"Z-Load Leaf Spring Corrected (mN)"],z_frameStiffDC)
    end
    data[!,"Displacement into Surface (nm)"] = dispInSurf

    cdepth = zeros(length(dispInSurf))
    for i in eachindex(cdepth)
        cdepth[i] = contactDepth(dispInSurf[i], data[i,"Z-Load Leaf Spring Corrected (mN)"], 
        data[i,"Z-Stiffness (N/m)"], epsilon)
    end
    data[!, "Contact Depth (nm)"] = cdepth

    carea = zeros(length(cdepth))
    dispInSurf = []
    for i in eachindex(carea)
        carea[i] = contactArea(cdepth[i], area2, area1, area05, area25, area0125, area006125, area0030625, area00153125, area00075625)
    end
    data[!, "Contact Area (nm²)"] = carea
    cdepth = []
    carea = []

    return data
end

function reducedModulus(acStiff, contArea)
    return 1e3 * 0.5 * pi^0.5 * acStiff / contArea^0.5
end

function modulusSample(redModul, samplePoisson, indenterModul, indenterPoisson)
    return (1-samplePoisson^2)/(1/redModul-((1-indenterPoisson^2)/indenterModul))
end

function reducedShearModulus(lateralStiffness, contactArea)
    return pi^0.5 * lateralStiffness/(contactArea^0.5*8)
end

function interfacialShear(reducedShear, indenterShearModul, samplePoisson, indenterPoisson)
    return (2-samplePoisson)/(1/reducedShear-((2-indenterPoisson)/indenterShearModul))
end

function hardness(lfCorrLoad, contArea)
    return 1e6 * lfCorrLoad / contArea
end

##Area per stiffnes sqr
function apce(contArea, acStiff)
    return 1E-18 * contArea/acStiff^2
end

## Constant E Area average(apce@drift-or unloading?) * z_stiffness^2 -> area calculated
function cnstEarea()

end

function harmonicPerStaticF(axis, data)
    fperf = zeros(length(data[:, axis*"-Force (μN)"]))
    for i in eachindex(fperf)
        fperf[i] = 1e-3*data[i, axis*"-Force (μN)"] / data[i, axis*"-Load Zeroed (mN)"]
    end
    data[!, axis*"-Harmonic Force / Quasi Static Load"] = fperf
    return data
end

function mechanicalProperties(data)
    redmod = zeros(length(data[:,"Contact Depth (nm)"]))
    apers = zeros(length(redmod))
    for i in eachindex(redmod)
        if(data[i, "Contact Area (nm²)"] > 0)
            redmod[i] = reducedModulus(data[i, "Z-Stiffness (N/m)"], data[i, "Contact Area (nm²)"])
            apers[i] = apce(data[i, "Contact Area (nm²)"],data[i, "Z-Stiffness (N/m)"])
        else
            redmod[i] = 0
            apers[i] = 0
        end
    end
    data[!, "Reduced Modulus (MPa)"] = redmod
    data[!, "Area per Stiffness Square (m⁴/N²)"] = apers
    apers=[]

    sammod = zeros(length(redmod))
    for i in eachindex(sammod)
        if(redmod[i] > 0)
            sammod[i] = modulusSample(redmod[i], poissonSample, youngIndenter, poissonIndenter)
        else
            sammod[i] = 0
        end
    end
    data[!, "Sample Elastic Modulus (MPa)"] = sammod
    redmod=[]

    hard = zeros(length(sammod))
    sammod=[]
    for i in eachindex(hard)
        if(data[i, "Contact Area (nm²)"] > 0)
            hard[i] = hardness(data[i,"Z-Load Leaf Spring Corrected (mN)"], data[i, "Contact Area (nm²)"])
        else
            hard[i] = 0
        end
    end
    data[!, "Hardness (GPa)"] = hard

    return data
end

function poissonMidlin(lateralStiff, verticalStiff)

end

function shearStrengthMindlin(data)

end

function staticFrictionMindlin(data)

end

## Apply Classic Oliver-Pharr method for elastic unloading part; estimate Modulus and Hardness for overall indentation 
function contactMechanicUnloading(data)

end
