import Interpolations: interpolate, gradient, Gridded, Linear
import Statistics: mean
include("hitTest.jl")
#=
Data Calculation Functions
#
#DataFrame/ or just arrays (16 inputs-max-) 
#-> calculations & new Calculated DataFrame (70 columns of calculation) 
=#


## ! Segment determination -via meta file or user defined.. data has blank rows
## in between segments but not indcates any name for segment

#=
To find range and voltage offset through the data ([1]th element)
=#
function dispOffset(dispArrayInput,range=3)
    return round(((dispArrayInput[1]-(range+1)*1e4)-5)/10)
end

function findRange(inputArray)
    return round((inputArray[1]/1e4)-1)
end

#= Range 3 is the range for Cerberus, but to make it more generic for any DCM system
  the function can be take range argument too.
Raw displacement data firstly will be taken into dispRange3 function which subtracts
the offset of voltage measurement at that range from all data
=#
function dispRange3(rawDisp,dispOffset,range=3)
    for i in 1:length(rawDisp) 
        rawDisp[i]=rawDisp[i]-5-dispOffset*10-(range+1)*1e4
    end
    return rawDisp
end



##### Displacement Offset (dispOffset) Values for Z, X and Y axes:
z_dispOff = 131
x_dispOff = 129
y_dispOff = 128
###############################

#=
After range & DAQ adjustment to Data subtract the residual 
offset voltage by taking the first data point zero
=#
function offsetAllVolt(arrayInp,indexOffset=0,offset=0)
    for i in 1:length(arrayInp)
        arrayInp[i]=arrayInp[i]-arrayInp[indexOffset]-offset
    end
    return arrayInp
end

#=
Finally, create the voltage data to actual displacement by
applying calibration values.
=#
function totalDisp(dispVolt, calValue, range=3)
    dispNano = convert(Array{Union{Float64,Missing}}, dispVolt)
    for i in 1:length(dispVolt)
        dispNano[i]=(dispVolt[i]*calValue[i])/2^range
    end
    return dispNano
end

##### Modcal Factor 
z_modCal = 1.15e-8
x_modCal = 1.15e-8
y_modCal = 1.15e-8

###### DC Calibrations #####

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

### Range 0 Gain --bitwise--
r0_gain = 1

### Range 1 Gain
r1_gain = 2

### Range 2 Gain
r2_gain = 4

### Range 3 Gain
r3_gain = 8

### Range 4 Gain
r4_gain = 16

### Range 5 Gain
r5_gain = 32

### Range 6 Gain
r6_gain = 64

### Range 7 Gain
r7_gain = 128

### DC Load Frame Stiffness
z_frameStiffDC = 2.2e5
x_frameStiffDC = 2.2e5
y_frameStiffDC = 2.2e5

###### AC Calibrations #####

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
x_dispCalAC2 = -0.531442
y_dispCalAC2 = -0.0013593

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

##### DC Displacement Calibration Values (calValue) for  Z, X and Y axes:
z_dispCalDC = z_dispCalDC0 + z_dispCalDC1 * (z_dispOff - 128) + z_dispCalDC2 * (z_dispOff - 128)^2
x_dispCalDC = x_dispCalDC0 + x_dispCalDC1 * (x_dispOff - 128) + x_dispCalDC2 * (x_dispOff - 128)^2
y_dispCalDC = y_dispCalDC0 + y_dispCalDC1 * (y_dispOff - 128) + y_dispCalDC2 * (y_dispOff - 128)^2

##### AC Displacement Calibration Values (calValue) for  Z, X and Y axes:
z_dispCalAC = z_dispCalAC0 + z_dispCalAC1 * (z_dispOff - 128) + z_dispCalAC2 * (z_dispOff - 128)^2
x_dispCalAC = x_dispCalAC0 + x_dispCalAC1 * (x_dispOff - 128) + x_dispCalAC2 * (x_dispOff - 128)^2
y_dispCalAC = y_dispCalAC0 + y_dispCalAC1 * (y_dispOff - 128) + y_dispCalAC2 * (y_dispOff - 128)^2

##### DC Load Calibration Values (calValue) for  Z, X and Y axes:
z_loadCalDC = z_loadCalDC0 + z_loadCalDC1 * (z_dispOff - 128) + z_loadCalDC2 * (z_dispOff - 128)^2
x_loadCalDC = x_loadCalDC0 + x_loadCalDC1 * (x_dispOff - 128) + x_loadCalDC2 * (x_dispOff - 128)^2
y_loadCalDC = y_loadCalDC0 + y_loadCalDC1 * (y_dispOff - 128) + y_loadCalDC2 * (y_dispOff - 128)^2

##### AC Load Calibration Values (calValue) for  Z, X and Y axes:
z_loadCalAC = z_loadCalAC0 + z_loadCalAC1 * (z_dispOff - 128) + z_loadCalAC2 * (z_dispOff - 128)^2
x_loadCalAC = x_loadCalAC0 + x_loadCalAC1 * (x_dispOff - 128) + x_loadCalAC2 * (x_dispOff - 128)^2
y_loadCalAC = y_loadCalAC0 + y_loadCalAC1 * (y_dispOff - 128) + y_loadCalAC2 * (y_dispOff - 128)^2

#=
Similar function for range of load
=#
function loadRange3(rawLoadArray,range=3)
    rangedLoad = convert(Array{Union{Missing,Float64}}, rawLoadArray)
    for i in 1:length(rawLoadArray)
        rangedLoad[i] = rawLoadArray[i]-(range-1)*10-5
    end
    return rangedLoad
end

#=
Zero Point Determination; live GUI need here!
after determination return 3 displacement 3 load and time arrays with zeroed form
=#
function zeroDetermination(dispNano, rangedLoad)
    return zeroPt(dispNano, rangedLoad)
end

#=
Convert Load data in voltage to mN
=#
function loadConvert(nulledLoad, loadCalDC)
    loadMn = convert(Array{Union{Missing,Float64}},nulledLoad)
    for i in 1:length(nulledLoad)
        loadMn[i]=nulledLoad[i]*loadCalDC/1000
    end
    return loadMn
end

#=
subtract the load on the leaf spring by calculating it via
slope of load-disp at approach *times displacement[i] 
=#
function leafspringCorrection(nulledLoad, nulledDisp, appInterval)
    fApproach = interpolate((nulledDisp[appInterval],),nulledLoad[appInterval],Gridded(Linear()))
    slope =  mean(only.(gradient.(Ref(fApproach),nulledDisp[appInterval])))
    lfcorrectedLoad = convert(Array{Union{Missing, Float64}},nulledLoad)
    for i in 1 : length(nulledLoad)
        lfcorrectedLoad[i] = nulledLoad[i]-slope*nulledDisp[i]
    end
    return lfcorrectedLoad
end

#=
Thermal drift correction to displacement:
subtract from the "zeroed" displacement (depth into surface) 
slope of that displacement-time at the peak hold *times zeroed time[i]
=#
function thermaldriftCorrection(nulledDisp, nulledTime, holdInterval)
    fHold = interpolate((nulledTime[holdInterval],), nulledLoad[holdInterval], Gridded(Linear()))
    slope = mean(only.(gradient.(Ref(fHold), nulledTime)))
    tdcorrectedDisp = convert(Array{Union{Missing,Float64}},nulledDisp)
    for i in 1:length(nulledDisp)
        tdcorrectedDisp[i] = nulledDisp[i] - slope * nulledTime[i]
    end
    return tdcorrectedDisp
end

#=
Convert Amplitude in Voltage to nm 
=#
function amplitudeNm(rawAmplitude, dispCalAC, range=3)
    ampNm = convert(Array{Union{Missing,Float64}},rawAmplitude)
    for i in 1:length(rawAmplitude)
        ampNm[i] = rawAmplitude[i] * dispCalAC / 2^range
    end
    return ampNm
end

#=
Convert AC -harmonic force into μN from voltage
=#
function acforceUn(rawforceAC, loadCalAC, modcal)
    forceUn = convert(Array{Union{Missing,Float64}}, rawforceAC)
    for i in 1:length(rawforceAC)
        forceUn[i] = rawforceAC[i] * loadCalAC * modcal
    end
    return forceUn
end

#=
Phase Correction
=#
function phaseCorrection(rawPhase, co=1)
    correctedPhase = convert(Array{Unioin{Missing,Float64}},rawPhase)
    for i in 1:length(rawPhase)
        correctedPhase[i]=rawPhase[i]*co- (-0.0844 * 125 - 0.0596)
    end
    return correctedPhase
end

#Cos Phase
cosPhase(correctedPhase) = cos(deg2rad(correctedPhase))

#Sin Phase
sinPhase(correctedPhase) = sin(deg2rad(correctedPhase))

#=
Harmonic Force / Harmonic Displacement (for getting stiffness from DCM)
=#
function forcePerAmp(ampNm, forceUn)
    fPerA = convert(Array{Union{Missing,Float64}},amplitudeNm)
    for i in 1:length(fPerA)
        fPerA[i] = forceUn[i] / ampNm[i] * 1000 #*1000 for mN/μN
    end
    return fPerA
end

function stiffnessAC(fperA, correctedPhase, kLoadFrame)
    acStiff = convert(Array{Union{Missing, Float64}}, fperA)
    for i in 1:length(acStiff)
        acStiff[i] = 1/((1/(fperA[i]*cosPhase(phase[i]))) - (1/kLoadFrame))
    end
    return acStiff
end

function dampingAC(fperA, correctedPhase, acStiffness, dampingCoeff)
    acDamping = convert(Array{Union{Missing, Float64}}, fperA)
    for i in 1:length(acDamping)
        acDamping[i] = (fperA[i]*sinPhase(correctedPhase[i])-2*pi*dampingCoeff*125)/acStiffness[i]
    end
    return acDamping
end

function loadPerStiffnesSqr(lfCorrectedLoad, acStiffness)
    lPerSq =  convert(Array{Union{Missing, Float64}}, lfCorrectedLoad)
    for i in 1:length(lPerSq)
        lPerSq[i] = 1E3 * lfCorrectedLoad[i]/acStiffness[i]^2
    end
    return lPerSq
end

function displacementIntoSurf(tdCorrectedDisp, lfCorrectedLoad, lfStiffness)
    dispInSur = convert(Array{Union{Missing, Float64}}, tdCorrectedDisp)
    for i in 1:length(dispInSur)
        dispInSur[i] = tdCorrectedDisp[i] - (lfCorrectedLoad[i] * 0.001 / acStiffness[i]) * 1E9
    end
    return dispInSur
end

function contactDepth(dispInSur, lfCorrectedLoad, acStiff, coEpsilon)
    contDepth = convert(Array{Union{Missing, Float64}}, dispInSur)
    for i in 1:length(contDepth)
        contDepth[i] = dispInSur[i] - 1E9 * coEpsilon * (0.001*lfCorrectedLoad[i] / acStiff[i])
    end
    return contDepth
end

function contactArea(contDepth, coQuad, coLin, co0_5, co0_25, co0_125, co_0_06125, co0_030625, co0_0153125, co0_00765625)
    contArea = convert(Array{Union{Missing, Float64}}, contDepth)
    for i in 1:length(contArea)
        contArea[i] = coQuad*contDepth[i]^2 + coLin*contDepth[i] + co0_5*contDepth[i]^0.5 +
        co0_25*contDepth[i]^0.25 + co0_125*contDepth[i]^0.125 + co_0_06125*contDepth[i]^0.06125 + co0_030625*contDepth[i]^0.030625 +
        co0_0153125*contDepth[i]^0.0153125 + co0_00765625*contDepth[i]^0.00765625 
    end
    return contArea
end

function reducedModulus(acStiff, contArea)
    redModul = convert(Array{Union{Missing, Float64}}, acStiff)
    for i in 1:length(redModul)
        redModul[i] = pi^0.5 * acStiff[i]/contArea[i]^0.5
    end
    return redModul
end

function modulusSample(redModul, samplePoisson, indenterModul, indenterPoisson)
    sampleModul = convert(Array{Union{Missing, Float64}}, redModul)
    for i in 1:length(sampleModul)
        sampleModul[i] = (1-samplePoisson^2)/(1/redModul[i]-((1-indenterPoisson^2)/indenterModul))
    end
    return sampleModul
end

function hardness(lfCorrLoad, contArea)
    hardness = convert(Array{Union{Missing, Float64}}, lfCorrLoad)
    for i in 1:length(hardness)
        hardness[i] = lfCorrLoad[i] / contArea[i]
    end
    return hardness
end

function apce(contArea, acStiff)
    apce = convert(Array{Union{Missing, Float64}}, contArea)
    for i in 1:length(apce)
        apce[i] = 1E-18 * contArea[i]/acStiff[i]^2
    end
    return apce
end

## Constant E Area ? average(apce@drift) * z_stiffness^2
function cnstEarea()

end

function friction(lateralHarmonicForce, quasiStaticForce)
    friction = convert(Array{Union{Missing, Float64}}, lateralHarmonicForce)
    for i in 1:length(friction)
        friction[i] = lateralHarmonicForce[i] / (quasiStaticForce[i]*1000)
    end
    return friction
end


