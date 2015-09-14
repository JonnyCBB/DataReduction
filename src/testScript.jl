# include("ReciprocalSpaceUtils.jl")
# include("ElementDatabase.jl")
# include("MtzdumpHandling.jl")
# include("SequenceFileParser.jl")
# include("UpdateAtomAndRefs.jl")
# using Gadfly
# using LsqFit

######### Inputs ##########
const xrayEnergy = 12.7 #Set X-ray Energy
const xrayWavelength =  12.4/xrayEnergy # Sort out this conversion at some point
const integrationFileLocation = "integration_scaling_files\\pointless.mtz"
# const integrationFileLocation = "integration_scaling_files\\test450images.mtz"
const sequenceFileLocation = "SequenceFiles\\2BN3fasta.txt"
# const sequenceFileLocation = "SequenceFiles\\4X4Vfasta.txt"
const useSeqFile = true #Choose whether to use a sequence file to get variance and B-factor estimates
const separateSymEquivs = false #Merge symmetry equivalents or keep them separate.
const sigIDiffTol = 0.1 #Tolerance level for difference between sigIpr and sigIsum
const numOfRefs = 20000 #Number of reflections to be used in data reduction analysis.

const intensityType = "Combined" #How to deal with Ipr and Isum
const numMtzColsFor1stRefLine = 9 #Number of columns in 1st MTZ Dump line for reflection information
const numMtzColsFor2ndand3rdRefLines = 4 #Number of columns in 2nd/3rd MTZ Dump line for reflection information
const estimateTotalIntensityFromPartialRef = true #Estimate the total intensity from partial information.
const additionalElements = ""

const minRefInResBin = 50 #choose minimum number of reflections in resolution bin.
const minRefPerImage = 3
################################################################################
#Section: Inputs - Extract sequence information
#-------------------------------------------------------------------------------
atomDict = getAtomicCompositon(sequenceFileLocation)
additionalElements!(atomDict, additionalElements)
#End Section: Inputs - Extract sequence information
################################################################################

################################################################################
#Section: Inputs - Extract reflection information
#-------------------------------------------------------------------------------
#This section implements the methods to extract the integrated intensity
#information using MTZ Dump.
mtzdumpOutput = runMtzdump(integrationFileLocation, numOfRefs)
spacegroup, unitcell, hklList, imageArray = parseMTZDumpOutput(mtzdumpOutput)
#End Section: Inputs - Extract reflection information
################################################################################

################################################################################
#Section: Create intermediate Parameters
#-------------------------------------------------------------------------------
#Here we create the parameters required to later update information about the
#atoms and the reflections. These include
#1) The scattering angles of the reflections
#2) The scattering factors of each element
scatteringAngles = getAllScatteringAngles(hklList)
elementDict = createElementDictionary()
#End Section: Create intermediate Parameters
################################################################################

################################################################################
#Section: Sort the resolution Bins
#-------------------------------------------------------------------------------
resbins = sortResolutionBins(hklList, minRefInResBin)
sortHKLIntoResBins!(resbins, hklList)
#End Section: Sort the resolution Bins
################################################################################

################################################################################
#Section: Update atom and reflection information
#-------------------------------------------------------------------------------
#In this section we use the information gathered from both the atomic
#composition and the reflection information to update various parameter values.
calcElementf0!(elementDict, scatteringAngles, xrayWavelength)
updateAtomDict!(atomDict, spacegroup)
f0SqrdDict = calcTotalf0Sqrd(atomDict, scatteringAngles, elementDict)
updateRefListAndImageArray!(hklList, imageArray, estimateTotalIntensityFromPartialRef)

function calcResbinMeanIntensity!(resbins::Vector{ResBin}, f0SqrdDict::Dict{Float64, Float64}, hklList::Dict{Vector{Int64}, Reflection})
    for resbin in resbins
        summedIntensity = 0.0
        for hkl in resbin.hkls
            summedIntensity += f0SqrdDict[hklList[hkl].scatteringAngle]
        end
        resbin.meanIntensity = summedIntensity/resbin.numOfRef
    end
end
calcResbinMeanIntensity!(resbins, f0SqrdDict, hklList)

bfactors = Vector{Float64}()
imageNumArray = Vector{Int64}()
imageNum = 0
for img in imageArray
    imageNum += 1
    intensityDict = Dict{Float64,Float64}()
    for hkl in keys(img.observationList)
        reflection = hklList[hkl]
        intensityDict[reflection.resolution] = img.observationList[hkl].intensity/reflection.epsilon
    end
    logIntensityRatio = Vector{Float64}()
    hSqrd = Vector{Float64}()
    for resbin in resbins
        summedIntensity = 0.0
        numRefs = 0
        for resolution in keys(intensityDict)
            if resbin.minRes > resolution >= resbin.maxRes
                numRefs += 1
                summedIntensity += intensityDict[resolution]
            end
        end
        meanIntensity = summedIntensity/numRefs
        if !isnan(meanIntensity) && meanIntensity > 0.0
            push!(logIntensityRatio, log(meanIntensity/resbin.meanIntensity))
            midres = (resbin.maxRes + resbin.minRes)/2
            push!(hSqrd, 1/midres^2)
        end
    end
    if length(logIntensityRatio) >= minRefPerImage
        model(x, p) = p[1] + p[2] * x
        fit = curve_fit(model, hSqrd, logIntensityRatio, [0.0, 0.0])
        push!(bfactors, -2*fit.param[2])
        push!(imageNumArray, imageNum)
    end
    # xvals = linspace(0, maximum(hSqrd), 100)
    # yvals = model(xvals, fit.param)
    # plot(
    # layer(x=hSqrd, y=logIntensityRatio, Geom.point),
    # layer(x=xvals, y=yvals, Geom.line)
    # )
end
model(x,p) = p[1] + p[2]*x
fit = curve_fit(model, imageNumArray, bfactors,[0.0,0.0])
fit.param[2]
#End Section: Update atom and reflection information
################################################################################
