#########################################################################
####################### using/import statements #########################
using PyCall
using Gadfly
using LsqFit
using Colors
using Distributions
using KernelDensity
import Gadfly.ElementOrFunction

# include("ReciprocalSpaceUtils.jl")
# include("ElementDatabase.jl")
# include("MtzdumpHandling.jl")
# include("SequenceFileParser.jl")
# include("UpdateAtomAndRefs.jl")
# include("FilteringUtils.jl")

######### Inputs ##########
const xrayEnergy = 12.7 #Set X-ray Energy
const xrayWavelength =  12.4/xrayEnergy # Sort out this conversion at some point
const integrationFileLocation = "integration_scaling_files\\pointless.mtz"
# const integrationFileLocation = "integration_scaling_files\\test450images.mtz"
const sequenceFileLocation = "SequenceFiles\\2BN3fasta.txt"
# const sequenceFileLocation = "SequenceFiles\\4X4Vfasta.txt"
const sfFileLocation = "integration_scaling_files\\test450images_scaled1.mtz"
const useSeqFile = true #Choose whether to use a sequence file to get variance and B-factor estimates
const separateSymEquivs = false #Merge symmetry equivalents or keep them separate.
const sigIDiffTol = 0.1 #Tolerance level for difference between sigIpr and sigIsum
const numOfRefs = 20000 #Number of reflections to be used in data reduction analysis.

const intensityType = "Combined" #How to deal with Ipr and Isum
const numMtzColsFor1stRefLine = 9 #Number of columns in 1st MTZ Dump line for reflection information
const numMtzColsFor2ndand3rdRefLines = 4 #Number of columns in 2nd/3rd MTZ Dump line for reflection information
const numMtzColsIntLineCTruncate = 6
const estimateTotalIntensityFromPartialRef = true #Estimate the total intensity from partial information.
const additionalElements = ""

const minRefInResBin = 50 #choose minimum number of reflections in resolution bin.
const minRefPerImage = 3
const displayPlots = false

const minFracCalc = 0.95
const applyBFacTof0 = true

const kdeStep = 0.0001
const keepPercentageScaleData = 0.9

const outputImageDir = "plots"

const processVarCoeff = 1.0
const observationVarCoeff = 1.0
################################################################################
#Section: Create plot directory
#-------------------------------------------------------------------------------
if !isempty(outputImageDir)
    if !isdir(outputImageDir)
        mkdir(outputImageDir)
    end
end
#End Section: Create plot directory
################################################################################

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
spacegroup, unitcell, hklList, imageArray = parseMosflmMTZDumpOutput(mtzdumpOutput)
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
calcResbinMeanIntensity!(resbins, f0SqrdDict, hklList)
changeInBfac, bGradSigma, bIntercept, bInterceptSigma, modalScale, sigmaScale = calcBandScaleParams(hklList, imageArray, resbins, outputImageDir, displayPlots)
#End Section: Update atom and reflection information
################################################################################


################################################################################
#Section: Inflate observation errors
#-------------------------------------------------------------------------------
#In this section we inflate the sigma values of the observed intensities
#according to their total calculated fraction values.
#Basically if the calculated intensity fraction is not close enough to 1
#then this means that the true observed intensity measurement has not
#been fully measured. Rather than estimating this true observed
#intensity we instead inflate the sigma value for the reflection. This
#basically means that we're increasing our uncertainty about the
#intensity measurement rather than trying to deterministically give an
#estimate of the true intensity.
tempFacDict, SFMultiplierDict = calcTempAndSFMultFactorDict(scatteringAngles, bIntercept, changeInBfac, xrayWavelength)
inflateObservedSigmas!(imageArray, hklList, changeInBfac, minFracCalc, applyBFacTof0)
#End Section: Inflate observation errors
################################################################################

################################################################################
#Section: Extract initial guess structure factor amplitudes
#-------------------------------------------------------------------------------
# Get initial amplitudes by method 1
# getInitialAmplitudes!(hklList, atomDict, scatteringAngles, elementDict, tempFacDict)

# Get initial amplitudes by method 2
getInitialAmplitudes!(hklList, f0SqrdDict, tempFacDict)

# Get initial amplitudes by method 3
# mtzDumpOutput = runMtzdump(sfFileLocation, 1200)
# refAmpDict, scaleFac = parseCTruncateMTZDumpOutput(mtzDumpOutput)
# getInitialAmplitudes!(hklList, refAmpDict, scaleFac)

#End Section: Extract initial guess structure factor amplitudes
################################################################################

########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
###########################################################################

################################################################################
#Section: Set up initial guess for the amplitudes
#-------------------------------------------------------------------------------
initialStateEstimate = getInitialState(hklList, f0SqrdDict)
#End Section: Extract initial guess structure factor amplitudes
################################################################################


################################################################################
#Section: Perform filtering
#-------------------------------------------------------------------------------
######## Create an index reference for the miller indices of the reflections.
hklIndexReference = createHKLIndexReferenceDict(hklList)

img = imageArray[100]
observationVector = Vector{Float64}(length(hklList))
processVarVector = Vector{Float64}(length(hklList))
observationVarVector = Vector{Float64}(length(hklList))
SFMultiplierVec = Vector{Float64}(length(hklList))
SFSigmaVec = Vector{Float64}(length(hklList))
observationIndices = Vector{Int64}()
scaleFactor = modalScale
hklCounter = 0
for hkl in keys(hklList)
    hklCounter += 1
    reflection = hklList[hkl]
    D = SFMultiplierDict[reflection.scatteringAngle]
    Σ = f0SqrdDict[reflection.scatteringAngle]
    σ = sqrt(abs(1 - D^2)*Σ)
    ############################################################################
    #NEED TO SORT THIS OUT. I'VE HAD TO USE A GAUSSIAN DISTRIBUTION INSTEAD OF A
    #RICIAN DISTRIBUTION.
    #processVarVector[hklCounter] = varRice(F, D, σ)
    processVarVector[hklCounter] = σ^2 * processVarCoeff
    SFMultiplierVec[hklCounter] = D
    SFSigmaVec[hklCounter] = σ
    if haskey(img.observationList, hkl)
        push!(observationIndices, hklCounter)
        observationVector[hklCounter] = img.observationList[hkl].intensity
        observationVarVector[hklCounter] = img.observationList[hkl].sigI
    else
        observationVector[hklCounter] = NaN
        observationVarVector[hklCounter] = (2 * scaleFactor * processFunction(reflection.amplitude, D, σ))^2 * processVarVector[hklCounter] * observationVarCoeff
    end
end
processFunction(amplitudes::Vector{Float64}) = processFunction(amplitudes, SFMultiplierVec, SFSigmaVec)
observationFunction(amplitudes::Vector{Float64}) = observationFunction(amplitudes, scaleFactor)



#End Section: Perform Filtering
################################################################################
