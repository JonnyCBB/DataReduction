#########################################################################
####################### using/import statements #########################
# using PyCall
# using Gadfly
# using LsqFit
# using Colors
# using Compose

# include("ReciprocalSpaceUtils.jl")
# include("ElementDatabase.jl")
# include("MtzdumpHandling.jl")
# include("SequenceFileParser.jl")
# include("UpdateAtomAndRefs.jl")

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
const displayBfacPlot = false

const minFracCalc = 0.95
const applyBFacTof0 = true

const outputImageDir = "plots"
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
calcResbinMeanIntensity!(resbins, f0SqrdDict, hklList)
changeInBfac, midBfac = calcBfactor(hklList, imageArray, resbins, outputImageDir, displayBfacPlot)
#End Section: Update atom and reflection information
################################################################################


########################################################################
#Section: Inflate observation errors
#-----------------------------------------------------------------------
#In this section we inflate the sigma values of the observed intensities
#according to their total calculated fraction values.
#Basically if the calculated intensity fraction is not close enough to 1
#then this means that the true observed intensity measurement has not
#been fully measured. Rather than estimating this true observed
#intensity we instead inflate the sigma value for the reflection. This
#basically means that we're increasing our uncertainty about the
#intensity measurement rather than trying to deterministically give an
#estimate of the true intensity.
tempFacDict, SFMultiplierDict = calcTempAndSFMultFactorDict(scatteringAngles, midBfac, changeInBfac, xrayWavelength)

calcD(ΔB::Float64, θ::Float64, λ::Float64) = exp(-2 * ΔB * (sin(deg2rad(θ))^2) / λ^2)

ΔB = changeInBfac
rotStart = imageArray[1].rotAngleStart
rotEnd = imageArray[end].rotAngleStop
ϕ = imageArray[1].rotAngleStop - rotStart
for diffractionImage in imageArray
    #Loop through each observation
    for hkl in keys(diffractionImage.observationList)
        centroidImageDiffErrFactor = 0.0
        reflection = hklList[hkl]
        refObservation = diffractionImage.observationList[hkl] # get the observation object
        imageNumAndInts = hcat(refObservation.imageNums, refObservation.imageIntensities)
        imageNumAndInts = sortrows(imageNumAndInts)
        imageNums = sort(refObservation.imageNums)
        numPartials = length(imageNums)
        firstImage = imageArray[imageNums[1]]
        lastImage = imageArray[imageNums[end]]

        ########################################################################
        #Mini Section: Calculate error factor - Partials observed on different images
        #-----------------------------------------------------------------------
        #In this section we have to take into account our uncertainty due to the
        #fact that partial observations of reflections are observed on different
        #images. The variances that we want to add are '(1 - Dⱼ^2) * Ip' Where Dⱼ
        #is the structure factor multiplier, j represents the number of images
        #away from the image on which the centroid of the reflection was found,
        #Ip represents the (partial) intensity of the reflection that was observed
        #on that image.
        #The definition of j may be confusing so here's an example: lets say we
        #observe a reflection on images 5, 6, 7, 8, 9, 10. Suppose the centroid
        #of the reflection was observed on image 8. Then when we're considering
        #the value Dⱼ for image 10 we get that j = 10 - 8 = 2.

        #We also have to consider the cases when the centroids are not observed
        #on any of the given images (i.e. for reflections that are only partly
        #observed). That's what the if statment is about below.
        if refObservation.rotCentroid < firstImage.rotAngleStart
            if imageNums[end] > 1
                numImagesAboveCentroid = Int(floor(abs(firstImage.rotAngleStop - imageArray[imageNums[end]].rotAngleStop)/ϕ))
                for i in 1:numImagesAboveCentroid
                    partialIntensity = imageNumAndInts[numPartials - numImagesAboveCentroid + i, 2]
                    centroidImageDiffErrFactor += abs(1 - calcD(i * ΔB, reflection.resolution, xrayWavelength)^2) * partialIntensity
                end
            end
        elseif refObservation.rotCentroid > lastImage.rotAngleStop
            if imageNums[1] < length(imageArray)
                numImagesBelowCentroid = Int(floor(abs(lastImage.rotAngleStart - imageArray[imageNums[1]].rotAngleStart)/ϕ))
                for i in 1:numImagesBelowCentroid
                    partialIntensity = imageNumAndInts[numImagesBelowCentroid + 1 - i]
                    centroidImageDiffErrFactor += abs(1 - calcD(i * ΔB, reflection.resolution, xrayWavelength)^2) * partialIntensity
                end
            end
        else
            if length(imageNums) > 1
                numImagesAboveCentroid = Int(floor(abs(refObservation.rotCentroid - imageArray[imageNums[end]].rotAngleStop)/ϕ))
                numImagesBelowCentroid = Int(floor(abs(refObservation.rotCentroid - imageArray[imageNums[1]].rotAngleStart)/ϕ))
                for i in 1:numImagesAboveCentroid
                    partialIntensity = imageNumAndInts[numPartials - numImagesAboveCentroid + i, 2]
                    centroidImageDiffErrFactor += abs(1 - calcD(i * ΔB, reflection.resolution, xrayWavelength)^2) * partialIntensity
                end
                for i in 1:numImagesBelowCentroid
                    partialIntensity = imageNumAndInts[numImagesBelowCentroid + 1 - i]
                    centroidImageDiffErrFactor += abs(1 - calcD(i * ΔB, reflection.resolution, xrayWavelength)^2) * partialIntensity
                end
            end
        end
        refObservation.imageIntensities = Array(Float64,0) #Clear the memory since we don't need this variable anymore
        #End Mini Section: Calculate error factor - Partials observed on different images
        ########################################################################

        ########################################################################
        #Mini Section: Calculate error factor - Incomplete fraction calculated.
        #-----------------------------------------------------------------------
        #In this section we consider our uncertainty of the intensity due to the
        #fact that the integration program tells us that the fraction of the
        #reflection calculated is less than 1.
        fracCalcErrFactor = 0.0
        if refObservation.fractionCalc < minFracCalc
            fracCalcErrFactor = 1.0 - refObservation.fractionCalc
        end
        #End Mini Section: Calculate error factor - Incomplete fraction calculated.
        ########################################################################

        ########################################################################
        #Mini Section: Inflate the uncertainty.
        #-----------------------------------------------------------------------
        #In this section we recalculate our uncertainties i.e. the measured sigmas.
        #Since it's the variances that add (not the sigmas), we need to square the
        #current sigma value before adding it. We can simply add the other factors.
        if applyBFacTof0
            refObservation.sigI = sqrt(refObservation.sigI^2 + centroidImageDiffErrFactor + fracCalcErrFactor * reflection.epsilon * f0SqrdDict[reflection.scatteringAngle] * tempFacDict[reflection.scatteringAngle])
        else
            refObservation.sigI = sqrt(refObservation.sigI^2 + centroidImageDiffErrFactor + fracCalcErrFactor * reflection.epsilon * f0SqrdDict[reflection.scatteringAngle])
        end
        diffractionImage.observationList[hkl] = refObservation # Update the reflection observation
        #End Mini Section: Inflate the uncertainty.
        ########################################################################
    end
end
#End Section: Inflate observation errors
########################################################################
