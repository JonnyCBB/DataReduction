################################################################################
#NEED TO ADD METHOD INFORMATION
################################################################################
function updateAtomDict!(atomDict::Dict{ASCIIString, Int64}, spacegroup::SpaceGroup)
    for element in keys(atomDict)
        atomDict[element] *= spacegroup.numSymOps
    end
end

################################################################################
#NEED TO ADD METHOD INFORMATION
################################################################################
function calcTotalf0Sqrd(atomDict::Dict{ASCIIString, Int64}, scatteringAngles::Vector{Float64}, elementDict::Dict{ASCIIString, Element})
    f0SqrdDict = Dict{Float64, Float64}()
    for scatAngle in scatteringAngles
        f0::Float64 = 0.0
        for atom in keys(atomDict)
            f0 += elementDict[atom].f0[scatAngle]^2 * atomDict[atom]
        end
        f0SqrdDict[scatAngle] = f0
    end
    return f0SqrdDict
end

################################################################################
#NEED TO ADD METHOD INFORMATION
################################################################################
function updateRefListAndImageArray!(hklList::Dict{Vector{Int64},Reflection}, imageArray::Vector{DiffractionImage}, estimatePartialIntensity::Bool=true)
    for hkl in keys(hklList)
        reflection = hklList[hkl]
        numOfExistingObs = length(reflection.observations)
        #Loop through each observation
        for obsNum = 1:numOfExistingObs
            ########################################################################
            #Section: Update observation information
            #-----------------------------------------------------------------------
            #Finish calculating the rotation centroid and the sigI of the
            #observation and update the observation info in the reflection list.
            refObservation = reflection.observations[obsNum] # get the observation object
            numPartials = length(refObservation.imageNums) # calculate the number of partial observations
            refObservation.rotCentroid = refObservation.rotCentroid/numPartials # divide by the rotCentroid number of partial observations to get the average.
            refObservation.sigI = sqrt(refObservation.sigI) # Square the sum of the variances to get the total standard deviation.
            reflection.observations[obsNum] = refObservation # return updated observation information for the reflection.
            #End Section: Update observation information
            ########################################################################

            ########################################################################
            #Section: Add observations to images - Fully observed reflections
            #-----------------------------------------------------------------------
            #In this section we search for the image on which the current
            #observation was "seen". This is determined by the image which contains
            #the rotation centroid of the reflection.
            foundCentroidImage = false #Check whether an image containing the centroid of the reflection has been found.
            for imageNum in refObservation.imageNums #loop through the images
                diffractionImage = imageArray[imageNum]
                if diffractionImage.rotAngleStart <= refObservation.rotCentroid < diffractionImage.rotAngleStop
                    diffractionImage.observationList[hkl] = refObservation
                    foundCentroidImage = true
                    break
                elseif numPartials == 1
                    if imageNum != 1 && imageNum != length(imageArray)
                        diffractionImage.observationList[hkl] = refObservation
                        foundCentroidImage = true
                        @printf("*************************WARNING**************************\n")
                        @printf("Observation %d of Reflection [%d,%d,%d] has been allocated to image %d which is not the correct image.\n",obsNum, hkl[1], hkl[2], hkl[3],imageNum)
                        @printf("Jonny you need to properly sort this out!!!\n")
                        @printf("The reflection centroid isn't actually on the image. This is just your crappy work around to continue work.\nSORT IT!!!\n\n")
                    end
                end
            end
            #End Section: Add observations to images - Fully observed reflections
            ########################################################################

            ########################################################################
            #Section: Add observations to images - Partially observed reflections
            #-----------------------------------------------------------------------
            #In this section we find the observations of reflections that were only
            #seen partially e.g. an observation where we only sampled the last half
            #of it's intensity on the first image. We then estimate the fraction of
            #of the reflection that was sampled and estimate the total intensity
            #from that fraction.
            #NOTE OF CAUTION:
            #The current fraction estimation is VERY crude. It erroneously uses the
            #phi angles in place of what should be true path lengths in reciprocal
            #space.
            if !foundCentroidImage #Check if the observation has been allocated to an image
                imageNums = sort(refObservation.imageNums) #Get the image numbers in order
                #Check whether the image was only partially observed on the first
                #image or whether it was partially observed on the first image.
                if imageNums[1] == 1
                    imageArray[1].observationList[hkl] = refObservation #Add reflection to the reflection list for the first image

                    #This is a crude (pretty rubbish) method to estimate the
                    #penetration of the reciprocal lattice sphere.
                    lastImageOfObs = imageArray[imageNums[end]]
                    phiImage = (lastImageOfObs.rotAngleStart + lastImageOfObs.rotAngleStop)/2
                    phiEdge = imageArray[1].rotAngleStart
                    spherePenetration = abs(phiImage - phiEdge)
                elseif imageNums[end] == length(imageArray)
                    imageArray[end].observationList[hkl] = refObservation #Add reflection to the reflection list for the last image

                    #This is a crude (pretty rubbish) method to estimate the
                    #penetration of the reciprocal lattice sphere.
                    firstImageOfObs = imageArray[imageNums[1]]
                    phiImage = (firstImageOfObs.rotAngleStart + firstImageOfObs.rotAngleStop)/2
                    phiEdge = imageArray[end].rotAngleStop
                    spherePenetration = abs(phiImage - phiEdge)
                else
                    #If the observation wasn't observed on either the first or last
                    #image and wasn't allocated any other image in between then
                    #there is a problem. It could be a problem with the input data
                    #or with the current program. Not sure which so we'll alert the
                    #user with a warning message.
                    println("************************WARNING************************")
                    @printf("Observation %d of reflection (%d,%d,%d) could not be allocated to an image\nPlease check that the rotation centroid of this reflection is a valid number.\nThe centroid calculated was %.3f.\nThe images checked to determine image for the centroid were:\n",obsNum, hkl[1], hkl[2], hkl[3], refObservation.rotCentroid)
                    for imageNum in refObservation.imageNums
                        @printf("image number: %d. Rotation Start and stop: %.3f deg - %.3f deg\n", imageNum, imageArray[imageNum].rotAngleStart, imageArray[imageNum].rotAngleStop)
                    end
                    @printf("Contact elspeth.garman@bioch.ox.ac.uk to sort out the MTZ Dump parser.\n\n")
                end

                #Only try to scale up intensity values if they are postive.
                if estimatePartialIntensity && refObservation.intensity > 0
                    sphereDiameter = 2*abs(refObservation.rotCentroid - phiImage) #estimate the diameter of the reciprocal lattice sphere.
                    spherePathFrac = spherePenetration/sphereDiameter #calculate partiality fraction
                    intensityFraction = volumeRatio(spherePathFrac) #estimate the fraction of the intensity that was measured
                    refObservation.intensity = refObservation.intensity/intensityFraction #estimate total intensity
                    refObservation.sigI = refObservation.sigI/intensityFraction #scale the sigma by the corresponding value
                    reflection.observations[obsNum] = refObservation # return updated observation information for the reflection.
                end
            end
            #End Section: Add observations to images - Partially observed reflections
            ########################################################################
        end
        hklList[hkl] = reflection # update the reflection in the reflection list
    end
end


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
function inflateObservedSigmas!(imageArray::Vector{DiffractionImage}, hklList::Dict{Vector{Int64},Reflection}, ΔB::Float64, minFracCalc::Float64=0.99, applyBFacTof0::Bool=true)
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
end
#End Section: Inflate observation errors
########################################################################
