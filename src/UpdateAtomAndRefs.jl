################################################################################
#NEED TO ADD METHOD INFORMATION
################################################################################
function updateAtomDict!(atomDict::Dict{ASCIIString, UInt32}, spacegroup::SpaceGroup)
    for element in keys(atomDict)
        atomDict[element] *= spacegroup.numSymOps
    end
end

################################################################################
#NEED TO ADD METHOD INFORMATION
################################################################################
function calcTotalf0Sqrd(atomDict::Dict{ASCIIString, UInt32}, scatteringAngles::Vector{Float32}, elementDict::Dict{ASCIIString, Element})
    f0SqrdDict = Dict{Float32, Float32}()
    for scatAngle in scatteringAngles
        f0::Float32 = 0.0
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
#This p(q) function is from Rossman et al. 1979. My volumeRatio function is an
#exact volume ratio for the partial reflection. However if you plot these
#functions it turns out that they pretty much overlap so they're pretty much
#equivalent for argument values between 0 - 1
p(q) = 3q^2 - 2q^3
function volumeRatio(penetrationDepth::Float32, sphereDiameter::Float32)
    d = sphereDiameter
    h = penetrationDepth
    ratio = 8 * h^2 * ((3/2)*d - h)/(4 * d^3)
    return ratio
end
volumeRatio(penetrationFraction::Float32) = volumeRatio(penetrationFraction, Float32(1.0))

################################################################################
#NEED TO ADD METHOD INFORMATION
################################################################################
function updateRefListAndImageArray!(hklList::Dict{Vector{Int16},Reflection},
    imageArray::Vector{DiffractionImage},
    estimatePartialIntensity::Bool=true,
    USE_LOW_FRAC_REF::Bool=false,
    USE_LOW_FRAC_CALC_IFF_SINGLE_OBS::Bool=true)
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
            refObservation.sigI = sqrt(refObservation.sigI) # Square root the sum of the variances to get the total standard deviation.
            reflection.observations[obsNum] = refObservation # return updated observation information for the reflection.
            #End Section: Update observation information
            ########################################################################

            addObsToImage = true # Assume we're initially happy to add the observation to an image.

            #Now we need to check whether we want to add the observation to any
            #of the images. We first check if the user actually wants to use
            #reflections where the calculated fraction is low. If not then check
            #that the calculated fraction of the current observation is below
            #the threshold. If both of these conditions are satisfied then we
            #need to check whether the user would like to use the observation if
            #it is the only observation for that reflection. If there is only
            #one observation then we will use it. Otherwise we'll only use it if
            #the current observation has a greater fraction calculated than all
            #of the other observations.
            if !USE_LOW_FRAC_REF && refObservation.fractionCalc < minFracCalc
                if USE_LOW_FRAC_CALC_IFF_SINGLE_OBS
                    if numOfExistingObs > 1
                        for otherObs in reflection.observations
                            if otherObs.fractionCalc > refObservation.fractionCalc
                                addObsToImage = false
                                break
                            end
                        end
                    end
                else
                    addObsToImage = false
                end
            end

            if addObsToImage
                ########################################################################
                #Section: Add observations to images - Fully observed reflections
                #-----------------------------------------------------------------------
                #In this section we search for the image on which the current
                #observation was "seen". This is determined by the image which
                #containsvthe rotation centroid of the reflection.
                foundCentroidImage = false # reflection initially hasn't been allocated to an image.
                for imageNum in refObservation.imageNums #loop through the images
                    if imageNum > length(imageArray)
                        println(hkl)
                        error("Image number exceeds image array bounds")
                    end
                    diffractionImage = imageArray[imageNum]
                    if diffractionImage.rotAngleStart <= refObservation.rotCentroid < diffractionImage.rotAngleStop
                        #Check to see if the current image already has an observation of a given reflection.
                        #If the reflection hasn't been observed on the diffraction image before then we can
                        #simply add the reflection to the list. If the reflection has already an observation on
                        #that image (e.g. a symmetry equivalent) then we can assume that they should have the same
                        #true intensity but they differ due to measurement error. Therefore the intensity of the
                        #reflection will be the average of the symmetry equivalents.
                        if !haskey(diffractionImage.observationList, hkl)
                            diffractionImage.observationList[hkl] = refObservation
                        else
                            duplicateObservation = diffractionImage.observationList[hkl]
                            #If the fraction of the reflection calculated by the integration software is above a certain
                            #amount then combine the data for the multiple observations by performing weigthed averages
                            #of the corresponding data.
                            if refObservation.fractionCalc >= minFracCalc
                                newIntensity = (duplicateObservation.fractionCalc * duplicateObservation.intensity + refObservation.fractionCalc * refObservation.intensity)/(duplicateObservation.fractionCalc + refObservation.fractionCalc)

                                newSigma = sqrt(duplicateObservation.fractionCalc^2 * duplicateObservation.sigI^2 + refObservation.fractionCalc^2 * refObservation.sigI^2/(duplicateObservation.fractionCalc^2 + refObservation.fractionCalc^2))

                                newImageNums = [duplicateObservation.imageNums; refObservation.imageNums]
                                newImageIntensities = [duplicateObservation.imageIntensities; refObservation.imageIntensities]
                                newFractionCalc = (duplicateObservation.fractionCalc + refObservation.fractionCalc)/2

                                newRotCentroid = (duplicateObservation.fractionCalc * duplicateObservation.rotCentroid + refObservation.fractionCalc * refObservation.rotCentroid)/(duplicateObservation.fractionCalc + refObservation.fractionCalc)

                                ####################################################
                                #THE MYISYM SHOULD BE SORTED OUT!!!
                                ####################################################
                                newMISYM = duplicateObservation.misym

                                #Create the new combined observation object.
                                diffractionImage.observationList[hkl] = ReflectionObservation(newRotCentroid, newFractionCalc, newMISYM, newIntensity, newSigma, newImageNums, newImageIntensities)
                            end
                        end
                        foundCentroidImage = true
                        break
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
                    allocatedToImage = false
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
                        allocatedToImage = true
                    elseif imageNums[end] == length(imageArray)
                        imageArray[end].observationList[hkl] = refObservation #Add reflection to the reflection list for the last image

                        #This is a crude (pretty rubbish) method to estimate the
                        #penetration of the reciprocal lattice sphere.
                        firstImageOfObs = imageArray[imageNums[1]]
                        phiImage = (firstImageOfObs.rotAngleStart + firstImageOfObs.rotAngleStop)/2
                        phiEdge = imageArray[end].rotAngleStop
                        spherePenetration = abs(phiImage - phiEdge)
                        allocatedToImage = true
                    else
                        #If the observation wasn't observed on either the first or last
                        #image and wasn't allocated any other image in between then
                        #there is a problem. It could be a problem with the input data
                        #or with the current program. Not sure which so we'll alert the
                        #user with a warning message.
                        #What we'll do here is just allocate it to the image
                        #corresponding to the calculated reflection centroid.
                        imageCounter = 0
                        for diffractionImage in imageArray
                            imageCounter += 1
                            ########################################################
                            #THIS BLOCK OF CODE IS DUPLICATED FROM THE SECTION ABOVE
                            #THIS SHOULD BE WRITTEN AS A SEPARATE FUNCTION.
                            if diffractionImage.rotAngleStart <= refObservation.rotCentroid < diffractionImage.rotAngleStop
                                #Check to see if the current image already has an observation of a given reflection.
                                #If the reflection hasn't been observed on the diffraction image before then we can
                                #simply add the reflection to the list. If the reflection has already an observation on
                                #that image (e.g. a symmetry equivalent) then we can assume that they should have the same
                                #true intensity but they differ due to measurement error. Therefore the intensity of the
                                #reflection will be the average of the symmetry equivalents.
                                if !haskey(diffractionImage.observationList, hkl)
                                    diffractionImage.observationList[hkl] = refObservation
                                    allocatedToImage = true

                                    println("************************WARNING************************")
                                    @printf("The rotation centroid for observation %d of reflection (%d,%d,%d) was calculated (by the integration software) to be outside of an image on which the observation was observed.\nPlease check that the rotation centroid of this reflection is a valid number.\nThe centroid calculated was %.2f. (Note: The estimated fraction of the observation that was measured is: %.2f).\nThe images for which the observation was measured were:\n",obsNum, hkl[1], hkl[2], hkl[3], refObservation.rotCentroid, refObservation.fractionCalc)
                                    for imageNum in refObservation.imageNums
                                        @printf("image number: %d. Rotation start and stop: %.2f deg - %.2f deg\n", imageNum, imageArray[imageNum].rotAngleStart, imageArray[imageNum].rotAngleStop)
                                    end
                                    @printf("The observation centroid was allocated to image: %d.\nRotation start and stop: %.2f deg - %.2f deg\n\n", imageCounter, imageArray[imageCounter].rotAngleStart, imageArray[imageCounter].rotAngleStop)
                                else
                                    duplicateObservation = diffractionImage.observationList[hkl]
                                    #If the fraction of the reflection calculated by
                                    #the integration software is above a certain
                                    #threshold then combine the data for the multiple
                                    #observations by performing weigthed averages of
                                    # the corresponding data.
                                    if refObservation.fractionCalc >= minFracCalc
                                        newIntensity = (duplicateObservation.fractionCalc * duplicateObservation.intensity + refObservation.fractionCalc * refObservation.intensity)/(duplicateObservation.fractionCalc + refObservation.fractionCalc)

                                        newSigma = sqrt(duplicateObservation.fractionCalc^2 * duplicateObservation.sigI^2 + refObservation.fractionCalc^2 * refObservation.sigI^2/(duplicateObservation.fractionCalc^2 + refObservation.fractionCalc^2))

                                        newImageNums = [duplicateObservation.imageNums; refObservation.imageNums]
                                        newImageIntensities = [duplicateObservation.imageIntensities; refObservation.imageIntensities]
                                        newFractionCalc = (duplicateObservation.fractionCalc + refObservation.fractionCalc)/2

                                        newRotCentroid = (duplicateObservation.fractionCalc * duplicateObservation.rotCentroid + refObservation.fractionCalc * refObservation.rotCentroid)/(duplicateObservation.fractionCalc + refObservation.fractionCalc)

                                        ####################################################
                                        #THE MYISYM SHOULD BE SORTED OUT!!!
                                        ####################################################
                                        newMISYM = duplicateObservation.misym

                                        #Create the new combined observation object.
                                        diffractionImage.observationList[hkl] = ReflectionObservation(newRotCentroid, newFractionCalc, newMISYM, newIntensity, newSigma, newImageNums, newImageIntensities)
                                        allocatedToImage = true

                                        println("************************WARNING************************")
                                        @printf("The rotation centroid for observation %d of reflection (%d,%d,%d) was calculated (by the integration software) to be outside of an image on which the observation was observed.\nPlease check that the rotation centroid of this reflection is a valid number.\nThe centroid calculated was %.2f. (Note: The estimated fraction of the observation that was measured is: %.2f).\nThe images for which the observation was measured were:\n",obsNum, hkl[1], hkl[2], hkl[3], refObservation.rotCentroid, refObservation.fractionCalc)
                                        for imageNum in refObservation.imageNums
                                            @printf("image number: %d. Rotation start and stop: %.2f deg - %.2f deg\n", imageNum, imageArray[imageNum].rotAngleStart, imageArray[imageNum].rotAngleStop)
                                        end
                                        @printf("The observation centroid was allocated to image: %d.\nRotation start and stop: %.2f deg - %.2f deg\n\n", imageCounter, imageArray[imageCounter].rotAngleStart, imageArray[imageCounter].rotAngleStop)
                                    end
                                end
                                break
                            end
                        end
                        if imageNums[end] < imageCounter
                            phiEdge = imageArray[imageNums[end]].rotAngleStop
                            phiImage = (imageArray[imageNums[1]].rotAngleStart + imageArray[imageNums[1]].rotAngleStop)/2
                            spherePenetration = abs(phiImage - phiEdge)
                        elseif imageNums[1] > imageCounter
                            phiEdge = imageArray[imageNums[1]].rotAngleStart
                            phiImage = (imageArray[imageNums[end]].rotAngleStart + imageArray[imageNums[end]].rotAngleStop)/2
                            spherePenetration = abs(phiImage - phiEdge)
                        else
                            errMsg = @sprintf("This means that the reflection was observed on the image which coincides with it's centroid, yet hasn't already been allocated an image.\nIf you're seeing this error message then something is VERY wrong with the code and you should contact elspeth.garman@bioch.ox.ac.uk ASAP this problem can be sorted.")
                            error(errMsg)
                        end
                    end

                    #If the observation can't be allocated to an image despite having a measured fraction above the threshold then something is wrong (I don't know what but something is) and it needs to be sorted.
                    if !allocatedToImage && refObservation.fractionCalc >= minFracCalc
                        errMsg = @sprintf("Observation %d of reflection (%d,%d,%d) couldn't be allocated to an image.\nThe centroid for this reflection was %.2f and the fraction calculated was %.2f.\nContact elspeth.garman@bioch.ox.ac.uk so this problem can be sorted.\n", obsNum, hkl[1], hkl[2], hkl[3], refObservation.rotCentroid, refObservation.fractionCalc)
                        error(errMsg)
                    end

                    #Only try to scale up intensity values if they are postive, have
                    #been allocated to an image and the fraction calculated is below
                    #the threshold.
                    if estimatePartialIntensity && refObservation.intensity > 0 && allocatedToImage
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
        end
        hklList[hkl] = reflection # update the reflection in the reflection list
    end
end


########################################################################
#Section: Inflate observation errors
#-----------------------------------------------------------------------
#In this section we inflate the sigma values of the observed intensities
#according to their total calculated fraction values (1) and also account for
#errors in making full reflection intensity estimates from partially observed
#data (2).
#(1) If the calculated intensity fraction is not close enough to 1
#then this means that the true observed intensity measurement has not
#been fully measured. Rather than estimating this true observed
#intensity we instead inflate the sigma value for the reflection. This
#means that we're increasing our uncertainty about the intensity measurement
#rather than trying to deterministically give an estimate of the true
#intensity.
#(2) If we combine partials we have to account for the fact that they were
#observed at different times i.e. on different images.
function inflateObservedSigmas!(imageArray::Vector{DiffractionImage}, hklList::Dict{Vector{Int16},Reflection}, ΔB::Float32, minFracCalc::Float32=Float32(0.99), applyBFacTof0::Bool=true)
    rotStart = imageArray[1].rotAngleStart # Phi angle at beginning of experiment
    rotEnd = imageArray[end].rotAngleStop # Phi angle at end of experiment.
    ϕ = imageArray[1].rotAngleStop - rotStart # Individual image oscillation
    for diffractionImage in imageArray
        #Loop through each observation
        for hkl in keys(diffractionImage.observationList)
            centroidImageDiffErrFactor::Float32 = 0.0
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
            fracCalcErrFactor::Float32 = 0.0
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

################################################################################
#Section: Get initial amplitudes - Method 1
#-----------------------------------------------------------------------
#This method doesn't require estimates from another program. This method takes
#longer to compute than method 2 and I'm not sure that the estimates are better.
#It adds the atomic scattering factors together and then multilplies by the
#temperature factor.
function getInitialAmplitudes!(hklList::Dict{Vector{Int16}, Reflection}, atomDict::Dict{ASCIIString, UInt32}, scatteringAngles::Vector{Float32}, elementDict::Dict{ASCIIString, Element}, tempFacDict::Dict{Float32, Float32})
    #First thing: calculate the initial amplitudes corrected by the initial
    #temperature factor.
    initialAmplitudeDict = Dict{Float64, Float64}()
    for scatAngle in scatteringAngles
        f0::Float32 = 0.0
        for atom in keys(atomDict)
            f0 += elementDict[atom].f0[scatAngle] * atomDict[atom]
        end
        initialAmplitudeDict[scatAngle] = f0 * tempFacDict[scatAngle]
    end

    #Now we loop through the list of reflections and assign the estimate of the
    #initial amplitude according to the scattering angle (i.e. reflections with
    #the same reflection angle will have the same initial amplitude as an
    #initial guess).
    for hkl in keys(hklList)
        reflection = hklList[hkl]
        reflection.amplitude = initialAmplitudeDict[reflection.scatteringAngle]
    end
end
#End Section: Get initial amplitudes - Method 1
################################################################################

################################################################################
#Section: Get initial amplitudes - Method 2
#-----------------------------------------------------------------------
#This method is really crude but should be quicker than the others. It basically
#assigns the initial amplitude according to the scattering angle of the
#reflection. It is calculated as the square root of the squared amplitudes
#multiplied by the temperature factor.
function getInitialAmplitudes!(hklList::Dict{Vector{Int16}, Reflection}, f0SqrdDict::Dict{Float32, Float32}, tempFacDict::Dict{Float32, Float32})
    for hkl in keys(hklList)
        reflection = hklList[hkl]
        reflection.amplitude = sqrt(f0SqrdDict[reflection.scatteringAngle]) * tempFacDict[reflection.scatteringAngle]
    end
end
#End Section: Get initial amplitudes - Method 2
################################################################################

################################################################################
#Section: Get initial amplitudes - Method 3
#-----------------------------------------------------------------------
#This method of extracting the initial amplitudes uses the amplitudes
#derived from the refAmpDict which is basically the amplitudes extracted
#from the Ctruncate MTZ file. If the reflection is missing from the
#Ctruncate output but is present in the integrated intensity file then
#we just use the square root of the observed intensity on the image.
#Otherwise if the intensity is negative then I just set the amplitude
#value to zero.
function getInitialAmplitudes!(hklList::Dict{Vector{Int16}, Reflection}, refAmpDict::Dict{Vector{Int16}, Vector{Float32}})
    for hkl in keys(hklList)
        maxIntensity::Float32 = -1000.0
        reflection = hklList[hkl]
        if haskey(refAmpDict, reflection.symEquivHKL)
            reflection.amplitude = refAmpDict[reflection.symEquivHKL][1]
            reflection.amplitudeSig = refAmpDict[reflection.symEquivHKL][2]
        else
            for observation in reflection.observations
                if observation.intensity > maxIntensity
                    maxIntensity = observation.intensity
                end
            end
            if maxIntensity > 0.0
                reflection.amplitude = sqrt(maxIntensity)
            else
                reflection.amplitude = 0.0
            end
        end
    end
end
#End Section: Get initial amplitudes - Method 3
################################################################################
