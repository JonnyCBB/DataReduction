#########################################################################
####################### using/import statements #########################
using PyCall
using Gadfly
using LsqFit
using Colors
using Distributions
using KernelDensity
using StateSpace
using DataFrames
import Gadfly.ElementOrFunction

# include("ReciprocalSpaceUtils.jl")
# include("ElementDatabase.jl")
# include("MtzdumpHandling.jl")
# include("SequenceFileParser.jl")
# include("UpdateAtomAndRefs.jl")
# include("FilteringUtils.jl")

######### Inputs ##########
const xrayEnergy = Float32(12.7) #Set X-ray Energy
const xrayWavelength =  Float32(12.4/xrayEnergy) # Sort out this conversion at some point
const integrationFileLocation = "integration_scaling_files\\pointless.mtz"
# const integrationFileLocation = "integration_scaling_files\\test450images.mtz"
const sequenceFileLocation = "SequenceFiles\\2BN3fasta.txt"
# const sequenceFileLocation = "SequenceFiles\\4X4Vfasta.txt"
const sfFileLocation = "integration_scaling_files\\test450images_scaled1.mtz"
const useSeqFile = true #Choose whether to use a sequence file to get variance and B-factor estimates
const separateSymEquivs = false #Merge symmetry equivalents or keep them separate.
const sigIDiffTol = Float32(0.1) #Tolerance level for difference between sigIpr and sigIsum
const numOfRefs = Int32(20000) #Number of reflections to be used in data reduction analysis.

const intensityType = "Combined" #How to deal with Ipr and Isum
const numMtzColsFor1stRefLine = UInt8(9) #Number of columns in 1st MTZ Dump line for reflection information
const numMtzColsFor2ndand3rdRefLines = UInt8(4) #Number of columns in 2nd/3rd MTZ Dump line for reflection information
const numMtzColsIntLineCTruncate = UInt8(6)
const estimateTotalIntensityFromPartialRef = true #Estimate the total intensity from partial information.
const additionalElements = "S 2"

const imageOscillation = Float32(0.1) #degrees of oscillation for each image.

const minRefInResBin = UInt16(50) #choose minimum number of reflections in resolution bin.
const minRefPerImage = UInt32(3)
const displayPlots = false

const minFracCalc = Float32(0.95)
const applyBFacTof0 = true

const kdeStep = Float32(0.0001)
const keepPercentageScaleData = Float32(0.9)

const outputImageDir = "plots"

const processVarCoeff = 1.0
const estimatedObservationVar = 1e20
const measurementVarCoeff = 1.0
const estMissObs = true

#Parameters for the Unscented Kalman Filter
const α = 1e-3
const β = 2.0
const κ = 0.0

const NUM_CYCLES = 5
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
spacegroup, unitcell, hklList, imageArray = parseMosflmMTZDumpOutput(mtzdumpOutput, imageOscillation)
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
# mtzDumpOutput = runMtzdump(sfFileLocation, Int32(1200))
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
########################################################################

################################################################################
#Section: Iteration section treating reflections independently
#-------------------------------------------------------------------------------
const NUM_IMAGES = length(imageArray)
scaleFactor = modalScale
hklCounter = 0
numPlotColours = 3
getColors = distinguishable_colors(numPlotColours, Color[LCHab(70, 60, 240)],
                                   transform=c -> deuteranopic(c, 0.5),
                                   lchoices=Float64[65, 70, 75, 80],
                                   cchoices=Float64[0, 50, 60, 70],
                                   hchoices=linspace(0, 330, 24))

for hkl in keys(hklList)
    hklCounter += 1
    if hklCounter == 2
        println("made it through a round :)")
        break
    end
    reflection = hklList[hkl]
    D = SFMultiplierDict[reflection.scatteringAngle]
    Σ = f0SqrdDict[reflection.scatteringAngle]
    σ = sqrt(abs(1.0 - D^2)*Σ)

    ############################ Get initial State #############################
    #NOTE: Setting the initial variance equal to the size of the initial state
    #seems to work well in simulations. There isn't any theory to support this
    #though so may need to be changed.
    initialGuess = MvNormal([Float64(reflection.amplitude)], [Float64(reflection.amplitude)])
    # println("Initial Guess State: ", mean(initialGuess))
    # println("Initial Guess Var: ", cov(initialGuess))
    ############################################################################

    ########################## Get observation info ############################
    #Initialise the observation vector and the variances
    observationVec = fill(NaN, NUM_IMAGES)
    observationVarVec = fill(estimatedObservationVar, NUM_IMAGES)

    #Insert the measured values of the intensities and their variances on the
    #images where the reflection was actually measured.
    #To prevent time going through every image we get the total number of
    #observations for a reflection and once we've recorded that many
    #observations we then break out of the loop.
    numRefObs = length(reflection.observations)
    for imgNum in 1:NUM_IMAGES
        diffImage = imageArray[imgNum]
        obsCount = 0
        if haskey(diffImage.observationList, hkl)
            obsCount += 1
            observationVec[imgNum] = diffImage.observationList[hkl].intensity
            observationVarVec[imgNum] = diffImage.observationList[hkl].sigI^2
            if obsCount == numRefObs
                break
            end
        end
    end
    ############################################################################

    #This value for m is only used to make sure the variable is in the required
    #scope. m is rewritten before it is actually used.
    m = AdditiveNonLinUKFSSM(processFunction, [σ^2]',
                               observationFunction, [observationVarVec[1]]')
    for iterNum in 1:NUM_CYCLES
        ########################################################################
        #Section: Perform Filtering
        #-----------------------------------------------------------------------
        y = observationVec'
        params = UKFParameters(α, β, κ)
        x_filtered = Array(AbstractMvNormal, size(y, 2) + 1)
        x_filtered[1] = initialGuess
        y_obs = zeros(y)
        loglik = 0.0
        for i in 1:size(y, 2)
            y_current = y[:, i]
            processFunction(state) = processFunction(mean(x_filtered[i])[1], D, σ)
            observationFunction(state) = observationFunction(state, scaleFactor)
            # if !isnan(y_current[1])
            #     println("reflection amplitude: ", reflection.amplitude)
            #     println("process Variance: ", σ^2)
            #     println("Image Num: ", i)
            #     println("observation is: ", y_current)
            #     println("observation variance is: ", observationVarVec[i])
            #     println()
            # end
            m = AdditiveNonLinUKFSSM(processFunction, [σ^2]',
                                       observationFunction, [observationVarVec[i]]')
            x_pred, sigma_points = StateSpace.predict(m, x_filtered[i], params)
            y_pred, P_xy = observe(m, x_pred, sigma_points, y_current)

            # Check for missing values in observation
            y_Boolean = isnan(y_current)
            if any(y_Boolean)
                if estMissObs
                    y_current, y_cov_mat = estimateMissingObs(m, x_pred, y_pred, y_current, y_Boolean)
                    x_filtered[i+1] = update(m, x_pred, sigma_points, y_current, y_cov_mat)
                    loglik += logpdf(observe(m, x_filtered[i+1], calcSigmaPoints(x_filtered[i+1], params), y_current)[1], y_current)
                else
                    x_filtered[i+1] = x_pred
                end
            else
                x_filtered[i+1] = update(m, x_pred, sigma_points, y_current)
                loglik += logpdf(observe(m, x_filtered[i+1], calcSigmaPoints(x_filtered[i+1], params), y_current)[1], y_current)
            end
            loglik += logpdf(x_pred, mean(x_filtered[i+1]))
            y_obs[:,i] = y_current
        end
        filtState = FilteredState(y_obs, x_filtered, loglik)
        ########################################################################
        #Mini Section: Plot Filtering Results
        #-----------------------------------------------------------------------
        df_fs = DataFrame(
            x = 0:NUM_IMAGES,
            y = vec(mean(filtState)),
            ymin = vec(mean(filtState)) - 2*sqrt(vec(cov(filtState))),
            ymax = vec(mean(filtState)) + 2*sqrt(vec(cov(filtState))),
            f = "Filtered values"
            )

        pltflt = plot(
        layer(df_fs, x=:x, y=:y, ymin=:ymin, ymax=:ymax, Geom.line, Geom.ribbon)
        )
        display(pltflt)
        #End Mini Section: Plot Filtering Results
        ########################################################################
        #End Section: Perform Filtering
        ########################################################################


        ########################################################################
        #Section: Perform Smoothing
        #-----------------------------------------------------------------------
        n = size(filtState.observations, 2)
        smooth_dist = Array(AbstractMvNormal, n)
        #If the final filtered amplitude value was negative then set it to zero.
        if mean(filtState.state[end])[1] >= 0.0
            smooth_dist[end] = filtState.state[end]
        else
            smooth_dist[end] = MvNormal([0.0], cov(filtState.state[end]))
        end
        loglik = logpdf(observe(m, smooth_dist[end], calcSigmaPoints(smooth_dist[end], params), filtState.observations[:, end])[1], filtState.observations[:, end])
        for i in (n - 1):-1:1
            sp = calcSigmaPoints(filtState.state[i+1], params)
            processFunction(state) = processFunction(mean(x_filtered[i+1])[1], D, σ)
            observationFunction(state) = observationFunction(state, scaleFactor)
            m = AdditiveNonLinUKFSSM(processFunction, [σ^2]',
                                       observationFunction, [observationVarVec[i]]')
            pred_state, cross_covariance = smoothedTimeUpdate(m, filtState.state[i+1], sp)
            smootherGain = cross_covariance * inv(cov(pred_state))
            x_smooth = mean(filtState.state[i+1]) + smootherGain * (mean(smooth_dist[i+1]) - mean(pred_state))
            P_smooth = cov(filtState.state[i+1]) + smootherGain * (cov(smooth_dist[i+1]) - cov(pred_state)) * smootherGain'
            ################################################################
            #set any negative amplitude values to zero.
            if x_smooth[1] < 0.0
                x_smooth[1] = 0.0
            end
            ################################################################
            smooth_dist[i] = MvNormal(x_smooth, P_smooth)
            loglik += logpdf(predictSmooth(m, smooth_dist[i], params), mean(smooth_dist[i+1]))
            if !any(isnan(filtState.observations[:, i]))
                loglik += logpdf(observe(m, smooth_dist[i], calcSigmaPoints(smooth_dist[i], params), filtState.observations[:, i])[1], filtState.observations[:, i])
            end
        end
        smoothedState = FilteredState(filtState.observations, smooth_dist, loglik)


        ########################################################################
        #Mini Section: Plot Smoothing Results
        #-----------------------------------------------------------------------
        df_ss = DataFrame(
            x = 1:NUM_IMAGES,
            y = vec(mean(smoothedState)),
            ymin = vec(mean(smoothedState)) - 2*sqrt(vec(cov(smoothedState))),
            ymax = vec(mean(smoothedState)) + 2*sqrt(vec(cov(smoothedState))),
            f = "Filtered values"
            )

        pltsmth = plot(
        #layer(x=1:NUM_IMAGES, y=zeros(NUM_IMAGES), Geom.line, Theme(default_color=colorant"black")),
        #layer(xintercept=randVec, Geom.vline, Theme(default_color=getColors[2], line_width=4px)),
        layer(df_ss, x=:x, y=:y, ymin=:ymin, ymax=:ymax, Geom.line, Geom.ribbon, Theme(line_width=4px))
        )
        display(pltsmth)
        #End Mini Section: Plot Smoothing Results
        ########################################################################
        #End Section: Perform Smoothing
        ########################################################################
    end
end

#End Section: Iteration section treating reflections independently
################################################################################

















# ################################################################################
# #Section: Set up initial guess for the amplitudes
# #-------------------------------------------------------------------------------
# initialStateEstimate = getInitialState(hklList, f0SqrdDict)
# #End Section: Extract initial guess structure factor amplitudes
# ################################################################################

# ################################################################################
# #Section: Perform filtering
# #-------------------------------------------------------------------------------
# α=1e-3
# β=2.0
# κ=0.0
# ukfParams = UKFParameters(α, β, κ)
# numImages = length(imageArray)
# filtered_states = Vector{AbstractMvNormal}(numImages)
# loglikFilt = 0.0
#
# ######## Create an index reference for the miller indices of the reflections.
# hklIndexReference = createHKLIndexReferenceDict(hklList)
#
# ################################################################################
# #Mini Section: Extract filtering parameters
# #-------------------------------------------------------------------------------
# imgNum = 1
# img = imageArray[imgNum]
# observationVector = Vector{Float64}(length(hklList))
# processVarVector = Vector{Float64}(length(hklList))
# observationVarVector = Vector{Float64}(length(hklList))
# SFMultiplierVec = Vector{Float64}(length(hklList))
# SFSigmaVec = Vector{Float64}(length(hklList))
# observationIndices = Vector{Int64}()
# scaleFactor = modalScale
# hklCounter = 0
# for hkl in keys(hklList)
#     hklCounter += 1
#     reflection = hklList[hkl]
#     D = SFMultiplierDict[reflection.scatteringAngle]
#     Σ = f0SqrdDict[reflection.scatteringAngle]
#     σ = sqrt(abs(1.0 - D^2)*Σ)
#     ############################################################################
#     #NEED TO SORT THIS OUT. I'VE HAD TO USE A GAUSSIAN DISTRIBUTION INSTEAD OF A
#     #RICIAN DISTRIBUTION.
#     #processVarVector[hklCounter] = varRice(F, D, σ)
#     if reflection.isCentric
#         processVarVector[hklCounter] = 2.0 * σ^2 * reflection.epsilon * processVarCoeff
#     else
#         processVarVector[hklCounter] = σ^2 * reflection.epsilon * processVarCoeff
#     end
#     SFMultiplierVec[hklCounter] = D
#     SFSigmaVec[hklCounter] = σ
#     if haskey(img.observationList, hkl)
#         push!(observationIndices, hklCounter)
#         observationVector[hklCounter] = img.observationList[hkl].intensity
#         observationVarVector[hklCounter] = img.observationList[hkl].sigI^2 * measurementVarCoeff
#     else
#         observationVector[hklCounter] = NaN
#         observationVarVector[hklCounter] = (2.0 * scaleFactor * processFunction(reflection.amplitude, D, σ))^2 * processVarVector[hklCounter] * observationVarCoeff
#     end
# end
# processFunction(amplitudes::Vector{Float64}) = processFunction(amplitudes, SFMultiplierVec, SFSigmaVec)
# observationFunction(amplitudes::Vector{Float64}) = observationFunction(amplitudes, scaleFactor)
#
# ukfStateModel = AdditiveNonLinUKFSSM(processFunction, diagm(processVarVector), observationFunction, diagm(observationVarVector))
# #End Mini Section: Extract filtering parameters
# ################################################################################
#
# ################################################################################
# #Mini Section: Peform Filtering
# #-------------------------------------------------------------------------------
# if imgNum == 1
#     amp_pred, sigma_points = predict(ukfStateModel, initialStateEstimate, ukfParams)
# else
#     amp_pred, sigma_points = predict(ukfStateModel, filtered_states[imgNum-1], ukfParams)
# end
# y_pred, P_xy = observe(ukfStateModel, amp_pred, sigma_points, observationVector)
# # Check for missing values in observation
# obs_Boolean = isnan(observationVector)
# if any(obs_Boolean)
#     if estMissObs
#         observationVector, obs_cov_mat = estimateMissingObs(ukfStateModel, amp_pred, y_pred, observationVector, obs_Boolean)
#         filtered_states[imgNum] = update(ukfStateModel, amp_pred, sigma_points, observationVector, obs_cov_mat)
#         ########################################################################
#         #NEED TO EXTRACT THE PARTS THAT ARE OBSERVED HERE.
#         ########################################################################
#         #loglikFilt += logpdf(observe(m, filtered_states[i], calcSigmaPoints(filtered_states[i], params), observationVector)[1], observationVector)
#     else
#         filtered_states[imgNum] = amp_pred
#     end
# else
#     filtered_states[imgNum] = update(ukfStateModel, amp_pred, sigma_points, observationVector)
# 	#loglikFilt += logpdf(observe(m, filtered_states[i], calcSigmaPoints(filtered_states[i], params), observationVector)[1], observationVector)
# end
# #loglikFilt += logpdf(amp_pred, mean(filtered_states[i]))
# #End Mini Section: Peform Filtering
# ################################################################################
#
# #End Section: Perform Filtering
# ################################################################################
#
# outfile = open("testOutputComp.txt", "w")
# writecsv(outfile, mean(filtered_states[1]))
# writecsv(outfile, cov(filtered_states[1]))
# close(outfile)
