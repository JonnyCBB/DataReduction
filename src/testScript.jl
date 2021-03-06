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

include("ReciprocalSpaceUtils.jl")
include("ElementDatabase.jl")
include("MtzdumpHandling.jl")
include("SequenceFileParser.jl")
include("UpdateAtomAndRefs.jl")
include("FilteringUtils.jl")
include("WriteOutputFiles.jl")

######### Inputs ##########
const xrayEnergy = Float32(12.7) #Set X-ray Energy
const xrayWavelength =  Float32(12.4/xrayEnergy) # Sort out this conversion at some point
const integrationFileLocation = "integration_scaling_files\\pointless.mtz"
# const integrationFileLocation = "integration_scaling_files\\test450images.mtz"
const sequenceFileLocation = "SequenceFiles\\2BN3fasta.txt"
# const sequenceFileLocation = "SequenceFiles\\4X4Vfasta.txt"
const sfFileLocation = "integration_scaling_files\\test450images_scaled1_withSeq.mtz"
const useSeqFile = true #Choose whether to use a sequence file to get variance and B-factor estimates
const separateSymEquivs = false #Merge symmetry equivalents or keep them separate.
const sigIDiffTol = Float32(0.1) #Tolerance level for difference between sigIpr and sigIsum
const numOfRefs = Int32(-1) #Number of reflections to be used in data reduction analysis.

const intensityType = "Combined" #How to deal with Ipr and Isum
const numMtzColsFor1stRefLine = UInt8(9) #Number of columns in 1st MTZ Dump line for reflection information
const numMtzColsFor2ndand3rdRefLines = UInt8(4) #Number of columns in 2nd/3rd MTZ Dump line for reflection information
const numMtzColsIntLineCTruncate = UInt8(6)
const estimateTotalIntensityFromPartialRef = true #Estimate the total intensity from partial information.
const additionalElements = ""

const imageOscillation = Float32(0.1) #degrees of oscillation for each image.

const minRefInResBin = UInt16(100) #choose minimum number of reflections in resolution bin.
const minRefPerImage = UInt32(3)
const displayPlots = false

const minFracCalc = Float32(0.95)
const USE_LOW_FRAC_REF = false
const USE_LOW_FRAC_CALC_IFF_SINGLE_OBS = false
const applyBFacTof0 = true

const kdeStep = Float32(0.0001)
const keepPercentageScaleData = Float32(0.95)

const outputImageDir = "plots3"

const processVarCoeff = 1e0
const estimatedObservationVar = 1e6
const measurementVarCoeff = 1e0
const estMissObs = true

#Parameters for the Unscented Kalman Filter
const α = 1e-3
const β = 2.0
const κ = 0.0

const NUM_CYCLES = 200
const MIN_CYCLE_NUM = 5
const MIN_SCALING_CYCLE_NUM = 3
const MAX_SCALING_CYCLES = 100
const INITIAL_VALUE_TOL = 1e0 # Need a more informative name for this.
const USE_WEAK_REF_PRIOR = false
const USE_BAYESIAN_EST_VALUES_AS_FINAL = false
const WEAK_AMP_THRESHOLD = 2.0
const NUM_STD_FOR_INTEGRATION = 3.0
const LOG_LIK_THRESHOLD = 1e0
const USE_RICE_ESTIMATE = false

const ANOMALOUS = false

const HKL_FILENAME = "ReflectionFile_Merged.hkl"
const F2MTZ_INPUT_FILENAME = "f2mtzInput.dat"
const PROJECT_NAME = "Filtering_Test"
const CRYSTAL_NAME = "Insulin_0259"
const DATASET_NAME = "Insulin_0259"
const MTZOUT_FILENAME = "Reduction3.mtz"
const SORTMTZ_INPUT_FILENAME ="SortmtzInput.txt"
################################################################################
#Assertions
@assert NUM_CYCLES > MAX_SCALING_CYCLES
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
updateRefListAndImageArray!(hklList, imageArray, estimateTotalIntensityFromPartialRef, USE_LOW_FRAC_REF, USE_LOW_FRAC_CALC_IFF_SINGLE_OBS)
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
# getInitialAmplitudes!(hklList, f0SqrdDict, tempFacDict)

# Get initial amplitudes by method 3
mtzDumpOutput = runMtzdump(sfFileLocation, Int32(-1))
refAmpDict, scaleFac = parseCTruncateMTZDumpOutput(mtzDumpOutput)
getInitialAmplitudes!(hklList, refAmpDict)

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
const NUM_REFLECTIONS = length(hklList)
refNumsToPlot = collect(1:300)
#scaleFactor = modalScale
scaleFactor = 1/40.5778 #Value I obtained from CTruncate
hklCounter = 0
numPlotColours = 3
getColors = distinguishable_colors(numPlotColours, Color[LCHab(70, 60, 240)],
                                   transform=c -> deuteranopic(c, 0.5),
                                   lchoices=Float64[65, 70, 75, 80],
                                   cchoices=Float64[0, 50, 60, 70],
                                   hchoices=linspace(0, 330, 24))
################################################################################
#Preallocate arrays to create a DataFrame for further statistical analysis.
allHKLs = collect(keys(hklList))
CtruncAmp = Vector{AbstractFloat}(NUM_REFLECTIONS)
CtruncSig = Vector{AbstractFloat}(NUM_REFLECTIONS)
MyAmp = Vector{AbstractFloat}(NUM_REFLECTIONS)
MySig = Vector{AbstractFloat}(NUM_REFLECTIONS)
relAmpDiff = Vector{AbstractFloat}(NUM_REFLECTIONS)
relSigDiff = Vector{AbstractFloat}(NUM_REFLECTIONS)
inBounds = Vector{ASCIIString}(NUM_REFLECTIONS)
resCol = Vector{AbstractFloat}(NUM_REFLECTIONS)
resbinsCol = Vector{ASCIIString}(NUM_REFLECTIONS)
weakRefCol = Vector{ASCIIString}(NUM_REFLECTIONS)
################################################################################
@printf("Beginning Filtering/Smoothing Cycles...\n")
@printf("=======================================\n")
for hkl in keys(hklList)
    hklCounter += 1
    @printf("Iteration %d of %d: Reflection [%d,%d,%d]\n",hklCounter, NUM_REFLECTIONS, hkl[1], hkl[2], hkl[3])
    if hklCounter == NUM_REFLECTIONS + 1
        println("Skipping the filtering loop.")
        break
    end
    reflection = hklList[hkl]
    D = SFMultiplierDict[reflection.scatteringAngle]
    Σ = f0SqrdDict[reflection.scatteringAngle]
    σ = sqrt(abs(1.0 - D^2)*Σ)
    scalingCycles = true

    ############################ Get initial State #############################
    #NOTE: Setting the initial variance equal to the size of the initial state
    #seems to work well in simulations. There isn't any theory to support this
    #though so may need to be changed.
    initialGuess = MvNormal([Float64(reflection.amplitude)], [Float64(reflection.amplitude)])
    intGuess = MvNormal([Float64(reflection.amplitude)], [Float64(reflection.amplitude)])
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
    imagesWithActualObs = Vector{Int64}()
    anyPositiveInts = false
    anyNegativeInts = false
    for imgNum in 1:NUM_IMAGES
        diffImage = imageArray[imgNum]
        obsCount = 0
        if haskey(diffImage.observationList, hkl)
            obsCount += 1
            push!(imagesWithActualObs, imgNum)
            if diffImage.observationList[hkl].intensity >= 0.0
                anyPositiveInts = true
            end
            if diffImage.observationList[hkl].intensity < 0.0
                anyNegativeInts = true
            end
            observationVec[imgNum] = diffImage.observationList[hkl].intensity
            observationVarVec[imgNum] = measurementVarCoeff * diffImage.observationList[hkl].sigI^2
            if obsCount == numRefObs
                break
            end
        end
    end

    #Remove any negative intensities if there are multiple observations and
    #some are postive.
    if length(imagesWithActualObs) > 1 && anyNegativeInts && anyPositiveInts
        for i in length(imagesWithActualObs):-1:1
            imgNum = imagesWithActualObs[i]
            if imageArray[imgNum].observationList[hkl].intensity < 0.0
                observationVec[imgNum] = NaN
                observationVarVec[imgNum] = estimatedObservationVar
                splice!(imagesWithActualObs, i)
            end
        end
    end
    ############################################################################

    #This value for m is only used to make sure the variable is in the required
    #scope. m is rewritten before it is actually used.
    oldInitState = MvNormal([0.0],[1.0])
    smoothedState = FilteredState(observationVec', Array(AbstractMvNormal,0), 0.0)
    m = AdditiveNonLinUKFSSM(processFunction, [σ^2]',
                               observationFunction, [observationVarVec[1]]')
    loglikVals = Vector{Float64}(NUM_CYCLES)
    totalIterNum = 0
    endScalingIter = 0
    for iterNum in 1:NUM_CYCLES
        if iterNum != 1
            newStateVec = 1/D * smoothedState.state[1].μ
            #newStateCov = D * cov(smoothedState.state[1]) * D + m.V + (oldInitState.μ - newStateVec) * (oldInitState.μ - newStateVec)'
            if iterNum != 2
                newStateCov = D * cov(smoothedState.state[1]) * D + m.V + (oldInitState.μ - newStateVec) * (oldInitState.μ - newStateVec)'
            else
                newStateCov = D * cov(smoothedState.state[1]) * D + m.V
            end
            #Calculate F/sigF to determine strength of amplitude
            amplitudeStrength = newStateVec[1]/sqrt(newStateCov[1])

            #If the amplitude is weak and the user is happy to use Bayesian
            #methods to get a "stronger/better" amplitude estimate then let's
            #do it: lets get a "better" amplitude estimate.
            if USE_WEAK_REF_PRIOR && amplitudeStrength < WEAK_AMP_THRESHOLD
                if reflection.isCentric
                    F = newStateVec[1]
                    varF = newStateCov[1]
                    expFun(x) = 2x/varF * sqrt(2/(π*Σ)) * exp(-(x^2 + F^2)/varF - x^2/(2*Σ)) * besseli(0,2*x*F/varF)
                    expFunNum(x) = x * expFun(x)
                    expFunDenom = quadgk(expFun, 0.0, F + NUM_STD_FOR_INTEGRATION * sqrt(varF))[1]
                    expFunPred(x) = expFunNum(x)/expFunDenom
                    newStateVec = [quadgk(expFunPred, 0.0, F + NUM_STD_FOR_INTEGRATION * sqrt(varF))[1]]
                else
                    F = newStateVec[1]
                    varF = newStateCov[1]
                    expFun(x) = 4*x*F/(Σ*varF) * exp(-(x^2 + F^2)/varF - x^2/Σ) * besseli(0,2*x*F/varF)
                    expFunNum(x) = x * expFun(x)
                    expFunDenom = quadgk(expFun, 0.0 , F + NUM_STD_FOR_INTEGRATION * sqrt(varF))[1]
                    expFunPred(x) = expFunNum(x)/expFunDenom
                    newStateVec = [quadgk(expFunPred, 0.0 , F + NUM_STD_FOR_INTEGRATION * sqrt(varF))[1]]
                end
            end
            initialGuess = MvNormal(newStateVec, newStateCov)
            #initialGuess = smoothedState.state[1]
        else
            amplitudeStrength = mean(initialGuess)[1]/sqrt(cov(initialGuess)[1])
        end
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

            if scalingCycles
                processVariance = abs(1.0 - D^2) * mean(initialGuess)[1]
            else
                processVariance = abs(1.0 - D^2) * mean(initialGuess)[1]^2
            end

            if reflection.isCentric
                σ = sqrt(2 * reflection.epsilon * processVariance)
            else
                σ = sqrt(reflection.epsilon * processVariance)
            end
            if amplitudeStrength < WEAK_AMP_THRESHOLD && !scalingCycles && USE_RICE_ESTIMATE
                processFunction(state) = processFunction(state, D, σ)
                processVariance = varRice(mean(x_filtered[i])[1], D, σ)
            else
                processFunction(state) = processFunction(state, D)
                processVariance = σ^2
            end

            observationFunction(state) = observationFunction(state, scaleFactor)
            m = AdditiveNonLinUKFSSM(processFunction, [processVarCoeff * processVariance]', observationFunction, [observationVarVec[i]]')
            x_pred, sigma_points = StateSpace.predict(m, x_filtered[i], params)
            y_pred, P_xy = observe(m, x_pred, sigma_points, y_current)
            # if !any(isnan(y_current))
            #     println("reflection amplitude: ", reflection.amplitude)
            #     println("process Variance: ", σ^2)
            #     println("Image Num: ", i)
            #     println("observation is: ", y_current)
            #     println("observation variance is: ", observationVarVec[i])
            #     println("predicted Intensity is: ", y_pred)
            # end

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
        oldInitState = filtState.state[1]
        ########################################################################
        #Mini Section: Plot Filtering Results
        #-----------------------------------------------------------------------
        # df_fs = DataFrame(
        #     x = 0:NUM_IMAGES,
        #     y = vec(mean(filtState)),
        #     ymin = vec(mean(filtState)) - 2*sqrt(vec(cov(filtState))),
        #     ymax = vec(mean(filtState)) + 2*sqrt(vec(cov(filtState))),
        #     f = "Filtered values"
        #     )
        #
        # pltflt = plot(
        # layer(xintercept=imagesWithActualObs, Geom.vline, Theme(default_color=getColors[2], line_width=4px)),
        # layer(df_fs, x=:x, y=:y, ymin=:ymin, ymax=:ymax, Geom.line, Geom.ribbon, Theme(line_width=4px))
        # )
        #display(pltflt)
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
            if amplitudeStrength < WEAK_AMP_THRESHOLD && !scalingCycles && USE_RICE_ESTIMATE
                processFunction(state) = processFunction(state, D, σ)
                processVariance = varRice(mean(filtState.state[i+1])[1], D, σ)
            else
                processFunction(state) = processFunction(state, D)
                processVariance = σ^2
            end
            observationFunction(state) = observationFunction(state, scaleFactor)
            m = AdditiveNonLinUKFSSM(processFunction, [processVarCoeff * processVariance]',
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
        loglikVals[iterNum] = smoothedState.loglik
        # println()
        # println("Interation: ", iterNum)
        # println("Mean estimate is: ", vec(mean(smoothedState))[1])
        # println("Sigma value is: ", sqrt(vec(cov(smoothedState))[1]))


        ########################################################################
        #Mini Section: Plot Smoothing Results
        #-----------------------------------------------------------------------
        # df_ss = DataFrame(
        #     x = 1:NUM_IMAGES,
        #     y = vec(mean(smoothedState)),
        #     ymin = vec(mean(smoothedState)) - 2*sqrt(vec(cov(smoothedState))),
        #     ymax = vec(mean(smoothedState)) + 2*sqrt(vec(cov(smoothedState))),
        #     f = "Smoothed values"
        #     )
        #
        # smthPltTitle = @sprintf("Reflection [%d,%d,%d]",hkl[1], hkl[2], hkl[3])
        # pltsmth = plot(
        # #layer(x=1:NUM_IMAGES, y=zeros(NUM_IMAGES), Geom.line, Theme(default_color=colorant"black")),
        # layer(xintercept=imagesWithActualObs, Geom.vline, Theme(default_color=getColors[2], line_width=2px)),
        # layer(df_ss, x=:x, y=:y, ymin=:ymin, ymax=:ymax, Geom.line, Geom.ribbon, Theme(line_width=4px)),
        # Guide.xlabel("Image Number"), Guide.ylabel("Amplitude"),
        # Guide.manual_color_key("Colour Key",["Smoothed estimate", "Observation Occurence"],[getColors[1],getColors[2]]),
        # Guide.title(smthPltTitle)
        # )
        #End Mini Section: Plot Smoothing Results
        ########################################################################
        #End Section: Perform Smoothing
        ########################################################################
        if iterNum > MIN_SCALING_CYCLE_NUM
            if scalingCycles
                if abs(loglikVals[iterNum] - loglikVals[iterNum - MIN_SCALING_CYCLE_NUM]) < LOG_LIK_THRESHOLD || iterNum > MAX_SCALING_CYCLES || abs(1/D * smoothedState.state[1].μ[1] - oldInitState.μ[1]) < INITIAL_VALUE_TOL
                    endScalingIter = iterNum
                    scalingCycles = false
                    if iterNum > MAX_SCALING_CYCLES && abs(loglikVals[iterNum] - loglikVals[iterNum - MIN_SCALING_CYCLE_NUM]) >= LOG_LIK_THRESHOLD
                        println("Scaling cycles didn't converge.") #Need more helpful error statement here.
                    end
                end
            elseif iterNum > MIN_CYCLE_NUM + endScalingIter
                if abs(loglikVals[iterNum] - loglikVals[iterNum - MIN_CYCLE_NUM]) < LOG_LIK_THRESHOLD || iterNum == NUM_CYCLES
                    totalIterNum = iterNum
                    ############################################################
                    #Sort data for DataFrame - Ctruncate values
                    CtruncAmp[hklCounter] = reflection.amplitude
                    CtruncSig[hklCounter] = reflection.amplitudeSig
                    ############################################################
                    reflection.amplitude = 1/D * smoothedState.state[1].μ[1]
                    amplitudeVariance = D * cov(smoothedState.state[1]) * D + m.V + (oldInitState.μ - newStateVec) * (oldInitState.μ - newStateVec)'
                    reflection.amplitudeSig = sqrt(amplitudeVariance[1])
                    if USE_BAYESIAN_EST_VALUES_AS_FINAL
                        if reflection.isCentric
                            F = reflection.amplitude
                            varF = reflection.amplitudeSig^2
                            expFun(x) = 2x/varF * sqrt(2/(π*Σ)) * exp(-(x^2 + F^2)/varF - x^2/(2*Σ)) * besseli(0,2*x*F/varF)
                            expFunNum(x) = x * expFun(x)
                            expFunDenom = quadgk(expFun, 0.0, F + NUM_STD_FOR_INTEGRATION * sqrt(varF))[1]
                            expFunPred(x) = expFunNum(x)/expFunDenom
                            reflection.amplitude = quadgk(expFunPred, 0.0, F + NUM_STD_FOR_INTEGRATION * sqrt(varF))[1]
                        else
                            F = reflection.amplitude
                            varF = reflection.amplitudeSig^2
                            expFun(x) = 4*x*F/(Σ*varF) * exp(-(x^2 + F^2)/varF - x^2/Σ) * besseli(0,2*x*F/varF)
                            expFunNum(x) = x * expFun(x)
                            expFunDenom = quadgk(expFun, 0.0 , F + NUM_STD_FOR_INTEGRATION * sqrt(varF))[1]
                            expFunPred(x) = expFunNum(x)/expFunDenom
                            reflection.amplitude = quadgk(expFunPred, 0.0 , F + NUM_STD_FOR_INTEGRATION * sqrt(varF))[1]
                        end
                    end
                    ############################################################
                    #Sort data for DataFrame -All other values
                    MyAmp[hklCounter] = reflection.amplitude
                    MySig[hklCounter] = reflection.amplitudeSig
                    relAmpDiff[hklCounter] = (MyAmp[hklCounter] - CtruncAmp[hklCounter])/MyAmp[hklCounter]
                    relSigDiff[hklCounter] = (MySig[hklCounter] - CtruncSig[hklCounter])/MySig[hklCounter]
                    if MyAmp[hklCounter] - 2*MySig[hklCounter] <= CtruncAmp[hklCounter] <=
                        MyAmp[hklCounter] + 2*MySig[hklCounter]
                        inBounds[hklCounter] = "true"
                    else
                        inBounds[hklCounter] = "false"
                    end
                    resCol[hklCounter] = reflection.resolution
                    findHKL(x) = x == hkl
                    for bin in resbins
                        if !isempty(find(findHKL, bin.hkls))
                            resbinsCol[hklCounter] = @sprintf("%.2f > x >= %.2f", bin.minRes, bin.maxRes)
                            break
                        end
                    end
                    if MyAmp[hklCounter]/MySig[hklCounter] < WEAK_AMP_THRESHOLD
                        weakRefCol[hklCounter] = "true"
                    else
                        weakRefCol[hklCounter] = "false"
                    end
                    df_ss = DataFrame(
                        x = 1:NUM_IMAGES,
                        y = vec(mean(smoothedState)),
                        ymin = vec(mean(smoothedState)) - 2*sqrt(vec(cov(smoothedState))),
                        ymax = vec(mean(smoothedState)) + 2*sqrt(vec(cov(smoothedState))),
                        f = "Smoothed values"
                        )

                    smthPltTitle = @sprintf("Reflection [%d,%d,%d]. Resolution: %.2fA",hkl[1], hkl[2], hkl[3], reflection.resolution)
                    pltsmth = plot(
                    layer(x=1:NUM_IMAGES, y=CtruncAmp[hklCounter]*ones(NUM_IMAGES), Geom.line, Theme(default_color=getColors[3])),
                    layer(xintercept=imagesWithActualObs, Geom.vline, Theme(default_color=getColors[2], line_width=2px)),
                    layer(df_ss, x=:x, y=:y, ymin=:ymin, ymax=:ymax, Geom.line, Geom.ribbon, Theme(line_width=4px)),
                    Guide.xlabel("Image Number"), Guide.ylabel("Amplitude"),
                    Guide.manual_color_key("Colour Key",["Smoothed estimate", "Observation Occurence"],[getColors[1],getColors[2]]),
                    Guide.title(smthPltTitle)
                    )
                    #display(pltsmth)
                    if any(hklCounter .== refNumsToPlot)
                        smthPltFilename = @sprintf("testplots3/SmoothedPlot_%d,%d,%d_res%.0f.pdf",hkl[1], hkl[2], hkl[3],reflection.resolution)
                        draw(PDF(smthPltFilename, 16cm, 9cm), pltsmth)
                    end
                    ############################################################
                    break
                end
            end
        end
    end
    loglikpltTitle = @sprintf("Log Likelihood VS Cycle Number. Ref [%d,%d,%d]. Resolution: %.2fA",hkl[1], hkl[2], hkl[3], reflection.resolution)
    pltloglik = plot(
    layer(x=1:totalIterNum, y=loglikVals[1:totalIterNum], Geom.line, Theme(line_width=4px)),
    layer(xintercept=[endScalingIter], Geom.vline, Theme(default_color=getColors[2], line_width=4px)),
    Guide.xlabel("Cycle Number"), Guide.ylabel("Log Likelihood"),
    Guide.manual_color_key("Colour Key",["Log Likelihood", "End of scaling cycles"],[getColors[1],getColors[2]]),
    Guide.title(loglikpltTitle)
    )
    #display(pltloglik)
    if any(hklCounter .== refNumsToPlot)
        loglikPltFilename = @sprintf("testplots3/LogLikplot_%d,%d,%d_res%.0f.pdf",hkl[1], hkl[2], hkl[3],reflection.resolution)
        draw(PDF(loglikPltFilename, 16cm, 9cm), pltloglik)
    end
end
@printf("End of filtering/smoothing cycles.\n\n")
#End Section: Iteration section treating reflections independently
################################################################################

########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################

################################################################################
#Section: Write Output files
#-------------------------------------------------------------------------------
@printf("Writing output files...\n")
writeOutputFiles(hklList, spacegroup, unitcell, HKL_FILENAME,
                 F2MTZ_INPUT_FILENAME, PROJECT_NAME, CRYSTAL_NAME, DATASET_NAME,
                 SORTMTZ_INPUT_FILENAME, MTZOUT_FILENAME)
@printf("Finished Writing output files.\n\n")
#End Section: Write Output files
################################################################################
df = DataFrame(HKL = allHKLs, CAmp = CtruncAmp, CSig = CtruncSig, MAmp = MyAmp,
               MSig = MySig, RelDiffA = relAmpDiff, RelDiffS = relSigDiff,
               WeakRef = weakRefCol, InBound = inBounds, Resolution = resCol, ResBins = resbinsCol
               )
if hklCounter >= NUM_REFLECTIONS
    writetable("OutputComparison3.csv", df)
end
@printf("Program run = SUCCESS :)\n")
