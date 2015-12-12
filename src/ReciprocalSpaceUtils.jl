
immutable Unitcell{T<:Float32}
    a::T
    b::T
    c::T
    α::T
    β::T
    γ::T
end

type ReflectionObservation
    rotCentroid::Float32
    fractionCalc::Float32
    misym::UInt16
    intensity::Float32
    sigI::Float32
    imageNums::Vector{UInt16}
    imageIntensities::Vector{Float32}
end

type DiffractionImage
    rotAngleStart::Float32
    rotAngleStop::Float32
    observationList::Dict{Vector{Int16}, ReflectionObservation}
end
DiffractionImage(startAng::Float32, stopAng::Float32) = DiffractionImage(startAng,
stopAng, Dict{Vector{Int16}, ReflectionObservation}())

immutable SpaceGroup
    symbol::ASCIIString
    number::UInt8
    numSymOps::UInt8
    symOps::Vector{Matrix{Float32}}
end

type Reflection
    hkl::Vector{Int16}
    symEquivHKL::Vector{Int16}
    resolution::Float32
    scatteringAngle::Float32
    epsilon::UInt8
    isCentric::Bool
    isSysAbs::Bool
    amplitude::Float32
    amplitudeSig::Float32
    observations::Vector{ReflectionObservation}

    function Reflection(hkl::Vector{Int16}, symEquivHKL::Vector{Int16},
                        spacegroup::SpaceGroup, unitcell::Unitcell, xrayWavelength::Float32)
        reflection = new()
        reflection.hkl = hkl
        reflection.symEquivHKL = symEquivHKL
        recipspaceunitcell = convertDirectToReciprocal(unitcell)
        reflection.resolution = calcResolution(recipspaceunitcell,hkl)
        reflection.scatteringAngle = calcScatteringAngle(reflection.resolution,
        xrayWavelength)
        determineHKLClass!(reflection, spacegroup, unitcell)
        reflection.amplitude = 0.0
        reflection.amplitudeSig = 0.0
        reflection.observations = Vector{ReflectionObservation}()
        return reflection
    end
end

type ResBin
    minRes::Float32
    maxRes::Float32
    numOfRef::UInt16
    meanIntensity::Float32
    hkls::Array{Vector{Int16},1}
end
ResBin(minRes::Float32, maxRes::Float32, numOfRef::UInt16) = ResBin(minRes, maxRes, numOfRef, Float32(0.0), Array(Vector{Int16},1))


function convertDirectToReciprocal(directUnitCell::Unitcell)
    a = directUnitCell.a
    b = directUnitCell.b
    c = directUnitCell.c
    α = deg2rad(directUnitCell.α)
    β = deg2rad(directUnitCell.β)
    γ = deg2rad(directUnitCell.γ)

    #Calculate volume of direct Unit cell
    V = a*b*c*sqrt(1 - (cos(α))^2 - (cos(β))^2 - (cos(γ))^2 + (2 * cos(α) * cos(β) * cos(γ)))

    #Calculate reciprocal lattice parameters
    aᵣ = (    b * c * sin(α))/V
    bᵣ = (a *     c * sin(β))/V
    cᵣ = (a * b *     sin(γ))/V
    αᵣ = rad2deg(acos((cos(β)  * cos(γ) - cos(α)) / (sin(β) * sin(γ))));
    βᵣ = rad2deg(acos((cos(α)  * cos(γ) - cos(β)) / (sin(α) * sin(γ))));
    γᵣ = rad2deg(acos((cos(α)  * cos(β) - cos(γ)) / (sin(α) * sin(β))));

    #Create reciprocal unit cell lattice
    reciprocalUnitcell = Unitcell(aᵣ, bᵣ, cᵣ, αᵣ, βᵣ, γᵣ)
    return reciprocalUnitcell
end

function calcResolution(reciprocalUnitcell::Unitcell, millerIndices::Vector{Int16})
    h = millerIndices[1]
    k = millerIndices[2]
    l = millerIndices[3]

    a = reciprocalUnitcell.a
    b = reciprocalUnitcell.b
    c = reciprocalUnitcell.c
    α = deg2rad(reciprocalUnitcell.α)
    β = deg2rad(reciprocalUnitcell.β)
    γ = deg2rad(reciprocalUnitcell.γ)

    resolution = 1/sqrt(h^2 * a^2 + k^2 * b^2 + l^2 * c^2 +
                            2 * h * k * a * b * cos(γ) +
                            2 * h * l * a * c * cos(β) +
                            2 * k * l * b * c * cos(α))

    return resolution
end

function calcScatteringAngle(resolution::Float32, xrayWavelength::Float32)
    if isinf(resolution)
        return 0.0
    else
        d = resolution
        λ = xrayWavelength
        θ = rad2deg(asin(λ/2d)) #calculate scattering angle from Bragg's Law
        return θ
    end
end

function getAllScatteringAngles(reflectionList::Dict{Vector{Int16},Reflection})
    scatteringAngles = Array(Float32,length(reflectionList))
    counter = 1
    for key in keys(reflectionList)
        scatteringAngles[counter] = reflectionList[key].scatteringAngle
        counter += 1
    end
    scatteringAngles = unique(scatteringAngles)
    return scatteringAngles
end

function determineHKLClass!(reflection::Reflection, spacegroup::SpaceGroup, unitcell::Unitcell)
    tolerance = 1e-3
    reflection.epsilon = 0
    reflection.isSysAbs = false
    reflection.isCentric = false
    for transformationMatrix in spacegroup.symOps
        symEquivReflection, phaseShift = applySymOpTrans(transformationMatrix, reflection.hkl, unitcell)
        if(symEquivReflection == reflection.hkl)
            if (cos(phaseShift) > (1 - tolerance))
                reflection.epsilon += 1
            else
                reflection.epsilon, reflection.isSysAbs = 0, true
                break
            end
        elseif(symEquivReflection == -reflection.hkl)
            reflection.isCentric = true
        end
    end
end

function applySymOpTrans(transformation::Array{Float32,2}, point::Vector{Int16}, unitcell::Unitcell)
    transformationMatrix = transformation[1:3,1:3]
    translation = transformation[1:3,4]
    if sort(translation)[end] > 1.0
        relativeTranslation = [translation[1]/unitcell.a, translation[2]/unitcell.b, translation[3]/unitcell.c]
    else
        relativeTranslation = translation
    end
    pointTransformation = map(x -> Int16(x), transformationMatrix * point)
    phaseShift = 2π*dot(relativeTranslation, point)
    return pointTransformation, phaseShift
end

function removeSysAbs!(reflectionList::Dict{Vector{Int16},Reflection})
    for key in keys(reflectionList)
        if reflectionList[key].isSysAbs
            delete!(reflectionList, key)
        end
    end
end

function sortResolutionBins(reflectionList::Dict{Vector{Int16},Reflection}, minRefPerBin::UInt16)
    resbins = ResBin[]
    counter = 1
    numOfRef = length(reflectionList)
    resolutionArray = Array(Float32,numOfRef)

    #First we need to sort the resolution values into increasing resolution order
    #(i.e. from biggest resolution number to smallest)
    for key in keys(reflectionList)
        resolutionArray[counter] = reflectionList[key].resolution
        counter += 1
    end
    resolutionArray = sort(resolutionArray, rev=true)

    #Now we can sort the resolution values
    previousRes = resolutionArray[1]
    numOfRefInResBinCounter = 0
    minRes = previousRes
    for i = 1:numOfRef
        if numOfRefInResBinCounter >= minRefPerBin && resolutionArray[i] != previousRes
            maxRes = previousRes
            resbin = ResBin(minRes, maxRes, UInt16(numOfRefInResBinCounter))
            push!(resbins, resbin)
            minRes = previousRes
            numOfRefInResBinCounter = 0
        end
        numOfRefInResBinCounter += 1
        previousRes = resolutionArray[i]
    end
    return resbins
end

function sortHKLIntoResBins!(resbins::Array{ResBin,1}, reflectionList::Dict{Vector{Int16},Reflection})
    minRes = resbins[1].minRes
    for resbin in resbins
        resbin.hkls = Array(Vector{Int16}, resbin.numOfRef)
        refCounter = 1
        for key in keys(reflectionList)
            reflection = reflectionList[key]
            if resbin.minRes > reflection.resolution ≥ resbin.maxRes || reflection.resolution == minRes
                resbin.hkls[refCounter] = reflection.hkl
                if refCounter == resbin.numOfRef
                    break
                else
                    refCounter += 1
                end
            elseif reflection.hkl == [Int16(0),Int16(0),Int16(0)] && resbin.minRes == Inf32
                resbin.hkls[refCounter] = reflection.hkl
                if refCounter == resbin.numOfRef
                    break
                else
                    refCounter += 1
                end
            end
        end
    end
end

function calcResbinMeanIntensity!(resbins::Vector{ResBin}, f0SqrdDict::Dict{Float32, Float32}, hklList::Dict{Vector{Int16}, Reflection})
    for resbin in resbins
        summedIntensity = 0.0
        for hkl in resbin.hkls
            summedIntensity += f0SqrdDict[hklList[hkl].scatteringAngle]
        end
        resbin.meanIntensity = summedIntensity/resbin.numOfRef
    end
end

################################################################################
#NEED TO COME BACK TO THIS LATER
################################################################################
function assignMeanIntensityToResBins!(resbins::Array{ResBin,1}, reflectionList::Dict{Vector{Int16},Reflection}, f0SqrdDict::Dict{Float32, Float32})
    binCounter = 1
    for resbin in resbins
        intensityArray = Array(Float64,resbin.numOfRef,1)
        refCounter = 0
        ########################################################################
        #Section: Calculate the average intensity in resolution shells
        #-----------------------------------------------------------------------
        #If there is was no sequence file/opted not to use a sequence file to
        #obtain the scattering factor information, then the average intensity
        #for a resolution shell will be calculated and used as the mean expected
        #intensity for the resolution shell. Note that this intensity is
        #corrected for the expected intensity factor (denoted epsilon).
        #Since there are many observations for a given reflection we use the
        #intensity from the observation that gave the biggest measured
        #intensity. The maximum is chosen because this will give the biggest
        #variance to reflect our biggest level of uncertainty since we have not
        #been given the atomic composition.
        if isempty(atomDict)
            for hkl in resbin.hkls
                refCounter += 1
                reflection = reflectionList[hkl]
                obsNumWithMaxIntensity = 1
                numRefObs = length(reflection.observations)
                if numRefObs > 1
                    for obsNum in 1:numRefObs
                        observation = reflection.observations[obsNum]
                        if observation.intensity > reflection.observations[obsNumWithMaxIntensity].intensity/reflection.epsilon
                            obsNumWithMaxIntensity = obsNum
                        end
                    end
                end
                intensityArray[refCounter] = reflection.observations[obsNumWithMaxIntensity].intensity
            end
        #-----------------------------------------------------------------------
        #If the sequence file has been given then the atomic scattering factors
        #will have been available. So we can use the sum of the squared atomic
        #scattering factors as an expected intensity estimate. These are stored
        #in the f0SrqdDict. Since it is resolution dependent the keys are the
        #corresponding scattering angles.
        else
            for hkl in resbin.hkls
                refCounter += 1
                reflection = reflectionList[hkl]
                intensityArray[refCounter] = f0SqrdDict[reflection.scatteringAngle]
            end
        end

        resbin.meanIntensity = mean(intensityArray)
        resbin.stdDev = std(intensityArray)
        resbins[binCounter] = resbin
        binCounter += 1
        #End Section: Calculate the average intensity in resolution shells
        ########################################################################
    end
end

function calcBandScaleParams(hklList::Dict{Vector{Int16},Reflection}, imageArray::Vector{DiffractionImage}, resbins::Vector{ResBin}, outputImageDir::ASCIIString="", displayPlots::Bool=false)
    bfactors = Vector{Float32}()
    bfactorSigmas = Vector{Float32}()
    scaleFactors = Vector{Float32}()
    scaleFactorSigmas = Vector{Float32}()
    imageNumArray = Vector{Int16}()
    imageNum::UInt16 = 0
    #loop through image array
    for img in imageArray
        imageNum += 1 # increment image number
        intensityDict = Dict{Float32,Float32}() #create intensity dict Dict{resolution, intensity/epsilon}
        #loop through all observations on the image
        for hkl in keys(img.observationList)
            reflection = hklList[hkl] # get reflection
            #extract epsilon corrected intensity
            intensityDict[reflection.resolution] = img.observationList[hkl].intensity/reflection.epsilon
        end
        logIntensityRatio = Vector{Float64}()
        hSqrd = Vector{Float64}()
        #now loop through each resolution bin
        for resbin in resbins
            summedIntensity::Float64 = 0.0
            numRefs::UInt16 = 0
            #loop through all resolutions in the list of observations
            for resolution in keys(intensityDict)
                if resbin.minRes > resolution >= resbin.maxRes
                    numRefs += 1
                    summedIntensity += intensityDict[resolution] # sum the intensities
                end
            end
            meanIntensity = summedIntensity/numRefs # calculate mean intensity for the resolution bin
            if !isnan(meanIntensity) && meanIntensity > 0.0
                push!(logIntensityRatio, log(meanIntensity/resbin.meanIntensity))
                midres = (resbin.maxRes + resbin.minRes)/2.0
                push!(hSqrd, 1.0/midres^2)
            end
        end
        #Check that there were enough data points for the image
        if length(logIntensityRatio) >= minRefPerImage
            model(x, p) = p[1] + p[2] * x #create linear function model
            fit = curve_fit(model, hSqrd, logIntensityRatio, [0.0, 0.0]) #fit straight line to the data
            paramSigmas = sqrt(diag(estimate_covar(fit))) # return the sigma values for the fitted values
            push!(bfactors, -2.0*fit.param[2])
            push!(bfactorSigmas, 2.0*paramSigmas[2])
            push!(scaleFactorSigmas, abs(exp(paramSigmas[1])))
            push!(scaleFactors, exp(fit.param[1]))
            push!(imageNumArray, imageNum)
        end
        # xvals = linspace(0, maximum(hSqrd), 100)
        # yvals = model(xvals, fit.param)
        # plot(
        # layer(x=hSqrd, y=logIntensityRatio, Geom.point),
        # layer(x=xvals, y=yvals, Geom.line)
        # )
    end
    #######################################################################################
    #######################################################################################
    #NEED TO SORT THIS OUTLIER IMPLEMENTATION
    #It works for the current dataset but I'm not yet sure how general this is. I would
    #really like to do some non-parametric regression but my brief play with it wasn't as
    #successful as I would've liked.

    ##### Remove troublesome data #####
    #For b factors we assume they come from a normal distribution. so we can calculate the
    #z-scores. We'll remove the ones that are more than 2 (i.e. more than 2 standard
    #deviations from the mean. This corresponds to keeping data that are 95% likely to be
    #observed if we repeatedly performed the experiment - the frequentist view).
    zScoresB = (bfactors - mean(bfactors))./std(bfactors)

    #For the scale factors We'll use a different method. We'll use a (made up) kernel density
    #method. Basically we fit a kernel density to the data and then decide that we'll only
    #keep a certain percentage of the data. By default this is set to 90% but the user can
    #decide this.
    kdeScaleFacsFull = kde(scaleFactors, 0:kdeStep:maximum(scaleFactors))
    cumProbScaleFull = cumsum(kdeScaleFacsFull.density)/sum(kdeScaleFacsFull.density)
    outlierThresholdValue = kdeScaleFacsFull.x[findfirst(cumProbScaleFull .> keepPercentageScaleData, true)]

    #Find indices where outliars appear in both the bfactor and scale factor arrays
    indicesOfOutliers = Vector{UInt16}()
    for i in 1:length(bfactors)
        if abs(zScoresB[i]) > 2 || scaleFactors[i] > outlierThresholdValue
            push!(indicesOfOutliers, i)
        end
    end

    ### Plot distributions before the removal of outliers so we can see the full
    #distributions.
    n = 4 # We need four colours for the different elements of all of the plots in total
    #Now we generate the distiguishable colours.
    getColors = distinguishable_colors(n, Color[LCHab(70, 60, 240)],
                   transform=c -> deuteranopic(c, 0.5),
                   lchoices=Float64[65, 70, 75, 80],
                   cchoices=Float64[0, 50, 60, 70],
                   hchoices=linspace(0, 330, 24))

    #Plot the B factor distribtion
    plotTitleB = @sprintf("B-factor distribution before outlier removal")
    bfactPltBefore = plot(
    layer(x=imageNumArray, y=bfactors, Geom.point, Theme(default_color=getColors[2])),
    Guide.xlabel("Image Number"), Guide.ylabel("B-factor"),
    Guide.manual_color_key("Colour Key",["Data"],[getColors[2]]),
    Guide.title(plotTitleB)
    )

    #Plot the scale factor distribution
    plotTitleScale = @sprintf("Scale factor distribution before outlier removal")
    scalefactPltBefore = plot(
    layer(x=imageNumArray, y=scaleFactors, Geom.point, Theme(default_color=getColors[2], highlight_width=0px)),
    Guide.xlabel("Image Number"), Guide.ylabel("Scale factor"),
    Guide.manual_color_key("Colour Key",["Data"],[getColors[2]]),
    Guide.title(plotTitleScale),
    )

    if isempty(outputImageDir)
        outputBPltBeforeLocation = @sprintf("BFac_Plot_Before_outlier_removal.pdf")
        outputScalePltBeforeLocation = @sprintf("ScaleFac_Plot_Before_outlier_removal.pdf")
    else
        outputBPltBeforeLocation = @sprintf("%s/BFac_Plot_Before_outlier_removal.pdf", outputImageDir)
        outputScalePltBeforeLocation = @sprintf("%s/ScaleFac_Plot_Before_outlier_removal.pdf", outputImageDir)
    end
    draw(PDF(outputBPltBeforeLocation, 16cm, 9cm), bfactPltBefore)
    draw(PDF(outputScalePltBeforeLocation, 16cm, 9cm), scalefactPltBefore)
    if displayPlots
        display(bfactPltBefore)
        display(scalefactPltBefore)
    end

    ### Remove the outliers
    for index in reverse(indicesOfOutliers)
        splice!(bfactors, index)
        splice!(scaleFactors, index)
        splice!(imageNumArray, index)
    end

    #######################################################################################
    #Fit a linear model to both the B factor values
    modelB(x,parameters) = parameters[2] + parameters[1]*x
    bFacFit = curve_fit(modelB, imageNumArray, bfactors, [0.0, mean(bfactors)])
    bfacGrad::Float32 = bFacFit.param[1]
    bParamSigmas = sqrt(diag(estimate_covar(bFacFit)))
    bGradSigma::Float32, bInterceptSigma::Float32 = bParamSigmas[1], bParamSigmas[2]
    bIntercept::Float32 = bFacFit.param[2]

    #Use kernel density estimation to generate the distribution of scale factors
    kdeScaleFacs = kde(scaleFactors, 0:kdeStep:maximum(scaleFactors))
    modalScale::Float32 = kdeScaleFacsFull.x[findlast(kdeScaleFacs.density, maximum(kdeScaleFacs.density))]
    meanScale::Float32 = mean(scaleFactors)
    sigmaScale::Float32 = std(scaleFactors)

    ########################################################################################
    #Plot the results
    xvals = linspace(0, length(imageArray), 100)
    yvalsB = modelB(xvals, bFacFit.param)

    ### plot B factor graph
    plotTitleB = @sprintf("B factors after outlier removal. Gradient: %.4f, Intercept: %.4f", bFacFit.param[1], bFacFit.param[2])
    bfactPlt = plot(
    layer(x=xvals, y=yvalsB, Geom.line, Theme(default_color=getColors[1])),
    layer(x=imageNumArray, y=bfactors, Geom.point, Theme(default_color=getColors[2])),
    Guide.xlabel("Image Number"), Guide.ylabel("B-factor"),
    Guide.manual_color_key("Colour Key",["Line of best fit", "Data"],[getColors[1],getColors[2]]),
    Guide.title(plotTitleB)
    )

    ### plot Scale factor graph
    plotTitleScale = @sprintf("Scale Factors after outlier removal. modal value: %.4f, mean value: %.4f", modalScale, meanScale)
    scalefactPlt = plot(
    layer(x=xvals, y=modalScale*ones(length(xvals)), Geom.line, Theme(default_color=getColors[3], line_width=2px)),
    layer(x=xvals, y=meanScale*ones(length(xvals)), Geom.line, Theme(default_color=getColors[4], line_width=2px)),
    layer(x=imageNumArray, y=scaleFactors, Geom.point, Theme(default_color=getColors[2])),
    Guide.xlabel("Image Number"), Guide.ylabel("Scale factor"),
    Guide.manual_color_key("Colour Key",["Modal Line", "Mean Line", "Data"],[getColors[3],getColors[4],getColors[2]]),
    Guide.title(plotTitleScale)
    )

    ###########################################################################
    #"Undo" the damage on the b factors (assuming linear evolution)
    correctedBfactors = bfactors - bfacGrad * imageNumArray
    ### plot B factor Histogram factor graph
    normBdist = fit(Normal, correctedBfactors)
    bvals = linspace(minimum(correctedBfactors), maximum(correctedBfactors), 100)
    bDistVals = pdf(normBdist, bvals)
    plotTitleBdist = @sprintf("Distribution of damage corrected B factors. Mean value: %.4f", mean(correctedBfactors))
    pltBDist = plot(
    layer(x=bvals, y=bDistVals, Geom.line, Theme(default_color=getColors[2], line_width=2px)),
    layer(x=correctedBfactors, Geom.histogram(bincount=100, density=true)),
    Guide.xlabel("B factor value"), Guide.ylabel("Density value"),
    Guide.manual_color_key("Colour Key",["Data", "Fitted Normal Distribution"],[getColors[1],getColors[2]]),
    Guide.title(plotTitleBdist)
    )

    # ### QQplot
    # normalLine(x) = x
    # zScoresCorrB = (correctedBfactors - mean(correctedBfactors))./std(correctedBfactors)
    # Gadfly.plot(qq::QQPair) = Gadfly.plot(
    # layer(x=qq.qx, y=qq.qy, Geom.point, Theme(highlight_width=0px)),
    # layer(normalLine,-3,3, Theme(default_color=getColors[2], line_width=2px)),
    # Guide.xlabel("Theoretical Normal quantiles"), Guide.ylabel("Corrected B factor quantiles"),
    # Guide.title("QQplot for the corrected B factors")
    # )
    # qqplot(x, y) = Gadfly.plot(qqbuild(x, y))
    # qqnorm(x) = qqplot(Normal(), x)
    # qqB = qqnorm(zScoresCorrB)

    ### plot Scale factor Histogram factor graph
    pltScaleDist = plot(
    layer(x=collect(kdeScaleFacs.x), y=kdeScaleFacs.density, Geom.line, Theme(default_color=getColors[2], line_width=2px)),
    layer(xintercept=[modalScale], Geom.vline, Theme(default_color=getColors[3], line_width=2px)),
    layer(xintercept=[meanScale], Geom.vline, Theme(default_color=getColors[4], line_width=2px)),
    layer(x=scaleFactors, Geom.histogram(bincount=100, density=true)),
    Guide.xlabel("Scale factor value"), Guide.ylabel("Density value"),
    Guide.manual_color_key("Colour Key",["Data", "Kernel Density Estimation","Mode","Mean"],[getColors[1],getColors[2], getColors[3], getColors[4]]),
    Guide.title("Distribution of Scale factors")
    )

    if isempty(outputImageDir)
        outputBPltAfterLocation = @sprintf("BFac_Plot_After_outlier_removal.pdf")
        outputScalePltAfterLocation = @sprintf("ScaleFac_Plot_After_outlier_removal.pdf")
        outputBDistLocation = @sprintf("BFac_Distribution.pdf")
        outputScaleDistLocation = @sprintf("ScaleFac_Distribution.pdf")
        # outputBQQpltLocation = @sprintf("BFac_QQplot.pdf")
    else
        outputBPltAfterLocation = @sprintf("%s/BFac_Plot_After_outlier_removal.pdf", outputImageDir)
        outputScalePltAfterLocation = @sprintf("%s/ScaleFac_Plot_After_outlier_removal.pdf", outputImageDir)
        outputBDistLocation = @sprintf("%s/BFac_Distribution.pdf", outputImageDir)
        outputScaleDistLocation = @sprintf("%s/ScaleFac_Distribution.pdf", outputImageDir)
        # outputBQQpltLocation = @sprintf("%s/BFac_QQplot.pdf", outputImageDir)
    end
    draw(PDF(outputBPltAfterLocation, 16cm, 9cm), bfactPlt)
    draw(PDF(outputScalePltAfterLocation, 16cm, 9cm), scalefactPlt)
    draw(PDF(outputBDistLocation, 16cm, 9cm), pltBDist)
    draw(PDF(outputScaleDistLocation, 16cm, 9cm), pltScaleDist)
    # draw(PDF(outputBQQpltLocation, 16cm, 9cm), qqB)

    if displayPlots
        display(bfactPlt)
        display(pltBDist)
        # display(qqB)
        display(scalefactPlt)
        display(pltScaleDist)
    end
    return bfacGrad, bGradSigma, bIntercept, bInterceptSigma, modalScale, sigmaScale
end

################################################################################
#Add method information
################################################################################
function calcTempAndSFMultFactorDict(scatteringAngles::Vector{Float32}, bFactor::Float32, bFacChange::Float32, wavelength::Float32)
    tempFacDict = Dict{Float32, Float32}()
    SFMultiplierDict = Dict{Float32, Float32}()
    for scatAngle in scatteringAngles
        B = bFactor
        θ = deg2rad(scatAngle)
        λ = wavelength
        ΔB = bFacChange
        tempFacDict[scatAngle] = exp(-2 * B * (sin(θ)^2) / λ^2)
        SFMultiplierDict[scatAngle] = exp(-2 * ΔB * (sin(θ)^2) / λ^2)
    end
    return tempFacDict, SFMultiplierDict
end

################################################################################
#Add method information
################################################################################
calcD(ΔB::Float32, θ::Float32, λ::Float32) = exp(-2 * ΔB * (sin(deg2rad(θ))^2) / λ^2)
