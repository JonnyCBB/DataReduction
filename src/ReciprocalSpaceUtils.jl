immutable Unitcell{T<:Float64}
    a::T
    b::T
    c::T
    α::T
    β::T
    γ::T
end

type ReflectionObservation
    rotCentroid::Float64
    fractionCalc::Float64
    misym::Int64
    intensity::Float64
    sigI::Float64
    imageNums::Array{Uint16}
end

type DiffractionImage
    rotAngleStart::Float32
    rotAngleStop::Float32
    observationList::Dict{Vector{Int64}, ReflectionObservation}
    scaleFactor::Float64
end
DiffractionImage(startAng::Float64, stopAng::Float64) = DiffractionImage(startAng,
stopAng, Dict{Vector{Int64}, ReflectionObservation}(), 0.0)

immutable SpaceGroup
    symbol::ASCIIString
    number::Uint8
    numSymOps::Uint8
    symOps::Vector{Matrix{Float64}}
end

type Reflection
    hkl::Vector{Int64}
    resolution::Float64
    scatteringAngle::Float64
    epsilon::Uint8
    isCentric::Bool
    isSysAbs::Bool
    amplitude::Float64
    observations::Vector{ReflectionObservation}

    function Reflection(hkl::Vector{Int64}, spacegroup::SpaceGroup,
        unitcell::Unitcell, xrayWavelength::Float64)
        reflection = new()
        reflection.hkl = hkl
        recipspaceunitcell = convertDirectToReciprocal(unitcell)
        reflection.resolution = calcResolution(recipspaceunitcell,hkl)
        reflection.scatteringAngle = calcScatteringAngle(reflection.resolution,
        xrayWavelength)
        determineHKLClass!(reflection, spacegroup, unitcell)
        reflection.amplitude = 0.0
        reflection.observations = Vector{ReflectionObservation}()
        return reflection
    end
end

type ResBin
    minRes::Float64
    maxRes::Float64
    numOfRef::Int64
    meanIntensity::Float64
    hkls::Array{Vector{Int64},1}
end
ResBin(minRes::Float64, maxRes::Float64, numOfRef::Int64) = ResBin(minRes, maxRes, numOfRef, 0.0, Array(Vector{Int64},1))


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

function calcResolution(reciprocalUnitcell::Unitcell, millerIndices::Vector{Int64})
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

function calcScatteringAngle(resolution::Float64, xrayWavelength::Float64)
    if isinf(resolution)
        return 0
    else
        d = resolution
        λ = xrayWavelength
        θ = rad2deg(asin(λ/2d)) #calculate scattering angle from Bragg's Law
        return θ
    end
end

function getAllScatteringAngles(reflectionList::Dict{Vector{Int64},Reflection})
    scatteringAngles = Array(Float64,length(reflectionList))
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

function applySymOpTrans(transformation::Array{Float64,2}, point::Vector{Int64}, unitcell::Unitcell)
    transformationMatrix = transformation[1:3,1:3]
    translation = transformation[1:3,4]
    if sort(translation)[end] > 1.0
        relativeTranslation = [translation[1]/unitcell.a, translation[2]/unitcell.b, translation[3]/unitcell.c]
    else
        relativeTranslation = translation
    end
    pointTransformation = transformationMatrix * point
    pointTransformation = convert(Array{Int64,1},round(pointTransformation))
    phaseShift = 2π*dot(relativeTranslation, point)
    return pointTransformation, phaseShift
end

function removeSysAbs!(reflectionList::Dict{Vector{Int64},Reflection})
    for key in keys(reflectionList)
        if reflectionList[key].isSysAbs
            delete!(reflectionList, key)
        end
    end
end

function sortResolutionBins(reflectionList::Dict{Vector{Int64},Reflection}, minRefPerBin::Int64)
    resbins = ResBin[]
    counter = 1
    numOfRef = length(reflectionList)
    resolutionArray = Array(Float64,numOfRef)

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
            resbin = ResBin(minRes, maxRes, numOfRefInResBinCounter)
            push!(resbins, resbin)
            minRes = previousRes
            numOfRefInResBinCounter = 0
        end
        numOfRefInResBinCounter += 1
        previousRes = resolutionArray[i]
    end
    return resbins
end

function sortHKLIntoResBins!(resbins::Array{ResBin,1}, reflectionList::Dict{Vector{Int64},Reflection})
    minRes = resbins[1].minRes
    for resbin in resbins
        resbin.hkls = Array(Vector{Int64}, resbin.numOfRef)
        refCounter = 1
        for key in keys(reflectionList)
            reflection = reflectionList[key]
            if reflection.hkl == [0,0,0] && resbin.minRes == Inf
                resbin.hkls[refCounter] = reflection.hkl
                if refCounter == resbin.numOfRef
                    break
                else
                    refCounter += 1
                end
            elseif resbin.minRes > reflection.resolution ≥ resbin.maxRes || reflection.resolution == minRes
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

################################################################################
#NEED TO COME BACK TO THIS LATER
################################################################################
function assignMeanIntensityToResBins!(resbins::Array{ResBin,1}, reflectionList::Dict{Vector{Int64},Reflection}, f0SqrdDict::Dict{Float64, Float64})
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
