"""
# Create a HKL Index refence

    createHKLIndexReferenceDict(hklList::Dict{Vector{Int64},Reflection})

This function takes in the dictionary of reflections as its input and returns a dictionary that tells the user which index in an array corresponds to a given reflection.
"""
function createHKLIndexReferenceDict(hklList::Dict{Vector{Int64},Reflection})
    hklIndexReference = Dict{Vector{Int64},Int64}()
    counter = 0
    for hkl in keys(hklList)
        counter += 1
        hklIndexReference[hkl] = counter
    end
    return hklIndexReference
end

"""
# Get initial filtering estimate

    getInitialState(hklList::Dict{Vector{Int64}, Reflection}, f0SqrdDict::Dict{Float64, Float64})

This function takes in as input a reflection list dictionary `hklList` and a squared structure factor amplitude dictionary `f0SqrdDict` and returns a multivariate normal distribution representing the estimate of the initial state of the structure factor amplitudes of the crystal.
"""
function getInitialState(hklList::Dict{Vector{Int64}, Reflection}, f0SqrdDict::Dict{Float64, Float64})
    initialAmplitudes = Vector{Float64}(length(hklList))
    initialVariance = Vector{Float64}(length(hklList))
    counter = 0
    for hkl in keys(hklList)
        counter += 1
        reflection = hklList[hkl]
        initialAmplitudes[counter] = reflection.amplitude
        initialVariance[counter] = f0SqrdDict[reflection.scatteringAngle] * reflection.epsilon
    end
    return MvNormal(initialAmplitudes, diagm(initialVariance))
end

"""
# Calculate the Laguerre Polynomial (with n = 1/2)

    L(x::Float64)

This function takes in as input a `Float64` type number, `x`, and returns the Laguerre polynomial function evaluated at `x`.
"""
L(x::Real) = exp(x/2) * ( (1-x)*besseli(0,-x/2) - x*besseli(1,-x/2) )

"""
# Calculate mean of the Rician Distribution

    meanRice(a::Float64, x::Float64)

This function takes the rician parameters `a` and `σ` (both of `Float64` type) and returns the mean of the Rice distribution.
"""
meanRice(a::Float64, σ::Float64) = σ * √(π/2) * L(- (a^2)/(2σ^2) )
meanRice(F::Float64, D::Float64, σ::Float64) = meanRice(D*F, σ)

"""
# Calculate variance of the Rician Distribution

    meanRice(a::Float64, x::Float64)

This function takes the rician parameters `a` and `σ` (both of `Float64` type) and returns the variance of the Rice distribution.
"""
varRice(a::Float64, σ::Float64) = 2σ^2 + a^2 - meanRice(a, σ)^2
varRice(F::Float64, D::Float64, σ::Float64) = varRice(D*F, σ)

"""
# The Process Function

    processFunction(amplitudes::Vector{Float64}, D::Float64, σ::Float64)

This is the process function for the Kalman Filter and describes how the mean of the amplitude distribution changes after irradiation with X-rays.
`amplitudes` is a vector of `Float64`'s that represent the amplitudes for each reflection at the previous time step.
`D` is a vector of `Float64`'s that are the corresponding structure factor multipliers (Luzzati 1952) for each reflection. These are defined in Read 1990.
`σ` is a vector of `Float64`'s that represent the standard deviations of the normally distributed structure factors (not the amplitudes).
"""
function processFunction(amplitudes::Vector{Float64}, D::Vector{Float64}, σ::Vector{Float64})
    newAmplitudes = Vector{Float64}(length(amplitudes))
    counter = 0
    for F in amplitudes
        counter += 1
        ########################################################################
        #NEED TO SORT THIS OUT. I'VE HAD TO USE A GAUSSIAN DISTRIBUTION INSTEAD
        #OF A RICIAN DISTRIBUTION.
        #newAmplitudes[counter] = meanRice(F, D[counter], σ[counter])
        newAmplitudes[counter] = D[counter] * F
    end
    return newAmplitudes
end
################################################################################
#NEED TO SORT THIS OUT. I'VE HAD TO USE A GAUSSIAN DISTRIBUTION INSTEAD OF A
#RICIAN DISTRIBUTION.
#processFunction(F::Float64, D::Float64, σ::Float64) = meanRice(F, D, σ)
processFunction(F::Float64, D::Float64, σ::Float64) = D * F

"""
# The Observation Function

    observationFunction(amplitudes::Vector{Float64}, K::Float64)

This is the observation function for the Kalman Filter and describes how the intensities (observations) are generated from the set of amplitudes.
`amplitudes` is a vector of `Float64`'s that represent the amplitudes for each reflection at the previous time step.
`K` is the scale factor for the current diffraction image.
"""
function observationFunction(amplitudes::Vector{Float64}, K::Float64)
    return K * amplitudes.^2
end
