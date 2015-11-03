"""
# Get initial filtering estimate

    getInitialState(hklList, f0SqrdDict)

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

L(x) = exp(x/2) * ( (1-x)*besseli(0,-x/2) - x*besseli(1,-x/2) )
meanRice(F, D, σ) = σ * √(π/2) * L(- ((D * F)^2)/(2σ^2) )
varRice(F, D, σ) = 2σ^2 + (D * F)^2 - meanRice(F, D, σ)^2
########## create Process function ##########
function processFunction(amplitudes, D, σ)
    newAmplitudes = Vector{Float64}(length(amplitudes))
    counter = 0
    for F in amplitudes
        counter += 1
        newAmplitudes[counter] = meanRice(F, D, σ)
    end
    return newAmplitudes
end

function observationFunction(amplitudes, K)
    predObservations = Vector{Float64}(length(amplitudes))
    counter = 0
    for F in amplitudes
        counter += 1
        predObservations[counter] = K * F^2
    end
    return predObservations
end
