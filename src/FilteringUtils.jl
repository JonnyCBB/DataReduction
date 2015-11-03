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
