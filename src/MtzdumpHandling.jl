#########################################################################
########################## Types Declarations ###########################
immutable MtzdumpParams
    inputFilename::ASCIIString
    nref::Int32
end

#########################################################################
############################## Functions ################################
"""
# Run MTZDump

    runMtzdump(numRef, mtzFile)

Takes the number of reflections to be output by MTZ Dump and the location of the MTZ file (`mtzFile`) that is to be read.
The function then uses these inputs to create an input file for MTZDump and runs the program to obtain the output.
The output of the program is then returned.

## Input Parameters

    numRef
Is an integer number referring to the number of reflections that is to be listed by MTZ Dump. If we want MTZ dump to return all of the reflections in the MTZ file then the input should be '-1'.

    mtzFile
This is ASCIIString type object that contains the path to the MTZ file which is to be read.

## Return Parameters
    mtzdumpOutput
This is an ASCIIString type object that contains the entire output of MTZDump.

## Example Usage
Here we will create a temporary input file with the name `mtzInput.txt` and we will tell the program MTZDump that we want 20 reflections.
Now we use the function runMtzdump on an MTZ file with the name "MyMTZFile.mtz"

    mtzOutput = runMtzdump("MyMTZFile.mtz", 20)

`mtzOutput` is now an ASCIIString type object containing the output from the MTZDump run.
"""
function runMtzdump(mtzFile::ASCIIString, numRef::Int32=Int32(20))
    inputParams = MtzdumpParams("mtzdumpinputs.txt", numRef) #Create MtzdumpParams object
    mtzdumpInputFile = open(inputParams.inputFilename,"w") #Create text file for writing
    write(mtzdumpInputFile, @sprintf("nref %d\r\n", inputParams.nref)) # write number of reflections line
    write(mtzdumpInputFile, "symmetry \r\n") # write line to return symmetries
    write(mtzdumpInputFile, "batch \r\n") # write line to return batch/image information
    write(mtzdumpInputFile, "end \r\n") # write "end" line
    close(mtzdumpInputFile) # close file

    #@pyimport doesn't automatically look for python modules in the current
    #directory so this statement checks to see if the current working directory
    #is one of the paths that will be searched. If it isn't then it adds the
    #current directory to the vector of directories that @pyimport searches.
    if !in(pwd(), PyVector(pyimport("sys")["path"]))
        unshift!(PyVector(pyimport("sys")["path"]), pwd())
    end
    #import python function to run MTZ dump. (This is horrible hack because I couldn't work out how to do it in Julia)
    @pyimport RunSystemCommand as runsys

    #Run MTZDump and return the output
    @printf("Running MTZ dump...\n")
    tic()
    mtzdumpOutput = runsys.run_system_command(@sprintf("mtzdump hklin %s < %s", mtzFile, inputParams.inputFilename))
    @printf("Finished running MTZ dump - ")
    toc()
    #If the input file exists (it should) then delete it because we don't need it anymore
    if isfile(inputParams.inputFilename)
        rm(inputParams.inputFilename)
    end

    return mtzdumpOutput # return MTZDump output
end

################################################################################
#NEED TO ADD METHOD INFORMATION
################################################################################
function combineObsIntensity(Ipr::Float32, Isum::Float32, LP::Float32)
    if abs(Ipr + Isum) < 0.001
        return Ipr
    else
        Imid = (Ipr + Isum)/2.0
        Iraw = Imid * LP
        Ipower = 3
        w = 1.0 / (1.0 + (Iraw/Imid)^Ipower)
        Icombined = w*Ipr + (1-w)*Isum
        return Icombined
    end
end

################################################################################
#NEED TO ADD METHOD INFORMATION
################################################################################
function parseMosflmMTZDumpOutput(mtzDumpOutput::ASCIIString, imageOsc::Float32=Float32(0.0), rotDiffTol::Float32=Float32(0.1))
    ################################################################################
    #Parameter types can be annotated if necessary to save memory if that
    #becomes an issue.
    ################################################################################
    hklList = Dict{Vector{Int16},Reflection}()

    batchNumber::UInt16 = 0
    numberOfImages::UInt16 = 0
    actualImageOsc::Float32 = 0.0
    rotStart::Float32 = 0.0
    rotFinish::Float32 = 0.0
    numSymOps::UInt8 = 0
    spaceGroupNumber::UInt8 = 0
    spaceGroupSymbol::ASCIIString = ""
    symmetryOps = Array(Matrix{Float32},1)
    symOpsLines::Bool = false
    symOpMatrix = Array(Float32,4,4)
    symOpMatrixRow::UInt8 = 1
    symOpNumber::UInt8 = 0
    searchCellDims::Bool = false
    searchReflections::Bool = false
    colNumH::UInt8 = 1
    colNumK::UInt8 = 2
    colNumL::UInt8 = 3
    colNumMIsym::UInt8 = 4
    colNumBatch::UInt8 = 5
    colNumIsum::UInt8 = 6
    colNumSigIsum::UInt8 = 7
    colNumIpr::UInt8 = 8
    colNumSigIpr::UInt8 = 9
    colNumFracCalc::UInt8 = 10
    colNumRot::UInt8 = 13
    colNumLP::UInt8 = 15
    refLine::UInt8 = 0
    hkl = Vector{Int16}(3)
    origHKL = Vector{Int16}(3)
    imageNumber::UInt16 = 0
    misymNum::UInt16 = 0
    Isum::Float32, sigIsum::Float32 = 0.0, 0.0
    Ipr::Float32, sigIpr::Float32 = 0.0, 0.0
    fractionCalc::Float32 = 0.0
    rot::Float32 = 0.0
    LP::Float32 = 0.0

    #Split the MTZ Dump output log by the newline character
    mtzdumpOutputLines = split(mtzDumpOutput, "\n")

    for line in mtzdumpOutputLines
        ############################################################################
        # Section: Parse number of symmetry Operations
        #---------------------------------------------------------------------------
        #In this section of code we extract the number of symmetry operators for
        #the space group. This allows us to use that number to preallocate an array
        #of the right size for all of the symmetry operators
        if contains(line, "Number of Symmetry Operations")
            numSymOps = parse(UInt8, split(line)[7])
            symmetryOps = Array(Matrix{Float16}, numSymOps)
        end
        #End of Section: Parse number of symmetry Operations
        ############################################################################


        ############################################################################
        # Section: Parse Space Group
        #---------------------------------------------------------------------------
        #In this section of code we want to extract information about the space
        #group. In the following "if...end" block we extract the space group number
        #and symbol from the corresponding lines in the MTZ dump output.
        if contains(line, "Space Group")
            spaceGroupLine = split(line)
            spaceGroupNumber = parse(UInt8, spaceGroupLine[5])
            spaceGroupSymbol = spaceGroupLine[6][2:end-1]
        end

        #Here we extract the symmtery operators (4x4 matrices containing the 3x3
        #rotation matrix and the 3x1 translation vector) for each symmetry operator
        #and store it in the corresponding element of the symmetry operator array.
        if search(line, r"Symmetry [0-9]* (-|X|Y|Z)") != 0:-1
            symOpsLines = true
            symOpNumber = parse(Int, split(line)[2])
        elseif symOpsLines == true
            matrixRowValues = split(line)
            symOpMatrix[symOpMatrixRow,:] = [parse(Float32,matrixRowValues[1]), parse(Float32,matrixRowValues[2]),
                                            parse(Float32,matrixRowValues[3]), parse(Float32,matrixRowValues[4])]
            if symOpMatrixRow == 4
                symmetryOps[symOpNumber] = symOpMatrix
                symOpsLines = false
                symOpMatrixRow = 1
                #If the current symmetry operator is the final one and we have now
                #added it to the symmetry operator array then we can now create the
                #space group object
                if symOpNumber == numSymOps
                    global spacegroup = SpaceGroup(spaceGroupSymbol, spaceGroupNumber,
                    numSymOps, symmetryOps)
                end
            else
                symOpMatrixRow += 1
            end
        end
        #End of Section: Parse Space Group
        ############################################################################

        ############################################################################
        # Section: Parse Unit Cell
        #---------------------------------------------------------------------------
        #In this section we extract the unit cell parameters for the crystal.
        #First we look for "Dataset ID" in the line. When we reach this line, we
        #know that the Unit cell parameters are going to be given so we set the
        #"searchCellDims" paramter to true so we start looking for the unit cell
        #lines
        if contains(line, "Dataset ID")
           searchCellDims = true
        end

        #When we reach the line cotaining the unit cell params we extract them.
        #Finally we create the Unitcell object
        if search(line, r"[0-9][0-9].[0-9][0-9][0-9][0-9]") != 0:-1 && searchCellDims == true
            unitcellDims = split(line)
            global unitcell = Unitcell(parse(Float32,unitcellDims[1]), parse(Float32,unitcellDims[2]),
            parse(Float32,unitcellDims[3]), parse(Float32,unitcellDims[4]),
            parse(Float32,unitcellDims[5]), parse(Float32,unitcellDims[6]))
            searchCellDims = false
        end
        #End of Section: Parse Unit Cell
        ############################################################################


        ############################################################################
        # Section: Parse batch/image phi angle information
        #---------------------------------------------------------------------------
        #In this section we extract the start and stop phi angle information for
        #each image. First we obtain the number of images so we can preallocate the
        #image array.
        if contains(line, "Number of Batches")
            numberOfImages = parse(UInt16, split(line)[6])
        end

        #Here we look for the line containing the start and stop phi angles. When we
        #find this line we extract the relevant angle information - this includes the
        #initial rotation angle, the final rotation angle and finally the oscillation
        #per image if it hasn't been supplied as an argument for the function.
        if contains(line, "Start & stop Phi angles (degrees)")
            batchNumber += 1
            if batchNumber == 1
                phiAngleInfoLine = split(line)
                rotStart = parse(Float32, phiAngleInfoLine[7])
            elseif batchNumber == numberOfImages
                phiAngleInfoLine = split(line)
                rotFinish = parse(Float32, phiAngleInfoLine[8])
            end
            #Calculate the actual rotation per image
            if actualImageOsc == 0.0
                phiAngleInfoLine = split(line)
                startAng, stopAng = parse(Float32, phiAngleInfoLine[7]), parse(Float32, phiAngleInfoLine[8])
                actualImageOsc = stopAng - startAng
            end
            #If the user hasn't defined the actual image oscillation then set it
            #equal to the actual image oscillation.
            if imageOsc == 0.0
                imageOsc = actualImageOsc
            end
        end
        #End of Section: Parse batch/image phi angle information
        ############################################################################


        #When we reach this line we know that there are no more column labels to
        #look for so we tell the program to stop looking for column labels.
        if contains(line, "No. of reflections used in FILE STATISTICS")
            fileStatsLines = false
        end
        #End of Section: Determine column numbers to obtain correct Reflection information
        ############################################################################


        ############################################################################
        # Section: Extract reflection data.
        #---------------------------------------------------------------------------
        #In this section we extract the data for each reflection.

        #If we come across this line then we tell the program to start looking for
        #lines containing information about reflections.
        if contains(line, "LIST OF REFLECTIONS")
            searchReflections = true
        end

        #Here is the meat of the code...
        if searchReflections == true && !isempty(strip(line)) # Check the line is non-empty and that we're searching for reflection info
            nonEmptyLine = split(line) #split the line by whitespace
            if isnumber(nonEmptyLine[1]) # Check the first non-whitespace string can be parsed as numeric (this only works for integers)
                if length(nonEmptyLine) == numMtzColsFor1stRefLine #Check that the line is of a given length, otherwise we'll run into trouble with the parser.
                    refLine = 1
                    hkl = [parse(Int16,nonEmptyLine[colNumH]), parse(Int16,nonEmptyLine[colNumK]), parse(Int16,nonEmptyLine[colNumL])]
                    misymNum = parse(UInt16, nonEmptyLine[colNumMIsym])
                    if separateSymEquivs #Check if we want to keep symmetry equivalents separate.
                        ############################################################
                        #Mini Section: Get original HKL indices
                        #-----------------------------------------------------------
                        #We use information from the M/ISYM column to determine the
                        #original HKL indices. I have done this according to the
                        #information given on the CCP4 MTZ Format page here:
                        #http://www.ccp4.ac.uk/html/mtzformat.html
                        Iplus = false
                        isym = misymNum
                        if isym > 256
                            isym = isym - 256
                        end
                        if !iseven(round(Int,isym))
                            isym += 1
                            Iplus = true
                        end
                        symopNum = round(Int, isym/2)
                        origHKL = map(x -> Int16(x),symmetryOps[symopNum][1:3,1:3] * hkl)
                        if !Iplus
                            origHKL = -origHKL
                        end
                        #End Mini Section: Get original HKL indices
                        ############################################################
                    else
                        #if we're happy to merge the data for symmetry equivalents
                        #then we don't have to do anything to the HKL indices.
                        origHKL = hkl
                    end

                    #If the HKL indices haven't been added to the reflection
                    #dictionary then we have to add this reflection to it.
                    if !haskey(hklList, origHKL)
                        hklList[origHKL] = Reflection(origHKL, hkl, spacegroup, unitcell, xrayWavelength)
                    end
                    #Extract some important reflection information.
                    trueImageNumber = parse(UInt16, nonEmptyLine[colNumBatch])
                    ############################################################
                    #*****
                    #!!!!!!!!!!!!!!!!!!!!JONNY SORT THIS OUT!!!!!!!!!!!!!!!!!!!!
                    #!!!!!!!!!!!!!!!!!!!!JONNY SORT THIS OUT!!!!!!!!!!!!!!!!!!!!
                    #!!!!!!!!!!!!!!!!!!!!JONNY SORT THIS OUT!!!!!!!!!!!!!!!!!!!!
                    #THIS CAN CAUSE PROBLEMS WITH THE INDEXING GOING OUT OF
                    #BOUNDS WHEN ALLOCATING REFLECTIONS IMAGES IN A DOWNSTREAM
                    #METHOD. THIS HAPPENS BECAUSE THE ERROR PROPAGATION OF THE
                    #DECIMAL CALCULATIONS BECOMES TOO HIGH AND CALCULATES AN
                    #IMAGE THAT IS OUT OF BOUNDS.
                    #
                    #A QUICK FIX WOULD BE TO INCLUDE AN ERROR STATEMENT AT THE
                    #END OF THIS FUNCTION TO CATCH THIS.
                    #The rounding is set to 3 decimal places in the line below.
                    #It would be more general if it was possible to determine
                    #the precision from the user defined input. A
                    ############################################################
                    imageNumber = Int16(ceil(round(actualImageOsc/imageOsc * trueImageNumber, 1)))
                    Isum, sigIsum = parse(Float32, nonEmptyLine[colNumIsum]), parse(Float32, nonEmptyLine[colNumSigIsum])
                    Ipr, sigIpr = parse(Float32, nonEmptyLine[colNumIpr]), parse(Float32, nonEmptyLine[colNumSigIpr])

                    #From the MTZ file from the MOSFLM output it seems that the
                    #sigma values for both the summed and profile fitted intensities
                    # are exactly the same. So I'm using the summed sigma intensity
                    #as the 'true' sigma (given that it doesn't matter which one I
                    #choose). However if the sigma's do deviate by a significant
                    #amount I haven't taken this into account so I have assigned
                    #included the following as a warning to the user just in case
                    #they differ a lot.
                    if sigIsum - sigIpr > sigIDiffTol
                        @printf("*****************WARNING*****************\nThe sigma of the profile fitted and summed intensities for reflection (%d, %d, %d) differ by a value greater than %0.2f\nSigIpr = %0.3f\nSigIsum = %0.3f.\nUsing the sigma of the summed intensity...\n\n", origHKL[1], origHKL[2], origHKL[3], sigIDiffTol, sigIpr, sigIsum)
                    end
                else
                    #If the number of columns is not the same as the expected ones then this throws an error because the parser will fail in that case.
                    error("The MTZ Dump output doesn't have %d columns for the reflection line.\nThis means the reflections in the file will not be parsed properly\nContact elspeth.garman@bioch.ox.ac.uk to sort out the MTZ Dump parser for your case.\n\n", numMtzColsFor1stRefLine)
                end
            elseif contains(nonEmptyLine[1], ".")
                #Check that the number of columns for the reflection information is correct.
                if length(nonEmptyLine) ≤ numMtzColsFor2ndand3rdRefLines
                    #Because information about a single reflection is stored on
                    #multiple lines in the MTZ Dump output we need to keep track of
                    #of the line number so we have to increment it.
                    refLine += 1
                    if refLine == 2
                        #Extract the relevant information in this line.
                        fractionCalc, rot = parse(Float32, nonEmptyLine[colNumFracCalc - numMtzColsFor1stRefLine]), parse(Float32, nonEmptyLine[colNumRot - numMtzColsFor1stRefLine])
                    elseif refLine == 3
                        #Extract the Lorentz-Polarisation correction factor from
                        #this line
                        LP = parse(Float32, nonEmptyLine[colNumLP - (numMtzColsFor1stRefLine + numMtzColsFor2ndand3rdRefLines)])

                        ############################################################
                        #Mini Section: Create observation object for reflection
                        #-----------------------------------------------------------
                        #In this mini section we use the information that we've
                        #extracted about the current reflection to update/create
                        #an observation object - i.e. an object that contains
                        #information about this particular observation of the
                        #reflection.

                        #Here we decide whether to use the summed, profile fiited or
                        #combined intensity.
                        if uppercase(intensityType) == "SUMMED"
                            intensity = Isum
                        elseif uppercase(intensityType) == "PROFILE"
                            intensity = Ipr
                        else
                            intensity = convert(Float32, combineObsIntensity(Ipr, Isum, LP))
                        end

                        #Here we need to check whether the current observation of
                        #the current reflection is a completely new observation in
                        #which case we create a new observation object, or whether
                        #this is a partial observation of an observation object that
                        #has already been created and we need to update that
                        #observation object.

                        #The easy case is when the observations vector for the
                        #reflection is empty. If it is that means that this is a new
                        #observation and so we create a new observation object and
                        #add it to the array.
                        if isempty(hklList[origHKL].observations)
                            push!(hklList[origHKL].observations, ReflectionObservation(rot, fractionCalc, misymNum, intensity, sigIsum^2, [imageNumber], [intensity]))
                        else
                            #If there are existing observation objects for the
                            #reflection then we need to loop through all
                            #observations to check whether this is a partial
                            #observation for an existing observation object. The
                            #criteria to decide if the partial observation is part
                            #of an existing observation object is: image at which
                            #the partial observation was observed is a consectuive
                            #image of an existing observation object AND it has
                            #the same M/ISYM number.
                            numOfExistingObs::UInt32 = length(hklList[origHKL].observations)
                            sameObservation::Bool = false
                            for obsNum = 1:numOfExistingObs #Loop through observations
                                refObservation = hklList[origHKL].observations[obsNum]
                                imageNumsOfObs = hklList[origHKL].observations[obsNum].imageNums
                                for obsImageNum in imageNumsOfObs #Loop through images
                                    if misymNum == refObservation.misym #Check that the M/ISYM number is the same
                                        if obsImageNum - 1 ≤ imageNumber ≤ obsImageNum + 1 #Check that the image is consecutive
                                            #If it is consecutive then update the corresponding observation information.
                                            sameObservation = true
                                            if isnan(rot) || isnan(intensity) || isnan(sigIsum) || isnan(fractionCalc)
                                                @printf("**************************WARNING**************************\n")
                                                @printf("Partial observation of reflection [%d,%d,%d] on image %d has NaN value.\n\n", origHKL[1], origHKL[2], origHKL[3], imageNumber)
                                            else
                                                refObservation.rotCentroid += rot
                                                refObservation.fractionCalc += fractionCalc
                                                refObservation.intensity += intensity
                                                refObservation.sigI += sigIsum^2
                                                push!(refObservation.imageNums, imageNumber)
                                                push!(refObservation.imageIntensities, intensity)
                                                hklList[origHKL].observations[obsNum] = refObservation
                                            end
                                            break
                                        end
                                    end
                                end
                                #First check if we've discovered that the current reflection record is a partial reflection whose observation object has
                                #already been created. If yes, then we don't need to check any more observations so we can break out of the loop. Otherwise
                                #if we've checked all images for all current reflection observations and found that the current image is not a consecutive
                                #image, then it's almost certainly a new observation of a reflection so we create a new ReflectionObservation object.
                                if sameObservation
                                    break
                                elseif obsNum == numOfExistingObs && !sameObservation
                                    push!(hklList[origHKL].observations, ReflectionObservation(rot, fractionCalc, misymNum, intensity, sigIsum^2, [imageNumber], [intensity]))
                                end
                            end
                        end
                        #End Mini Section: Create observation object for reflection
                        ############################################################
                    end
                else
                    error("The MTZ Dump output has more than %d columns for the reflection line.\nThis means the reflections in the file will not be parsed properly\n Contact elspeth.garman@bioch.ox.ac.uk to sort out the MTZ Dump parser for your case.", numMtzColsFor2ndand3rdRefLines)
                end
            elseif contains(line, "<B>")
                searchReflections = false
            end
        end
        #End of Section: Extract reflection data.
        ############################################################################
    end

    #Calculate number of total images (effectively these are discrete time intervals)
    #with the oscillation provided by the user.
    numEffectiveImages = Int((rotFinish - rotStart)/imageOsc)
    #Check to make sure that all reflection observations have been allocated within
    #the allowed range of images
    for hkl in keys(hklList)
        obsCounter = 0
        for refObservation in hklList[hkl].observations
            obsCounter += 1
            for imgNums in refObservation.imageNums
                if imgNums > numEffectiveImages
                    errMsg = @sprintf("Observation %d for reflection (%d, %d, %d) has been allocated to image %d.\nThis shouldn't happen because the maximum number of images according to the user input is %d.\nIf you're seeing this warning then you'll need to contact Professor Elspeth Garman: elspeth.garman@bioch.ox.ac.uk to sort this out. Sorry for the inconvenience. ", obsCounter, origHKL[1], origHKL[2], origHKL[3], imgNums, numEffectiveImages)
                    error(errMsg)
                end
            end
        end
    end

    #Finally add diffraction images to an array.
    imageArray = Vector{DiffractionImage}(numEffectiveImages)
    for imgNum = 1:numEffectiveImages
        oscStart = Float32(rotStart + (imgNum - 1) * imageOsc)
        imageArray[imgNum] = DiffractionImage(oscStart, oscStart+imageOsc)
    end

    return spacegroup, unitcell, hklList, imageArray
end

################################################################################
#NEED TO ADD METHOD INFORMATION
################################################################################
function parseCTruncateMTZDumpOutput(mtzDumpOutput::ASCIIString)
    #Split the MTZ Dump output log by the newline character
    mtzdumpOutputLines = split(mtzDumpOutput, "\n")

    searchReflections::Bool = false
    firstRef::Bool = true
    hkl = Vector{Int16}(3)
    refLine::UInt8 = 0
    refAmpDict = Dict{Vector{Int16},Vector{Float32}}()
    scaleFac::Float32 = 0.0
    amplitude::Float32 = 0.0
    sigAmp::Float32 = 0.0
    for line in mtzdumpOutputLines
        if contains(line, "LIST OF REFLECTIONS")
            searchReflections = true
        end

        if searchReflections == true && !isempty(strip(line))
            nonEmptyLine = split(line) #split the line by whitespace
            if isnumber(nonEmptyLine[1]) # Check the first non-whitespace string can be parsed as numeric (this only works for integers)
                if length(nonEmptyLine) == numMtzColsFor1stRefLine #Check that the line is of a given length, otherwise we'll run into trouble with the parser.
                    hkl = [parse(Int16,nonEmptyLine[1]), parse(Int16,nonEmptyLine[2]), parse(Int16,nonEmptyLine[3])]
                    amplitude = parse(Float32,nonEmptyLine[4])
                    ampSig = parse(Float32,nonEmptyLine[5])
                    refAmpDict[hkl] = [amplitude, ampSig]
                else
                    #If the number of columns is not the same as the expected ones then this throws an error because the parser will fail in that case.
                    errMsg = @sprintf("The MTZ Dump output doesn't have %d columns for the reflection line.\nThis means the reflections in the file will not be parsed properly\nContact elspeth.garman@bioch.ox.ac.uk to sort out the MTZ Dump parser for your case.\n\n", numMtzColsFor1stRefLine)
                    error(errMsg)
                end
            elseif contains(nonEmptyLine[1], ".") || contains(nonEmptyLine[1], "?")
                if firstRef
                    if length(nonEmptyLine) == numMtzColsIntLineCTruncate #Check that the line is of a given length, otherwise we'll run into trouble with the parser.
                        refLine += 1
                        intensity = parse(Float32, nonEmptyLine[4])
                        scaleFac = amplitude / sqrt(intensity)
                        firstRef = false
                    else
                        @printf("hkl: [%d, %d, %d]\n", hkl[1], hkl[2], hkl[3])
                        print(line)
                        #If the number of columns is not the same as the expected ones then this throws an error because the parser will fail in that case.
                        errMsg = @sprintf("The MTZ Dump output doesn't have %d columns for the reflection line.\nThis means the reflections in the file will not be parsed properly\nContact elspeth.garman@bioch.ox.ac.uk to sort out the MTZ Dump parser for your case.\n\n", numMtzColsIntLineCTruncate)
                        error(errMsg)
                    end
                end
            end
        end
    end
    return refAmpDict, scaleFac
end
