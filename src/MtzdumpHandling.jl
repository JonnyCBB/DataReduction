#########################################################################
####################### using/import statements #########################
using PyCall

#########################################################################
########################## Types Declarations ###########################
immutable MtzdumpParams
    inputFilename::ASCIIString
    nref::Int64
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
function runMtzdump(mtzFile::ASCIIString, numRef::Int64=20)
    inputParams = MtzdumpParams("mtzdumpinputs.txt", numRef) #Create MtzdumpParams object
    mtzdumpInputFile = open(inputParams.inputFilename,"w") #Create text file for writing
    write(mtzdumpInputFile, @sprintf("nref %d\r\n", inputParams.nref)) # write number of reflections line
    write(mtzdumpInputFile, "symmetry \r\n") # write line to return symmetries
    write(mtzdumpInputFile, "batch \r\n") # write line to return batch/image information
    write(mtzdumpInputFile, "end \r\n") # write "end" line
    close(mtzdumpInputFile) # close file

    #import python function to run MTZ dump. (This is horrible hack because I couldn't work out how to do it in Julia)
    @pyimport getMtzdumpOutput as getmtz

    #Run MTZDump and return the output
    mtzdumpOutput = getmtz.runMtzdump(@sprintf("mtzdump hklin %s < %s", mtzFile, inputParams.inputFilename))

    #If the input file exists (it should) then delete it because we don't need it anymore
    if isfile(inputParams.inputFilename)
        rm(inputParams.inputFilename)
    end

    return mtzdumpOutput # return MTZDump output
end

################################################################################
#NEED TO ADD METHOD INFORMATION
################################################################################
function combineObsIntensity(Ipr::Float64, Isum::Float64, LP::Float64)
    Imid = (Ipr + Isum)/2
    Iraw = Imid * LP
    Ipower = 3
    w = 1 / (1 + (Iraw/Imid)^Ipower)
    Icombined = w*Ipr + (1-w)*Isum
    return Icombined
end

################################################################################
#NEED TO ADD METHOD INFORMATION
################################################################################
#This p(q) function is from Rossman et al. 1979. My volumeRatio function is an
#exact volume ratio for the partial reflection. However if you plot these
#functions it turns out that they pretty much overlap so they're pretty much
#equivalent for argument values between 0 - 1
p(q) = 3q^2 - 2q^3
function volumeRatio(penetrationDepth::Float64, sphereDiameter::Float64)
    d = sphereDiameter
    h = penetrationDepth
    ratio = 8 * h^2 * ((3/2)*d - h)/(4 * d^3)
    return ratio
end
volumeRatio(penetrationFraction::Float64) = volumeRatio(penetrationFraction, 1.0)

################################################################################
#NEED TO ADD METHOD INFORMATION
################################################################################
function parseMTZDumpOutput(mtzDumpOutput::ASCIIString, rotDiffTol::Float64=0.1)
    ################################################################################
    #Parameter types can be annotated if necessary to save memory if that
    #becomes an issue.
    ################################################################################
    hklList = Dict{Vector{Int64},Reflection}()

    batchNumber = 0
    numSymOps = 0
    spaceGroupNumber = 0
    spaceGroupSymbol = ""
    symmetryOps = Array(Matrix{Float16})
    symOpsLines = false
    symOpMatrix = [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]
    symOpMatrixRow = 1
    symOpNumber = 0
    searchCellDims = false
    searchReflections = false
    fileStatsLines = false
    colNumH = 0
    colNumK = 0
    colNumL = 0
    colNumMIsym = 0
    colNumBatch = 0
    colNumIsum = 0
    colNumSigIsum = 0
    colNumIpr = 0
    colNumSigIpr = 0
    colNumFracCalc = 0
    colNumRot = 0
    colNumLP = 0
    refLine = 0
    origHKL = [0, 0, 0]
    imageNumber = 0
    misymNum = 0
    Isum, sigIsum = 0.0, 0.0
    Ipr, sigIpr = 0.0, 0.0
    fractionCalc = 0.0
    rot = 0.0
    LP = 0.0

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
            numSymOps = parse(Int, split(line)[7])
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
            spaceGroupNumber = parse(Int, spaceGroupLine[5])
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
            symOpMatrix[symOpMatrixRow,:] = [parse(Float64,matrixRowValues[1]), parse(Float64,matrixRowValues[2]),
                                             parse(Float64,matrixRowValues[3]), parse(Float64,matrixRowValues[4])]
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
            global unitcell = Unitcell(parse(Float64,unitcellDims[1]), parse(Float64,unitcellDims[2]),
            parse(Float64,unitcellDims[3]), parse(Float64,unitcellDims[4]),
            parse(Float64,unitcellDims[5]), parse(Float64,unitcellDims[6]))
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
            numberOfImages = parse(Int, split(line)[6])
            global imageArray = Array(DiffractionImage,numberOfImages)
        end

        #Here we look for the line containing the start and stop phi angles. When we
        #find this line we extract the relevant angle information and create a
        #diffraction image object in the corresponding element of the image array.
        if contains(line, "Start & stop Phi angles (degrees)")
            batchNumber += 1
            phiAngleInfoLine = split(line)
            startAng, stopAng = parse(Float64, phiAngleInfoLine[7]), parse(Float64, phiAngleInfoLine[8])
            imageArray[batchNumber] = DiffractionImage(startAng, stopAng)
        end
        #End of Section: Parse batch/image phi angle information
        ############################################################################


        ############################################################################
        # Section: Determine column numbers to obtain correct Reflection information
        #---------------------------------------------------------------------------
        #In this section we extract the column numbers for which the information
        #corresponding for the column label can be found in the MTZ file.
        #This information can be found in the section titled "OVERALL FILE
        #STATISTICS" in the MTZ Dump output and so if we see this line then we set
        #the fileStatsLines (boolean) variable to 'true' to tell the program that
        #we should start looking for the column label lines.
        if contains(line, "OVERALL FILE STATISTICS")
            fileStatsLines = true
        end

        #When the program is ready to look for the lines containing the column
        #numbers for the corresponding column labels it checks to make sure that the
        #first non whitespace character is a number. If so then we know that this
        #line should tell us the column number for the column label. The column
        #label is always the last non-whitespace block, whereas the corresponding
        #column number is always the first.
        if !isempty(strip(line)) && fileStatsLines # Check that the line isn't empty and that we should be looking for column labels
            nonEmptyLine = split(line) #split the split by whitespace so we have an array containing Substring objects
            if ismatch(r"[0-9]",convert(ASCIIString, nonEmptyLine[1])) # Check that the first character of the line is numeric
                columnLabel = convert(ASCIIString, nonEmptyLine[end]) #Extract the column label
                #The following code checks the type of the column label and assigns
                #the corresponding column number
                if "H" == columnLabel
                    colNumH = parse(Int, nonEmptyLine[1])
                elseif "K" == columnLabel
                    colNumK = parse(Int, nonEmptyLine[1])
                elseif "L" == columnLabel
                    colNumL = parse(Int, nonEmptyLine[1])
                elseif "M/ISYM" == columnLabel
                    colNumMIsym = parse(Int, nonEmptyLine[1])
                elseif "BATCH" == columnLabel
                    colNumBatch = parse(Int, nonEmptyLine[1])
                elseif "I" == columnLabel
                    colNumIsum = parse(Int, nonEmptyLine[1])
                elseif "SIGI" == columnLabel
                    colNumSigIsum = parse(Int, nonEmptyLine[1])
                elseif "IPR" == columnLabel
                    colNumIpr = parse(Int, nonEmptyLine[1])
                elseif "SIGIPR" == columnLabel
                    colNumSigIpr = parse(Int, nonEmptyLine[1])
                elseif "FRACTIONCALC" == columnLabel
                    colNumFracCalc = parse(Int, nonEmptyLine[1])
                elseif "ROT" == columnLabel
                    colNumRot = parse(Int, nonEmptyLine[1])
                elseif "LP" == columnLabel
                    colNumLP = parse(Int, nonEmptyLine[1])
                end
            end
        end

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
                    hkl = [parse(Int,nonEmptyLine[colNumH]), parse(Int,nonEmptyLine[colNumK]), parse(Int,nonEmptyLine[colNumL])]
                    misymNum = parse(Float64, nonEmptyLine[colNumMIsym])
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
                        origHKL = round(Int,symmetryOps[symopNum][1:3,1:3] * hkl)
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
                        hklList[origHKL] = Reflection(origHKL, spacegroup, unitcell, xrayWavelength)
                    end
                    #Extract some important reflection information.
                    imageNumber = parse(Int, nonEmptyLine[colNumBatch])
                    Isum, sigIsum = parse(Float64, nonEmptyLine[colNumIsum]), parse(Float64, nonEmptyLine[colNumSigIsum])
                    Ipr, sigIpr = parse(Float64, nonEmptyLine[colNumIpr]), parse(Float64, nonEmptyLine[colNumSigIpr])

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
                if length(nonEmptyLine) <= numMtzColsFor2ndand3rdRefLines
                    #Because information about a single reflection is stored on
                    #multiple lines in the MTZ Dump output we need to keep track of
                    #of the line number so we have to increment it.
                    refLine += 1
                    if refLine == 2
                        #Extract the relevant information in this line.
                        fractionCalc, rot = parse(Float64, nonEmptyLine[colNumFracCalc - numMtzColsFor1stRefLine]), parse(Float64, nonEmptyLine[colNumRot - numMtzColsFor1stRefLine])
                    elseif refLine == 3
                        #Extract the Lorentz-Polarisation correction factor from
                        #this line
                        LP = parse(Float64, nonEmptyLine[colNumLP - (numMtzColsFor1stRefLine + numMtzColsFor2ndand3rdRefLines)])

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
                        if intensityType == "Summed"
                            intensity = Isum
                        elseif intensityType == "Profile"
                            intensity = Ipr
                        else
                            intensity = combineObsIntensity(Ipr, Isum, LP)
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
                            push!(hklList[origHKL].observations, ReflectionObservation(rot, fractionCalc, misymNum, intensity, sigIsum^2, [imageNumber]))
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
                            numOfExistingObs = length(hklList[origHKL].observations)
                            sameObservation = false
                            for obsNum = 1:numOfExistingObs #Loop through observations
                                refObservation = hklList[origHKL].observations[obsNum]
                                imagesOfObs = hklList[origHKL].observations[obsNum].imageNums
                                for image in imagesOfObs #Loop through images
                                    if misymNum == refObservation.misym #Check that the M/ISYM number is the same
                                        if imageNumber == image + 1 || imageNumber == image - 1 #Check is the image is consecutive
                                            #If it is consecutive then update the corresponding observation information.
                                            sameObservation = true
                                            refObservation.rotCentroid += rot
                                            refObservation.fractionCalc += fractionCalc
                                            refObservation.intensity += intensity
                                            refObservation.sigI += sigIsum^2
                                            push!(refObservation.imageNums, imageNumber)
                                            hklList[origHKL].observations[obsNum] = refObservation
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
                                    push!(hklList[origHKL].observations, ReflectionObservation(rot, fractionCalc, misymNum, intensity, sigIsum^2, [imageNumber]))
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

    return spacegroup, unitcell, hklList, imageArray
end
