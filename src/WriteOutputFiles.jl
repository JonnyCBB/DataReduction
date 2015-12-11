function writeOutputFiles(hklList::Dict{Vector{Int16},Reflection}, spacegroup::SpaceGroup,
                  unitcell::Unitcell, HKL_FILENAME::ASCIIString, F2MTZ_INPUT_FILENAME::ASCIIString,
                  PROJECT_NAME::ASCIIString, CRYSTAL_NAME::ASCIIString, DATASET_NAME::ASCIIString,
                  SORTMTZ_INPUT_FILENAME::ASCIIString, MTZOUT_FILENAME::ASCIIString)
    outputFilename = runF2MTZ(hklList, spacegroup, unitcell, HKL_FILENAME, F2MTZ_INPUT_FILENAME,
                      PROJECT_NAME, CRYSTAL_NAME, DATASET_NAME)
    runSORTMTZ(outputFilename, SORTMTZ_INPUT_FILENAME, MTZOUT_FILENAME)
end

function runF2MTZ(hklList::Dict{Vector{Int16},Reflection}, spacegroup::SpaceGroup,
                  unitcell::Unitcell, HKL_FILENAME::ASCIIString,
                  F2MTZ_INPUT_FILENAME::ASCIIString, PROJECT_NAME::ASCIIString,
                  CRYSTAL_NAME::ASCIIString, DATASET_NAME::ASCIIString)
    ############################################################################
    #Section: Write HKL file with reflection information
    #---------------------------------------------------------------------------
    #THIS STILL NEEDS TO BE CORRECTED FOR THE SITUATION WHEN THE PROGRAM IS
    #EXTENDED TO HANDLE ANOMALOUS DATA.
    hklFile = open(HKL_FILENAME, "w")
    write(hklFile, @sprintf("Number of Reflections = %d\n", length(hklList)))
    if ANOMALOUS
        write(hklFile, @sprintf("Anomalous data = %s\n","FALSe"))
    else
        write(hklFile, @sprintf("Anomalous data = %s\n","TRUE"))
    end
    write(hklFile, @sprintf("Column Headers: H K L F SigF DANO SIGDANO F(+) SIGF(+) F(-) SIGF(-)\n"))
    for hkl in keys(hklList)
        reflection = hklList[hkl]
        write(hklFile, @sprintf("%5d", hkl[1]))
        write(hklFile, @sprintf("%5d", hkl[2]))
        write(hklFile, @sprintf("%5d", hkl[3]))
        write(hklFile, @sprintf("%10.3f", reflection.amplitude))
        write(hklFile, @sprintf("%10.3f", reflection.amplitudeSig))
        write(hklFile, @sprintf("%10.3f", 0.0))
        write(hklFile, @sprintf("%10.3f", 0.0))
        write(hklFile, @sprintf("%10.3f", reflection.amplitude))
        write(hklFile, @sprintf("%10.3f", reflection.amplitudeSig))
        write(hklFile, @sprintf("%10.3f", reflection.amplitude))
        write(hklFile, @sprintf("%10.3f\n", reflection.amplitudeSig))
    end
    close(hklFile)
    #End Section: Write HKL file with reflection information
    #---------------------------------------------------------------------------

    ############################################################################
    #Section: Write f2mtz input file
    #---------------------------------------------------------------------------
    f2mtzInputFile = open(F2MTZ_INPUT_FILENAME, "w")
    write(f2mtzInputFile, @sprintf("SKIP 3\n"))
    write(f2mtzInputFile, @sprintf("TITLE HKL to MTZ\n"))
    write(f2mtzInputFile, @sprintf("NAME PROJECT %s CRYSTAL %s DATASET %s\n", PROJECT_NAME, CRYSTAL_NAME, DATASET_NAME))
    write(f2mtzInputFile, @sprintf("CELL %.2f %.2f %.2f %.2f %.2f %.2f\n", unitcell.a, unitcell.b, unitcell.c, unitcell.α, unitcell.β, unitcell.γ))
    write(f2mtzInputFile, @sprintf("SYMMETRY %d\n", spacegroup.number))
    write(f2mtzInputFile, @sprintf("LABOUT H K L F_mean SIGF_mean DANO SIGDANO F(+) SIGF(+) F(-) SIGF(-)\n"))
    write(f2mtzInputFile, @sprintf("CTYPOUT H H H F Q D Q G L G L\n"))
    write(f2mtzInputFile, @sprintf("FORMAT '(3F5.0,8F10.3)'\n"))
    write(f2mtzInputFile, @sprintf("END\n"))
    close(f2mtzInputFile)
    #End Section: Write f2mtz input file
    ############################################################################

    ############################################################################
    #Section: Run F2MTZ
    #---------------------------------------------------------------------------
    #import python function to run F2MTZ. (This is horrible hack because I
    #couldn't work out how to do it in Julia)
    @pyimport RunSystemCommand as runsys

    #Create temporary filename for the unsorted mtz file
    f2mtzOutFilename = @sprintf("%s_unsorted.mtz", split(MTZOUT_FILENAME,".")[1])

    #Run f2mtz
    runsys.run_system_command(@sprintf("f2mtz HKLIN %s HKLOUT %s < %s", HKL_FILENAME, f2mtzOutFilename, F2MTZ_INPUT_FILENAME))
    #End Section: Run F2MTZ
    ############################################################################

    #Delete f2mtz input files
    if isfile(HKL_FILENAME)
        rm(HKL_FILENAME)
    end
    if isfile(F2MTZ_INPUT_FILENAME)
        rm(F2MTZ_INPUT_FILENAME)
    end

    return f2mtzOutFilename #return the filename used for the output MTZ file.
end

function runSORTMTZ(f2mtzOutFilename::ASCIIString, SORTMTZ_INPUT_FILENAME::ASCIIString,
                    MTZOUT_FILENAME::ASCIIString)
    #Create input file for SORTMTZ
    sortmtzInputFile = open(SORTMTZ_INPUT_FILENAME, "w")
    write(sortmtzInputFile, "H K L\n")
    close(sortmtzInputFile)

    #import python function to run SORTMTZ. (This is horrible hack because I
    #couldn't work out how to do it in Julia)
    @pyimport RunSystemCommand as runsys

    #Run SORTMTZ
    runsys.run_system_command(@sprintf("sortmtz HKLIN %s HKLOUT %s < %s", f2mtzOutFilename, MTZOUT_FILENAME, SORTMTZ_INPUT_FILENAME))

    #Delete unused files
    if isfile(SORTMTZ_INPUT_FILENAME)
        rm(SORTMTZ_INPUT_FILENAME)
    end

    if isfile(f2mtzOutFilename)
        rm(f2mtzOutFilename)
    end
end
