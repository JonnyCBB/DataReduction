################################################################################
#Section: Inputs - Extract sequence information
#-------------------------------------------------------------------------------

immutable Residue
    identifier::ASCIIString
    carbons::Int64
    hydrogens::Int64
    nitrogens::Int64
    oxygens::Int64
    sulphurs::Int64
    phosphoruses::Int64
    seleniums::Int64
    molecularWeight::Float64
end

type ResidueType
    resType::ASCIIString
end
################################################################################
#Section: Residue composition
#-------------------------------------------------------------------------------
#In this section we define the composition of each residue type, whether it
#is protein, DNA or RNA. We also define an identifier for the residue.
#-------------------------------------------------------------------------------
const ALANINE_ID = "ALA" #Alanine identifier
const ALANINE_MW = 89.09 #Alanine molecular weight (g/mol)
#Alanine no. of atoms per residue.
const ALANINE_C = 3
const ALANINE_H = 7
const ALANINE_N = 1
const ALANINE_O = 2
const ALANINE_S = 0
const ALANINE_P = 0
const ALANINE_SE = 0
#-------------------------------------------------------------------------------
const ARGININE_ID = "ARG"# Arginine identifier.
const ARGININE_MW = 174.20 # Arginine molecular weight (g/mol).
# Arginine no. of atoms per residue.
const ARGININE_C = 6
const ARGININE_H = 15
const ARGININE_N = 4
const ARGININE_O = 1
const ARGININE_S = 0
const ARGININE_P = 0
const ARGININE_SE = 0
#-------------------------------------------------------------------------------
const ASPARAGINE_ID = "ASN" # Asparagine identifier.
const ASPARAGINE_MW = 132.12 # Asparagine molecular weight (g/mol).
# Asparagine no. of atoms per residue.
const ASPARAGINE_C = 4
const ASPARAGINE_H = 8
const ASPARAGINE_N = 2
const ASPARAGINE_O = 2
const ASPARAGINE_S = 0
const ASPARAGINE_P = 0
const ASPARAGINE_SE = 0
#-------------------------------------------------------------------------------
const ASPARTATE_ID = "ASP" # Aspartate identifier.
const ASPARTATE_MW = 133.10 # Aspartate molecular weight (g/mol).
# Aspartate no. of atoms per residue.
const ASPARTATE_C = 4
const ASPARTATE_H = 5
const ASPARTATE_N = 1
const ASPARTATE_O = 3
const ASPARTATE_S = 0
const ASPARTATE_P = 0
const ASPARTATE_SE = 0
#-------------------------------------------------------------------------------
const CYSTEINE_ID = "CYS" # Cysteine identifier.
const CYSTEINE_MW = 121.16 # Cysteine molecular weight (g/mol).
# Cysteine no. of atoms per residue.
const CYSTEINE_C = 3
const CYSTEINE_H = 5
const CYSTEINE_N = 1
const CYSTEINE_O = 1
const CYSTEINE_S = 1
const CYSTEINE_P = 0
const CYSTEINE_SE = 0
#-------------------------------------------------------------------------------
const GLUTAMINE_ID = "GLN" # Glutamine identifier.
const GLUTAMINE_MW = 146.15 # Glutamine molecular weight (g/mol).
# Glutamine no. of atoms per residue.
const GLUTAMINE_C = 5
const GLUTAMINE_H = 8
const GLUTAMINE_N = 2
const GLUTAMINE_O = 2
const GLUTAMINE_S = 0
const GLUTAMINE_P = 0
const GLUTAMINE_SE = 0
#-------------------------------------------------------------------------------
const GLUTAMATE_ID = "GLU" # Glutamate identifier.
const GLUTAMATE_MW = 147.13 # Glutamate molecular weight (g/mol).
# Glutamate no. of atoms per residue.
const GLUTAMATE_C = 5
const GLUTAMATE_H = 7
const GLUTAMATE_N = 1
const GLUTAMATE_O = 3
const GLUTAMATE_S = 0
const GLUTAMATE_P = 0
const GLUTAMATE_SE = 0
#-------------------------------------------------------------------------------
const GLYCINE_ID = "GLY" # Glycine identifier.
const GLYCINE_MW = 75.07 # Glycine molecular weight (g/mol).
# Glycine no. of atoms per residue.
const GLYCINE_C = 2
const GLYCINE_H = 3
const GLYCINE_N = 1
const GLYCINE_O = 1
const GLYCINE_S = 0
const GLYCINE_P = 0
const GLYCINE_SE = 0
#-------------------------------------------------------------------------------
const HISTIDINE_ID = "HIS" # Histidine identifier.
const HISTIDINE_MW = 155.16 # Histidine molecular weight (g/mol).
# Histidine no. of atoms per residue.
const HISTIDINE_C = 6
const HISTIDINE_H = 7
const HISTIDINE_N = 3
const HISTIDINE_O = 1
const HISTIDINE_S = 0
const HISTIDINE_P = 0
const HISTIDINE_SE = 0
#-------------------------------------------------------------------------------
const ISOLEUCINE_ID = "ILE" # Isoleucine identifier
const ISOLEUCINE_MW = 131.17 # Isoleucine molecular weight (g/mol).
# Isoleucine no. of atoms per residue.
const ISOLEUCINE_C = 6
const ISOLEUCINE_H = 11
const ISOLEUCINE_N = 1
const ISOLEUCINE_O = 1
const ISOLEUCINE_S = 0
const ISOLEUCINE_P = 0
const ISOLEUCINE_SE = 0
#-------------------------------------------------------------------------------
const LEUCINE_ID = "LEU" # Leucine identifier.
const LEUCINE_MW = 131.17 # Leucine molecular weight (g/mol).
# Leucine no. of atoms per residue.
const LEUCINE_C = 6
const LEUCINE_H = 11
const LEUCINE_N = 1
const LEUCINE_O = 1
const LEUCINE_S = 0
const LEUCINE_P = 0
const LEUCINE_SE = 0
#-------------------------------------------------------------------------------
const LYSINE_ID = "LYS" # Lysine identifier.
const LYSINE_MW = 146.19 # Lysine molecular weight (g/mol).
# Lysine no. of atoms per residue.
const LYSINE_C = 6
const LYSINE_H = 12
const LYSINE_N = 2
const LYSINE_O = 1
const LYSINE_S = 0
const LYSINE_P = 0
const LYSINE_SE = 0
#-------------------------------------------------------------------------------
const METHIONINE_ID = "MET" # Methionine identifier.
const METHIONINE_MW = 149.21 # Methionine molecular weight (g/mol).
# Methionine no. of atoms per residue.
const METHIONINE_C = 5
const METHIONINE_H = 9
const METHIONINE_N = 1
const METHIONINE_O = 1
const METHIONINE_S = 1
const METHIONINE_P = 0
const METHIONINE_SE = 0
#-------------------------------------------------------------------------------
const SELENOMETHIONINE_ID = "MSE" # Selenomethionine identifier.
const SELENOMETHIONINE_MW = 196.11 # Selenomethionine molecular weight (g/mol).
# Methionine no. of atoms per residue.
const SELENOMETHIONINE_C  = 5
const SELENOMETHIONINE_H = 9
const SELENOMETHIONINE_N = 1
const SELENOMETHIONINE_O = 1
const SELENOMETHIONINE_S = 0
const SELENOMETHIONINE_P = 0
const SELENOMETHIONINE_SE = 1
#-------------------------------------------------------------------------------
const PHENYLALANINE_ID = "PHE" # Phenylalanine identifier.
const PHENYLALANINE_MW = 165.19 # Phenylalanine molecular weight (g/mol).
# Phenylalanine no. of atoms per residue.
const PHENYLALANINE_C = 9
const PHENYLALANINE_H = 9
const PHENYLALANINE_N = 1
const PHENYLALANINE_O = 1
const PHENYLALANINE_S = 0
const PHENYLALANINE_P = 0
const PHENYLALANINE_SE = 0
#-------------------------------------------------------------------------------
const PROLINE_ID = "PRO" # Proline identifier.
const PROLINE_MW = 115.13 # Proline molecular weight (g/mol).
# Proline no. of atoms per residue.
const PROLINE_C = 5
const PROLINE_H = 7
const PROLINE_N = 1
const PROLINE_O = 1
const PROLINE_S = 0
const PROLINE_P = 0
const PROLINE_SE = 0
#-------------------------------------------------------------------------------
const SERINE_ID = "SER" # Serine identifier.
const SERINE_MW = 105.09 # Serine molecular weight (g/mol).
# Serine no. of atoms per residue.
const SERINE_C = 3
const SERINE_H = 5
const SERINE_N = 1
const SERINE_O = 2
const SERINE_S = 0
const SERINE_P = 0
const SERINE_SE = 0
#-------------------------------------------------------------------------------
const THREONINE_ID = "THR" # Threonine identifier.
const THREONINE_MW = 119.12 # Threonine molecular weight (g/mol).
  # Threonine no. of atoms per residue.
const THREONINE_C = 4
const THREONINE_H = 7
const THREONINE_N = 1
const THREONINE_O = 2
const THREONINE_S = 0
const THREONINE_P = 0
const THREONINE_SE = 0
#-------------------------------------------------------------------------------
const TRYPTOPHAN_ID = "TRP" # Tryptophan identifier.
const TRYPTOPHAN_MW = 204.23 # Tryptophan molecular weight (g/mol).
# Tryptophan no. of atoms per residue.
const TRYPTOPHAN_C = 11
const TRYPTOPHAN_H = 10
const TRYPTOPHAN_N = 2
const TRYPTOPHAN_O = 1
const TRYPTOPHAN_S = 0
const TRYPTOPHAN_P = 0
const TRYPTOPHAN_SE = 0
#-------------------------------------------------------------------------------
const TYROSINE_ID = "TYR" # Tyrosine identifier.
const TYROSINE_MW = 181.19 # Tyrosine molecular weight (g/mol).
# Tyrosine no. of atoms per residue.
const TYROSINE_C = 9
const TYROSINE_H = 9
const TYROSINE_N = 1
const TYROSINE_O = 2
const TYROSINE_S = 0
const TYROSINE_P = 0
const TYROSINE_SE = 0
#------------------------------------------------------------------------------#
const VALINE_ID = "VAL" # Valine identifier.
const VALINE_MW = 117.15 # Valine molecular weight (g/mol).
# Valine no. of atoms per residue.
const VALINE_C = 5
const VALINE_H = 9
const VALINE_N = 1
const VALINE_O = 1
const VALINE_S = 0
const VALINE_P = 0
const VALINE_SE = 0
#-------------------------------------------------------------------------------
const AMP_ID = "  A" # AMP identifier.
const AMP_MW = 347.2 # AMP molecular weight (g/mol).
# AMP no. of atoms per residue.
const AMP_C = 10
const AMP_H = 12
const AMP_N = 5
const AMP_O = 6
const AMP_S = 0
const AMP_P = 1
const AMP_SE = 0
#-------------------------------------------------------------------------------
const UMP_ID = "  U" # UMP identifier.
const UMP_MW = 324.2 # UMP molecular weight (g/mol).
# UMP no. of atoms per residue.
const UMP_C = 9
const UMP_H = 11
const UMP_N = 2
const UMP_O = 8
const UMP_S = 0
const UMP_P = 1
const UMP_SE = 0
#-------------------------------------------------------------------------------
const GMP_ID = "  G" # GMP identifier.
const GMP_MW = 363.2 # GMP molecular weight (g/mol).
# GMP no. of atoms per residue.
const GMP_C = 10
const GMP_H = 12
const GMP_N = 5
const GMP_O = 7
const GMP_S = 0
const GMP_P = 1
const GMP_SE = 0
#-------------------------------------------------------------------------------
const CMP_ID = "  C" # CMP identifier.
const CMP_MW = 323.2 # CMP molecular weight (g/mol).
# CMP no. of atoms per residue.
const CMP_C = 9
const CMP_H = 12
const CMP_N = 3
const CMP_O = 7
const CMP_S = 0
const CMP_P = 1
const CMP_SE = 0
#-------------------------------------------------------------------------------
const DAMP_ID = " DA" # dAMP identifier.
const DAMP_MW = 331.2 # dAMP molecular weight (g/mol).
  # dAMP no. of atoms per residue.
const DAMP_C = 10
const DAMP_H = 12
const DAMP_N = 5
const DAMP_O = 5
const DAMP_S = 0
const DAMP_P = 1
const DAMP_SE = 0
#-------------------------------------------------------------------------------
const DTMP_ID = " DT" # dTMP identifier.
const DTMP_MW = 322.2 # dTMP molecular weight (g/mol).
# dTMP no. of atoms per residue.
const DTMP_C = 10
const DTMP_H = 11
const DTMP_N = 2
const DTMP_O = 7
const DTMP_S = 0
const DTMP_P = 1
const DTMP_SE = 0
#-------------------------------------------------------------------------------
const DGMP_ID = " DG" # dGMP identifier.
const DGMP_MW = 347.2 # dGMP molecular weight (g/mol).
# dGMP no. of atoms per residue.
const DGMP_C = 10
const DGMP_H = 12
const DGMP_N = 5
const DGMP_O = 6
const DGMP_S = 0
const DGMP_P = 1
const DGMP_SE = 0
#-------------------------------------------------------------------------------
const DCMP_ID = " DC" # dCMP identifier.
const DCMP_MW = 307.2 # dCMP molecular weight (g/mol).
# dCMP no. of atoms per residue.
const DCMP_C = 9
const DCMP_H = 12
const DCMP_N = 3
const DCMP_O = 6
const DCMP_S = 0
const DCMP_P = 1
const DCMP_SE = 0
#-------------------------------------------------------------------------------
#End Section: Residue composition
################################################################################


################################################################################
# Section: Create residue dictionary
#-------------------------------------------------------------------------------
#In this section we create a dictionary containing the residue information. The
#keys for the dictionary will be the three letter ID for the residue.
function createResidueDict()
    residueDict = Dict{ASCIIString, Residue}() #Create the empty dictionary.

    residueDict[ALANINE_ID] = Residue(ALANINE_ID, ALANINE_C, ALANINE_H, ALANINE_N,
    ALANINE_O, ALANINE_S, ALANINE_P, ALANINE_SE, ALANINE_MW)

    residueDict[ARGININE_ID] = Residue(ARGININE_ID, ARGININE_C, ARGININE_H, ARGININE_N,
    ARGININE_O, ARGININE_S, ARGININE_P, ARGININE_SE, ARGININE_MW)

    residueDict[ASPARAGINE_ID] = Residue(ASPARAGINE_ID, ASPARAGINE_C, ASPARAGINE_H, ASPARAGINE_N,
    ASPARAGINE_O, ASPARAGINE_S, ASPARAGINE_P, ASPARAGINE_SE, ASPARAGINE_MW)

    residueDict[ASPARTATE_ID] = Residue(ASPARTATE_ID, ASPARTATE_C, ASPARTATE_H, ASPARTATE_N,
    ASPARTATE_O, ASPARTATE_S, ASPARTATE_P, ASPARTATE_SE, ASPARTATE_MW)

    residueDict[CYSTEINE_ID] = Residue(CYSTEINE_ID, CYSTEINE_C, CYSTEINE_H, CYSTEINE_N,
    CYSTEINE_O, CYSTEINE_S, CYSTEINE_P, CYSTEINE_SE, CYSTEINE_MW)

    residueDict[GLUTAMINE_ID] = Residue(GLUTAMINE_ID, GLUTAMINE_C, GLUTAMINE_H, GLUTAMINE_N,
    GLUTAMINE_O, GLUTAMINE_S, GLUTAMINE_P, GLUTAMINE_SE, GLUTAMINE_MW)

    residueDict[GLUTAMATE_ID] = Residue(GLUTAMATE_ID, GLUTAMATE_C, GLUTAMATE_H, GLUTAMATE_N,
    GLUTAMATE_O, GLUTAMATE_S, GLUTAMATE_P, GLUTAMATE_SE, GLUTAMATE_MW)

    residueDict[GLYCINE_ID] = Residue(GLYCINE_ID, GLYCINE_C, GLYCINE_H, GLYCINE_N,
    GLYCINE_O, GLYCINE_S, GLYCINE_P, GLYCINE_SE, GLYCINE_MW)

    residueDict[HISTIDINE_ID] = Residue(HISTIDINE_ID, HISTIDINE_C, HISTIDINE_H, HISTIDINE_N,
    HISTIDINE_O, HISTIDINE_S, HISTIDINE_P, HISTIDINE_SE, HISTIDINE_MW)

    residueDict[ISOLEUCINE_ID] = Residue(ISOLEUCINE_ID, ISOLEUCINE_C, ISOLEUCINE_H, ISOLEUCINE_N,
    ISOLEUCINE_O, ISOLEUCINE_S, ISOLEUCINE_P, ISOLEUCINE_SE, ISOLEUCINE_MW)

    residueDict[LEUCINE_ID] = Residue(LEUCINE_ID, LEUCINE_C, LEUCINE_H, LEUCINE_N,
    LEUCINE_O, LEUCINE_S, LEUCINE_P, LEUCINE_SE, LEUCINE_MW)

    residueDict[LYSINE_ID] = Residue(LYSINE_ID, LYSINE_C, LYSINE_H, LYSINE_N,
    LYSINE_O, LYSINE_S, LYSINE_P, LYSINE_SE, LYSINE_MW)

    residueDict[METHIONINE_ID] = Residue(METHIONINE_ID, METHIONINE_C, METHIONINE_H, METHIONINE_N,
    METHIONINE_O, METHIONINE_S, METHIONINE_P, METHIONINE_SE, METHIONINE_MW)

    residueDict[SELENOMETHIONINE_ID] = Residue(SELENOMETHIONINE_ID, SELENOMETHIONINE_C, SELENOMETHIONINE_H, SELENOMETHIONINE_N,
    SELENOMETHIONINE_O, SELENOMETHIONINE_S, SELENOMETHIONINE_P, SELENOMETHIONINE_SE, SELENOMETHIONINE_MW)

    residueDict[PHENYLALANINE_ID] = Residue(PHENYLALANINE_ID, PHENYLALANINE_C, PHENYLALANINE_H, PHENYLALANINE_N,
    PHENYLALANINE_O, PHENYLALANINE_S, PHENYLALANINE_P, PHENYLALANINE_SE, PHENYLALANINE_MW)

    residueDict[PROLINE_ID] = Residue(PROLINE_ID, PROLINE_C, PROLINE_H, PROLINE_N,
    PROLINE_O, PROLINE_S, PROLINE_P, PROLINE_SE, PROLINE_MW)

    residueDict[SERINE_ID] = Residue(SERINE_ID, SERINE_C, SERINE_H, SERINE_N,
    SERINE_O, SERINE_S, SERINE_P, SERINE_SE, SERINE_MW)

    residueDict[THREONINE_ID] = Residue(THREONINE_ID, THREONINE_C, THREONINE_H, THREONINE_N,
    THREONINE_O, THREONINE_S, THREONINE_P, THREONINE_SE, THREONINE_MW)

    residueDict[TRYPTOPHAN_ID] = Residue(TRYPTOPHAN_ID, TRYPTOPHAN_C, TRYPTOPHAN_H, TRYPTOPHAN_N,
    TRYPTOPHAN_O, TRYPTOPHAN_S, TRYPTOPHAN_P, TRYPTOPHAN_SE, TRYPTOPHAN_MW)

    residueDict[TYROSINE_ID] = Residue(TYROSINE_ID, TYROSINE_C, TYROSINE_H, TYROSINE_N,
    TYROSINE_O, TYROSINE_S, TYROSINE_P, TYROSINE_SE, TYROSINE_MW)

    residueDict[VALINE_ID] = Residue(VALINE_ID, VALINE_C, VALINE_H, VALINE_N,
    VALINE_O, VALINE_S, VALINE_P, VALINE_SE, VALINE_MW)

    residueDict[AMP_ID] = Residue(AMP_ID, AMP_C, AMP_H, AMP_N,
    AMP_O, AMP_S, AMP_P, AMP_SE, AMP_MW)

    residueDict[UMP_ID] = Residue(UMP_ID, UMP_C, UMP_H, UMP_N,
    UMP_O, UMP_S, UMP_P, UMP_SE, UMP_MW)

    residueDict[GMP_ID] = Residue(GMP_ID, GMP_C, GMP_H, GMP_N,
    GMP_O, GMP_S, GMP_P, GMP_SE, GMP_MW)

    residueDict[CMP_ID] = Residue(CMP_ID, CMP_C, CMP_H, CMP_N,
    CMP_O, CMP_S, CMP_P, CMP_SE, CMP_MW)

    residueDict[DAMP_ID] = Residue(DAMP_ID, DAMP_C, DAMP_H, DAMP_N,
    DAMP_O, DAMP_S, DAMP_P, DAMP_SE, DAMP_MW)

    residueDict[DTMP_ID] = Residue(DTMP_ID, DTMP_C, DTMP_H, DTMP_N,
    DTMP_O, DTMP_S, DTMP_P, DTMP_SE, DTMP_MW)

    residueDict[DGMP_ID] = Residue(DGMP_ID, DGMP_C, DGMP_H, DGMP_N,
    DGMP_O, DGMP_S, DGMP_P, DGMP_SE, DGMP_MW)

    residueDict[DCMP_ID] = Residue(DCMP_ID, DCMP_C, DCMP_H, DCMP_N,
    DCMP_O, DCMP_S, DCMP_P, DCMP_SE, DCMP_MW)

    return residueDict
end
# End Section: Create residue dictionary
################################################################################

################################################################################
#Section: Convert single letter ID to three letter ID
#-------------------------------------------------------------------------------
function resIDLetterConversion(singleLetterResID::Char, residueType::ResidueType)
    id = uppercase(singleLetterResID)
    if residueType.resType == "AminoAcid"
        if id == 'A'
            return ALANINE_ID
        elseif id == 'B'
            return ASPARAGINE_ID
        elseif id == 'C'
            return CYSTEINE_ID
        elseif id == 'D'
            return ASPARTATE_ID
        elseif id == 'E'
            return GLUTAMATE_ID
        elseif id == 'F'
            return PHENYLALANINE_ID
        elseif id == 'G'
            return GLYCINE_ID
        elseif id == 'H'
            return HISTIDINE_ID
        elseif id == 'I'
            return ISOLEUCINE_ID
        elseif id == 'J'
            return LEUCINE_ID
        elseif id == 'K'
            return LYSINE_ID
        elseif id == 'L'
            return LEUCINE_ID
        elseif id == 'M'
            return METHIONINE_ID
        elseif id == 'N'
            return ASPARAGINE_ID
        elseif id == 'O'
            @printf("***********************WARNING***********************\n")
            @printf("Amino acid residue id: '%s' unrecognised!\n",id)
            @printf("If you are hoping this will be interpreted as Pyrrolysine please contact us - ")
            @printf("E-mail: elspeth.garman@bioch.ox.ac.uk\n\n")
            return nothing
        elseif id == 'P'
            return PROLINE_ID
        elseif id == 'Q'
            return GLUTAMINE_ID
        elseif id == 'R'
            return ARGININE_ID
        elseif id == 'S'
            return SERINE_ID
        elseif id == 'T'
            return THREONINE_ID
        elseif id == 'U'
            @printf("***********************WARNING***********************\n")
            @printf("Amino acid residue id: '%s' unrecognised!\n",id)
            @printf("If you are hoping this will be interpreted as Selenocysteine please contact us - ")
            @printf("E-mail: elspeth.garman@bioch.ox.ac.uk\n\n")
            return nothing
        elseif id == 'V'
            return VALINE_ID
        elseif id == 'W'
            return TRYPTOPHAN_ID
        elseif id == 'X'
            return ALANINE_ID
        elseif id == 'Y'
            return TYROSINE_ID
        elseif id == 'Z'
            return GLUTAMATE_ID
        elseif id == '-'
            @printf("***********************WARNING***********************\n")
            @printf("Amino acid residue id: '%s' Does not exist", id)
            return nothing
        elseif id == '*'
            @printf("***********************WARNING***********************\n")
            @printf("Amino acid residue id: '%s' Does not exist", id)
            return nothing
        end
    elseif residueType.resType == "DNA"
        if id == 'A'
            return DAMP_ID
        elseif id == 'B'
            return DTMP_ID
        elseif id == 'C'
            return DCMP_ID
        elseif id == 'D'
            return DGMP_ID
        elseif id == 'E'
            @printf("***********************WARNING***********************\n")
            @printf("Nucleic acid residue id: '%s' Does not exist", id)
            return nothing
        elseif id == 'F'
            @printf("***********************WARNING***********************\n")
            @printf("Nucleic acid residue id: '%s' Does not exist", id)
            return nothing
        elseif id == 'G'
            return DGMP_ID
        elseif id == 'H'
            return DTMP_ID
        elseif id == 'I'
            @printf("***********************WARNING***********************\n")
            @printf("Nucleic acid residue id: '%s' Does not exist", id)
            return nothing
        elseif id == 'J'
            @printf("***********************WARNING***********************\n")
            @printf("Nucleic acid residue id: '%s' Does not exist", id)
            return nothing
        elseif id == 'K'
            return DGMP_ID
        elseif id == 'L'
            @printf("***********************WARNING***********************\n")
            @printf("Nucleic acid residue id: '%s' Does not exist", id)
            return nothing
        elseif id == 'M'
            return DAMP_ID
        elseif id == 'N'
            return DCMP_ID
        elseif id == 'O'
            @printf("***********************WARNING***********************\n")
            @printf("Nucleic acid residue id: '%s' Does not exist", id)
            return nothing
        elseif id == 'P'
            @printf("***********************WARNING***********************\n")
            @printf("Nucleic acid residue id: '%s' Does not exist", id)
            return nothing
        elseif id == 'Q'
            @printf("***********************WARNING***********************\n")
            @printf("Nucleic acid residue id: '%s' Does not exist", id)
            return nothing
        elseif id == 'R'
            return DGMP_ID
        elseif id == 'S'
            return DCMP_ID
        elseif id == 'T'
            return DTMP_ID
        elseif id == 'U'
            @printf("***********************WARNING***********************\n")
            @printf("DNA residue id: %s Does not exist", id)
            return nothing
        elseif id == 'V'
            return DAMP_ID
        elseif id == 'W'
            return DTMP_ID
        elseif id == 'X'
            @printf("***********************WARNING***********************\n")
            @printf("DNA residue id: %s Does not exist", id)
            return nothing
        elseif id == 'Y'
            return DCMP_ID
        elseif id == 'Z'
            @printf("***********************WARNING***********************\n")
            @printf("DNA residue id: %s Does not exist", id)
            return nothing
        elseif id == '-'
            @printf("***********************WARNING***********************\n")
            @printf("DNA residue id: %s Does not exist", id)
            return nothing
        elseif id == '*'
            @printf("***********************WARNING***********************\n")
            @printf("DNA residue id: %s Does not exist", id)
            return nothing
        end
    elseif residueType.resType == "RNA"
        if id == 'A'
            return AMP_ID
        elseif id == 'B'
            return UMP_ID
        elseif id == 'C'
            return CMP_ID
        elseif id == 'D'
            return GMP_ID
        elseif id == 'E'
            @printf("***********************WARNING***********************\n")
            @printf("Nucleic acid residue id: '%s' unrecognised!\n",id)
            return nothing
        elseif id == 'F'
            @printf("***********************WARNING***********************\n")
            @printf("Nucleic acid residue id: '%s' unrecognised!\n",id)
            return nothing
        elseif id == 'G'
            return GMP_ID
        elseif id == 'H'
            return UMP_ID
        elseif id == 'I'
            @printf("***********************WARNING***********************\n")
            @printf("Nucleic acid residue id: '%s' unrecognised!\n",id)
            return nothing
        elseif id == 'J'
            @printf("***********************WARNING***********************\n")
            @printf("Nucleic acid residue id: '%s' unrecognised!\n",id)
            return nothing
        elseif id == 'K'
            return GMP_ID
        elseif id == 'L'
            @printf("***********************WARNING***********************\n")
            @printf("Nucleic acid residue id: '%s' unrecognised!\n",id)
            return nothing
        elseif id == 'M'
            return AMP_ID
        elseif id == 'N'
            return CMP_ID
        elseif id == 'O'
            @printf("***********************WARNING***********************\n")
            @printf("Nucleic acid residue id: '%s' unrecognised!\n",id)
            return nothing
        elseif id == 'P'
            @printf("***********************WARNING***********************\n")
            @printf("Nucleic acid residue id: '%s' unrecognised!\n",id)
            return nothing
        elseif id == 'Q'
            @printf("***********************WARNING***********************\n")
            @printf("Nucleic acid residue id: '%s' unrecognised!\n",id)
            return nothing
        elseif id == 'R'
            return GMP_ID
        elseif id == 'S'
            return CMP_ID
        elseif id == 'T'
            @printf("***********************WARNING***********************\n")
            @printf("RNA residue id: %s Does not exist", id)
            return nothing
        elseif id == 'U'
            return UMP_ID
        elseif id == 'V'
            return AMP_ID
        elseif id == 'W'
            return UMP_ID
        elseif id == 'X'
            @printf("***********************WARNING***********************\n")
            @printf("RNA residue id: %s Does not exist", id)
            return nothing
        elseif id == 'Y'
            return CMP_ID
        elseif id == 'Z'
            @printf("***********************WARNING***********************\n")
            @printf("RNA residue id: %s Does not exist", id)
            return nothing
        elseif id == '-'
            @printf("***********************WARNING***********************\n")
            @printf("RNA residue id: %s Does not exist", id)
            return nothing
        elseif id == '*'
            @printf("***********************WARNING***********************\n")
            @printf("RNA residue id: %s Does not exist", id)
            return nothing
        else
            @printf("***********************WARNING***********************\n")
            @printf("Nucleic acid residue id: '%s' unrecognised!\n",id)
            return nothing
        end
    end
end
#End Section: Convert single letter ID to three letter ID
################################################################################
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

function getAtomOccurence!(atomOccurenceDict::Dict{ASCIIString, Int64}, residueID::Char, residueType::ResidueType,
    residues::Dict{ASCIIString, Residue})
    threeLetterID = resIDLetterConversion(residueID, residueType)
    residue = residues[threeLetterID]
    atomOccurenceDict["H"] += residue.hydrogens
    atomOccurenceDict["C"] += residue.carbons
    atomOccurenceDict["N"] += residue.nitrogens
    atomOccurenceDict["O"] += residue.oxygens
    atomOccurenceDict["P"] += residue.phosphoruses
    atomOccurenceDict["S"] += residue.sulphurs
    atomOccurenceDict["H"] += residue.seleniums
end

function parseCommentLine!(residueType::ResidueType, line::ASCIIString)
    if contains(line, "ISDNA")
        residueType.resType = "DNA"
    elseif contains(line, "ISRNA")
        residueType.resType = "RNA"
    else
        residueType.resType = "AminoAcid"
    end
end

function parseSequenceLine!(atomOccurenceDict::Dict{ASCIIString, Int64}, residueType::ResidueType, line::ASCIIString,
    residues::Dict{ASCIIString, Residue})
    for character in line
        if character != '\n'
            getAtomOccurence!(atomOccurenceDict, character, residueType, residues)
        end
    end
end

function parseSeqFileLine!(atomOccurenceDict::Dict{ASCIIString, Int64}, residueType::ResidueType, line::ASCIIString,
    residues::Dict{ASCIIString, Residue})
    if !isempty(strip(line))
        firstCharacter = line[1]
        if firstCharacter == '>'
            parseCommentLine!(residueType, line)
        else
            parseSequenceLine!(atomOccurenceDict, residueType, line, residues)
        end
    end
end

function getAtomicCompositon(sequenceFileLocation::ASCIIString)
    atomOccurenceDict = Dict{ASCIIString, Int64}()
    atomOccurenceDict["H"] = 0
    atomOccurenceDict["C"] = 0
    atomOccurenceDict["N"] = 0
    atomOccurenceDict["O"] = 0
    atomOccurenceDict["S"] = 0
    atomOccurenceDict["P"] = 0
    atomOccurenceDict["SE"] = 0
    residueType = ResidueType("AminoAcid")
    residues = createResidueDict()
    if isfile(sequenceFileLocation)
        open(sequenceFileLocation) do seqFile
            for line in eachline(seqFile)
                parseSeqFileLine!(atomOccurenceDict, residueType, line, residues)
            end
        end
    else
        @printf("*************************WARNING**************************\nThe sequence file path given: '%s' can't be found.\nPlease check the given path.\nMeanwhile the program will proceed without the file...\n\n", sequenceFileLocation)
        error("Currently there is no implementation to work without a sequence file. I need to get to this.")
    end
    return atomOccurenceDict
end

function additionalElements!(atomDict::Dict{ASCIIString, Int64}, additionalElements::ASCIIString)
    composition = split(additionalElements)
    if iseven(length(composition))
        for inputs = 1:2:length(composition)
            element = uppercase(composition[inputs])
            elementNums = parse(Int, composition[inputs+1])
            if haskey(atomDict, element)
                atomDict[element] += elementNums
            else
                atomDict[element] = elementNums
            end
        end
    else
        error("There should be an even number of inputs for the 'additionalElements' line.")
    end
end
#End Section: Inputs - Extract sequence information
################################################################################
