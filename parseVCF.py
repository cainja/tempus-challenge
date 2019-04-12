from cyvcf2 import VCF
import pandas as pd

def loadVCF(filepath):
    """ Get an iterator from cyvcf2 library

    Parameters:
    filepath (str or os.path object): location of VCF file

    Returns:
    generator: iterator for each variants from file

    """

    variant_iterator = VCF(filepath)
    return variant_iterator

def parseVCFtoDataFrame(filepath,
                        keys = ['CHROM', 'POS', 'REF' ,'ALT'],
                        info_keys = ['TYPE', 'DP']):
    """ Grab selected information from VCF to put into a more easily
        manipulated pandas dataframe

    Parameters:
    filepath (str or os.path object): location of VCF file
    keys (list): names of columns (str) from VCF file, defaults to ['CHROM', 'POS', 'REF' ,'ALT']
            (Options are CHROM, POS, ID, REF, ALT, QUAL, FILTER)
    info_keys (list): names of site-level annotations - for this project TYPE and DP are defaults.

    Returns:
    Dataframe: pandas dataframe with specified columns and site level annotations
    """

    # Build list of lists to put into a dataframe
    variant_iterator = loadVCF(filepath)
    row_container = []
    VCFDataFrame = pd.DataFrame(columns=keys+info_keys)
    for variant in variant_iterator:

        #keys and info_keys are in different nested locations
        row = [getattr(variant, key) for key in keys]
        row += [variant.INFO[info_key] for info_key in info_keys]
        row_container.append(row)


    VCFDataFrame = pd.DataFrame(row_container, columns=keys+info_keys)

    return VCFDataFrame

def prioritizeLabel(VCFDataFrame,
                    priority = ['complex', 'ins', 'del', 'mnp', 'snp'], inplace=True):
    """ Picks single label out of possible labels in TYPE according to priority list

    Parameters:
    VCFDataFrame (pd.DataFrame): table containing parsed VCF data
    priority (list): list of priorities (str) in descending order. Defaults to
        ['complex', 'ins', 'del', 'mnp', 'snp'].
    inplace (bool): option to update VCFDataframe

    Returns:
    DataFrame: updated dataframe based on prioritized variant types

    """
    selected_type = []
    selected_alt = []

    for altType, alt in zip(VCFDataFrame['TYPE'], VCFDataFrame['ALT']):
        altType = altType.split(',')
        assert (len(altType)==len(alt))

        #most cases where only 1 variant in locus
        if len(altType) == 1:
            selected_type.append(altType[0])
            selected_alt.append(alt[0])
            continue

        #pick one possible variant if multiple
        for type in priority:
            if type in altType:
                idx = altType.index(type)
                selected_alt.append(alt[idx])
                selected_type.append(altType[idx])
                break
        else: raise Exception('/'.join(altType) + ' were not found in priority.')

    #either edit VCFDataFrame or make a copy based on inplace
    if inplace:
        VCFDataFrame['VAR_TYPE']=selected_type
        VCFDataFrame['VAR']=selected_alt
    else:
        copy = VCFDataFrame.copy()
        copy['VAR_TYPE']=selected_type
        copy['VAR']=selected_alt
        return copy
