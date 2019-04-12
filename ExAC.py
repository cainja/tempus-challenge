import requests

def callBulkMethod(VCFDataFrame):
    """ Function to interact with exacs API to get variant array of information from
    database

    Parameters:
    VCFDataFrame: Starting dataframe with variants.
                  Needs columns 'CHROM', 'POS', 'REF', 'ALT_Priority'
    url (str): This is the API that we want to interact with (only 2 options)
            Either /rest/bulk/variant or /rest/bulk/variant/variant

    Return:
    json object from API call
    """

    url= 'http://exac.hms.harvard.edu/rest/bulk/variant/variant'

    #need list of alternates to look up in API in the form 'CHROM-POS-REF-VAR'
    variants = ['-'.join([v['CHROM'], str(v['POS']), v['REF'], v['VAR']]) \
                for idx , v in VCFDataFrame.iterrows()]
    data = requests.post(url, json = variants)

    try:
        return data.json()
    except JSONDecodeError as e:
        raise ValueError('Variant list was not compatible with exacs API.')

def getFreq(VCFDataFrame, json, title = 'FREQ_ExAC',inplace = True):
    """ Function to parse json from ExAC API to get allelle frequency
    Returns updated DataFrame with new column of frequencies

    Parameters:
    VCFDataFrame: Starting dataframe. Needs columns 'CHROM', 'POS', 'REF', 'ALT_Priority'
    title (str): name for dataframe columns
    json (json): result from ExAC API call on /rest/bulk/variant/variant
    inplace (bool): whether or not to do the transformation in place

    Returns:
    Dataframe: pandas dataframe with new column containing frequency of mutation
    """
    freqs = []
    for key, item in json.items():
        try:
            freqs.append(item['allele_freq'])
        except KeyError as e:
            freqs.append('.') #General empty value for VCF files incase we want to annotate

    assert (len(freqs)==len(VCFDataFrame.index)), 'JSON size does not correspond to Dataframe'

    if inplace:
        VCFDataFrame[title] = freqs
    else:
        copy = VCFDataFrame.copy()
        copy[title]=freqs
        return copy

def getConsequence(VCFDataFrame, json, title = 'TYPE_vep', inplace = True):
    """ Function to parse json from ExAC API to get major consequence from vep annotation
    Returns updated DataFrame with new column of frequencies

    Parameters:
    VCFDataFrame: Starting dataframe. Needs columns 'CHROM', 'POS', 'REF', 'ALT_Priority'
    json (json): result from ExAC API call on /rest/bulk/variant/variant
    title (str): name for dataframe columns
    inplace (bool): whether or not to do the transformation in place

    Returns:
    Dataframe: pandas dataframe with new column containing frequency of mutation
    """
    consequences = []
    for key, item in json.items():
        try:
            consequences.append(item['vep_annotations'][0]['major_consequence'])
        except (KeyError, IndexError) as e:
            consequences.append('.') #General empty value for VCF files incase we want to annotate

    assert (len(consequences)==len(VCFDataFrame.index)), 'JSON size does not correspond to Dataframe'

    if inplace:
        VCFDataFrame[title] = consequences
    else:
        copy = VCFDataFrame.copy()
        copy[title]=consequences
        return copy
