from cyvcf2 import VCF, Writer

def annotate(filepath, VCFDataFrame):
    """Function to write specific calculated and API values into desired VCF
    This is a very explicit function.

    Parameters:
    filepath: File path to desired .vcf to annotate
    VCFDataFrame: Dataframe with values that we want to annotate

    Returns:
    updates .vcf file at filepath
    """

    vcf = VCF(filepath)

    #This is hardcoded as it is curated
    list_of_annotations = [
        {'ID': 'VAR',
        'Description': "Selected Variant based on prioritization of \
        (1) 'complex', (2) 'ins', (3) 'del', (4) 'mnp', (5) 'snp'",
        'Type':'String','Number': '1'},
        {'ID': 'VAR_TYPE',
        'Description':"Annotated variant type based on prioritization of \
        (1) 'complex', (2) 'ins', (3) 'del', (4) 'mnp', (5) 'snp'",
        'Type':'String','Number': '1'},
        {'ID': 'VAR_COUNT',
        'Description':"Count of times selected variant was observed",
        'Type':'Interger','Number': '1'},
        {'ID': 'VAR_FRAC',
        'Description':"Fraction of total reads that the variant was observed",
        'Type':'Float','Number': '1'},
        {'ID': 'FREQ_ExAC',
        'Description':"Allele frequency of prioritized alt according to ExAC",
        'Type':'Float','Number': '1'},
        {'ID': 'TYPE_vep',
        'Description':"vep annotation of major consequence of variant",
        'Type':'String','Number': '1'}
    ]

    for annotation in list_of_annotations:
        vcf.add_info_to_header(annotation)

    w = Writer('Annotated_{}'.format(filepath), vcf)

    for i, variant in enumerate(vcf):
        for annotation in list_of_annotations:
            variant.INFO[annotation['ID']] = str(VCFDataFrame.at[i, annotation['ID']])
        w.write_record(variant)

    w.close(); vcf.close()
