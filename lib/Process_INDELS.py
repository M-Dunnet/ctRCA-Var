## ctRCA module imports
from lib import INDEL_Classifier as indel
from lib import ReferenceSet

def process_dels(params, file):
    '''
    Function to process deletion data from BAM files.
    '''
    ## Construct refset for the test sample
    print(f"\tLooking for deletion variants")
    del_refset = ReferenceSet.Del_ReferenceSet(
        file_path=file,
        bed_path=params['paths']['bedfile'],
        genome_path=params['paths']['genome'],
        collapse_umi=True
    ).reference_set

    del_scores = indel.INDEL_Analyzer(
        test_file=del_refset,
        target_indels=params['paths']['del_targets'],
        control_data=params['paths']['del_controls']
    )
    
    return del_scores.variant_calls

def process_ins(params, file):
    '''
    Function to process deletion data from BAM files.
    '''
    ## Construct refset for the test sample
    print(f"\tLooking for insertion variants")
    ins_refset = ReferenceSet.Ins_ReferenceSet(
        file_path=file,
        bed_path=params['paths']['bedfile'],
        genome_path=params['paths']['genome'],
        collapse_umi=True
    ).reference_set

    ins_scores = indel.INDEL_Analyzer(
        test_file=ins_refset,
        target_indels=params['paths']['ins_targets'],
        control_data=params['paths']['ins_controls']
    )

    return ins_scores.variant_calls
