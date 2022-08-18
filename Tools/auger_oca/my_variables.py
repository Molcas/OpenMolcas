import init as it

# my_variables.py

# Collecting data from init()

def my_variables(input_file):
    file_based,folder_based,fink_projection,conjdys,RAES,NAES,NAES_T,NAES_S,print_direct_dys =it.init1()

    OCA_atom,OCA_center,OCA_line,benergy,totalSymmetry,symmetry,nbasf,nash,nmo,cmotab,tdmtab,ncmo,nbasft,nasht,nosht,\
    comtaboff,comtbasoff,comtbasoff2,comtcmoff,comtnashoff,cmob,cmoa,tdmab,dyson=it.init2(input_file)

    hd5_file,basis_id_hd5,n_elements_desym,n_elements,element_desym,element= it.init3(nbasft,symmetry)

    return file_based,folder_based,fink_projection,conjdys,RAES,NAES,NAES_T,NAES_S,print_direct_dys,\
    OCA_atom,OCA_center,OCA_line,benergy,totalSymmetry,symmetry,nbasf,nash,nmo,cmotab,tdmtab,ncmo,nbasft,\
    nasht,nosht,comtaboff,comtbasoff,comtbasoff2,comtcmoff,comtnashoff,cmob,cmoa,tdmab,dyson,\
    hd5_file,basis_id_hd5,n_elements_desym,n_elements,element_desym,element
