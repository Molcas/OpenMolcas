import numpy as np
import os

# oca_integrals.py

def oca_integrals(OCA_atom):
    # The file OCA.dat has 21 lines of integrals for each element. 
    # Available elements (currently) are: C, O, N, Ne
    # Reads 5 columns: I,J,L,M,G
    #  I,J: 1   2    3    4
    #       2s  2pz  2px  2py
    oca=list()
    nele_oca=21
    oca_dir = str(os.path.dirname(os.path.abspath(__file__)))+'/OCA.dat'  # give the full real path of OCA.dat
    with open(oca_dir,"r") as oc :
        for curline in oc:
            if curline.startswith("#"):
                pass
            else:
                line3=[float(elem) for elem in curline.split()]
                oca.append(line3)
    temp_oca = [x for x in oca if x != []] # 
    oca=np.array(temp_oca)
    #Carbon set
    oca_C=oca[:nele_oca]
    #Oxygen set
    oca_O=oca[nele_oca:nele_oca*2]
    #N set
    oca_N=oca[nele_oca*2:nele_oca*3]
    #Ne set
    oca_NE=oca[nele_oca*3:nele_oca*4]
    #Cl set
    oca_CL=oca[nele_oca*4:nele_oca*5]

    # Define One Center Integral according to the OCA_atom variable
    if OCA_atom=='C':
        OCI=oca_C
    elif OCA_atom=='O':
        OCI=oca_O
    elif OCA_atom=='N':
        OCI=oca_N
    elif OCA_atom=='NE':
        OCI=oca_NE
    elif OCA_atom=='CL':
        OCI=oca_CL

    return OCI

def elmij(OCA_atom,OCA_c,c,i,j,l,m):
    OCI = oca_integrals(OCA_atom)
    ver = 0.0
    if c ==OCA_c:
        for ez in OCI:
            if i == '2s':
                ii=1
            elif i == '2pz':
                ii=2
            elif i == '2px':
                ii=3
            elif i == '2py':
                ii=4
            else:
                break
            if j == '2s':
                jj=1
            elif j == '2pz':
                jj=2
            elif j == '2px':
                jj=3
            elif j == '2py':
                jj=4
            else:
                break

            if all([ii,jj,l,m] == ez[:4]) :
                ver = ez[4]
    return ver
# print('Eml',elmij('C 1s', 'C 1s','2py','2px',2,-2)) # this test should give: Eml 0.006919730033
# ------------------------------------------
