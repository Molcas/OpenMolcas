!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2022, Ignacio Fdez. Galvan                             *
!***********************************************************************

module RunFile_data

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

! lw       : Label width
! DS_cache : In-memory cache of real scalar values
! IS_cache : In-memory cache of integer scalar values

integer(kind=iwp), parameter :: IDRun = 34676777, lw = 16, nHdrSz = 128, nToc = 1024, nTocCA = 32, nTocDA = 256, nTocDS = 64, &
                                nTocIA = 128, nTocIS = 128, NulPtr = -1, VNRun = 4096
integer(kind=iwp), parameter :: icWr = 1, icRd = 2, &
                                rcOK = 0, rcNotFound = 1, rcWrongType = 2, &
                                sNotUsed = 0, sRegularField = 1, sSpecialField = 2, &
                                TypUnk = 0, TypInt = 1, TypDbl = 2, TypStr = 3, TypLgl = 4

type DS_cache_item
  real(kind=wp) :: val = Zero
  character(len=lw) :: lab = ''
end type DS_cache_item

type IS_cache_item
  integer(kind=iwp) :: val = 0
  character(len=lw) :: lab = ''
end type IS_cache_item

type RunHdr_type
  integer(kind=iwp) :: ID, Ver, Next, Items, DaLab, DaPtr, DaLen, DaMaxLen, DaTyp
end type RunHdr_type

type Toc_item
  character(len=lw) :: Lab = 'Empty'
  integer(kind=iwp) :: Ptr = NulPtr, Len = 0, MaxLen = 0, Typ = TypUnk
end type Toc_item

integer(kind=iwp) :: i_run_CA_used(nTocCA), i_run_DA_used(nTocDA), i_run_DS_used(nTocDS), i_run_IA_used(nTocIA), &
                     i_run_IS_used(nTocIS), num_DS_init, num_IS_init
character(len=8) :: RunName, RnNmStk(4)
type(DS_cache_item) :: DS_cache(nTocDS)
type(IS_cache_item) :: IS_cache(nTocIS)
type(RunHdr_type) :: RunHdr
type(Toc_item) :: Toc(nToc)

!***********************************************************************
!                            RUNFILE LABELS                            *
!***********************************************************************

!> List of known character array labels:
!>
!> - '``DFT functional``'     Name of the functional used for the KS-DFT calculation.
!> - '``Irreps``'             Names of the irreducible representations.
!> - '``Relax Method``'       Name of the method used for geometry optimizations.
!> - '``Seward Title``'       The title of the calculation as specified in module SEWARD.
!> - '``Slapaf Info 3``'      Misc. information for module SLAPAF.
!> - '``Unique Atom Names``'  List of the names of the symmetry unique atoms.
!> - '``Unique Basis Names``' List of the basis function names.
!> - '``MkNemo.lMole``'       The labels of molecules as specified in mknemo module.
!> - '``MkNemo.lCluster``'    The labels of clusters as specified in mknemo module.
!> - '``MkNemo.lEnergy``'     The labels of energies as specified in mknemo module.
!> - '``Frag_Type``'          EFP fragment labels
!> - '``ABC``'                EFP atom labels

character(len=lw), parameter :: LabelsCA(nTocCA) = [ &
                                'DFT functional  ','Irreps          ','Relax Method    ','Seward Title    ', & !  1- 4
                                'Slapaf Info 3   ','Unique Atom Name','Unique Basis Nam','LP_L            ', & !  5- 8
                                'MkNemo.lMole    ','MkNemo.lCluster ','MkNemo.lEnergy  ','Symbol ZMAT     ', & !  9-12
                                'Tinker Name     ','ESPF Filename   ','ChDisp          ','cmass           ', & ! 13-16
                                'BirthCertificate','LastEnergyMethod','MMO Labels      ','MCLR Root       ', & ! 17-20
                                'Frag_Type       ','ABC             ','Un_cen Names    ','cDmp            ', & ! 21-24
                                'dc: cDmp        ','SymmetryCInfo   ','SewardXTitle    ','Align_Weights   ', & ! 25-28
                                'Quad_c          ','                ','                ','                ']   ! 29-32

!> List of known real array labels:
!>
!> - '``Analytic Hessian``'         Analytic Hessian.
!> - '``Center of Charge``'         Nuclear center of charge.
!> - '``Center of Mass``'           Nuclear center of mass.
!> - '``CMO_ab``'
!> - '``D1ao``'                     One particle density matrix, AO basis.
!> - '``D1ao_ab``'
!> - '``D1aoVar``'                  Generalized one particle density matrix, AO basis.
!> - '``D1av``'                     Average one particle density matrix, AO basis.
!> - '``D1mo``'                     One particle density matrix, MO basis.
!> - '``D1sao``'                    One particle spin density matrix, AO basis.
!> - '``D1activeao``'               One particle density matrix, AO basis, active orbitals
!> - '``D2av``'                     Average two particle density matrix for the active space, AO basis.
!> - '``dExcdRa``'                  The potential of the exchange-correlation functional.
!> - '``DLAO``'
!> - '``DLMO``'
!> - '``Effective nuclear charge``' Effective nuclear charge for each unique atom.
!> - '``FockO_ab``'
!> - '``FockOcc``'                  Generalized Fock matrix, AO basis.
!> - '``GeoNew``'                   Next guess for Cartesian coordinates for the unique atoms.
!> - '``GeoNewPC``'                 Next guess for Cartesian coordinates for the unique point charges.
!> - '``GRAD``'                     Gradient with respect to nuclear displacements, for all unique atoms.
!> - '``Hess``'
!> - '``HF-forces``'                Hellmann--Feynman forces.
!> - '``Last orbitals``'            Last set of orbital computed.
!> - '``LCMO``'
!> - '``MEP-Coor``'                 List of nuclear coordinates along a minimum energy path, for unique atoms.
!> - '``MEP-Energies``'             List of energies along a minimum energy path.
!> - '``MEP-Grad``'                 List of nuclear gradients along a minimum energy path.
!> - '``MP2 restart``'              Information for restarting direct MP2 calculations.
!> - '``Mulliken Charge``'          Mulliken population charges for each unique atom.
!> - '``NEMO TPC``'
!> - '``Nuclear charge``'           Actual nuclear charge for each unique atom.
!> - '``OrbE``'                     SCF orbital energies.
!> - '``OrbE_ab``'
!> - '``P2MO``'
!> - '``PCM Charges``'              Charges for each tessera in the PCM model.
!> - '``PCM Info``'                 Misc. information needed for the PCM model.
!> - '``PLMO``'
!> - '``RASSCF orbitals``'          Last orbitals generated by the RASSCF module.
!> - '``Reaction field``'           Misc. information needed for the Kirkwood model.
!> - '``SCFInfoR``'                 Misc. information needed by the SCF module.
!> - '``SCF orbitals``'             Last orbitals generated by the SCF module.
!> - '``Slapaf Info 2``'            Misc. information needed by the SLAPAF module.
!> - '``Unique Coordinates``'       Cartesian coordinates for the symmetry unique atoms.
!> - '``Last energies``'            Energies for all roots in the last calculation.
!> - '``Dipole moment``'            The last computed dipole moment.
!> - '``MkNemo.vDisp``'             The displacements matrix as specified in the mknemo module.
!> - '``MkNemo.tqCluster``'         The transformation matrix for clusters as specified in the mknemo module.
!> - '``MkNemo.Energies``'          The energies of super-system and clusters as specified in the mknemo module.
!> - '``RASSCF OrbE``'              RASSCF orbital energies.
!> - '``GRD1``'                     MR-CISD gradient state1.
!> - '``GRD2``'                     MR-CISD gradient state2.
!> - '``NADC``'                     MR-CISD NADC vector state1/state2.
!> - '``MR-CISD energy``'           MR-CISD energies state1,state2.
!> - '``NOSEHOOVER``'               Extra-degrees of fredom needed to generated canonical ensemble.
!> - '``T-Matrix``'                 T-Matrix associated with PCO.
!> - '``rInt0``'                    Stored constrained values.
!> - '``Weights``'                  Weights used for alignment and hypersphere constraint.
!> - '``MEP-Lengths``'              Lengths of the MEP steps.
!> - '``MEP-Curvatures``'           Curvatures of the MEP steps.
!> - '``D1ao-``'                    Antisymmetric transition density matrix, in AO.
!> - '``P2MOT``'                    "Fake" two-body density needed for MC-PDFT gradient calculations.
!> - '``Guessorb energies``'
!> - '``Pseudo Coordinates``'
!> - '``Initial Coordinates``'
!> - '``Isotopes``'                 Atom masses in a.u.
!> - '``EPF_Coors``'                EFP fragment coordinates.
!> - '``Last Dipole Moments``'
!> - '``Un_cen Effective Charge``'
!> - '``Un_cen Coordinates``'

character(len=lw), parameter :: LabelsDA(nTocDA) = [ &
                                'Analytic Hessian','Center of Charge','Center of Mass  ','CMO_ab          ', & !   1-  4
                                'D1ao            ','D1ao_ab         ','D1aoVar         ','D1av            ', & !   5-  8
                                'D1mo            ','D1sao           ','D2av            ','dExcdRa         ', & !   9- 12
                                'DLAO            ','DLMO            ','Effective nuclea','FockO_ab        ', & !  13- 16
                                'FockOcc         ','GeoNew          ','GeoNewPC        ','GRAD            ', & !  17- 20
                                'Hess            ','HF-forces       ','Last orbitals   ','LCMO            ', & !  21- 24
                                'MEP-Coor        ','MEP-Energies    ','MEP-Grad        ','MP2 restart     ', & !  25- 28
                                'Mulliken Charge ','NEMO TPC        ','Nuclear charge  ','OrbE            ', & !  29- 32
                                'OrbE_ab         ','P2MO            ','PCM Charges     ','PCM Info        ', & !  33- 36
                                'PLMO            ','RASSCF orbitals ','Reaction field  ','SCFInfoR        ', & !  37- 40
                                'SCF orbitals    ','Slapaf Info 2   ','Unique Coordinat','Vxc_ref         ', & !  41- 44
                                'PotNuc00        ','h1_raw          ','h1    XX        ','HEFF            ', & !  45- 48
                                'PotNucXX        ','Quad_r          ','RCTFLD          ','RFrInfo         ', & !  49- 52
                                'DKH_Info        ','Real_Info       ','Last orbitals_ab','SCFInfoI_ab     ', & !  53- 56
                                'SCFInfoR_ab     ','Transverse      ','SM              ','LP_Coor         ', & !  57- 60
                                'LP_Q            ','DFT_TwoEl       ','Unit Cell Vector','SCF orbitals_ab ', & !  61- 64
                                'Guessorb        ','Guessorb energie','Last energies   ','LoProp Dens 0   ', & !  65- 68
                                'LoProp Dens 1   ','LoProp Dens 2   ','LoProp Dens 3   ','LoProp Dens 4   ', & !  69- 72
                                'LoProp Dens 5   ','LoProp Dens 6   ','LoProp Integrals','MpProp Orb Ener ', & !  73- 76
                                'LoProp H0       ','Dipole moment   ','RICD_Info       ','BMtrx           ', & !  77- 80
                                'CList           ','DList           ','MkNemo.vDisp    ','MkNemo.tqCluster', & !  81- 84
                                'MkNemo.Energies ','MMHessian       ','Bfn Coordinates ','Pseudo Coordinat', & !  85- 88
                                'Pseudo Charge   ','RASSCF OrbE     ','Ref_Geom        ','LoProp Charge   ', & !  89- 92
                                'Initial Coordina','Grad State1     ','Grad State2     ','NADC            ', & !  93- 96
                                'MR-CISD energy  ','Saddle          ','Reaction Vector ','IRC-Coor        ', & !  97-100
                                'IRC-Energies    ','IRC-Grad        ','MM Grad         ','Velocities      ', & ! 101-104
                                'FC-Matrix       ','umass           ','ESO_SINGLE      ','UMATR_SINGLE    ', & ! 105-108
                                'UMATI_SINGLE    ','ANGM_SINGLE     ','TanVec          ','Nuc Potential   ', & ! 109-112
                                'RF CASSCF Vector','Cholesky BkmThr ','NOSEHOOVER      ','T-Matrix        ', & ! 113-116
                                'rInt0           ','Weights         ','MEP-Lengths     ','MEP-Curvatures  ', & ! 117-120
                                'Hss_X           ','Hss_Q           ','KtB             ','BMxOld          ', & ! 121-124
                                'TROld           ','qInt            ','dqInt           ','Fragment_Fock   ', & ! 125-128
                                'RAmatrixV       ','IAmatrixV       ','AllCIP          ','AllCIPP         ', & ! 129-132
                                'VenergyP        ','K               ','MMO Coords      ','MMO Grad        ', & ! 133-136
                                'Hss_upd         ','TR              ','D1ao-           ','ESFS_SINGLE     ', & ! 117-140
                                'LA Fact         ','primitives      ','Isotopes        ','P2AO            ', & ! 141-144
                                'State Overlaps  ','EFP_Coors       ','DIP1_SINGLE     ','P2MOT           ', & ! 145-148
                                'ONTOPO          ','ONTOPT          ','OE_OT           ','TEG_OT          ', & ! 149-152
                                'FI_V            ','FA_V            ','FOCK_PDFT       ','AMFI_SINGLE     ', & ! 153-156
                                'HAMSOR_SINGLE   ','HAMSOI_SINGLE   ','Last Dipole Mome','Un_cen Effective', & ! 157-160
                                'Un_cen Coordinat','                ','                ','                ', & ! 161-164
                                '                ','                ','                ','Proj_Coord      ', & ! 165-168
                                'd1activeao      ','Keep_Coord      ','PCMSph          ','PCMTess         ', & ! 169-172
                                'Vert            ','Centr           ','SSph            ','PCMDM           ', & ! 173-176
                                'EF_Centers      ','OAM_Center      ','OMQ_Center      ','DMS_Centers     ', & ! 177-180
                                'Wel_Info        ','AMP_Center      ','RP_Centers      ','XF              ', & ! 181-184
                                'rDmp            ','rDmp:A          ','rDmp:S          ','D1saoVar        ', & ! 185-188
                                'ESFS_SINGLEAU   ','ESO_LOW         ','SFS_HAM         ','SFS_OVLP        ', & ! 189-192
                                'FocMS           ','MSPDFTD5        ','MSPDFTD6        ','TwoEIntegral    ', & ! 193-196
                                'D1MOt           ','D1INTER         ','P2INTER         ','D1AO_MS         ', & ! 297-200
                                'D1SAO_MS        ','MS_FINAL_ROT    ','F1MS            ','F2MS            ', & ! 201-204
                                'FxyMS           ','SH_Ovlp_Save    ','Old_Phase       ','<rhoB|VnucA>    ', & ! 205-208
                                '                ','                ','                ','                ', & ! 209-212
                                '                ','                ','                ','                ', & ! 213-216
                                '                ','                ','                ','                ', & ! 217-220
                                '                ','                ','                ','                ', & ! 221-224
                                '                ','                ','                ','                ', & ! 225-228
                                '                ','                ','                ','                ', & ! 229-232
                                '                ','                ','                ','                ', & ! 233-236
                                '                ','                ','                ','                ', & ! 237-240
                                '                ','                ','                ','                ', & ! 241-244
                                '                ','                ','                ','                ', & ! 245-248
                                '                ','                ','                ','                ', & ! 249-252
                                '                ','                ','                ','                ']   ! 253-256

!> List of known real scalar labels:
!>
!> - '``CASDFT energy``'             Energy for the last CASDFT calculation.
!> - '``CASPT2 energy``'             Energy for the last CASPT2 calculation.
!> - '``CASSCF energy``'             Energy for the last CASSCF calculation.
!> - '``Ener_ab``'
!> - '``KSDFT energy``'              Energy for the last KS-DFT calculation.
!> - '``Last energy``'               Last energy computed.
!> - '``PC Self Energy``'            Self energy for point charges.
!> - '``PotNuc``'                    Nuclear repusion energy.
!> - '``RF Self Energy``'            Self energy in the Kirkwood model.
!> - '``SCF energy``'                Energy for the last SCF calculation.
!> - '``EThr``'                      Energy convergence threshold.
!> - '``Thrs``'
!> - '``UHF energy``'
!> - '``DFT exch coeff``'            Scaling factor for exchange terms of a density functional.
!> - '``DFT corr coeff``'            Scaling factor for correlation terms of a density functional.
!> - '``Cholesky Threshold``'
!> - '``Total Nuclear Charge``'
!> - '``Numerical Gradient rDelta``'
!> - '``Total Charge``'              total number of electrons.

character(len=lw), parameter :: LabelsDS(nTocDS) = [ &
                                'CASDFT energy   ','CASPT2 energy   ','CASSCF energy   ','Ener_ab         ', & !  1- 4
                                'KSDFT energy    ','Last energy     ','PC Self Energy  ','PotNuc          ', & !  5- 8
                                'RF Self Energy  ','SCF energy      ','Thrs            ','UHF energy      ', & !  9-12
                                'E_0_NN          ','W_or_el         ','W_or_Inf        ','EThr            ', & ! 11-16
                                'Cholesky Thresho','Total Nuclear Ch','Numerical Gradie','MpProp Energy   ', & ! 17-20
                                'UHFSPIN         ','S delete thr    ','T delete thr    ','MD_Etot0        ', & ! 21-24
                                'MD_Time         ','LDF Accuracy    ','NAD dft energy  ','GradLim         ', & ! 25-28
                                'Average energy  ','Timestep        ','MD_Etot         ','Max error       ', & ! 29-32
                                'Total Charge    ','DFT exch coeff  ','DFT corr coeff  ','Value_l         ', & ! 33-36
                                'R_WF_HMC        ','                ','                ','                ', & ! 37-40
                                '                ','                ','                ','                ', & ! 41-44
                                '                ','                ','                ','                ', & ! 45-48
                                '                ','                ','                ','                ', & ! 49-52
                                '                ','                ','                ','                ', & ! 53-56
                                '                ','                ','                ','                ', & ! 57-60
                                '                ','                ','                ','                ']   ! 61-64

!> List of known integer array labels:
!>
!> - '``Center Index``'
!> - '``Ctr Index Prim``'       Idem with primitive basis set.
!> - '``nAsh``'                 The number of active orbitals per irreducible representation.
!> - '``nBas``'                 The number of basis functions per irreducible representation.
!> - '``nDel``'                 The number of deleted orbitals per irreducible representation.
!> - '``nFro``'                 The number of frozen orbitals per irreducible representation, i.e. orbitals that are not optimized.
!> - '``nIsh``'                 The number of inactive orbitals per irreducible representation.
!> - '``nIsh beta``'
!> - '``nOrb``'                 The total number of orbitals per irreducible representation.
!> - '``Orbital Type``'
!> - '``Slapaf Info 1``'        Misc. information for module SLAPAF.
!> - '``Symmetry operations``'  The symmetry operations of the point group.
!> - '``Non valence orbitals``' The total number of non valence orbitals per irreducible representation.
!> - '``MkNemo.hDisp``'         The hash matrix for displacements as specified in the mknemo module.
!> - '``NumCho``'               Number of Cholesky vectors.
!> - '``nFroPT``'               Number of Frozen for PT.
!> - '``nDelPT``'               Number of Deleted for PT.

character(len=lw), parameter :: LabelsIA(nTocIA) = [ &
                                'Center Index    ','nAsh            ','nBas            ','nDel            ', & !   1-  4
                                'nFro            ','nIsh            ','nIsh beta       ','nOrb            ', & !   5-  8
                                'Orbital Type    ','Slapaf Info 1   ','Symmetry operati','nIsh_ab         ', & !   9- 12
                                'nStab           ','Quad_i          ','RFcInfo         ','RFiInfo         ', & !  13- 16
                                'RFlInfo         ','SCFInfoI        ','Misc            ','SewIInfo        ', & !  17- 20
                                'SCFInfoI_ab     ','icDmp           ','Symmetry Info   ','Sizes           ', & !  21- 24
                                'IndS            ','LP_A            ','NumCho          ','nFroPT          ', & !  25- 28
                                'nDelPT          ','BasType         ','Spread of Coord.','Unit Cell Atoms ', & !  29- 32
                                'iSOShl          ','Non valence orbi','LoProp nInts    ','LoProp iSyLbl   ', & !  33- 36
                                'nDel_go         ','nBas_Prim       ','IsMM            ','Atom -> Basis   ', & !  37- 40
                                'Logical_Info    ','SCF nOcc        ','SCF nOcc_ab     ','iAOtSO          ', & !  41- 44
                                'iSOInf          ','AuxShell        ','nVec_RI         ','MkNemo.hDisp    ', & !  45- 48
                                'Index ZMAT      ','NAT ZMAT        ','nDisp           ','DegDisp         ', & !  49- 52
                                'LBList          ','Ctr Index Prim  ','MLTP_SINGLE     ','JBNUM_SINGLE    ', & !  53- 56
                                'LROOT_SINGLE    ','GeoInfo         ','Cholesky BkmDim ','Cholesky BkmVec ', & !  57- 60
                                'Atom Types      ','LA Def          ','Basis IDs       ','Desym Basis IDs ', & !  61- 64
                                'primitive ids   ','Root Mapping    ','Fermion IDs     ','IsMM Atoms      ', & !  65- 68
                                'Un_cen Charge   ','PCM_N           ','PCMiSph         ','NVert           ', & !  69- 72
                                'IntSph          ','NewSph          ','XMolnr          ','XEle            ', & !  73- 76
                                'iDmp            ','iDmp:S          ','NSTAT_SINGLE    ','cmsNACstates    ', & !  77- 80
                                'NACstatesOpt    ','                ','                ','                ', & !  81- 84
                                '                ','                ','                ','                ', & !  85- 88
                                '                ','                ','                ','                ', & !  89- 92
                                '                ','                ','                ','                ', & !  93- 96
                                '                ','                ','                ','                ', & !  97-100
                                '                ','                ','                ','                ', & ! 101-104
                                '                ','                ','                ','                ', & ! 105-108
                                '                ','                ','                ','                ', & ! 109-112
                                '                ','                ','                ','                ', & ! 113-116
                                '                ','                ','                ','                ', & ! 117-120
                                '                ','                ','                ','                ', & ! 121-124
                                '                ','                ','                ','                ']   ! 125-128

!> List of known integer scalar labels:
!>
!> - '``Multiplicity``'               The spin multiplicity of the last SCF or RASSCF calculation.
!> - '``nMEP``'                       Number of points on the minimum energy path.
!> - '``No of Internal coordinates``' The number of internal coordinates for the molecule that is allowed within the given point
!>                                    group.
!> - '``nSym``'                       The number of irreducible representations of the molecule.
!> - '``PCM info length``'            Length of the block containing misc. info for the PCM model.
!> - '``Relax CASSCF root``'          Signals which root to perform geometry optimization for in a state average CASSCF geometry
!>                                    optimization.
!> - '``SA ready``'                   Signals that SA wavefunction is ready for gradient calculations.
!> - '``System BitSwitch``'           A bit switch controlling various functions. Will be replaced!
!> - '``Unique atoms``'
!> - '``nActel``'                     The number of active electrons in CASSCF calculation.
!> - '``MkNemo.nMole``'               The number of molecules as specified in the mknemo module.
!> - '``nLambda``'                    The number of constraints in the PCO.
!> - '``DNG``'                        Force numerical gradients.
!> - '``HessIter``'                   Last iteration where the analytical Hessian was computed.
!> - '``TS Search``'                  Flag to mark if a TS search has been activated.
!> - '``CHCCLarge``'                  Segmentation of VOs in CHCC.
!> - '``Seed``'                       The seed number for random number generator used in surface hoping.
!> - '``Rotational Symmetry Number``'
!> - '``mp2prpt``'                    True(=1) if mbpt2 was run with prpt
!> - '``LDF Status``'                 Initialized or not
!> - '``DF Mode``'                    Local (1) or non-local (0) DF
!> - '``LDF Constraint``'             Constraint type for LDF
!> - '``OptimType``'                  Optimization type in hyper
!> - '``STSYM``'                      symmetry of the CAS root(s)
!> - '``nCoordFiles``'                number of xyz-files in gateway
!> - '``CHCCLarge``'                  Segmentation of VOs in CHCC
!> - '``Invert constraints``'
!> - '``Keep old gradient``'
!> - '``embpot``'                     Flag whether an embedding potential is present
!> - '``EFP``'                        Flag Effective fragment potentials
!> - '``Coor_Type``'                  EFP fragment coordinate format
!> - '``nEFP_Coor``'                  Associated number of coordinates per fragment
!> - '``Relax Original root``'
!> - '``NCONF``'                      For MS-PDFT gradient

character(len=lw), parameter :: LabelsIS(nTocIS) = [ &
                                'Multiplicity    ','nMEP            ','No of Internal c','nSym            ', & !   1-  4
                                'PCM info length ','Relax CASSCF roo','System BitSwitch','Unique atoms    ', & !   5-  8
                                'LP_nCenter      ','ChoIni          ','Unit Cell NAtoms','Cholesky Reorder', & !   9- 12
                                'ChoVec Address  ','SA ready        ','NumGradRoot     ','Number of roots ', & !  13- 16
                                'LoProp Restart  ','MpProp nOcOb    ','Highest Mltpl   ','nActel          ', & !  17- 20
                                'Run_Mode        ','Grad ready      ','ISPIN           ','SCF mode        ', & !  21- 24
                                'MkNemo.nMole    ','N ZMAT          ','Bfn Atoms       ','FMM             ', & !  25- 28
                                'Pseudo atoms    ','nChDisp         ','iOff_Iter       ','Columbus        ', & !  29- 32
                                'ColGradMode     ','IRC             ','MaxHops         ','nRasHole        ', & !  33- 36
                                'nRasElec        ','Rotational Symme','Saddle Iter     ','iMass           ', & !  37- 40
                                'mp2prpt         ','NJOB_SINGLE     ','MXJOB_SINGLE    ','NSS_SINGLE      ', & !  41- 44
                                'NSTATE_SINGLE   ','LDF Status      ','DF Mode         ','LDF Constraint  ', & !  45- 48
                                'OptimType       ','STSYM           ','RF CASSCF root  ','RF0CASSCF root  ', & !  49- 52
                                'nCoordFiles     ','nLambda         ','DNG             ','HessIter        ', & !  53- 56
                                'CHCCLarge       ','TS Search       ','Number of Hops  ','hopped          ', & !  57- 60
                                'Invert constrain','Keep old gradien','embpot          ','nPrim           ', & !  61- 64
                                'Seed            ','Track Done      ','MaxHopsTully    ','EFP             ', & !  65- 68
                                'nEFP_fragments  ','Coor_Type       ','nEFP_Coor       ','Relax Original r', & !  69- 72
                                'Unique centers  ','nXF             ','CSPF            ','NCONF           ', & !  73- 76
                                'SH RASSI run    ','isCMSNAC        ','isMECIMSPD      ','CalcNAC_Opt     ', & !  77- 80
                                'MECI_via_SLAPAF ','                ','                ','                ', & !  81- 84
                                '                ','                ','                ','                ', & !  85- 88
                                '                ','                ','                ','                ', & !  89- 92
                                '                ','                ','                ','                ', & !  93- 96
                                '                ','                ','                ','                ', & !  97-100
                                '                ','                ','                ','                ', & ! 101-104
                                '                ','                ','                ','                ', & ! 105-108
                                '                ','                ','                ','                ', & ! 109-112
                                '                ','                ','                ','                ', & ! 113-116
                                '                ','                ','                ','                ', & ! 117-120
                                '                ','                ','                ','                ', & ! 121-124
                                '                ','                ','                ','                ']   ! 125-128

public :: Arr2RunHdr, DS_cache, i_run_CA_used, i_run_DA_used, i_run_DS_used, i_run_IA_used, i_run_IS_used, icRd, icWr, IDRun, &
          IS_cache, LabelsCA, LabelsDA, LabelsDS, LabelsIA, LabelsIS, lw, nHdrSz, nToc, nTocCA, nTocDA, nTocDS, nTocIA, nTocIS, &
          NulPtr, num_DS_init, num_IS_init, rcNotFound, rcOK, rcWrongType, RnNmStk, RunHdr, RunHdr2Arr, RunName, sNotUsed, &
          sRegularField, sSpecialField, Toc, TypDbl, TypInt, TypLgl, TypStr, TypUnk, VNRun

contains

subroutine RunHdr2Arr(Arr)

  integer(kind=iwp), intent(out) :: Arr(nHdrSz)

  Arr(1) = RunHdr%ID
  Arr(2) = RunHdr%Ver
  Arr(3) = RunHdr%Next
  Arr(4) = RunHdr%Items
  Arr(5) = RunHdr%DaLab
  Arr(6) = RunHdr%DaPtr
  Arr(7) = RunHdr%DaLen
  Arr(8) = RunHdr%DaMaxLen
  Arr(9) = RunHdr%DaTyp
  Arr(10:) = 0

end subroutine RunHdr2Arr

subroutine Arr2RunHdr(Arr)

  integer(kind=iwp), intent(in) :: Arr(nHdrSz)

  RunHdr%ID = Arr(1)
  RunHdr%Ver = Arr(2)
  RunHdr%Next = Arr(3)
  RunHdr%Items = Arr(4)
  RunHdr%DaLab = Arr(5)
  RunHdr%DaPtr = Arr(6)
  RunHdr%DaLen = Arr(7)
  RunHdr%DaMaxLen = Arr(8)
  RunHdr%DaTyp = Arr(9)

end subroutine Arr2RunHdr

end module RunFile_data
