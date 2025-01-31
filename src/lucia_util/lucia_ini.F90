!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Lucia_Ini()

use rasscf_lucia, only: Sigma_on_disk, ini_h0
use lucia_data, only: ECORE
use lucia_data, only: NGAS, NGSSH, NCISPC, NCMBSPC, ICMBSPC, IGSOCCX, LCMBSPC
use lucia_data, only: IPRCIX, IPRORB, IPRDIA, IPRXT, IPROCC, IPRDEN, IPRRSP, IPRPRO, IPRCC, IPRNCIV
use lucia_data, only: INTIMP, ENVIRO, MXCIV, ICISTR, MXINKA, THRES_E, INOCALC, ISAVE_EXP, IEXPAND, LCSBLK, NOMOFL, IDENSI, &
                      IUSE_PH, IADVICE, ITRACI, ITRACI_CR, ITRACI_CN, IDIAG, MXP1, MXP2, MXQ, MAXIT, IRESTR, INCORE, MOCAA, MOCAB, &
                      IUSE_PA, NOINT, ICJKAIB, INIREF, IRESTRF, IPERT, NPERT, IAPRREF, IAPRZER, IEXTKOP, IC1DSC, IH0SPC, NPTSPC, &
                      IRFROOT, NH0EXSPC, INIDEG, XLAMBDA, IFINMO, E_THRE, C_THRE, E_CONV, C_CONV, ICLSSEL, IPTEKT, IH0ROOT, IRST2, &
                      ISKIPEI, NPROP, ITRAPRP, IRESPONS, NRESP, N_AVE_OP, MXITLE, IROOTHOMING, IPTFOCK, ITRA_FI, ITRA_IN, MULSPC, &
                      IFMULSPC, LPAT, NCNV_RT, I_RE_MS2_SPACE, I_RE_MS2_VALUE, ISIMSYM, NOCSF, IE0AVEX, IEXCSYM, NEXCSTATE, &
                      NPSSPC, IH0EXSPC, IH0INSPC, IOCPTSPC, ISEQCI, IXYZSYM, MAXORD_OP, NPSSH, PROPER, RESP_OP, RESP_W, AVE_OP, &
                      NSEQCI, CSEQCI
use lucia_data, only: MS2, MULTS, IREFSM, NROOT, PSSIGN, INTSEL, PLSIGN, IDC, IROOT
use lucia_data, only: I_ELIMINATE_GAS, N_ELIMINATED_GAS, N_2ELIMINATED_GAS, I2ELIMINATED_IN_GAS, IELIMINATED_IN_GAS
use lucia_data, only: irat
use lucia_data, only: EXTSPC, PNTGRP, NIRREP, NSMCMP, NSMOB, MAXML, MAXL, NACTEL, MXR4TP, MXER4, MNRS1R, MNRS10, MXRS3R, MXRS30, &
                      INTSPC, MNRS1RE, MXRS3RE, MNRS1ZE, MXRS3ZE, NDELSH, NINASH, NRS0SH, NRS4SH, NRSSH
use lucia_data, only: IPART
use lucia_data, only: NMOS_ENV, NAOS_ENV
use spinfo, only: NSYM_MOLCAS, NACTEL_MOLCAS, MS2_MOLCAS, ISPIN_MOLCAS, LSYM_MOLCAS, NROOTS_MOLCAS, NGAS_MOLCAS, THRE_MOLCAS, &
                      ITMAX_MOLCAS, INOCALC_MOLCAS, ISAVE_EXP_MOLCAS, IEXPAND_MOLCAS, IPT2_MOLCAS, I_ELIMINATE_GAS_MOLCAS, &
                      N_ELIMINATED_GAS_MOLCAS, N_2ELIMINATED_GAS_MOLCAS, IPRCI_MOLCAS, POTNUC_MOLCAS, I2ELIMINATED_IN_GAS_MOLCAS, &
                      IELIMINATED_IN_GAS_MOLCAS, IGSOCCX_MOLCAS, ISPEED, NBAS_MOLCAS, NGSSH_MOLCAS, NISH_MOLCAS, NORB_MOLCAS
use Definitions, only: RtoI

implicit none
integer, parameter :: MXPKW = 125
integer isetkw(MXPKW)
integer IEXPERT, NERROR, NWARN, I, IRREP, IGAS, J, NMISS, IR4TP, ICISPC, JPTSPC, IUSED, IDOPERT, JCMBSPC, JSEQCI, ICOMP, IPROP, &
        IAVE, IRESP, IPSSPC
real*8 ECORE_ENV

! ================================================
!  Some initialization to avoid compiler warnings
! ================================================
!MSCOMB_CC = 0
!I_USE_NEWCCP = 0
! =======================
!  Some initial settings
! =======================
INI_H0 = 1
! Flag for compatibility with normal MOLCAS input format

IEXPERT = 0
NERROR = 0
NWARN = 0
EXTSPC = 0
! No cc as default
! Start out with normal integrals
!I_USE_SIMTRH = 0
! Initialize flag to be used by densi_master to check whether the
! sigma-vector is on disk.
Sigma_on_disk = .false.

! ===================
!  Initialize isetkw
! ===================
do i=1,mxpkw
  isetkw(i) = 0
end do

! =========================
!  Point group of orbitals
! =========================
pntgrp = 1

! ==============================
!  Number of irreps of orbitals
! ==============================
nirrep = nsym_molcas
nsmcmp = nsym_molcas
nsmob = nsym_molcas
maxml = -1
maxl = -1

! ============================
!  Number of active electrons
! ============================
nactel = nactel_Molcas

! =================
!  Inactive shells
! =================
do irrep=1,nirrep
  ninash(irrep) = nish_molcas(irrep)
end do

! ==================
!  Secondary shells
! ==================

! ===========================
!  Two times spin projection
! ===========================
ms2 = ms2_molcas

! ===================
!  Spin multiplicity
! ===================
mults = iSpin_molcas

! ====================
!  Reference symmetry
! ====================
irefsm = lsym_molcas

! =======
!  Roots
! =======
nroot = nroots_molcas
do i=1,nroot
  !iroot(i) = iroot_molcas(i)
  iroot(i) = 0
end do

! ============
!  Enviroment
! ============
intimp = 6
enviro(1:6) = 'RASSCF'

! ====================
! 27: Ms combinations
! ====================
if ((MULTS == 1) .and. (ISPEED(1) == 1)) then
  ISETKW(27) = 1
  PSSIGN = 1.0d0
else
  ISETKW(27) = 0
end if

! ==========================================================
! 33 : Largest allowed number of Vectors in diagonalization
! ==========================================================
MXCIV = 2*NROOT

! ==========================
!  Storage mode for vectors
! ==========================
icistr = 2
! =================================
!  Tell Lucia to use CSF-expansion
! =================================
!noscf = 0

! =============================================
! 39 : Define dimension of resolution matrices
! =============================================
MXINKA = 150

! ==========================================================
!  Generalized active space concept invoked, orbital spaces
! ==========================================================
ngas = ngas_molcas
ngssh = 0
do irrep=1,nirrep
  do igas=1,ngas
    ngssh(irrep,igas) = ngssh_molcas(igas,irrep)
  end do
end do

! ==================================================
!  Generalized active space occupation restrictions
! ==================================================
ncispc = 1
do i=1,ngas
  do j=1,2
    igsoccx(i,j,1) = igsoccx_molcas(i,j)
  end do
end do

! ==========================
!  Energy convergence of CI
! ==========================
thres_e = thre_molcas

! ==========================================
!  Sequency of calculations to wavefunction
! ==========================================
ncmbspc = ncispc
nseqci(1) = 1
cseqci(1,1) = 'CI'
iseqci(1,1) = itmax_molcas

! =======================================================
!  No calculation, save CI vector info, expand CI vector
! =======================================================
INOCALC = INOCALC_MOLCAS
ISAVE_EXP = ISAVE_EXP_MOLCAS
IEXPAND = IEXPAND_MOLCAS

! ================================================
!  Length of smallest block og C an Sigma vectors
! ================================================
lcsblk = 3000000

! ===================
! 66 : No MO-AO file
! ===================

! No MO-AO file
NOMOFL = 1

! =================================
!  Calculation of density matrices
! =================================
idensi = 2

! ==================================================================
!  Isimsym: Treat all symmetryblocks with given type simultaneously
! ==================================================================

! ==============================================
!  Particle hole simplifications, default is no
! ==============================================
IUSE_PH = 0

! ============================================
! 87 : Allow the sigma routine to take advice
! ============================================
!
IADVICE = 1

! ============================================================
!  Transform CI vectors to alternative orbital representation
! ============================================================
itraci = 1
itraci_cr = 'REST'
itraci_cn = ''
if (ipt2_molcas == 1) then
  itraci_cn = 'CANO'
else
  itraci_cn = 'NATU'
end if

! ==============================
!  HEXS - Highly excited states
!  DEXS - Doubly excited states
! ==============================
I_ELIMINATE_GAS = I_ELIMINATE_GAS_MOLCAS
N_ELIMINATED_GAS = N_ELIMINATED_GAS_MOLCAS
N_2ELIMINATED_GAS = N_2ELIMINATED_GAS_MOLCAS
!write(6,*) I_ELIMINATE_GAS_MOLCAS,N_ELIMINATED_GAS_MOLCAS
do IGAS=1,N_ELIMINATED_GAS
  IELIMINATED_IN_GAS(IGAS) = IELIMINATED_IN_GAS_MOLCAS(IGAS)
end do
do IGAS=1,N_2ELIMINATED_GAS
  I2ELIMINATED_IN_GAS(IGAS) = I2ELIMINATED_IN_GAS_MOLCAS(IGAS)
end do

! =============
!  Printlevels
! =============
iprcix = iprci_molcas-2
iprorb = 0
iprdia = 0
iprxt = 0
iprocc = 0
iprden = (iprci_molcas-2)/2
iprrsp = 0
iprpro = 0
iprcc = 0

! ==============
!  Set defaults
! ==============

!***********************************************************************
!                                                                      *
! Part 2 : Insert defaults for missing optional keywords               *
!          and print error messages for missing mandatory keywords     *
!                                                                      *
!***********************************************************************

NMISS = 0

! 8 : Core orbitals, only of interest if EXTSPC /= 0

call ISETVC(NRS0SH,0,NIRREP)

! 9 : RAS 1 orbitals

call ISETVC(NRSSH(1,1),0,NIRREP)

! 10 : RAS 2 orbitals

call ISETVC(NRSSH(1,2),0,NIRREP)

! 11 : RAS 3 orbitals

call ISETVC(NRSSH(1,3),0,NIRREP)

! 13 : Secondary space

MXR4TP = 1
do IR4TP=1,MXR4TP
  call ISETVC(NRS4SH(1,IR4TP),0,NIRREP)
end do
MXER4 = 0

! 14 : occupation restrictions for Reference space

MNRS1R = MNRS10
MXRS3R = MXRS30

! 15 : Selection of active configurations

! Standard is no selection
INTSEL = 0

! 20 : Diagonalization routine

! Standard is currently MICDV*
IDIAG = 1

! 21 : Explicit Hamiltonian

! Default is no explicit Hamiltonian
MXP1 = 0
MXP2 = 0
MXQ = 0

! 22 : Largest allowed number of CI iterations per root

! Default is 20 ( not active I expect )
MAXIT = 20

! 23 : Restart option

! Default is no explicit Hamiltonian
IRESTR = 0

! 25 : INCORE option for integrals

if (EXTSPC == 0) then
  INCORE = 1
else
  INCORE = 0
end if

! 26 : DELETEd shells

! If CAS + Active have been set or RAS + Ras3 have been set,
! obtain for MOLCAS Interface from number of basis functions
if (INTSPC /= 1) call ISETVC(NDELSH,0,NIRREP)

! 27 : Ms combinations

if (ISETKW(27) == 0) then
  PSSIGN = 0.0d0
  ISETKW(27) = 2
else if (MS2 /= 0) then
  write(6,*) ' Spin combinations only allowed with MS2 = 0'
  write(6,*) ' Your value of MS2 = ',MS2,' differs from zero'
  write(6,*) ' LUCIA will neglect your nice suggestion'
  write(6,*) ' to use spin combinations'
  PSSIGN = 0.0d0
  ISETKW(27) = 2
end if

! 28 : Ml combinations

PLSIGN = 0.0d0

if (PSSIGN == 0.0d0) then
  IDC = 1
else if (PSSIGN /= 0.0d0) then
  IDC = 2
end if

!write(6,* ) ' TEST readin IDC = ',IDC

! =======================================================================
! 44 : Use Minimal operatioon count method for alpha-alpha and beta-beta
! =======================================================================
! TEST OF PERFORMANCE
MOCAA = ISPEED(3)
MOCAB = ISPEED(4)

IUSE_PA = 0

! 33 : Number of Ci vectors in subspace

if ((IDIAG == 2) .and. (MXCIV > 2)) then
  MXCIV = 2
  NWARN = NWARN+1
  write(6,*) ' Warning : You have specified TERACI'
  write(6,*) '           I allow myself to set MXCIV = 2'
  write(6,*)
  write(6,*) '                   Best Wishes'
  write(6,*) '                      Lucia'
end if

! 34 : CI storage mode

! 35 : Employ CSF expansion ?

! Default is no ( only possibility at the moment )
! CSF expansion must only be used when two vectors are stored in CORE

! 37 : Avoid any readin of integrals ( useful for obtaining
!      size of CI expansion etc.

NOINT = 0

! 38 : Dump integrals in formatted form : Default is No

!IDMPIN = 0

! 40 : Use CJKAKB intermediate matrices in alpha-beta loop,
!      Default is  YES !!!!!

ICJKAIB = 1

!  41 : Initial CI in reference space, default is : No

INIREF = 0

!  42 : Restart with CI in reference space

IRESTRF = 0

! Core energy : Default is 0 / MOLCAS : Value read in !

ECORE = 0.0d0

! Use perturbation theory for zero order space . def is no !

if (ISETKW(47) == 0) then
  IPERT = 0
  NPERT = 0
  ISETKW(47) = 2
  ! Else ensure that a CI in reference space is performed
else
  INIREF = 1
end if

! 48 : Approximate Hamiltonian in reference space : NO !!

if (ISETKW(48) == 0) then
  IAPRREF = 0
  MNRS1RE = MNRS1R
  MXRS3RE = MXRS3R
  ISETKW(48) = 2
end if

! 49 : Approximate Hamiltonian in zero order space : NO !!

if (ISETKW(49) == 0) then
  IAPRZER = 0
  MNRS1ZE = MNRS10
  MXRS3ZE = MXRS30
  ISETKW(49) = 2
end if

! 50 : GAS shells must be defined

! 52 : Combination of gasspaces : Default is just to take each  space
!      By itself

if (ISETKW(52) == 0) then
  NCMBSPC = NCISPC
  do ICISPC=1,NCISPC
    LCMBSPC(ICISPC) = 1
    ICMBSPC(1,ICISPC) = ICISPC
  end do
  ISETKW(52) = 2
end if

! 53 : Convergence threshold for CI

! 54 : General sequencer : default is just straight sequence
!      of CI with default number of iterations

! 55 : EKT calculation : Default is no

if (ISETKW(55) == 0) then
  IEXTKOP = 0
  ISETKW(55) = 2
end if

!. 56 : Default Machine : Good old BIM machine

! 57 : Allow first order correction to be saved on DISC
!     (For vector free calculations )
!     Default is : NO !!
if (ISETKW(57) == 0) then
  IC1DSC = 0
  ISETKW(57) = 2
end if

! 58 : Restrictions on interactions of perturbation

! Default is : no
if (ISETKW(58) == 0) then
  IH0SPC = 0
  ISETKW(58) = 2
end if

! 59 : Type of perturbation in subspaces spaces

! Default is specified by IPART from keyword PERTU
if (ISETKW(59) == 0) then
  ISETKW(59) = 2
  if (IH0SPC /= 0) then
    do JPTSPC=1,NPTSPC
      IH0INSPC(JPTSPC) = IPART
    end do
  end if
end if

! 60 : Reference Root, default is NROOT

! Should be less or equal to NROOT
if (ISETKW(60) == 1) then
  if (IRFROOT > NROOT) then
    write(6,*) ' Reference root (RFROOT) larger'
    write(6,*) ' than total number of roots (NROOT)'
    write(6,*) ' CHANGE NROOT or RFROOT'
    NMISS = NMISS+1
  end if
end if

if (ISETKW(60) == 0) then
  ISETKW(60) = 2
  IRFROOT = NROOT
end if

! 61 : H0EX : Orbital spaces in which exaxt Hamiltonian is used
!      No default

! Is H0EX required ( Has H0FRM = 5 been used )
IUSED = 0
if (ISETKW(59) == 1) then
  IUSED = 0
  do JPTSPC=1,NPTSPC
    if (IH0INSPC(JPTSPC) == 5) IUSED = 1
  end do
end if
if ((IUSED == 0) .and. (ISETKW(61) == 0)) then
  ! No exact spaces included and none have been defined !
  NH0EXSPC = 0
  IH0EXSPC(1) = -1
end if
if ((IUSED == 1) .and. (ISETKW(61) == 0)) then
  ! Needed, but not supplied
  write(6,*) ' You have specified that zero order operator'
  write(6,*) ' Include exact Hamilton operator in subspace'
  write(6,*) ' You should then also supply Keyword H0EX'
  NMISS = NMISS+1
end if

! If perturbation theory will be invoked be sure that the
! form of perturbation theory has been specified through
! KEYWORD PERTU ( number 47 as you maybe know )
IDOPERT = 0
do JCMBSPC=1,NCMBSPC
  do JSEQCI=1,NSEQCI(JCMBSPC)
    if (ISEQCI(JSEQCI,JCMBSPC) == -5) IDOPERT = 1
  end do
end do

if ((IDOPERT == 1) .and. (IPERT == 0)) then
  write(6,*) ' Perturbation theory will be used'
  write(6,*) ' Please specify form through PERTU keyword'
  NMISS = NMISS+1
end if

! 62 : Default Handling of degenrences of initial CI vectors
!      Default is : No action

if (ISETKW(62) == 0) then
  INIDEG = 0
  ISETKW(62) = 2
end if

! 63 : Use F + Lambda(H-F) as operator instead of H
!      Default is : No i.e Lambda = 1

if (ISETKW(63) == 0) then
  XLAMBDA = 1.0d0
  ISETKW(63) = 2
end if

! 64 : Smallest block in batch of C and sigma
!      Default is zero

! 66 : NO MO file : Default is no access to MO-AO file

! 68 : Type of natural orbitals, default is natural orbitals

if (ISETKW(68) == 0) then
  ISETKW(68) = 2
  IFINMO = 1
end if

! 69 : Default Threshold for individual energy correction = 0.0

if (ISETKW(69) == 0) then
  E_THRE = 0.0d0
  ISETKW(69) = 2
end if

! 70 : Default Threshold for wave individual function corrections = 0.0

if (ISETKW(70) == 0) then
  C_THRE = 0.0d0
  ISETKW(70) = 2
end if

! 71 : Default Threshold for total energy corrections = 0.0

if (ISETKW(71) == 0) then
  E_CONV = 0.0d0
  ISETKW(71) = 2
end if

! 72 : Default Threshold for total wave function correction = 0.0

if (ISETKW(72) == 0) then
  C_CONV = 0.0d0
  ISETKW(72) = 2
end if

! 73 : Perform Class selection : Default if Yes if TERACI is used

if (ISETKW(73) == 0) then
  if (IDIAG == 1) then
    ICLSSEL = 0
  else if (IDIAG == 2) then
    ICLSSEL = 1
  end if
  ISETKW(73) = 2
end if

! 74 : Calculation of density matrices : Default is
!       calculation of one-body density

! If IDENSI was set to zero and properties were requested
! overwrite input to obtain 1-el matrix
if ((IDENSI == 0) .and. (ISETKW(80) == 1)) then
  write(6,*) ' You have specified calculation of'
  write(6,*) ' one-electron properties, and this'
  write(6,*) ' requires the calculation of the'
  write(6,*) ' one-electron density.'
  write(6,*)
  write(6,*) ' You have, however, specified IDENSI=0'
  write(6,*) ' corresponding  to no densities'
  write(6,*)
  write(6,*) ' I will allow myself to modify your'
  write(6,*) ' input to allow calculation of the'
  write(6,*) ' one-electron densities, so property'
  write(6,*) ' calculation can proceed as planned'
  write(6,*)
  write(6,*) ' Lucia'
  ! and do it
  IDENSI = 1
end if
! If CC is performed, one- and two- particle densities are
! used in current simple-minded implementation.
! Two-electron density also needed for MCSCF

! 75 : Perturbation expansion of EKT, default is no

if (ISETKW(75) == 0) then
  IPTEKT = 0
  ISETKW(75) = 2
end if

! 76 : Root for zero order operator, default is NROOT

! Should be less or equal to NROOT
if (ISETKW(76) == 1) then
  if (IH0ROOT > NROOT) then
    write(6,*) ' Zero order operator root (H0ROOT) larger'
    write(6,*) ' than total number of roots (NROOT)'
    write(6,*) ' CHANGE NROOT or H0ROOT'
    NMISS = NMISS+1
  end if
end if
if (ISETKW(76) == 0) then
  ISETKW(76) = 2
  IH0ROOT = NROOT
end if

! 77 : NO restart from previous vectors in calculation 2
!      Deafault is NO NO, ie. restart in calc 2

if (ISETKW(77) == 0) then
  IRST2 = 1
  ISETKW(77) = 2
end if

! 78 : skip initial energy evaluations - if possible

if (ISETKW(78) == 0) then
  ISKIPEI = 1
  ISETKW(78) = 2
end if

! 79 : Symmetry of x,y,z - needed for property calculations

if (ISETKW(79) == 0) then
  ! Problematic if Properties should be calculated
  if ((ISETKW(80) == 1) .or. (ISETKW(81) == 1) .or. (ISETKW(82) == 1)) then
    write(6,*) ' Symmetry of X,Y,Z has not been given'
    write(6,*) ' You have to specify this for property calc'
    write(6,*) ' Please add KEYWORD XYZSYM'
    NMISS = NMISS+1
    ISETKW(79) = -1
  else
    ! Is not needed, just supply zeroes
    do ICOMP=1,3
      IXYZSYM(ICOMP) = 0
    end do
    ISETKW(79) = 2
  end if
end if

! 80 : Property calculation, default is no

if (ISETKW(80) == 0) then
  NPROP = 0
  ISETKW(80) = 2
end if

! 81 : Transition properties, default is no

if (ISETKW(81) == 0) then
  ITRAPRP = 0
  ISETKW(81) = 2
end if

! 82 : Response properties, default is no

if (ISETKW(82) == 0) then
  IRESPONS = 0
  ISETKW(82) = 2
  NRESP = 0
  N_AVE_OP = 0
end if
! Properties should be defined if transition properties are
! invoked
if ((ITRAPRP /= 0) .and. (NPROP == 0)) then
  write(6,*) ' You have specified transition property calculation'
  write(6,*) ' (keyword TRAPRP) but no property labels have been supplied'
  write(6,*) '(Keyword PROPER). Transition densities will be obtained'
end if

! 83 : Max number of iterations in linear equations

if (ISETKW(83) == 0) then
  MXITLE = 20
  ISETKW(83) = 2
end if

! 85 : Root homing, default is no

if (ISETKW(85) == 0) then
  IROOTHOMING = 0
  ISETKW(85) = 2
end if

! 90 : Perturbation expansion of Fock matrix : default is no

if (ISETKW(90) == 0) then
  IPTFOCK = 0
  ISETKW(90) = 2
end if

! 91 : Print final CI vectors : default is no

if (ISETKW(91) == 0) then
  IPRNCIV = 0
  ISETKW(91) = 2
end if

! 93 : End Calculation with integral transformation

if (ISETKW(93) == 0) then
  ITRA_FI = 0
  ISETKW(93) = 2
end if

! 94 : Initialize Calculation with integral transformation

if (ISETKW(94) == 0) then
  ITRA_IN = 0
  ISETKW(94) = 2
end if
! Requires access to MO-AO file
if ((ITRA_IN == 1) .and. (NOMOFL == 1)) then
  write(6,*) ' Integral transformation required,'
  write(6,*) ' but no mo-ao file accessible'
  write(6,*) ' REMOVE KEWORD NOMOFL'
  ISETKW(94) = -1
  NERROR = NERROR+1
end if

! 95 : Multispace optimization in each run, default is no

if (ISETKW(95) == 0) then
  MULSPC = 0
  IFMULSPC = 0
  LPAT = 0
  ISETKW(95) = 2
end if

! Expert mode ( neglect mistyped keywords ) : default is no expert

if (ISETKW(97) == 0) then
  IEXPERT = 0
  ISETKW(97) = 2
end if

! Number of roots to be converged : default is total number of roots

if (ISETKW(98) == 0) then
  NCNV_RT = NROOT
  ISETKW(98) = 2
end if

! 100 : Do quantum dot calculation, default is no

!if (ISETKW(100) == 0) then
!  IDODQ = 0
!  ISETKW(100) = 2
!end if

! 101: Restrict MS2 at some intermediate level : default is no way

if (ISETKW(101) == 0) then
  I_RE_MS2_SPACE = 0
  I_RE_MS2_VALUE = 0
  ISETKW(101) = 2
end if

! 103 : Treat all TT blocks with given types simultaneously : Default is no

ISIMSYM = 0

! Thresholds only active in connection with IDIAG = 2,
! Check and maybe issue a warning
if (IDIAG == 2) then
  ! Check to ensure that zero or two thresholds were  set,
  if (ISETKW(69) /= ISETKW(70)) then
    write(6,*) ' Only a single threshold (E_THRE or C_THRE)'
    write(6,*) ' on individual determinants given.'
    write(6,*) ' One of the thresholds vanishes therefore and'
    write(6,*) ' all determinants will therefore be included'
    write(6,*)
    write(6,*) '                   Warns'
    write(6,*)
    write(6,*) '                   LUCIA'
  end if
else
  ! Good old diagonalization, thresholds not active
  if ((ISETKW(69) == 1) .or. (ISETKW(70) == 1)) then
    write(6,*) ' Thresholds on selection of individual coefficients'
    write(6,*) ' are only active in connection with keyword TERACI'
    write(6,*)
    write(6,*) '                   Warns'
    write(6,*)
    write(6,*) '                   LUCIA'
  end if
end if

if (NMISS /= 0) then
  write(6,'(1X,A,I9)') ' Number of missing required keyword ',NMISS
  write(6,*) ' You have wounded me I give up'
  write(6,*)
  write(6,*)
  write(6,*)
  write(6,*) '     An expert is a man who has made all the mistakes,'
  write(6,*) '     which can be made, in a very narrow field'
  write(6,*)
  write(6,*) '                                      Niels Bohr'
  if (IEXPERT == 0) then
    !stop
    call SYSABENDMSG('lucia_util/lucia_ini','Input error','')
  else
    write(6,*) ' Processing continues (EXPERT mode )'
  end if
end if
! Open one-electron file to obtain core energy and
! Number of MO's and AO's
if ((NOINT == 0) .and. (ENVIRO(1:4) /= 'NONE')) then
  ecore_env = potnuc_molcas
  ! Manual setting of naos_env and nmos_env
  nmos_env(1) = nbas_molcas(1)
  nmos_env(2) = nbas_molcas(2)
  nmos_env(3) = nbas_molcas(4)
  nmos_env(4) = nbas_molcas(3)
  nmos_env(5) = nbas_molcas(7)
  nmos_env(6) = nbas_molcas(6)
  nmos_env(7) = nbas_molcas(8)
  nmos_env(8) = nbas_molcas(5)
  naos_env(1) = norb_molcas(1)
  naos_env(2) = norb_molcas(2)
  naos_env(3) = norb_molcas(4)
  naos_env(4) = norb_molcas(3)
  naos_env(5) = norb_molcas(7)
  naos_env(6) = norb_molcas(6)
  naos_env(7) = norb_molcas(8)
  naos_env(8) = norb_molcas(5)
  ECORE = ECORE_ENV
else
  write(6,*) ' GETOBS and CHK_ORBDIM not called'
  ECORE = 0.0d0
end if
! Check number of orbitals and insert occupations for ALL/REST

!***********************************************************
!                                                          *
! Part 3 : Print input                                     *
!                                                          *
!***********************************************************

! Machine in use

! Type of reference state
if (IPRCIX >= 100) then
  write(6,*)
  write(6,*) '*************************************'
  write(6,*) '*  Symmetry and spin of CI vectors  *'
  write(6,*) '*************************************'
  write(6,*)
  ! Point group
  write(6,'(1X,A)') '     Point group ............ D2H'
  ! Spatial symmetry
  write(6,'(1X,A,I1)') '     Spatial symmetry ....... ',IREFSM
  ! Spin
  write(6,'(1X,A,I2)') '     2 times spinprojection  ',MS2
  ! Number of active electrons
  write(6,'(1X,A,I2)') '     Active electrons .....  ',NACTEL
  write(6,*)
  write(6,*) '*********************************************'
  write(6,*) '*  Shell spaces and occupation constraints  *'
  write(6,*) '*********************************************'
  write(6,*)

  if (XLAMBDA /= 1.0d0) then
    write(6,*)
    write(6,'(A,F13.8)') ' Modified operator H(l) = l*F + l*(H-F) used with l =',XLAMBDA
    !if (IUSEH0P == 0) then
    write(6,'(A)') ' Zero-order operator without projection used'
    !else
    !  write(6,'(A)') ' Zero-order operator with projection used'
    !end if
    if (IRESTR == 0) then
      write(6,*) ' Notice : This madness starts  in second calculation'
    else
      write(6,*) ' You have specified a calculation with modified'
      write(6,*) ' Hamiltonian (the LAMBDA option) and RESTART'
      write(6,*) ' so this is what I will do'
      write(6,*)
      write(6,*) '   1:) Perform CI in space 1 to obtain Hamiltonian'
      write(6,*) '       (no RESTART in this space )'
      write(6,*) '   2:) CI calculation in space 2  with'
      write(6,*) '       modified Hamiltonian and RESTART from LU21'
      write(6,*) ' Space 2 should therefore correspond to the'
      write(6,*) ' restarted calculation'
    end if
  end if

  write(6,*)
  write(6,*) '***********'
  write(6,*) '*  Roots  *'
  write(6,*) '***********'
  write(6,*)
  write(6,'(1X,A,I3)') '     Number of roots to be included  ',NROOT
  write(6,'(1X,A,(20I3))') '     Roots to be obtained ',(IROOT(I),I=1,NROOT)
  write(6,'(1X,A,I3)') '     Number of roots to be converged ',NCNV_RT

  write(6,*)
  write(6,*) '**************************'
  write(6,*) '*  Run time definitions  *'
  write(6,*) '**************************'
  write(6,*)
  ! Program environment
  write(6,'(A,A6)') '      Program environment... ',ENVIRO

  if (NOINT == 1) then
    ! Integral import
    write(6,'(1X,A)') '     No integrals will be read in'
  else if (NOINT == 0) then
    ! Integral storage
    if (INCORE == 1) write(6,'(1X,A)') '     All integrals stored in core'
  end if
  write(6,*)
  ! ( END IF for NOINT
  ! CSF or SD expansion
  if (NOCSF == 0) then
    write(6,'(1X,A)') "     CI optimization performed with CSF's"
  else
    write(6,'(1X,A)') "     CI optimization performed with SD's"
  end if
  ! Ms,Ml combinations
  if (ISETKW(27) == 1) write(6,'(1X,A,F8.3)') '     Spin combinations used with sign ',PSSIGN
  if (ISETKW(28) == 1) write(6,'(1X,A,F8.3)') '     ML   combinations used with sign ',PLSIGN
  ! Initial approximation to vectors
  write(6,*)
  if ((IRESTR == 1) .and. (IRESTRF == 0)) then
    write(6,'(1X,A)') '     Restarted calculation'
  else if (IRESTRF == 1) then
    write(6,'(1X,A)') '     Restarted calculation from REFERENCE space expansion'
  else
    if (MXP1 /= 0) then
      write(6,'(1X,A)') '     Initial vectors obtained from explicit Hamiltonian'
    else if (MXP1 == 0) then
      write(6,'(1X,A)') '     Initial vectors obtained from diagonal'
    end if
  end if
  ! Handling of degenerencies of initial vectors
  if (INIDEG == 1) then
    write(6,'(1X,A)') '     Symmetric combination of degenerate initial vectors'
  else if (INIDEG == -1) then
    write(6,'(1X,A)') '     Antiymmetric combination of degenerate initial vectors'
  else if (INIDEG == 0) then
    write(6,'(1X,A)') '     No combination of degenerate initial vectors'
  end if
  ! No calculation, save CI vector information and expand CI vector
  if (INOCALC == 1) write(6,'(1X,A)') ' No calculation will be performed'
  if (ISAVE_EXP == 1) write(6,'(1X,A)') ' Save CI vector information'
  if (IEXPAND == 1) write(6,'(1X,A)') ' Expand shorter CI vector in longer one'
  ! Ms,Ml combinations
  !if (ISETKW(27) == 1) write(6,'(1X,A,F8.3)') '     Spin combinations used with sign ',PSSIGN
  !if (ISETKW(28) == 1) write(6,'(1X,A,F8.3)') '     ML   combinations used with sign ',PLSIGN
  ! CI storage mode
  write(6,*)

  write(6,*) '     3 symmetry blocks will be held in core'

  if (LCSBLK /= 0) write(6,'(A,I10)') '      Smallest allowed size of sigma- and C-batch ',LCSBLK
  write(6,'(1X,A,I4)') '     Dimension of block of resolution strings ',MXINKA
  !if (IUSE_PH == 1) then
  !  WRITE(6,'(1X,A)') '     Particle-hole separation used'
  !else
  write(6,'(1X,A)') '      Particle-hole separation not used'
  !end if

  if (IADVICE == 1) write(6,'(1X,A)') '     Advice routine call to optimize sigma generation'

  !if (IUSE_PA == 1) then
  !  write(6,'(1X,A)') '     Strings divided into active and passive parts'
  !else
  write(6,'(1X,A)') '     Strings not divided into active and passive parts'
  !end if
  !if (ISIMSYM == 1) write(6,'(1X,A)') '     ALl TTS blocks with given types treated in sigma'
  !if (IUSE_HW == 1) write(6,*) ' Hardwired routines in use'

  write(6,*)
  if (IDENSI == 0) then
    write(6,'(1X,A)') '     No calculation of density matrices'
  else if (IDENSI == 1) then
    write(6,'(1X,A)') '     One-body density matrix calculated'
  else if (IDENSI == 2) then
    write(6,'(1X,A)') '     One- and two-body density matrices  calculated'
  end if
  write(6,*)
  !if (MOCAA /= 0) WRITE(6,*) '     MOC method used for alpha-alpha+beta-beta loop'
  !if (MOCAB /= 0) WRITE(6,*) '     MOC method used for alpha-beta loop'

  ! Diagonalization information
  write(6,'(1X,A)') '     CI diagonalization :'
  write(6,'(1X,A)') '     ===================='
  ! Subspace Hamiltinian
  if (MXP1+MXP2+MXQ == 0) then
    write(6,'(1X,A)') '        No subspace Hamiltonian'
  else
    write(6,'(1X,A,3I4)') '        Dimensions of subspace Hamiltonian ',MXP1,MXP2,MXQ
  end if
  ! Diagonalizer
  if ((IDIAG == 1) .and. (ICISTR == 1)) then
    write(6,'(1X,A)') '        Diagonalizer : MINDV4'
  else if ((IDIAG == 1) .and. (ICISTR >= 2)) then
    write(6,'(1X,A)') '        Diagonalizer : MICDV*'
  else if (IDIAG == 2) then
    write(6,'(1X,A)') '        Diagonalizer : PICO*'
  end if
  write(6,'(1X,A)') '        Simple diagonal used as preconditioner'
  ! Root homing
  if (IROOTHOMING == 1) then
    write(6,'(1X,A)') '        Root homing will be used'
  else
    write(6,'(1X,A)') '        No root homing'
  end if
  ! No restart in CI calc 2
  if (IRST2 == 0) write(6,'(1X,A)') '        No restart from previous vectors in second calc'
  if (ISKIPEI == 1) then
    write(6,'(1X,A)') '        Initial energy evaluations skipped after first calc'
    write(6,'(1X,A)') '        (Only active in connection with TERACI )'
  end if
  ! Number of iterations
  !write(6,'(1X,A,I2)') '        Allowed number of iterations    ',MAXIT
  ! Number of CI vectors in subspace
  write(6,'(1X,A,I2)') '        Allowed Dimension of CI subspace',MXCIV

  write(6,'(1X,A,ES12.5)') '        Convergence threshold for energy',THRES_E
  !. Multispace (multigrid info )
  !if (MULSPC == 1) then
  !  write(6,'(1X,A,I3)') '        Multispace method in use from space ',IFMULSPC
  !  write(6,*) '        Pattern'
  !  call IWRTMA(IPAT,1,LPAT,1,LPAT)
  !else
  write(6,'(1X,A)') '        No multispace method in use'
  !end if

  write(6,*)
  if (IDIAG == 2) then
    write(6,'(1X,A,ES12.5)') '        Individual second order energy threshold',E_THRE
    write(6,'(1X,A,ES12.5)') '        Individual first order wavefunction threshold',C_THRE
    if (ICLSSEL == 1) then
      write(6,*)
      write(6,'(1X,A)') '         Class selection will be performed :'
      write(6,'(1X,A)') '         ==================================='
      write(6,'(1X,A,ES12.5)') '          Total second order energy threshold',E_CONV
      write(6,'(1X,A,ES12.5)') '          Total first order wavefunction threshold',C_CONV
    else
      write(6,'(1X,A)') '            No class selection in iterative procedure'
    end if
  end if

  if (IPERT /= 0) then
    write(6,'(1X,A)') '     Perturbation calculation'
    write(6,'(1X,A)') '     ======================='
    write(6,*) '        Root Choosen as zero order state ',IRFROOT
    write(6,*) '        Root used for zero order operator ',IH0ROOT
    !OLD if (MPORENP == 1) then
    !OLD   write(6,*) '        Moller Plesset partitioning'
    !OLD else if (MPORENP == 2) then
    !OLD   write(6,*) '        Epstein-Nesbet partitioning'
    !OLD else if  (MPORENP == 0) then
    !OLD   write(6,*) '        One-body Hamiltonian readin'
    !OLD end if
    if (IE0AVEX == 1) then
      write(6,*) '        Expectation value of H0 used as zero order energy'
    else if (IE0AVEX == 2) then
      write(6,*) '        Exact energy of reference used as zero order energy'
    end if
    write(6,*) '        Correction vectors obtained through  order ',NPERT
    if (IH0SPC == 0) then
      write(6,*) '        No restrictions on perturbation interactions'
    else
      write(6,*) '        Perturbation restricted to interactions in subspaces'
    end if

    if (IH0SPC /= 0) then
      write(6,*)
      write(6,*) '        Number of perturbation subspaces ',NPTSPC
      write(6,*)
      write(6,*) '        ========================'
      write(6,*) '        Perturbation subspaces :'
      write(6,*) '        ========================'
      do JPTSPC=1,NPTSPC
      !OLD write(6,'(A)') ' ======================================================'
        write(6,'(A)')
        write(6,'(7X,A)') '         Min. occ    Max. occ'
        write(6,'(7X,A)') '         ========    ========'
        do IGAS=1,NGAS
          write(6,'(7X,A,I2,3X,I3,9X,I3)') '   GAS',IGAS,IOCPTSPC(1,IGAS,JPTSPC),IOCPTSPC(2,IGAS,JPTSPC)
        end do
      end do

      write(6,*)
      write(6,'(7X,A)') ' ======================================='
      write(6,'(7X,A)') ' Approximate Hamiltonian in CI subspaces'
      write(6,'(7X,A)') ' ======================================='
      write(6,'(7X,A)')
      write(6,'(7X,A)') '    Subspace          H(apr)'
      write(6,'(7X,A)') '  ============================='
      write(6,'(7X,A)')
      do JPTSPC=1,NPTSPC
        if (IH0INSPC(JPTSPC) == 1) then
          write(6,'(12X,I3,8X,A)') JPTSPC,' Diagonal Fock operator'
        else if (IH0INSPC(JPTSPC) == 2) then
          write(6,'(12X,I3,8X,A)') JPTSPC,' Epstein-Nesbet operator'
        else if (IH0INSPC(JPTSPC) == 3) then
          write(6,'(12X,I3,8X,A)') JPTSPC,' Nondiagonal Fock operator'
        else if (IH0INSPC(JPTSPC) == 4) then
          write(6,'(12X,I3,8X,A)') JPTSPC,' Complete Hamiltonian'
        else if (IH0INSPC(JPTSPC) == 5) then
          write(6,'(12X,I3,8X,A)') JPTSPC,' Mix of Fock and Exact operator'
        end if
      end do
      if (ISETKW(61) > 0) then
        write(6,*)
        write(6,'(7X,A)') ' Orbital subspaces where exact Hamiltonian is used :'
        write(6,'(7X,A)') '===================================================='
        write(6,*)
        write(6,'(10X,10(2X,I3))') (IH0EXSPC(I),I=1,NH0EXSPC)
        write(6,*)
      end if

    end if
  end if

  if (NPROP == 0) then
    write(6,*)
    !write(6,*) '     No calculation of properties'
  else
    write(6,'(7X,A,I3)') ' Number of properties to be calculated',NPROP
    write(6,*)
    write(6,'(9X,A)') ' Properties :'
    write(6,'(9X,A)') ' ============'
    do IPROP=1,NPROP
      write(6,'(16X,A)') PROPER(IPROP)
    end do

    !if (IRELAX == 0) then
    write(6,'(7X,A)') ' No use of relaxed densities'
    !else
    !  write(6,'(7X,A)') ' Relaxed densities used for property evaluation'
    !  write(6,'(7X,A)') ' (implemented only for pert)'
    !end if
  end if

  if ((IEXTKOP == 0) .and. (IPTEKT == 0)) then
    !write(6,'(5X,A)') " No extended Koopmans' calculations"
  else if (IEXTKOP /= 0) then
    write(6,'(5X,A)') " Extended Koopmans' calculations"
  else if (IPTEKT /= 0) then
    write(6,'(5X,A)') ' Perturbation expansion of EKT equations'
  end if

  if (IPTFOCK == 1) then
    write(6,*) ' Perturbation expansion of Fock matrix'
  else
    !write(6,*) 'No  Perturbation expansion of Fock matrix'
  end if

  if (ITRAPRP == 0) then
    !write(6,*)
    !write(6,'(5X,A)') ' No transition properties will be calculated'
  else
    write(6,*)
    write(6,'(5X,A)') ' Transition properties will be calculated'
    write(6,*) ' Symmetry of additional states :',IEXCSYM
    write(6,*) ' Number   of additional states :',NEXCSTATE
    write(6,*)
  end if

  if (IRESPONS /= 0) then
    write(6,*)
    write(6,*) '**************************'
    write(6,*) '*  Response Calculation  *'
    write(6,*) '**************************'
    write(6,*)
    write(6,*) ' CI-Response will be called after each CI calculation'
    write(6,*) ' Root used for response calculations (RFROOT) ',IRFROOT
    write(6,*)
    !write(6,*) ' Number of A-operators : ', N_AVE_OP
    write(6,*) ' Labels of A-operators'
    write(6,*) ' ====================='
    write(6,*)
    do IAVE=1,N_AVE_OP
      write(6,'(1X, 6X,A)') AVE_OP(IAVE)
    end do
    write(6,*)
    !write(6,*) ' Number of response calculations ', NRESP
    write(6,*) ' Perturbations :'
    write(6,*) ' ==============='
    write(6,*)
    write(6,*) ' Calc  Op1    Op2    Mxord1     Mxord2    Freq'
    do IRESP=1,NRESP
      write(6,'(1X,I2,2X,A,A,3X,I4,3X,I4,2X,F12.7)') IRESP,RESP_OP(1,IRESP),RESP_OP(2,IRESP),MAXORD_OP(1,IRESP), &
                                                     MAXORD_OP(2,IRESP),RESP_W(IRESP)
    end do
  end if

  !if (NOMOFL == 0) then
  write(6,*)
  write(6,'(7X,A)') ' Final orbitals :'
  write(6,'(7X,A)') ' ================'
  write(6,*)
  if (IFINMO == 1) then
    write(6,'(10X,A)') ' Natural orbitals'
  else if (IFINMO == 2) then
    write(6,'(10X,A)') ' Canonical orbitals'
  else if (IFINMO == 3) then
    write(6,'(10X,A)') ' Pseudo-natural orbitals'
    write(6,'(10X,A)') ' (Density matrix diagonalized in orbital subspaces )'
  else if (IFINMO == 4) then
    write(6,'(10X,A)') ' Pseudo-canonical orbitals'
    write(6,'(10X,A)') ' (FI+FA  diagonalized in orbital subspaces )'
  else if (IFINMO == 5) then
    write(6,'(10X,A)') ' Pseudo-natural-canonical orbitals (sic)'
    write(6,'(10X,A)') ' (Pseudo natural orbitals are first obtained'
    write(6,'(10X,A)') '  by diagonalizing density matrix in orbital subpspaces.'
    write(6,'(10X,A)') '  FI+FA is transformed to this basis, and the transformed'
    write(6,'(10X,A)') '  matrix is block diagonalized)'
    write(6,*)
    write(6,'(10X,A)') ' Orbital spaces in which transformed FIFA is diagonalized'
    write(6,'(10X,A)') ' ========================================================'
    do IPSSPC=1,NPSSPC
      write(6,'(A,I2,A,10I4,6X,2I6)') '     SPACE',IPSSPC,'          ',(NPSSH(IRREP,IPSSPC),IRREP=1,NIRREP)
    end do
  end if
  !end if
  ! Transformation of CI vectors
  if (ITRACI == 0) then
    !write(6,'(5X,A)') ' No transformation of CI vectors'
  else
    write(6,'(5X,A)') ' CI vectors transformed in each run'
    write(6,'(7X,A,A)') ' Complete or restricted rotations :',ITRACI_CR
    write(6,'(7X,A,A)') ' Type of Final orbitals           :',ITRACI_CN
  end if

  ! Integral Transformations

  !if (ITRA_FI == 1) write(6,*) "      Integrals transformed to final MO's"
  if (ITRA_IN == 1) write(6,*) "      Integrals transformed to initial  MO's"

  ! Print levels

  write(6,*)
  write(6,'(1X,A,ES18.9)') '      Core energy : ',ECORE

  !if (IDMPIN == 1) then
  !  write(6,'(1X,A)')
  !  write(6,*) '      Integrals written in formatted form (ES22.15)'
  !  write(6,*) '      on file 90'
  !end if

  write(6,*) ' IPART before leaving READIN = ',IPART
end if

! ================================================
!  Set ratio beteeen real and integer word length
! ================================================
irat = RtoI

! =============================================
!  Find largest unused vector for use in Lucia
! =============================================

! Allocate memory for lucia
! Set up string info

call lucia()

end subroutine Lucia_Ini
