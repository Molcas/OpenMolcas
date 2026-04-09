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

!#define _DEBUGPRINT_
subroutine Lucia_Ini()

use lucia_data, only: ECORE, ENVIRO, I2ELIMINATED_IN_GAS, I_ELIMINATE_GAS, IADVICE, ICISTR, ICJKAIB, ICMBSPC, IDC, IDIAG, &
                      IELIMINATED_IN_GAS, IGSOCCX, IH0INSPC, IH0SPC, ini_h0, IPART, IPRCIX, IPRDEN, IREFSM, IRESTR, ISIMSYM, &
                      LCMBSPC, LCSBLK, MOCAA, MS2, MULTS, MXINKA, N_2ELIMINATED_GAS, N_ELIMINATED_GAS, NACTEL, NCISPC, NCMBSPC, &
                      NGAS, NGSSH, NIRREP, NOINT, NPTSPC, NROOT, NSMOB, PSSIGN, Sigma_on_disk
use spinfo, only: DoComb, I2ELIMINATED_IN_GAS_MOLCAS, I_ELIMINATE_GAS_MOLCAS, IELIMINATED_IN_GAS_MOLCAS, IGSOCCX_MOLCAS, &
                  IPRCI_MOLCAS, ISPEED, ISPIN_MOLCAS, ITMAX_MOLCAS, LSYM_MOLCAS, MS2_MOLCAS, N_2ELIMINATED_GAS_MOLCAS, &
                  N_ELIMINATED_GAS_MOLCAS, NACTEL_MOLCAS, NGAS_MOLCAS, NGSSH_MOLCAS, NROOTS_MOLCAS, NSYM_MOLCAS, POTNUC_MOLCAS
#ifdef _DEBUGPRINT_
use lucia_data, only: NOCSF
use spinfo, only: IEXPAND_MOLCAS, INOCALC_MOLCAS, IPT2_MOLCAS, ISAVE_EXP_MOLCAS, THRE_MOLCAS
#endif
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), parameter :: MXPKW = 125
integer(kind=iwp) :: I, IDENSI, IDOPERT, IEXPERT, IRREP, isetkw(MXPKW), IUSED, MXCIV, NMISS, NWARN
real(kind=wp) :: ECORE_ENV
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ICLSSEL, IEXPAND, IFINMO, INOCALC, IRST2, ISKIPEI, ISAVE_EXP
real(kind=wp) :: PLSIGN, THRES_E
character(len=4) :: ITRACI_CN, ITRACI_CR
#endif

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
NWARN = 0
#ifdef _DEBUGPRINT_
ICLSSEL = 0
IFINMO = 0
IRST2 = 0
ISKIPEI = 0
#endif
! No cc as default
! Start out with normal integrals
!I_USE_SIMTRH = 0
! Initialize flag to be used by densi_master to check whether the
! sigma-vector is on disk.
Sigma_on_disk = .false.

! ===================
!  Initialize isetkw
! ===================
isetkw(:) = 0

! =========================
!  Point group of orbitals
! =========================
!PNTGRP = 1

! ==============================
!  Number of irreps of orbitals
! ==============================
nirrep = nsym_molcas
nsmob = nsym_molcas

! ============================
!  Number of active electrons
! ============================
nactel = nactel_Molcas

! =================
!  Inactive shells
! =================
!ninash(1:nirrep) = nish_molcas(1:nirrep)

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
!iroot(1:nroot) = iroot_molcas(1:nroot)
!iroot(1:nroot) = 0

! ============
!  Enviroment
! ============
enviro(1:6) = 'RASSCF'

! ====================
! 27: Ms combinations
! ====================
if ((MULTS == 1) .and. (ISPEED(1) == 1)) then
  ISETKW(27) = 1
  PSSIGN = One
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
  ngssh(irrep,1:ngas) = ngssh_molcas(1:ngas,irrep)
end do

! ==================================================
!  Generalized active space occupation restrictions
! ==================================================
ncispc = 1
igsoccx(1:ngas,:,1) = igsoccx_molcas(1:ngas,:)

#ifdef _DEBUGPRINT_
! ==========================
!  Energy convergence of CI
! ==========================
thres_e = thre_molcas
#endif

! ==========================================
!  Sequence of calculations to wavefunction
! ==========================================
ncmbspc = ncispc

#ifdef _DEBUGPRINT_
! =======================================================
!  No calculation, save CI vector info, expand CI vector
! =======================================================
INOCALC = INOCALC_MOLCAS
ISAVE_EXP = ISAVE_EXP_MOLCAS
IEXPAND = IEXPAND_MOLCAS
#endif

! ================================================
!  Length of smallest block og C an Sigma vectors
! ================================================
lcsblk = 3000000

! =================================
!  Calculation of density matrices
! =================================
idensi = 2

! ============================================
! 87 : Allow the sigma routine to take advice
! ============================================
!
IADVICE = 1

! ============================================================
!  Transform CI vectors to alternative orbital representation
! ============================================================
#ifdef _DEBUGPRINT_
itraci_cr = 'REST'
if (ipt2_molcas == 1) then
  itraci_cn = 'CANO'
else
  itraci_cn = 'NATU'
end if
#endif

! ==============================
!  HEXS - Highly excited states
!  DEXS - Doubly excited states
! ==============================
I_ELIMINATE_GAS = I_ELIMINATE_GAS_MOLCAS
N_ELIMINATED_GAS = N_ELIMINATED_GAS_MOLCAS
N_2ELIMINATED_GAS = N_2ELIMINATED_GAS_MOLCAS
!write(u6,*) I_ELIMINATE_GAS_MOLCAS,N_ELIMINATED_GAS_MOLCAS
IELIMINATED_IN_GAS(1:N_ELIMINATED_GAS) = IELIMINATED_IN_GAS_MOLCAS(1:N_ELIMINATED_GAS)
I2ELIMINATED_IN_GAS(1:N_2ELIMINATED_GAS) = I2ELIMINATED_IN_GAS_MOLCAS(1:N_2ELIMINATED_GAS)

! =============
!  Printlevels
! =============
iprcix = iprci_molcas-2
iprden = (iprci_molcas-2)/2

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

! 9 : RAS 1 orbitals

! 10 : RAS 2 orbitals

! 11 : RAS 3 orbitals

! 13 : Secondary space

! 14 : occupation restrictions for Reference space

! 15 : Selection of active configurations

! Standard is no selection

! 20 : Diagonalization routine

! Standard is currently MICDV*
IDIAG = 1

! 21 : Explicit Hamiltonian

! Default is no explicit Hamiltonian

! 22 : Largest allowed number of CI iterations per root

! Default is 20 (not active I expect)

! 23 : Restart option

! Default is no explicit Hamiltonian
IRESTR = 0

! 26 : DELETEd shells

! If CAS + Active have been set or RAS + Ras3 have been set,
! obtain for MOLCAS Interface from number of basis functions
!NDELSH(1:NIRREP) = 0

! 27 : Ms combinations

if (ISETKW(27) == 0) then
  PSSIGN = Zero
  ISETKW(27) = 2
else if (MS2 /= 0) then
  write(u6,*) ' Spin combinations only allowed with MS2 = 0'
  write(u6,*) ' Your value of MS2 = ',MS2,' differs from zero'
  write(u6,*) ' LUCIA will neglect your nice suggestion'
  write(u6,*) ' to use spin combinations'
  PSSIGN = Zero
  ISETKW(27) = 2
end if

! 28 : Ml combinations

#ifdef _DEBUGPRINT_
PLSIGN = Zero
#endif

if (PSSIGN == Zero) then
  DoComb = .false.
  IDC = 1
else
  DoComb = .true.
  IDC = 2
end if

!write(u6,*) ' TEST readin IDC = ',IDC

! ======================================================================
! 44 : Use Minimal operation count method for alpha-alpha and beta-beta
! ======================================================================
! TEST OF PERFORMANCE
MOCAA = ISPEED(3)
!MOCAB = ISPEED(4)

! 33 : Number of Ci vectors in subspace

if ((IDIAG == 2) .and. (MXCIV > 2)) then
  MXCIV = 2
  NWARN = NWARN+1
  write(u6,*) ' Warning : You have specified TERACI'
  write(u6,*) '           I allow myself to set MXCIV = 2'
  write(u6,*)
  write(u6,*) '                   Best Wishes'
  write(u6,*) '                      Lucia'
end if

! 34 : CI storage mode

! 35 : Employ CSF expansion ?

! Default is no (only possibility at the moment)
! CSF expansion must only be used when two vectors are stored in CORE

! 37 : Avoid any readin of integrals (useful for obtaining
!      size of CI expansion etc.)

NOINT = 0

! 38 : Dump integrals in formatted form : Default is No

!IDMPIN = 0

! 40 : Use CJKAKB intermediate matrices in alpha-beta loop,
!      Default is  YES !!!!!

ICJKAIB = 1

!  41 : Initial CI in reference space, default is : No

!  42 : Restart with CI in reference space

! Core energy : Default is 0 / MOLCAS : Value read in !

ECORE = Zero

! Use perturbation theory for zero order space . def is no !

if (ISETKW(47) == 0) ISETKW(47) = 2

! 48 : Approximate Hamiltonian in reference space : NO !!

if (ISETKW(48) == 0) ISETKW(48) = 2

! 49 : Approximate Hamiltonian in zero order space : NO !!

if (ISETKW(49) == 0) ISETKW(49) = 2

! 50 : GAS shells must be defined

! 52 : Combination of gasspaces : Default is just to take each  space
!      By itself

if (ISETKW(52) == 0) then
  NCMBSPC = NCISPC
  LCMBSPC(1:NCISPC) = 1
  ICMBSPC(1,1:NCISPC) = [(i,i=1,NCISPC)]
  ISETKW(52) = 2
end if

! 53 : Convergence threshold for CI

! 54 : General sequencer : default is just straight sequence
!      of CI with default number of iterations

! 55 : EKT calculation : Default is no

if (ISETKW(55) == 0) ISETKW(55) = 2

! 56 : Default Machine : Good old BIM machine

! 57 : Allow first order correction to be saved on DISC
!     (For vector free calculations)
!     Default is : NO !!
if (ISETKW(57) == 0) ISETKW(57) = 2

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
  if (IH0SPC /= 0) IH0INSPC(1:NPTSPC) = IPART
end if

! 60 : Reference Root, default is NROOT

! Should be less or equal to NROOT
!if (ISETKW(60) == 1) then
!  if (IRFROOT > NROOT) then
!    write(u6,*) ' Reference root (RFROOT) larger'
!    write(u6,*) ' than total number of roots (NROOT)'
!    write(u6,*) ' CHANGE NROOT or RFROOT'
!    NMISS = NMISS+1
!  end if
!end if

if (ISETKW(60) == 0) ISETKW(60) = 2

! 61 : H0EX : Orbital spaces in which exaxt Hamiltonian is used
!      No default

! Is H0EX required (Has H0FRM = 5 been used)
IUSED = 0
if (ISETKW(59) == 1) then
  IUSED = 0
  if (any(IH0INSPC(1:NPTSPC) == 5)) IUSED = 1
end if
!if ((IUSED == 0) .and. (ISETKW(61) == 0)) then
!  ! No exact spaces included and none have been defined !
!  NH0EXSPC = 0
!end if
if ((IUSED == 1) .and. (ISETKW(61) == 0)) then
  ! Needed, but not supplied
  write(u6,*) ' You have specified that zero order operator'
  write(u6,*) ' Include exact Hamilton operator in subspace'
  write(u6,*) ' You should then also supply Keyword H0EX'
  NMISS = NMISS+1
end if

! If perturbation theory will be invoked be sure that the
! form of perturbation theory has been specified through
! KEYWORD PERTU (number 47 as you maybe know)
IDOPERT = 0
!do JCMBSPC=1,NCMBSPC
!  do JSEQCI=1,NSEQCI(JCMBSPC)
!    if (ISEQCI(JSEQCI,JCMBSPC) == -5) IDOPERT = 1
!  end do
!end do
if (itmax_molcas == -5) IDOPERT = 1

if (IDOPERT == 1) then
  write(u6,*) ' Perturbation theory will be used'
  write(u6,*) ' Please specify form through PERTU keyword'
  NMISS = NMISS+1
end if

! 62 : Default Handling of degenrences of initial CI vectors
!      Default is : No action

if (ISETKW(62) == 0) ISETKW(62) = 2

! 63 : Use F + Lambda(H-F) as operator instead of H
!      Default is : No i.e Lambda = 1

if (ISETKW(63) == 0) ISETKW(63) = 2

! 64 : Smallest block in batch of C and sigma
!      Default is zero

! 66 : NO MO file : Default is no access to MO-AO file

! 68 : Type of natural orbitals, default is natural orbitals

if (ISETKW(68) == 0) then
  ISETKW(68) = 2
# ifdef _DEBUGPRINT_
  IFINMO = 1
# endif
end if

! 69 : Default Threshold for individual energy correction = 0.0

if (ISETKW(69) == 0) ISETKW(69) = 2

! 70 : Default Threshold for wave individual function corrections = 0.0

if (ISETKW(70) == 0) ISETKW(70) = 2

! 71 : Default Threshold for total energy corrections = 0.0

if (ISETKW(71) == 0) ISETKW(71) = 2

! 72 : Default Threshold for total wave function correction = 0.0

if (ISETKW(72) == 0) ISETKW(72) = 2

! 73 : Perform Class selection : Default if Yes if TERACI is used

if (ISETKW(73) == 0) then
# ifdef _DEBUGPRINT_
  if (IDIAG == 1) then
    ICLSSEL = 0
  else if (IDIAG == 2) then
    ICLSSEL = 1
  end if
# endif
  ISETKW(73) = 2
end if

! 74 : Calculation of density matrices : Default is
!      calculation of one-body density

! If IDENSI was set to zero and properties were requested
! overwrite input to obtain 1-el matrix
if ((IDENSI == 0) .and. (ISETKW(80) == 1)) then
  write(u6,*) ' You have specified calculation of'
  write(u6,*) ' one-electron properties, and this'
  write(u6,*) ' requires the calculation of the'
  write(u6,*) ' one-electron density.'
  write(u6,*)
  write(u6,*) ' You have, however, specified IDENSI=0'
  write(u6,*) ' corresponding  to no densities'
  write(u6,*)
  write(u6,*) ' I will allow myself to modify your'
  write(u6,*) ' input to allow calculation of the'
  write(u6,*) ' one-electron densities, so property'
  write(u6,*) ' calculation can proceed as planned'
  write(u6,*)
  write(u6,*) ' Lucia'
  ! and do it
  IDENSI = 1
end if
! If CC is performed, one- and two- particle densities are
! used in current simple-minded implementation.
! Two-electron density also needed for MCSCF

! 75 : Perturbation expansion of EKT, default is no

if (ISETKW(75) == 0) then
  ISETKW(75) = 2
end if

! 76 : Root for zero order operator, default is NROOT

! Should be less or equal to NROOT
!if (ISETKW(76) == 1) then
!  if (IH0ROOT > NROOT) then
!    write(u6,*) ' Zero order operator root (H0ROOT) larger'
!    write(u6,*) ' than total number of roots (NROOT)'
!    write(u6,*) ' CHANGE NROOT or H0ROOT'
!    NMISS = NMISS+1
!  end if
!end if
if (ISETKW(76) == 0) ISETKW(76) = 2

! 77 : NO restart from previous vectors in calculation 2
!      Deafault is NO NO, ie. restart in calc 2

if (ISETKW(77) == 0) then
# ifdef _DEBUGPRINT_
  IRST2 = 1
# endif
  ISETKW(77) = 2
end if

! 78 : skip initial energy evaluations - if possible

if (ISETKW(78) == 0) then
# ifdef _DEBUGPRINT_
  ISKIPEI = 1
# endif
  ISETKW(78) = 2
end if

! 79 : Symmetry of x,y,z - needed for property calculations

if (ISETKW(79) == 0) then
  ! Problematic if Properties should be calculated
  if ((ISETKW(80) == 1) .or. (ISETKW(81) == 1) .or. (ISETKW(82) == 1)) then
    write(u6,*) ' Symmetry of X,Y,Z has not been given'
    write(u6,*) ' You have to specify this for property calc'
    write(u6,*) ' Please add KEYWORD XYZSYM'
    NMISS = NMISS+1
    ISETKW(79) = -1
  else
    ISETKW(79) = 2
  end if
end if

! 80 : Property calculation, default is no

if (ISETKW(80) == 0) ISETKW(80) = 2

! 81 : Transition properties, default is no

if (ISETKW(81) == 0) ISETKW(81) = 2

! 82 : Response properties, default is no

if (ISETKW(82) == 0) ISETKW(82) = 2
! Properties should be defined if transition properties are
! invoked
!if ((ITRAPRP /= 0) .and. (NPROP == 0)) then
!  write(u6,*) ' You have specified transition property calculation'
!  write(u6,*) ' (keyword TRAPRP) but no property labels have been supplied'
!  write(u6,*) '(Keyword PROPER). Transition densities will be obtained'
!end if

! 83 : Max number of iterations in linear equations

if (ISETKW(83) == 0) ISETKW(83) = 2

! 85 : Root homing, default is no

if (ISETKW(85) == 0) ISETKW(85) = 2

! 90 : Perturbation expansion of Fock matrix : default is no

if (ISETKW(90) == 0) ISETKW(90) = 2

! 91 : Print final CI vectors : default is no

if (ISETKW(91) == 0) ISETKW(91) = 2

! 93 : End Calculation with integral transformation

if (ISETKW(93) == 0) ISETKW(93) = 2

! 94 : Initialize Calculation with integral transformation

if (ISETKW(94) == 0) ISETKW(94) = 2
! Requires access to MO-AO file
!if ((ITRA_IN == 1) .and. (NOMOFL == 1)) then
!  write(u6,*) ' Integral transformation required,'
!  write(u6,*) ' but no mo-ao file accessible'
!  write(u6,*) ' REMOVE KEWORD NOMOFL'
!  ISETKW(94) = -1
!  NERROR = NERROR+1
!end if

! 95 : Multispace optimization in each run, default is no

if (ISETKW(95) == 0) ISETKW(95) = 2

! Expert mode (neglect mistyped keywords) : default is no expert

if (ISETKW(97) == 0) then
  IEXPERT = 0
  ISETKW(97) = 2
end if

! Number of roots to be converged : default is total number of roots

if (ISETKW(98) == 0) ISETKW(98) = 2

! 100 : Do quantum dot calculation, default is no

!if (ISETKW(100) == 0) then
!  IDODQ = 0
!  ISETKW(100) = 2
!end if

! 101: Restrict MS2 at some intermediate level : default is no way

if (ISETKW(101) == 0) ISETKW(101) = 2

! 103 : Treat all TT blocks with given types simultaneously : Default is no

ISIMSYM = 0

! Thresholds only active in connection with IDIAG = 2,
! Check and maybe issue a warning
if (IDIAG == 2) then
  ! Check to ensure that zero or two thresholds were  set,
  if (ISETKW(69) /= ISETKW(70)) then
    write(u6,*) ' Only a single threshold (E_THRE or C_THRE)'
    write(u6,*) ' on individual determinants given.'
    write(u6,*) ' One of the thresholds vanishes therefore and'
    write(u6,*) ' all determinants will therefore be included'
    write(u6,*)
    write(u6,*) '                   Warns'
    write(u6,*)
    write(u6,*) '                   LUCIA'
  end if
else
  ! Good old diagonalization, thresholds not active
  if ((ISETKW(69) == 1) .or. (ISETKW(70) == 1)) then
    write(u6,*) ' Thresholds on selection of individual coefficients'
    write(u6,*) ' are only active in connection with keyword TERACI'
    write(u6,*)
    write(u6,*) '                   Warns'
    write(u6,*)
    write(u6,*) '                   LUCIA'
  end if
end if

if (NMISS /= 0) then
  write(u6,'(1X,A,I9)') ' Number of missing required keyword ',NMISS
  write(u6,*) ' You have wounded me I give up'
  write(u6,*)
  write(u6,*)
  write(u6,*)
  write(u6,*) '     An expert is a man who has made all the mistakes,'
  write(u6,*) '     which can be made, in a very narrow field'
  write(u6,*)
  write(u6,*) '                                      Niels Bohr'
  if (IEXPERT == 0) then
    !stop
    call SYSABENDMSG('lucia_util/lucia_ini','Input error','')
  else
    write(u6,*) ' Processing continues (EXPERT mode)'
  end if
end if
! Open one-electron file to obtain core energy and
! Number of MO's and AO's
if ((NOINT == 0) .and. (ENVIRO(1:4) /= 'NONE')) then
  ecore_env = potnuc_molcas
  ECORE = ECORE_ENV
else
  write(u6,*) ' GETOBS and CHK_ORBDIM not called'
  ECORE = Zero
end if
! Check number of orbitals and insert occupations for ALL/REST

#ifdef _DEBUGPRINT_
!***********************************************************
!                                                          *
! Part 3 : Print input                                     *
!                                                          *
!***********************************************************

write(u6,*)
write(u6,*) '*************************************'
write(u6,*) '*  Symmetry and spin of CI vectors  *'
write(u6,*) '*************************************'
write(u6,*)
! Point group
write(u6,'(1X,A)') '     Point group ............ D2H'
! Spatial symmetry
write(u6,'(1X,A,I1)') '     Spatial symmetry ....... ',IREFSM
! Spin
write(u6,'(1X,A,I2)') '     2 times spinprojection  ',MS2
! Number of active electrons
write(u6,'(1X,A,I2)') '     Active electrons .....  ',NACTEL
write(u6,*)
write(u6,*) '*********************************************'
write(u6,*) '*  Shell spaces and occupation constraints  *'
write(u6,*) '*********************************************'
write(u6,*)

!if (XLAMBDA /= One) then
!  write(u6,*)
!  write(u6,'(A,F13.8)') ' Modified operator H(l) = l*F + l*(H-F) used with l =',XLAMBDA
!  !if (IUSEH0P == 0) then
!  write(u6,'(A)') ' Zero-order operator without projection used'
!  !else
!  !  write(u6,'(A)') ' Zero-order operator with projection used'
!  !end if
!  if (IRESTR == 0) then
!    write(u6,*) ' Notice : This madness starts  in second calculation'
!  else
!    write(u6,*) ' You have specified a calculation with modified'
!    write(u6,*) ' Hamiltonian (the LAMBDA option) and RESTART'
!    write(u6,*) ' so this is what I will do'
!    write(u6,*)
!    write(u6,*) '   1:) Perform CI in space 1 to obtain Hamiltonian'
!    write(u6,*) '       (no RESTART in this space)'
!    write(u6,*) '   2:) CI calculation in space 2  with'
!    write(u6,*) '       modified Hamiltonian and RESTART from LU21'
!    write(u6,*) ' Space 2 should therefore correspond to the'
!    write(u6,*) ' restarted calculation'
!  end if
!end if

write(u6,*)
write(u6,*) '***********'
write(u6,*) '*  Roots  *'
write(u6,*) '***********'
write(u6,*)
write(u6,'(1X,A,I3)') '     Number of roots to be included  ',NROOT
!write(u6,'(1X,A,(20I3))') '     Roots to be obtained ',(IROOT(I),I=1,NROOT)
!write(u6,'(1X,A,I3)') '     Number of roots to be converged ',NCNV_RT

write(u6,*)
write(u6,*) '**************************'
write(u6,*) '*  Run time definitions  *'
write(u6,*) '**************************'
write(u6,*)
! Program environment
write(u6,'(A,A6)') '      Program environment... ',ENVIRO

if (NOINT == 1) then
  ! Integral import
  write(u6,'(1X,A)') '     No integrals will be read in'
else if (NOINT == 0) then
  ! Integral storage
  write(u6,'(1X,A)') '     All integrals stored in core'
end if
write(u6,*)
! END IF for NOINT
! CSF or SD expansion
if (NOCSF == 0) then
  write(u6,'(1X,A)') "     CI optimization performed with CSF's"
else
  write(u6,'(1X,A)') "     CI optimization performed with SD's"
end if
! Ms,Ml combinations
if (ISETKW(27) == 1) write(u6,'(1X,A,F8.3)') '     Spin combinations used with sign ',PSSIGN
if (ISETKW(28) == 1) write(u6,'(1X,A,F8.3)') '     ML   combinations used with sign ',PLSIGN
! Initial approximation to vectors
write(u6,*)
!if ((IRESTR == 1) .and. (IRESTRF == 0)) then
if (IRESTR == 1) then
  write(u6,'(1X,A)') '     Restarted calculation'
!else if (IRESTRF == 1) then
!  write(u6,'(1X,A)') '     Restarted calculation from REFERENCE space expansion'
else
  !if (MXP1 /= 0) then
  !  write(u6,'(1X,A)') '     Initial vectors obtained from explicit Hamiltonian'
  !else if (MXP1 == 0) then
  write(u6,'(1X,A)') '     Initial vectors obtained from diagonal'
  !end if
end if
! Handling of degenerencies of initial vectors
!if (INIDEG == 1) then
!  write(u6,'(1X,A)') '     Symmetric combination of degenerate initial vectors'
!else if (INIDEG == -1) then
!  write(u6,'(1X,A)') '     Antiymmetric combination of degenerate initial vectors'
!else if (INIDEG == 0) then
!  write(u6,'(1X,A)') '     No combination of degenerate initial vectors'
!end if
! No calculation, save CI vector information and expand CI vector
if (INOCALC == 1) write(u6,'(1X,A)') ' No calculation will be performed'
if (ISAVE_EXP == 1) write(u6,'(1X,A)') ' Save CI vector information'
if (IEXPAND == 1) write(u6,'(1X,A)') ' Expand shorter CI vector in longer one'
! Ms,Ml combinations
!if (ISETKW(27) == 1) write(u6,'(1X,A,F8.3)') '     Spin combinations used with sign ',PSSIGN
!if (ISETKW(28) == 1) write(u6,'(1X,A,F8.3)') '     ML   combinations used with sign ',PLSIGN
! CI storage mode
write(u6,*)

write(u6,*) '     3 symmetry blocks will be held in core'

if (LCSBLK /= 0) write(u6,'(A,I10)') '      Smallest allowed size of sigma- and C-batch ',LCSBLK
write(u6,'(1X,A,I4)') '     Dimension of block of resolution strings ',MXINKA
!if (IUSE_PH == 1) then
!  WRITE(u6,'(1X,A)') '     Particle-hole separation used'
!else
write(u6,'(1X,A)') '      Particle-hole separation not used'
!end if

if (IADVICE == 1) write(u6,'(1X,A)') '     Advice routine call to optimize sigma generation'

!if (IUSE_PA == 1) then
!  write(u6,'(1X,A)') '     Strings divided into active and passive parts'
!else
write(u6,'(1X,A)') '     Strings not divided into active and passive parts'
!end if
!if (ISIMSYM == 1) write(u6,'(1X,A)') '     ALl TTS blocks with given types treated in sigma'
!if (IUSE_HW == 1) write(u6,*) ' Hardwired routines in use'

write(u6,*)
if (IDENSI == 0) then
  write(u6,'(1X,A)') '     No calculation of density matrices'
else if (IDENSI == 1) then
  write(u6,'(1X,A)') '     One-body density matrix calculated'
else if (IDENSI == 2) then
  write(u6,'(1X,A)') '     One- and two-body density matrices  calculated'
end if
write(u6,*)
!if (MOCAA /= 0) write(u6,*) '     MOC method used for alpha-alpha+beta-beta loop'
!if (MOCAB /= 0) write(u6,*) '     MOC method used for alpha-beta loop'

! Diagonalization information
write(u6,'(1X,A)') '     CI diagonalization :'
write(u6,'(1X,A)') '     ===================='
! Subspace Hamiltinian
!if (MXP1+MXP2+MXQ == 0) then
write(u6,'(1X,A)') '        No subspace Hamiltonian'
!else
!  write(u6,'(1X,A,3I4)') '        Dimensions of subspace Hamiltonian ',MXP1,MXP2,MXQ
!end if
! Diagonalizer
if ((IDIAG == 1) .and. (ICISTR == 1)) then
  write(u6,'(1X,A)') '        Diagonalizer : MINDV4'
else if ((IDIAG == 1) .and. (ICISTR >= 2)) then
  write(u6,'(1X,A)') '        Diagonalizer : MICDV*'
else if (IDIAG == 2) then
  write(u6,'(1X,A)') '        Diagonalizer : PICO*'
end if
write(u6,'(1X,A)') '        Simple diagonal used as preconditioner'
! Root homing
!if (IROOTHOMING == 1) then
!  write(u6,'(1X,A)') '        Root homing will be used'
!else
write(u6,'(1X,A)') '        No root homing'
!end if
! No restart in CI calc 2
if (IRST2 == 0) write(u6,'(1X,A)') '        No restart from previous vectors in second calc'
if (ISKIPEI == 1) then
  write(u6,'(1X,A)') '        Initial energy evaluations skipped after first calc'
  write(u6,'(1X,A)') '        (Only active in connection with TERACI)'
end if
! Number of iterations
!write(u6,'(1X,A,I2)') '        Allowed number of iterations    ',MAXIT
! Number of CI vectors in subspace
write(u6,'(1X,A,I2)') '        Allowed Dimension of CI subspace',MXCIV

write(u6,'(1X,A,ES12.5)') '        Convergence threshold for energy',THRES_E
! Multispace (multigrid info)
!if (MULSPC == 1) then
!  write(u6,'(1X,A,I3)') '        Multispace method in use from space ',IFMULSPC
!  write(u6,*) '        Pattern'
!  call IWRTMA(IPAT,1,LPAT,1,LPAT)
!else
write(u6,'(1X,A)') '        No multispace method in use'
!end if

write(u6,*)
if (IDIAG == 2) then
  !write(u6,'(1X,A,ES12.5)') '        Individual second order energy threshold',E_THRE
  !write(u6,'(1X,A,ES12.5)') '        Individual first order wavefunction threshold',C_THRE
  if (ICLSSEL == 1) then
    write(u6,*)
    write(u6,'(1X,A)') '         Class selection will be performed :'
    write(u6,'(1X,A)') '         ==================================='
    !write(u6,'(1X,A,ES12.5)') '          Total second order energy threshold',E_CONV
    !write(u6,'(1X,A,ES12.5)') '          Total first order wavefunction threshold',C_CONV
  else
    write(u6,'(1X,A)') '            No class selection in iterative procedure'
  end if
end if

!if (IPERT /= 0) then
!  write(u6,'(1X,A)') '     Perturbation calculation'
!  write(u6,'(1X,A)') '     ======================='
!  write(u6,*) '        Root Choosen as zero order state ',IRFROOT
!  !write(u6,*) '        Root used for zero order operator ',IH0ROOT
!  !OLD if (MPORENP == 1) then
!  !OLD   write(u6,*) '        Moller Plesset partitioning'
!  !OLD else if (MPORENP == 2) then
!  !OLD   write(u6,*) '        Epstein-Nesbet partitioning'
!  !OLD else if  (MPORENP == 0) then
!  !OLD   write(u6,*) '        One-body Hamiltonian readin'
!  !OLD end if
!  !if (IE0AVEX == 1) then
!  !  write(u6,*) '        Expectation value of H0 used as zero order energy'
!  !else if (IE0AVEX == 2) then
!  !  write(u6,*) '        Exact energy of reference used as zero order energy'
!  !end if
!  write(u6,*) '        Correction vectors obtained through  order ',NPERT
!  if (IH0SPC == 0) then
!    write(u6,*) '        No restrictions on perturbation interactions'
!  else
!    write(u6,*) '        Perturbation restricted to interactions in subspaces'
!  end if
!
!  if (IH0SPC /= 0) then
!    write(u6,*)
!    write(u6,*) '        Number of perturbation subspaces ',NPTSPC
!    write(u6,*)
!    write(u6,*) '        ========================'
!    write(u6,*) '        Perturbation subspaces :'
!    write(u6,*) '        ========================'
!    do JPTSPC=1,NPTSPC
!      !OLD write(u6,'(A)') ' ======================================================'
!      write(u6,'(A)')
!      write(u6,'(7X,A)') '         Min. occ    Max. occ'
!      write(u6,'(7X,A)') '         ========    ========'
!      do IGAS=1,NGAS
!        write(u6,'(7X,A,I2,3X,I3,9X,I3)') '   GAS',IGAS,IOCPTSPC(1,IGAS,JPTSPC),IOCPTSPC(2,IGAS,JPTSPC)
!      end do
!    end do

!    write(u6,*)
!    write(u6,'(7X,A)') ' ======================================='
!    write(u6,'(7X,A)') ' Approximate Hamiltonian in CI subspaces'
!    write(u6,'(7X,A)') ' ======================================='
!    write(u6,'(7X,A)')
!    write(u6,'(7X,A)') '    Subspace          H(apr)'
!    write(u6,'(7X,A)') '  ============================='
!    write(u6,'(7X,A)')
!    do JPTSPC=1,NPTSPC
!      if (IH0INSPC(JPTSPC) == 1) then
!        write(u6,'(12X,I3,8X,A)') JPTSPC,' Diagonal Fock operator'
!      else if (IH0INSPC(JPTSPC) == 2) then
!        write(u6,'(12X,I3,8X,A)') JPTSPC,' Epstein-Nesbet operator'
!      else if (IH0INSPC(JPTSPC) == 3) then
!        write(u6,'(12X,I3,8X,A)') JPTSPC,' Nondiagonal Fock operator'
!      else if (IH0INSPC(JPTSPC) == 4) then
!        write(u6,'(12X,I3,8X,A)') JPTSPC,' Complete Hamiltonian'
!      else if (IH0INSPC(JPTSPC) == 5) then
!        write(u6,'(12X,I3,8X,A)') JPTSPC,' Mix of Fock and Exact operator'
!      end if
!    end do
!    !if (ISETKW(61) > 0) then
!    !  write(u6,*)
!    !  write(u6,'(7X,A)') ' Orbital subspaces where exact Hamiltonian is used :'
!    !  write(u6,'(7X,A)') '===================================================='
!    !  write(u6,*)
!    !  write(u6,'(10X,10(2X,I3))') (IH0EXSPC(I),I=1,NH0EXSPC)
!    !  write(u6,*)
!    !end if
!
!  end if
!end if

!if (NPROP == 0) then
write(u6,*)
!write(u6,*) '     No calculation of properties'
!else
!  write(u6,'(7X,A,I3)') ' Number of properties to be calculated',NPROP
!  write(u6,*)
!  write(u6,'(9X,A)') ' Properties :'
!  write(u6,'(9X,A)') ' ============'
!  do IPROP=1,NPROP
!    write(u6,'(16X,A)') PROPER(IPROP)
!  end do
!
!  !if (IRELAX == 0) then
!  write(u6,'(7X,A)') ' No use of relaxed densities'
!  !else
!  !  write(u6,'(7X,A)') ' Relaxed densities used for property evaluation'
!  !  write(u6,'(7X,A)') ' (implemented only for pert)'
!  !end if
!end if

!if ((IEXTKOP == 0) .and. (IPTEKT == 0)) then
!  !write(u6,'(5X,A)') " No extended Koopmans' calculations"
!else if (IEXTKOP /= 0) then
!  write(u6,'(5X,A)') " Extended Koopmans' calculations"
!else if (IPTEKT /= 0) then
!  write(u6,'(5X,A)') ' Perturbation expansion of EKT equations'
!end if

!if (IPTFOCK == 1) then
!  write(u6,*) ' Perturbation expansion of Fock matrix'
!else
!  !write(u6,*) 'No  Perturbation expansion of Fock matrix'
!end if

!if (ITRAPRP == 0) then
!  !write(u6,*)
!  !write(u6,'(5X,A)') ' No transition properties will be calculated'
!else
!  write(u6,*)
!  write(u6,'(5X,A)') ' Transition properties will be calculated'
!  !write(u6,*) ' Symmetry of additional states :',IEXCSYM
!  write(u6,*) ' Number   of additional states :',NEXCSTATE
!  write(u6,*)
!end if

!if (IRESPONS /= 0) then
!  write(u6,*)
!  write(u6,*) '**************************'
!  write(u6,*) '*  Response Calculation  *'
!  write(u6,*) '**************************'
!  write(u6,*)
!  write(u6,*) ' CI-Response will be called after each CI calculation'
!  write(u6,*) ' Root used for response calculations (RFROOT) ',IRFROOT
!  write(u6,*)
!  !write(u6,*) ' Number of A-operators : ', N_AVE_OP
!  !write(u6,*) ' Labels of A-operators'
!  !write(u6,*) ' ====================='
!  !write(u6,*)
!  !do IAVE=1,N_AVE_OP
!  !  write(u6,'(1X, 6X,A)') AVE_OP(IAVE)
!  !end do
!  !write(u6,*)
!  !write(u6,*) ' Number of response calculations ', NRESP
!  write(u6,*) ' Perturbations :'
!  write(u6,*) ' ==============='
!  write(u6,*)
!  write(u6,*) ' Calc  Op1    Op2    Mxord1     Mxord2    Freq'
!  do IRESP=1,NRESP
!    write(u6,'(1X,I2,2X,A,A,3X,I4,3X,I4,2X,F12.7)') IRESP,RESP_OP(1,IRESP),RESP_OP(2,IRESP),MAXORD_OP(1,IRESP), &
!                                                    MAXORD_OP(2,IRESP),RESP_W(IRESP)
!  end do
!end if

!if (NOMOFL == 0) then
write(u6,*)
write(u6,'(7X,A)') ' Final orbitals :'
write(u6,'(7X,A)') ' ================'
write(u6,*)
if (IFINMO == 1) then
  write(u6,'(10X,A)') ' Natural orbitals'
else if (IFINMO == 2) then
  write(u6,'(10X,A)') ' Canonical orbitals'
else if (IFINMO == 3) then
  write(u6,'(10X,A)') ' Pseudo-natural orbitals'
  write(u6,'(10X,A)') ' (Density matrix diagonalized in orbital subspaces)'
else if (IFINMO == 4) then
  write(u6,'(10X,A)') ' Pseudo-canonical orbitals'
  write(u6,'(10X,A)') ' (FI+FA  diagonalized in orbital subspaces)'
else if (IFINMO == 5) then
  write(u6,'(10X,A)') ' Pseudo-natural-canonical orbitals (sic)'
  write(u6,'(10X,A)') ' (Pseudo natural orbitals are first obtained'
  write(u6,'(10X,A)') '  by diagonalizing density matrix in orbital subpspaces.'
  write(u6,'(10X,A)') '  FI+FA is transformed to this basis, and the transformed'
  write(u6,'(10X,A)') '  matrix is block diagonalized)'
  write(u6,*)
  !write(u6,'(10X,A)') ' Orbital spaces in which transformed FIFA is diagonalized'
  !write(u6,'(10X,A)') ' ========================================================'
  !do IPSSPC=1,NPSSPC
  !  write(u6,'(A,I2,A,10I4,6X,2I6)') '     SPACE',IPSSPC,'          ',(NPSSH(IRREP,IPSSPC),IRREP=1,NIRREP)
  !end do
end if
!end if
! Transformation of CI vectors
!if (ITRACI == 0) then
!  write(u6,'(5X,A)') ' No transformation of CI vectors'
!else
write(u6,'(5X,A)') ' CI vectors transformed in each run'
write(u6,'(7X,A,A)') ' Complete or restricted rotations :',ITRACI_CR
write(u6,'(7X,A,A)') ' Type of Final orbitals           :',ITRACI_CN
!end if

! Integral Transformations

!if (ITRA_FI == 1) write(u6,*) "      Integrals transformed to final MO's"
!if (ITRA_IN == 1) write(u6,*) "      Integrals transformed to initial MO's"

! Print levels

write(u6,*)
write(u6,'(1X,A,ES18.9)') '      Core energy : ',ECORE

!if (IDMPIN == 1) then
!  write(u6,'(1X,A)')
!  write(u6,*) '      Integrals written in formatted form (ES22.15)'
!  write(u6,*) '      on file 90'
!end if

write(u6,*) ' IPART before leaving READIN = ',IPART
#endif

! =============================================
!  Find largest unused vector for use in Lucia
! =============================================

! Allocate memory for lucia
! Set up string info

call lucia()

end subroutine Lucia_Ini
