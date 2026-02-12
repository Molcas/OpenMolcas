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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!***********************************************************************

subroutine ChkInp()
!***********************************************************************
!                                                                      *
!     Check the input for obvious errors or violation of limits        *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

#ifdef _DMRG_
use qcmaquis_interface_cfg
#endif
use linalg_mod, only: abort_
use gas_data, only: iDoGAS, nGAS, iGSOCCX, nGSSH
use rasscf_global, only: ITMAX, lRoots, MAXIT, MAXJT, NAC, NIN, nRoots, OutFmt1, OutFmt2, PreThr, ProThr, ThFact, ThrE, ThrEn, &
                         ThrSX, ThrTE, iRoot
use output_ras, only: LF
use general_data, only: NTOT, NACTEL, NHOLE1, NRS1T, NELEC3, NRS3T, NRS2T, NSYM, ISPIN, STSYM, NSEL, NALTER, INVEC, NASH, NBAS, &
                        NDEL, NFRO, NISH, NORB, NRS1, NRS2, NRS3, NSSH, MALTER
use Molcas, only: MxAct, MxBas, MxGAS, MxIna, MxOrb, MxRoot
use RASDim, only: MxCIIt, MxIter, MxSXIt

implicit none
#include "warnings.h"
integer :: ierr, ierr1, ierr2
integer :: i, iSym, iAlter
integer :: iB0, iA0, iC0
integer :: c_orbs_per_GAS(nGAS)

!c_orbs_per_GAS: These are the cumulated number of spin orbitals per GAS space.

!----------------------------------------------------------------------*
! Local print level (if any)
IERR = 0
if (NTOT > mxOrb) then
  write(LF,*)
  write(LF,*) '************ ERROR ******************'
  call WarningMessage(2,'Too many orbitals.')
  write(LF,'(1X,A,I8)') 'Too many orbitals NTOT=',NTOT
  write(LF,'(1X,A,I8)') 'Limit is MXORB=',MXORB
  write(LF,*) '*************************************'
  IERR = 1
end if
if (NAC > mxAct) then
  write(LF,*)
  write(LF,*) '*************** ERROR ***************'
  call WarningMessage(2,'Too many active orbitals.')
  write(LF,'(1X,A,I8)') 'Too many active orbitals NAC=',NAC
  write(LF,'(1X,A,I8)') 'Limit is MXACT=',MXACT
  write(LF,*) '*************************************'
  IERR = 1
end if
if (NIN > mxIna) then
  write(LF,*)
  write(LF,*) '*************** ERROR ***************'
  call WarningMessage(2,'Too many inactive orbitals.')
  write(LF,'(1X,A,I8)') 'Too many inactive orbitals NIN=',NIN
  write(LF,'(1X,A,I8)') 'Limit is MXINA=',MXINA
  write(LF,*) '*************************************'
  IERR = 1
end if
if (NACTEL > 2*NAC) then
  write(LF,*)
  write(LF,*) '********************* ERROR **********************'
  call WarningMessage(2,'Too many active electrons.')
  write(LF,'(1X,A,I6)') 'Too many active electrons NACTEL=',NACTEL
  write(LF,'(1X,A,I6)') 'Cannot be more than 2*Nr of active orbitals=',2*NAC
  write(LF,*) '**************************************************'
  IERR = 1
end if
if (NHOLE1 > 2*NRS1T) then
  write(LF,*)
  write(LF,*) '******************** WARNING *********************'
  call WarningMessage(1,'Too many holes in Ras1.')
  write(LF,'(1X,A,I6)') 'You allow too many holes in Ras1 NHOLE1=',NHOLE1
  write(LF,'(1X,A,I6)') 'Cannot be more than 2*Nr of Ras1 orbitals=',2*NRS1T
  NHOLE1 = 2*NRS1T
  write(LF,'(1X,A,I6)') 'NHOLE1 has been reset to ',NHOLE1
  write(LF,*) '**************************************************'
end if
if (NELEC3 > 2*NRS3T) then
  write(LF,*)
  write(LF,*) '******************** WARNING *********************'
  call WarningMessage(1,'Too many electrons in Ras3.')
  write(LF,'(1X,A,I6)') 'You allow too many electrons in Ras3 NELEC3=',NELEC3
  write(LF,'(1X,A,I6)') 'Cannot be more than 2*Nr of Ras3 orbitals=',2*NRS3T
  NELEC3 = 2*NRS3T
  write(LF,'(1X,A,I6)') 'NELEC3 has been reset to ',NELEC3
  write(LF,*) '**************************************************'
end if
if (.not. iDoGas) then !(DM)
  if (NACTEL > 2*(NRS1T+NRS2T)+NELEC3) then
    write(LF,*)
    write(LF,*) '********************* ERROR **********************'
    call WarningMessage(2,'Too many active electrons.')
    write(LF,'(1X,A,I8)') 'Too many active electrons NACTEL=',NACTEL
    write(LF,'(1X,A,I8)') '(Incompatible with RAS restrictions).'
    write(LF,*) '**************************************************'
    IERR = 1
  end if
  if (NACTEL < 2*NRS1T-NHOLE1) then
    write(LF,*)
    write(LF,*) '********************* ERROR **********************'
    call WarningMessage(2,'Too few active electrons.')
    write(LF,'(1X,A,I8)') 'Too few active electrons NACTEL=',NACTEL
    write(LF,'(1X,A,I8)') '(Incompatible with RAS restrictions).'
    write(LF,*) '**************************************************'
    IERR = 1
  end if
  if ((NHOLE1 == 0) .and. (NRS1T > 0)) then
    write(LF,*)
    write(LF,*) '******************* WARNING *******************'
    call WarningMessage(1,'You allow no holes in Ras1')
    write(LF,*) 'You allow no holes in Ras1.                    '
    write(LF,*) 'This may be deliberate, but may give numerical '
    write(LF,*) 'problems in SXCTL section.'
    write(LF,*) '***********************************************'
  end if
  if ((NELEC3 == 0) .and. (NRS3T > 0)) then
    write(LF,*)
    write(LF,*) '******************* WARNING *******************'
    call WarningMessage(1,'You allow no electrons in Ras3')
    write(LF,*) 'You allow no electrons in Ras3.'
    write(LF,*) 'This may be deliberate, but may give numerical '
    write(LF,*) 'problems in SXCTL section.'
    write(LF,*) '***********************************************'
  end if
! for GAS
else
  if (NGAS > mxGAS) then !(SJS)
    write(LF,*)
    call WarningMessage(2,'GASSCF: too many GAS spaces.')
    write(LF,*) ' ************* ERROR **************'
    write(LF,*) ' GASSCF: Too many GAS spaces. '
    write(LF,*) ' Can be increased up to 16 by changing'
    write(LF,*) ' mxGAS in the Molcas module'
    write(LF,*) ' **********************************'
    call Quit(_RC_INPUT_ERROR_)
  end if
  if (nactel /= igsoccx(ngas,2)) then
    write(lf,*)
    write(lf,*) '**************** ERROR *************************'
    write(lf,*) 'nactel not match occupation'
    write(lf,*) 'nactel=',nactel,'igsoccx:',igsoccx(ngas,2)
    write(lf,*) '************************************************'
  end if

  c_orbs_per_GAS = cumsum(sum(ngssh(:nGAS,:nSym),2)*2)
  if (any(c_orbs_per_GAS < igsoccx(:nGAS,1))) then
    write(lf,*)
    write(lf,*) 'In at least one GAS space, the minimum required '
    write(lf,*) 'particle number by GAS constraints '
    write(lf,*) 'is larger than the particle number '
    write(lf,*) 'allowed by the Pauli principle.'
    write(lf,*)
    call abort_('GASSCF: Pauli forbidden.')
  end if
  ! Conceptionally this should not be a problem, but the code
  ! assumes it to be not the case.
  if (any(c_orbs_per_GAS < igsoccx(:nGAS,2))) then
    write(lf,*)
    write(lf,*) 'In at least one GAS space, the maximum allowed '
    write(lf,*) 'particle number by GAS constraints '
    write(lf,*) 'is larger than the particle number '
    write(lf,*) 'allowed by the Pauli principle.'
    write(lf,*)
    call abort_('GASSCF: Pauli forbidden.')
  end if
end if
if ((NSYM /= 1) .and. (NSYM /= 2) .and. (NSYM /= 4) .and. (NSYM /= 8)) then
  write(LF,*)
  call WarningMessage(2,'Wrong nr of symmetries.')
  write(LF,*) '*************** ERROR ***************'
  write(LF,'(1X,A,I8)') 'Nr of symmetries NSYM=',NSYM
  write(LF,*) ' Only possible values are 1,2,4 or 8.'
  write(LF,*) '*************************************'
  IERR = 1
end if
if (IERR == 1) call Quit(_RC_INPUT_ERROR_)

IERR1 = 0
do iSym=1,nSym
  if (nBas(iSym) < 0) IERR1 = 1
  if (nFro(iSym) < 0) IERR1 = 1
  if (nDel(iSym) < 0) IERR1 = 1
  if (nOrb(iSym) < 0) IERR1 = 1
  if (nISh(iSym) < 0) IERR1 = 1
  if (nASh(iSym) < 0) IERR1 = 1
  if (nSSh(iSym) < 0) IERR1 = 1
  if (nRS1(iSym) < 0) IERR1 = 1
  if (nRS2(iSym) < 0) IERR1 = 1
  if (nRS3(iSym) < 0) IERR1 = 1
end do
IERR2 = 0
do iSym=1,nSym
  if (nBas(iSym) > mxBas) IERR2 = 1
  if (nFro(iSym) > mxBas) IERR2 = 1
  if (nDel(iSym) > mxBas) IERR2 = 1
  if (nOrb(iSym) > mxBas) IERR2 = 1
  if (nISh(iSym) > mxBas) IERR2 = 1
  if (nASh(iSym) > mxBas) IERR2 = 1
  if (nSSh(iSym) > mxBas) IERR2 = 1
  if (nRS1(iSym) > mxBas) IERR2 = 1
  if (nRS2(iSym) > mxBas) IERR2 = 1
  if (nRS3(iSym) > mxBas) IERR2 = 1
end do
if (IERR1+IERR2 > 0) then
  write(LF,*)
  write(LF,*) '****************** ERROR *******************'
  call WarningMessage(2,'Erroneous nr of orbitals.')
  write(LF,*) 'Inappropriate nr of orbitals. One or more of'
  write(LF,*) 'these orbital counts is wrong or too large.'
  if (IERR1 > 0) write(LF,*) ' Negative values.'
  if (IERR2 > 0) write(LF,*) ' Extremely large values.'
  write(LF,'(1X,A,8I4)') '   All orbitals:',(NBAS(ISYM),ISYM=1,NSYM)
  write(LF,'(1X,A,8I4)') '         Frozen:',(NFRO(ISYM),ISYM=1,NSYM)
  write(LF,'(1X,A,8I4)') '       Inactive:',(NISH(ISYM),ISYM=1,NSYM)
  write(LF,'(1X,A,8I4)') '         Active:',(NASH(ISYM),ISYM=1,NSYM)
  write(LF,'(1X,A,8I4)') '          RAS-1:',(NRS1(ISYM),ISYM=1,NSYM)
  write(LF,'(1X,A,8I4)') '          RAS-2:',(NRS2(ISYM),ISYM=1,NSYM)
  write(LF,'(1X,A,8I4)') '          RAS-3:',(NRS3(ISYM),ISYM=1,NSYM)
  write(LF,'(1X,A,8I4)') '      Secondary:',(NSSH(ISYM),ISYM=1,NSYM)
  write(LF,'(1X,A,8I4)') '        Deleted:',(NDEL(ISYM),ISYM=1,NSYM)
  write(LF,'(1X,A,8I4)') 'Basis functions:',(NBAS(ISYM),ISYM=1,NSYM)
  if (IERR1 > 0) write(LF,*) ' Have you used a too small basis set?'
  write(LF,*) '********************************************'
  call Quit(_RC_INPUT_ERROR_)
end if

IERR = 0
IB0 = ISPIN-1
IA0 = (NACTEL-IB0)/2
IC0 = (NAC-IA0-IB0)
if ((2*IA0+IB0) /= NACTEL) IERR = 1
if (IA0 < 0) IERR = 1
if (IB0 < 0) IERR = 1
if (IC0 < 0) IERR = 1
if (IERR == 1) then
  write(LF,*)
  write(LF,*) '******************** ERROR *********************'
  call WarningMessage(2,'No such wave function.')
  write(LF,*) 'The following combined specifications are wrong.'
  write(LF,'(1X,A,8I4)') 'Nr of active electrons NACTEL=',NACTEL
  write(LF,'(1X,A,8I4)') 'Nr of active orbitals     NAC=',NAC
  write(LF,'(1X,A,8I4)') 'Spin degeneracy         ISPIN=',ISPIN
  write(LF,*) 'There can be no such wave function.'
  write(LF,*) '************************************************'
  call Quit(_RC_INPUT_ERROR_)
end if

if (MAXIT > mxIter) then
  write(LF,*)
  write(LF,*) '*************** WARNING ****************'
  call WarningMessage(1,'Too many macro-iterations.')
  write(LF,*) 'Too many macro-iterations requested.'
  MAXIT = MXITER
  write(LF,'(1X,A,I8)') 'Reset to maximum, new MAXIT=',MAXIT
  write(LF,*) '****************************************'
end if
if (MAXJT > (mxCiIt-2)) then
  write(LF,*)
  write(LF,*) '*************** WARNING ****************'
  call WarningMessage(1,'Too many CI-iterations.')
  write(LF,*) 'Too many CI-iterations requested.'
  MAXJT = mxCiIt-2
  write(LF,'(1X,A,I8)') 'Reset to maximum, new MAXJT=',MAXJT
  write(LF,*) '****************************************'
end if
if (ITMAX > MXSXIT) then
  write(LF,*)
  write(LF,*) '*************** WARNING ****************'
  call WarningMessage(1,'Too many SX-iterations.')
  write(LF,*) 'Too many SX-iterations requested.'
  ITMAX = MXSXIT
  write(LF,'(1X,A,I8)') 'Reset to maximum, new ITMAX=',ITMAX
  write(LF,*) '****************************************'
end if

THRE = max(0.0d0,THRE)
THRTE = max(0.0d0,THRTE)
THRSX = max(0.0d0,THRSX)
THREN = max(0.0d0,THREN)
THFACT = max(0.0d0,THFACT)

IERR = 0
if (NROOTS > mxRoot) IERR = 1
if (LROOTS > mxRoot) IERR = 1
if (IERR == 1) then
  write(LF,*)
  write(LF,*) '***************** ERROR *****************'
  call WarningMessage(2,'Max roots exceeded.')
  write(LF,'(1X,A,I6)') 'Input Error: Max roots exceeded.',mxRoot
  write(LF,'(1X,A,I4)') 'Nr of CI roots        LROOTS=',LROOTS
  write(LF,'(1X,A,I4)') 'Nr of optimized roots NROOTS=',NROOTS
  write(LF,*) '************************************************'
  call Quit(_RC_INPUT_ERROR_)
end if
do i=1,NROOTS
  if ((IROOT(i) < 0) .or. (IROOT(i) > LROOTS)) IERR = 1
end do
if (IERR == 1) then
  write(LF,*)
  write(LF,*) '***************** ERROR *****************'
  call WarningMessage(2,'Root specifications are wrong.')
  write(LF,*) 'CHKINP Error: Root specifications are wrong.'
  write(LF,'(1X,A,I4)') 'Nr of CI roots        LROOTS=',LROOTS
  write(LF,'(1X,A,I4)') 'Nr of optimized roots NROOTS=',NROOTS
  write(LF,*) '************************************************'
  call Quit(_RC_INPUT_ERROR_)
end if

IERR = 0
if ((NSYM /= 1) .and. (NSYM /= 2) .and. (NSYM /= 4) .and. (NSYM /= 8)) IERR = 1
if (STSYM > NSYM) IERR = 1
if (IERR == 1) then
  write(LF,*)
  write(LF,*) '***************** ERROR *****************'
  call WarningMessage(2,'Wrong symmetry.')
  write(LF,'(1X,A,I8)') 'CHKINP Error: Wrong symmetry.'
  write(LF,'(1X,A,I8)') 'State symmetry   STSYM=',STSYM
  write(LF,'(1X,A,I8)') 'Point group order NSYM=',NSYM
  write(LF,*) '************************************************'
  call Quit(_RC_INPUT_ERROR_)
end if

if (NSEL < LROOTS+1) then
  write(LF,*)
  write(LF,*) '***************** WARNING ***************'
  call WarningMessage(1,'Too small explicit Hamiltonian.')
  write(LF,*) 'CHKINP Warning: Too small explicit Hamiltonian.'
  write(LF,'(1X,A,I8)') 'Nr of CI roots LROOTS=',LROOTS
  write(LF,'(1X,A,I8)') 'You requested NSEL=',NSEL
  NSEL = LROOTS+1
  write(LF,'(1X,A,I8)') 'It has been reset to NSEL=',NSEL
  write(LF,*) '************************************************'
end if

IERR = 0
if ((NHOLE1 < 0) .and. (NRS1T /= 0)) IERR = 1
if ((NELEC3 < 0) .and. (NRS3T /= 0)) IERR = 1
if ((NACTEL < 0) .and. (NRS2T /= 0)) IERR = 1
if (IERR == 1) then
  write(LF,*)
  write(LF,*) '***************** ERROR *****************'
  call WarningMessage(2,'Error in RAS specification.')
  write(LF,*) 'Error in RAS specification.'
  write(LF,'(1X,A,I8)') 'Max holes in Ras1,     NHOLE1=',NHOLE1
  write(LF,'(1X,A,I8)') 'Max electrons in Ras3, NELEC3=',NELEC3
  write(LF,'(1X,A,I8)') 'Nr of active electrons NACTEL=',NACTEL
  write(LF,*) '*****************************************'
  call Quit(_RC_INPUT_ERROR_)
end if

!if (IPT2 == 1) then
!  if ((NHOLE1 /= 0) .or. (NELEC3 /= 0)) then
!    write(LF,*)
!    write(LF,*) '******************* WARNING *******************'
!    call WarningMessage(1,'"Quasi-canonical" is ignored.')
!    write(LF,*) 'You requested quasicanonical orbitals, but this'
!    write(LF,*) 'is not possible with a true RASSCF calculation.'
!    write(LF,*) 'Your request will be ignored.'
!    write(LF,*) '***********************************************'
!    IPT2 = 0
!  end if
!end if

!GG Sep 03 Check ALTEr
if (NAlter > 0) then
  do iAlter=1,NAlter
    if ((MAlter(iAlter,1) < 1) .or. (MAlter(iAlter,1) > NSym)) then
      write(LF,*)
      write(LF,*) '***************** ERROR *****************'
      call WarningMessage(1,'MAlter input is wrong.')
      write(LF,*) 'Wrong symmetry specie in pair ',iAlter
      write(LF,*) '*****************************************'
      call Quit(_RC_INPUT_ERROR_)
    end if
    if ((MAlter(iAlter,2) < 1) .or. (MAlter(iAlter,3) < 1)) then
      write(LF,*)
      write(LF,*) '***************** ERROR *****************'
      call WarningMessage(1,'MAlter input is wrong.')
      write(LF,*) 'Wrong orbital to exchange in pair ',iAlter
      write(LF,*) '*****************************************'
      call Quit(_RC_INPUT_ERROR_)
    end if
    if ((MAlter(iAlter,2) > nBas(MAlter(iAlter,1))) .or. (MAlter(iAlter,3) > nBas(MAlter(iAlter,1)))) then
      write(LF,*)
      write(LF,*) '***************** ERROR *****************'
      call WarningMessage(1,'MAlter input is wrong.')
      write(LF,*) 'Wrong orbital to exchange in pair ',iAlter
      write(LF,*) '*****************************************'
      call Quit(_RC_INPUT_ERROR_)
    end if
  end do
end if

!BOR  Check INVEC
if ((INVEC < 0) .or. (INVEC > 6)) then
  write(LF,*)
  write(LF,*) '************* ERROR ***************'
  ! This should be impossible:...
  call WarningMessage(2,'Keyword for start orbitals is missing.')
  write(LF,*) 'Keyword for start orbitals missing.'
  write(LF,*) 'Use either CORE, LUMORB, or JOBIPH.'
  write(LF,*) '***********************************'
  call Quit(_RC_INPUT_ERROR_)
end if

! PAM Krapperup Nov 05: Orbital print format.
! First question: Which orbital spaces are eligible for printing?
! No user selection, so fall back on default choice.
if (OutFmt1 == 'DEFAULT ') OutFmt1 = 'FEW     '
! Second question: How should they be printed?
if (OutFmt2 == 'DEFAULT ') then
! No user selection, so fall back on default choice.
  if (NTOT < 256) then
    OutFmt2 = 'FULL    '
  else
    OutFmt2 = 'COMPACT '
  end if
end if
! Third: has the user provided input for energy/occupation thresholds?
! A negative PROTHR shows no user value was given in input.
if (PROTHR < 0.0d0) then
  if (OutFmt1 == 'ALL     ') then
    PROTHR = 0.0d0
    PRETHR = 1.0d100
  else
! Else, format is FEW or NOCORE (or NOTHING, but then nothing is printed)
    PROTHR = 0.0d0
    PRETHR = 0.15d0
  end if
end if
!----------------------------------------------------------------------*

contains

pure function cumsum(X) result(res)
  integer, intent(in) :: X(:)
  integer :: res(size(X))
  integer :: i
  if (size(X) > 0) res(1) = X(1)
  do i=2,size(res)
    res(i) = res(i-1)+X(i)
  end do
end function cumsum

end subroutine ChkInp
