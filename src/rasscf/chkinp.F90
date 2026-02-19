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

use linalg_mod, only: abort_
use gas_data, only: iDoGAS, iGSOCCX, nGAS, nGSSH
use rasscf_global, only: iRoot, ITMAX, lRoots, MAXIT, MAXJT, NAC, NIN, nRoots, OutFmt1, OutFmt2, PreThr, ProThr, ThFact, ThrE, &
                         ThrEn, ThrSX, ThrTE
use general_data, only: INVEC, ISPIN, MALTER, NACTEL, NALTER, NASH, NBAS, NDEL, NELEC3, NFRO, NHOLE1, NISH, NORB, NRS1, NRS1T, &
                        NRS2, NRS2T, NRS3, NRS3T, NSEL, NSSH, NSYM, NTOT, STSYM
use Molcas, only: MxAct, MxBas, MxGAS, MxIna, MxOrb, MxRoot
use RASDim, only: MxCIIt, MxIter, MxSXIt
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: c_orbs_per_GAS(nGAS), iA0, iAlter, iB0, iC0, ierr, ierr1, ierr2, iSym
#include "warnings.h"

!c_orbs_per_GAS: These are the cumulated number of spin orbitals per GAS space.

!----------------------------------------------------------------------*
! Local print level (if any)
IERR = 0
if (NTOT > mxOrb) then
  write(u6,*)
  write(u6,*) '************ ERROR ******************'
  call WarningMessage(2,'Too many orbitals.')
  write(u6,'(1X,A,I8)') 'Too many orbitals NTOT=',NTOT
  write(u6,'(1X,A,I8)') 'Limit is MXORB=',MXORB
  write(u6,*) '*************************************'
  IERR = 1
end if
if (NAC > mxAct) then
  write(u6,*)
  write(u6,*) '*************** ERROR ***************'
  call WarningMessage(2,'Too many active orbitals.')
  write(u6,'(1X,A,I8)') 'Too many active orbitals NAC=',NAC
  write(u6,'(1X,A,I8)') 'Limit is MXACT=',MXACT
  write(u6,*) '*************************************'
  IERR = 1
end if
if (NIN > mxIna) then
  write(u6,*)
  write(u6,*) '*************** ERROR ***************'
  call WarningMessage(2,'Too many inactive orbitals.')
  write(u6,'(1X,A,I8)') 'Too many inactive orbitals NIN=',NIN
  write(u6,'(1X,A,I8)') 'Limit is MXINA=',MXINA
  write(u6,*) '*************************************'
  IERR = 1
end if
if (NACTEL > 2*NAC) then
  write(u6,*)
  write(u6,*) '********************* ERROR **********************'
  call WarningMessage(2,'Too many active electrons.')
  write(u6,'(1X,A,I6)') 'Too many active electrons NACTEL=',NACTEL
  write(u6,'(1X,A,I6)') 'Cannot be more than 2*Nr of active orbitals=',2*NAC
  write(u6,*) '**************************************************'
  IERR = 1
end if
if (NHOLE1 > 2*NRS1T) then
  write(u6,*)
  write(u6,*) '******************** WARNING *********************'
  call WarningMessage(1,'Too many holes in Ras1.')
  write(u6,'(1X,A,I6)') 'You allow too many holes in Ras1 NHOLE1=',NHOLE1
  write(u6,'(1X,A,I6)') 'Cannot be more than 2*Nr of Ras1 orbitals=',2*NRS1T
  NHOLE1 = 2*NRS1T
  write(u6,'(1X,A,I6)') 'NHOLE1 has been reset to ',NHOLE1
  write(u6,*) '**************************************************'
end if
if (NELEC3 > 2*NRS3T) then
  write(u6,*)
  write(u6,*) '******************** WARNING *********************'
  call WarningMessage(1,'Too many electrons in Ras3.')
  write(u6,'(1X,A,I6)') 'You allow too many electrons in Ras3 NELEC3=',NELEC3
  write(u6,'(1X,A,I6)') 'Cannot be more than 2*Nr of Ras3 orbitals=',2*NRS3T
  NELEC3 = 2*NRS3T
  write(u6,'(1X,A,I6)') 'NELEC3 has been reset to ',NELEC3
  write(u6,*) '**************************************************'
end if
if (.not. iDoGas) then !(DM)
  if (NACTEL > 2*(NRS1T+NRS2T)+NELEC3) then
    write(u6,*)
    write(u6,*) '********************* ERROR **********************'
    call WarningMessage(2,'Too many active electrons.')
    write(u6,'(1X,A,I8)') 'Too many active electrons NACTEL=',NACTEL
    write(u6,'(1X,A,I8)') '(Incompatible with RAS restrictions).'
    write(u6,*) '**************************************************'
    IERR = 1
  end if
  if (NACTEL < 2*NRS1T-NHOLE1) then
    write(u6,*)
    write(u6,*) '********************* ERROR **********************'
    call WarningMessage(2,'Too few active electrons.')
    write(u6,'(1X,A,I8)') 'Too few active electrons NACTEL=',NACTEL
    write(u6,'(1X,A,I8)') '(Incompatible with RAS restrictions).'
    write(u6,*) '**************************************************'
    IERR = 1
  end if
  if ((NHOLE1 == 0) .and. (NRS1T > 0)) then
    write(u6,*)
    write(u6,*) '******************* WARNING *******************'
    call WarningMessage(1,'You allow no holes in Ras1')
    write(u6,*) 'You allow no holes in Ras1.                    '
    write(u6,*) 'This may be deliberate, but may give numerical '
    write(u6,*) 'problems in SXCTL section.'
    write(u6,*) '***********************************************'
  end if
  if ((NELEC3 == 0) .and. (NRS3T > 0)) then
    write(u6,*)
    write(u6,*) '******************* WARNING *******************'
    call WarningMessage(1,'You allow no electrons in Ras3')
    write(u6,*) 'You allow no electrons in Ras3.'
    write(u6,*) 'This may be deliberate, but may give numerical '
    write(u6,*) 'problems in SXCTL section.'
    write(u6,*) '***********************************************'
  end if
! for GAS
else
  if (NGAS > mxGAS) then !(SJS)
    write(u6,*)
    call WarningMessage(2,'GASSCF: too many GAS spaces.')
    write(u6,*) ' ************* ERROR **************'
    write(u6,*) ' GASSCF: Too many GAS spaces. '
    write(u6,*) ' Can be increased up to 16 by changing'
    write(u6,*) ' mxGAS in the Molcas module'
    write(u6,*) ' **********************************'
    call Quit(_RC_INPUT_ERROR_)
  end if
  if (nactel /= igsoccx(ngas,2)) then
    write(u6,*)
    write(u6,*) '**************** ERROR *************************'
    write(u6,*) 'nactel not match occupation'
    write(u6,*) 'nactel=',nactel,'igsoccx:',igsoccx(ngas,2)
    write(u6,*) '************************************************'
  end if

  c_orbs_per_GAS = cumsum(sum(ngssh(:nGAS,:nSym),2)*2)
  if (any(c_orbs_per_GAS < igsoccx(:nGAS,1))) then
    write(u6,*)
    write(u6,*) 'In at least one GAS space, the minimum required '
    write(u6,*) 'particle number by GAS constraints '
    write(u6,*) 'is larger than the particle number '
    write(u6,*) 'allowed by the Pauli principle.'
    write(u6,*)
    call abort_('GASSCF: Pauli forbidden.')
  end if
  ! Conceptionally this should not be a problem, but the code
  ! assumes it to be not the case.
  if (any(c_orbs_per_GAS < igsoccx(:nGAS,2))) then
    write(u6,*)
    write(u6,*) 'In at least one GAS space, the maximum allowed '
    write(u6,*) 'particle number by GAS constraints '
    write(u6,*) 'is larger than the particle number '
    write(u6,*) 'allowed by the Pauli principle.'
    write(u6,*)
    call abort_('GASSCF: Pauli forbidden.')
  end if
end if
if ((NSYM /= 1) .and. (NSYM /= 2) .and. (NSYM /= 4) .and. (NSYM /= 8)) then
  write(u6,*)
  call WarningMessage(2,'Wrong nr of symmetries.')
  write(u6,*) '*************** ERROR ***************'
  write(u6,'(1X,A,I8)') 'Nr of symmetries NSYM=',NSYM
  write(u6,*) ' Only possible values are 1,2,4 or 8.'
  write(u6,*) '*************************************'
  IERR = 1
end if
if (IERR == 1) call Quit(_RC_INPUT_ERROR_)

IERR1 = 0
if (any(nBas(1:nSym) < 0)) IERR1 = 1
if (any(nFro(1:nSym) < 0)) IERR1 = 1
if (any(nDel(1:nSym) < 0)) IERR1 = 1
if (any(nOrb(1:nSym) < 0)) IERR1 = 1
if (any(nISh(1:nSym) < 0)) IERR1 = 1
if (any(nASh(1:nSym) < 0)) IERR1 = 1
if (any(nSSh(1:nSym) < 0)) IERR1 = 1
if (any(nRS1(1:nSym) < 0)) IERR1 = 1
if (any(nRS2(1:nSym) < 0)) IERR1 = 1
if (any(nRS3(1:nSym) < 0)) IERR1 = 1
IERR2 = 0
if (any(nBas(1:nSym) > mxBas)) IERR2 = 1
if (any(nFro(1:nSym) > mxBas)) IERR2 = 1
if (any(nDel(1:nSym) > mxBas)) IERR2 = 1
if (any(nOrb(1:nSym) > mxBas)) IERR2 = 1
if (any(nISh(1:nSym) > mxBas)) IERR2 = 1
if (any(nASh(1:nSym) > mxBas)) IERR2 = 1
if (any(nSSh(1:nSym) > mxBas)) IERR2 = 1
if (any(nRS1(1:nSym) > mxBas)) IERR2 = 1
if (any(nRS2(1:nSym) > mxBas)) IERR2 = 1
if (any(nRS3(1:nSym) > mxBas)) IERR2 = 1
if (IERR1+IERR2 > 0) then
  write(u6,*)
  write(u6,*) '****************** ERROR *******************'
  call WarningMessage(2,'Erroneous nr of orbitals.')
  write(u6,*) 'Inappropriate nr of orbitals. One or more of'
  write(u6,*) 'these orbital counts is wrong or too large.'
  if (IERR1 > 0) write(u6,*) ' Negative values.'
  if (IERR2 > 0) write(u6,*) ' Extremely large values.'
  write(u6,'(1X,A,8I4)') '   All orbitals:',(NBAS(ISYM),ISYM=1,NSYM)
  write(u6,'(1X,A,8I4)') '         Frozen:',(NFRO(ISYM),ISYM=1,NSYM)
  write(u6,'(1X,A,8I4)') '       Inactive:',(NISH(ISYM),ISYM=1,NSYM)
  write(u6,'(1X,A,8I4)') '         Active:',(NASH(ISYM),ISYM=1,NSYM)
  write(u6,'(1X,A,8I4)') '          RAS-1:',(NRS1(ISYM),ISYM=1,NSYM)
  write(u6,'(1X,A,8I4)') '          RAS-2:',(NRS2(ISYM),ISYM=1,NSYM)
  write(u6,'(1X,A,8I4)') '          RAS-3:',(NRS3(ISYM),ISYM=1,NSYM)
  write(u6,'(1X,A,8I4)') '      Secondary:',(NSSH(ISYM),ISYM=1,NSYM)
  write(u6,'(1X,A,8I4)') '        Deleted:',(NDEL(ISYM),ISYM=1,NSYM)
  write(u6,'(1X,A,8I4)') 'Basis functions:',(NBAS(ISYM),ISYM=1,NSYM)
  if (IERR1 > 0) write(u6,*) ' Have you used a too small basis set?'
  write(u6,*) '********************************************'
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
  write(u6,*)
  write(u6,*) '******************** ERROR *********************'
  call WarningMessage(2,'No such wave function.')
  write(u6,*) 'The following combined specifications are wrong.'
  write(u6,'(1X,A,8I4)') 'Nr of active electrons NACTEL=',NACTEL
  write(u6,'(1X,A,8I4)') 'Nr of active orbitals     NAC=',NAC
  write(u6,'(1X,A,8I4)') 'Spin degeneracy         ISPIN=',ISPIN
  write(u6,*) 'There can be no such wave function.'
  write(u6,*) '************************************************'
  call Quit(_RC_INPUT_ERROR_)
end if

if (MAXIT > mxIter) then
  write(u6,*)
  write(u6,*) '*************** WARNING ****************'
  call WarningMessage(1,'Too many macro-iterations.')
  write(u6,*) 'Too many macro-iterations requested.'
  MAXIT = MXITER
  write(u6,'(1X,A,I8)') 'Reset to maximum, new MAXIT=',MAXIT
  write(u6,*) '****************************************'
end if
if (MAXJT > (mxCiIt-2)) then
  write(u6,*)
  write(u6,*) '*************** WARNING ****************'
  call WarningMessage(1,'Too many CI-iterations.')
  write(u6,*) 'Too many CI-iterations requested.'
  MAXJT = mxCiIt-2
  write(u6,'(1X,A,I8)') 'Reset to maximum, new MAXJT=',MAXJT
  write(u6,*) '****************************************'
end if
if (ITMAX > MXSXIT) then
  write(u6,*)
  write(u6,*) '*************** WARNING ****************'
  call WarningMessage(1,'Too many SX-iterations.')
  write(u6,*) 'Too many SX-iterations requested.'
  ITMAX = MXSXIT
  write(u6,'(1X,A,I8)') 'Reset to maximum, new ITMAX=',ITMAX
  write(u6,*) '****************************************'
end if

THRE = max(Zero,THRE)
THRTE = max(Zero,THRTE)
THRSX = max(Zero,THRSX)
THREN = max(Zero,THREN)
THFACT = max(Zero,THFACT)

IERR = 0
if (NROOTS > mxRoot) IERR = 1
if (LROOTS > mxRoot) IERR = 1
if (IERR == 1) then
  write(u6,*)
  write(u6,*) '***************** ERROR *****************'
  call WarningMessage(2,'Max roots exceeded.')
  write(u6,'(1X,A,I6)') 'Input Error: Max roots exceeded.',mxRoot
  write(u6,'(1X,A,I4)') 'Nr of CI roots        LROOTS=',LROOTS
  write(u6,'(1X,A,I4)') 'Nr of optimized roots NROOTS=',NROOTS
  write(u6,*) '************************************************'
  call Quit(_RC_INPUT_ERROR_)
end if
if (any(IROOT(1:NROOTS) < 0) .or. any(IROOT(1:NROOTS) > LROOTS)) IERR = 1
if (IERR == 1) then
  write(u6,*)
  write(u6,*) '***************** ERROR *****************'
  call WarningMessage(2,'Root specifications are wrong.')
  write(u6,*) 'CHKINP Error: Root specifications are wrong.'
  write(u6,'(1X,A,I4)') 'Nr of CI roots        LROOTS=',LROOTS
  write(u6,'(1X,A,I4)') 'Nr of optimized roots NROOTS=',NROOTS
  write(u6,*) '************************************************'
  call Quit(_RC_INPUT_ERROR_)
end if

IERR = 0
if ((NSYM /= 1) .and. (NSYM /= 2) .and. (NSYM /= 4) .and. (NSYM /= 8)) IERR = 1
if (STSYM > NSYM) IERR = 1
if (IERR == 1) then
  write(u6,*)
  write(u6,*) '***************** ERROR *****************'
  call WarningMessage(2,'Wrong symmetry.')
  write(u6,'(1X,A,I8)') 'CHKINP Error: Wrong symmetry.'
  write(u6,'(1X,A,I8)') 'State symmetry   STSYM=',STSYM
  write(u6,'(1X,A,I8)') 'Point group order NSYM=',NSYM
  write(u6,*) '************************************************'
  call Quit(_RC_INPUT_ERROR_)
end if

if (NSEL < LROOTS+1) then
  write(u6,*)
  write(u6,*) '***************** WARNING ***************'
  call WarningMessage(1,'Too small explicit Hamiltonian.')
  write(u6,*) 'CHKINP Warning: Too small explicit Hamiltonian.'
  write(u6,'(1X,A,I8)') 'Nr of CI roots LROOTS=',LROOTS
  write(u6,'(1X,A,I8)') 'You requested NSEL=',NSEL
  NSEL = LROOTS+1
  write(u6,'(1X,A,I8)') 'It has been reset to NSEL=',NSEL
  write(u6,*) '************************************************'
end if

IERR = 0
if ((NHOLE1 < 0) .and. (NRS1T /= 0)) IERR = 1
if ((NELEC3 < 0) .and. (NRS3T /= 0)) IERR = 1
if ((NACTEL < 0) .and. (NRS2T /= 0)) IERR = 1
if (IERR == 1) then
  write(u6,*)
  write(u6,*) '***************** ERROR *****************'
  call WarningMessage(2,'Error in RAS specification.')
  write(u6,*) 'Error in RAS specification.'
  write(u6,'(1X,A,I8)') 'Max holes in Ras1,     NHOLE1=',NHOLE1
  write(u6,'(1X,A,I8)') 'Max electrons in Ras3, NELEC3=',NELEC3
  write(u6,'(1X,A,I8)') 'Nr of active electrons NACTEL=',NACTEL
  write(u6,*) '*****************************************'
  call Quit(_RC_INPUT_ERROR_)
end if

!if (IPT2 == 1) then
!  if ((NHOLE1 /= 0) .or. (NELEC3 /= 0)) then
!    write(u6,*)
!    write(u6,*) '******************* WARNING *******************'
!    call WarningMessage(1,'"Quasi-canonical" is ignored.')
!    write(u6,*) 'You requested quasicanonical orbitals, but this'
!    write(u6,*) 'is not possible with a true RASSCF calculation.'
!    write(u6,*) 'Your request will be ignored.'
!    write(u6,*) '***********************************************'
!    IPT2 = 0
!  end if
!end if

!GG Sep 03 Check ALTEr
if (NAlter > 0) then
  do iAlter=1,NAlter
    if ((MAlter(iAlter,1) < 1) .or. (MAlter(iAlter,1) > NSym)) then
      write(u6,*)
      write(u6,*) '***************** ERROR *****************'
      call WarningMessage(1,'MAlter input is wrong.')
      write(u6,*) 'Wrong symmetry specie in pair ',iAlter
      write(u6,*) '*****************************************'
      call Quit(_RC_INPUT_ERROR_)
    end if
    if ((MAlter(iAlter,2) < 1) .or. (MAlter(iAlter,3) < 1)) then
      write(u6,*)
      write(u6,*) '***************** ERROR *****************'
      call WarningMessage(1,'MAlter input is wrong.')
      write(u6,*) 'Wrong orbital to exchange in pair ',iAlter
      write(u6,*) '*****************************************'
      call Quit(_RC_INPUT_ERROR_)
    end if
    if ((MAlter(iAlter,2) > nBas(MAlter(iAlter,1))) .or. (MAlter(iAlter,3) > nBas(MAlter(iAlter,1)))) then
      write(u6,*)
      write(u6,*) '***************** ERROR *****************'
      call WarningMessage(1,'MAlter input is wrong.')
      write(u6,*) 'Wrong orbital to exchange in pair ',iAlter
      write(u6,*) '*****************************************'
      call Quit(_RC_INPUT_ERROR_)
    end if
  end do
end if

!BOR  Check INVEC
if ((INVEC < 0) .or. (INVEC > 6)) then
  write(u6,*)
  write(u6,*) '************* ERROR ***************'
  ! This should be impossible:...
  call WarningMessage(2,'Keyword for start orbitals is missing.')
  write(u6,*) 'Keyword for start orbitals missing.'
  write(u6,*) 'Use either CORE, LUMORB, or JOBIPH.'
  write(u6,*) '***********************************'
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
if (PROTHR < Zero) then
  if (OutFmt1 == 'ALL     ') then
    PROTHR = Zero
    PRETHR = 1.0e100_wp
  else
! Else, format is FEW or NOCORE (or NOTHING, but then nothing is printed)
    PROTHR = Zero
    PRETHR = 0.15_wp
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
