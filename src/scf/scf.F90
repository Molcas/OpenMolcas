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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2003-2005, Valera Veryazov                             *
!               2017, Roland Lindh                                     *
!***********************************************************************

subroutine SCF(ireturn)
!***********************************************************************
!                                                                      *
!     purpose: perform RHF and UHF calculations                        *
!                                                                      *
!***********************************************************************

use Interfaces_SCF, only: OccDef
use OFembed, only: Do_OFemb
use InfSCF, only: Atom, AufB, BName, BType, CMO, DSCF, HDiag, iStatPrn, KSDFT, mOV, nBB, nCore, nD, nDisc, nnB, OccNo, OnlyProp
use SCFFiles, only: LuInp
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: iReturn
integer(kind=iwp) :: iTerm, LthH, LUOrb, MemLow, MemSew
real(kind=wp) :: SIntTh, TCPU1, TCPU2, TWALL1, TWALL2
logical(kind=iwp) :: FstItr, Semi_Direct
integer(kind=iwp), external :: isStructure

#include "warnings.h"

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

call CWTime(TCPU1,TWall1)
call SCF_Init()
iTerm = 0

call OpnFls_SCF()
call ReadIn_SCF(SIntTh)
LuOrb = LuInp

Semi_Direct = (DSCF .and. ((nDisc /= 0) .or. ((nCore /= 0) .and. (nDisc == 0))))
if (Semi_Direct) then
  call mma_MaxDBLE(MemSew)
  MemLow = min(MemSew/2,1024*1024)
  MemSew = max(MemSew/10,MemLow)
  call xsetmem_ints(MemSew)
end if

call Init_SCF()
!                                                                      *
!***********************************************************************
!                                                                      *
! Pick up the starting orbitals and set the occupation numbers.

LuOrb = LuInp
call SOrb(LuOrb,SIntTh,iTerm)
call OccDef(OccNo,nnB,nD,CMO,nBB)

call mma_deallocate(HDiag)
if (Aufb) then
  lthH = nBB*nD
else
  lthH = mOV
end if
call mma_allocate(HDiag,lthH,Label='HDiag')
!                                                                      *
!***********************************************************************
!                                                                      *
call WrInp_SCF(SIntTh)

call Cre_SCFWfn()

FstItr = .true.

if (.not. OnlyProp) call WfCtl_SCF(iTerm,KSDFT,FstItr,SIntTh)

call FinalSCF()
call Free_TLists()
call mma_deallocate(BType)
call mma_deallocate(Atom)
call mma_deallocate(BName)

call CWTime(TCPU2,TWall2)

call GMFree()
call ClsFls_SCF()
if (Semi_Direct) call xRlsMem_Ints()

! Call MolDen Interface

if (nD == 1) then
  call Molden_Interface(nD-1,'SCFORB','MD_SCF')
else
  call Molden_Interface(nD-1,'UHFORB','MD_SCF')
end if
if (iStatPRN > 0) call FastIO('STATUS')

! Everything has to come to an end...

iReturn = iTerm

if (Do_OFemb) then
  if (isStructure() == 1) then
    if (iReturn /= _RC_ALL_IS_WELL_) call WarningMessage(1,'SCF: non-zero return code.')
    iReturn = _RC_CONTINUE_LOOP_
    call Check_FThaw(iReturn)
  end if
end if

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
return

end subroutine SCF
