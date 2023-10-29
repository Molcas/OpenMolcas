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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_PTS_WrRst(irc,NVT,l_NVT)
!
! Thomas Bondo Pedersen, April 2010.
!
! Purpose: Write restart files (parallel two-step algorithm).

use Cholesky, only: InfRed, InfVec, LuCho, LuMap, LuRed, LuRst, MaxVec, nnBstR, NumCho, nSym
use Cholesky_procedures, only: Cho_X_GetIP_InfVec
#ifdef _DEBUGPRINT_
use Cholesky, only: LuPri
use stdalloc, only: mma_allocate, mma_deallocate
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: l_NVT
integer(kind=iwp), intent(inout) :: NVT(l_NVT)
integer(kind=iwp) :: iAdr, iSym, iV
integer(kind=iwp), pointer :: InfVcT(:,:,:)
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: myNumCho(8)
integer(kind=iwp), allocatable :: IDV(:)
character(len=*), parameter :: SecNam = 'Cho_PTS_WrRst'
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
! Init return code
irc = 0

#ifdef _DEBUGPRINT_
! check that NumCho agrees with distribution
if (l_NVT < nSym) then
  irc = -1
  return
end if
do iSym=1,nSym

  call mma_allocate(IDV,NVT(iSym),Label='IDV')
  myNumCho(iSym) = 0
  call Cho_P_Distrib_Vec(1,NVT(iSym),IDV,myNumCho(iSym))
  call mma_deallocate(IDV)
  if (NumCho(iSym) /= myNumCho(iSym)) then
    write(LuPri,*) SecNam,': NumCho discrepancy in sym. ',iSym
    write(LuPri,*) '  NumCho=',NumCho(iSym)
    write(LuPri,*) 'myNumCho=',myNumCho(iSym)
    write(LuPri,*) '     NVT=',NVT(iSym)
    irc = 1
  end if
end do
myNumCho(1:nSym) = NumCho(1:nSym)
call Cho_GAIGOp(myNumCho,nSym,'+')
do iSym=1,nSym
  if (myNumCho(iSym) /= NVT(iSym)) then
    write(LuPri,*) SecNam,': NumCho discrepancy in sym. ',iSym
    write(LuPri,*) 'Sum of NumCho across nodes=',myNumCho(iSym)
    write(LuPri,*) '                       NVT=',NVT(iSym)
    irc = 2
  end if
end do
if (irc /= 0) return
#endif

! Erase existing files from first step
! (Keep map and vector files, though)
if (LuRst < 1) then
  LuRst = 7
  call DAName_MF_WA(LuRst,'CHORST')
end if
call DAEras(LuRst)
LuRst = 0
if (LuRed < 1) then
  LuRed = 7
  call DAName_MF_WA(LuRst,'CHRED')
end if
call DAEras(LuRed)
LuRed = 0
if (LuMap > 0) then
  call DAClos(LuMap)
  LuMap = 0
end if
do iSym=1,nSym
  if (LuCho(iSym) > 0) then
    call DAClos(LuCho(iSym))
    LuCho(iSym) = 0
  end if
end do

! Re-open files
call Cho_OpenVR(1,2)

! Set InfRed corresponding to only one reduced set
! (All vectors are now stored in 1st reduced set)
InfRed(1) = 0
! Set InfVec data (only disk addresses need updating)
call Cho_X_GetIP_InfVec(InfVcT)
do iSym=1,nSym
  ! InfVec(iV,1,iSym): parent diagonal in rs1 of vector iV of
  ! symmetry iSym
  InfVec(1:NVT(iSym),1,iSym) = InfVcT(1:NVT(iSym),1,iSym)
  InfVec(NVT(iSym)+1:MaxVec,1,iSym) = 0
  ! InfVec(iV,2,iSym): reduced set of vector iV of
  ! symmetry iSym (always rs1 in this implementation)
  InfVec(1:NVT(iSym),2,iSym) = 1
  InfVec(NVT(iSym)+1:MaxVec,2,iSym) = 0
  ! InfVec(iV,3,iSym): disk address of vector iV of
  ! symmetry iSym (always rs1 in this implementation)
  ! Note: this information is for the vectors on this
  ! node ONLY!!
  iAdr = 0
  do iV=1,NumCho(iSym)
    InfVec(iV,3,iSym) = iAdr
    iAdr = iAdr+nnBstR(iSym,1)
  end do
  InfVec(NumCho(iSym)+1:MaxVec,3,iSym) = -1
  ! InfVec(iV,4,iSym): WA disk address of vector iV of
  ! symmetry iSym (redundant? Oh yes, but used for statistics...)
  iAdr = 0
  do iV=1,NVT(iSym)
    InfVec(iV,4,iSym) = iAdr
    iAdr = iAdr+nnBstR(iSym,1)
  end do
  InfVec(NVT(iSym)+1:MaxVec,4,iSym) = -1
  ! InfVec(iV,5,iSym): is defined automatically by Cho_X_Init
  ! No need to set it here - use a dummy value
  InfVec(:,5,iSym) = -1
end do

! Make sure that all vector info is written
call iSwap(nSym,NVT,1,NumCho,1)

! Write restart file
call Cho_WrRstC(1)

! Write reduced set
call Cho_PutRed(1,1)

! Swap back
call iSwap(nSym,NVT,1,NumCho,1)

end subroutine Cho_PTS_WrRst
