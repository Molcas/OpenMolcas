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

subroutine Freezer(EAll,nFre,nFro,nFro1,nOcc,nBas,nSym,LocPrt)

use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: EAll(*)
integer(kind=iwp), intent(in) :: nFre, nSym, nFro(nSym), nOcc(nSym), nBas(nSym)
integer(kind=iwp), intent(out) :: nFro1(nSym)
logical(kind=iwp), intent(in) :: LocPrt
integer(kind=iwp) :: iCount, iFre, ipEOcc, iOcc(8), ipPivot, ipPoint, iSym, jOcc, jSym, kAll, kOcc, lEOcc, lPoint, NumFre
real(kind=wp) :: xMin
integer(kind=iwp), external :: Cho_iRange
character(len=7), parameter :: SecNam = 'Freezer'
#include "WrkSpc.fh"

! For nSym=1, simply transfer nFre to nFro1.
! Else initialize nFro1 array.
! ------------------------------------------

if ((nSym < 1) .or. (nSym > 8)) then
  write(u6,*) SecNam,': illegal nSym = ',nSym
  call SysAbendMsg(SecNam,'illegal nSym',' ')
else if (nSym == 1) then
  nFro1(1) = nFre
  return
else
  call Cho_iZero(nFro1,nSym)
end if

! Set up array of active occupied orbital energies.
! -------------------------------------------------

lPoint = nFre

iOcc(1) = 0
lEOcc = nOcc(1)
do iSym=2,nSym
  iOcc(iSym) = lEOcc
  lEOcc = lEOcc+nOcc(iSym)
end do

call GetMem('ScrOcc','Allo','Real',ipEOcc,lEOcc)
call GetMem('Pivot','Allo','Inte',ipPivot,lEOcc)
call GetMem('Point','Allo','Inte',ipPoint,lPoint)

iCount = 1
do iSym=1,nSym
  kAll = iCount+nFro(iSym)
  kOcc = ipEOcc+iOcc(iSym)
  call dCopy_(nOcc(iSym),EAll(kAll),1,Work(kOcc),1)
  iCount = iCount+nBas(iSym)
end do

! Find pointers to lowest nFre occupied orbital energies.
! -------------------------------------------------------

xMin = -1.0e15_wp
NumFre = nFre
call dScal_(lEOcc,-One,Work(ipEOcc),1) ! DiaMax finds MAX values
call CD_DiaMax(Work(ipEOcc),lEOcc,iWork(ipPivot),iWork(ipPoint),NumFre,xMin)
if (NumFre /= nFre) then
  write(u6,*) SecNam,': an error occurred in CD_DiaMax!'
  write(u6,*) 'NumFre = ',NumFre,' != ',nFre,' = nFre'
  call SysAbendMsg(SecNam,'CD_DiaMax failure',' ')
end if

! Set up nFro1 array.
! -------------------

do iFre=1,nFre
  iSym = Cho_iRange(iWork(ipPoint-1+iFre),iOcc,nSym,.false.)
  nFro1(iSym) = nFro1(iSym)+1
end do

! If requested, print.
! --------------------

if (LocPrt) then
  write(u6,'(/,3X,A,A,A)') 'Output from ',SecNam,':'
  write(u6,'(1X,A,I5,A)') 'The',nFre,' lowest occupied orbitals have been frozen.'
  write(u6,'(1X,A)') 'List of frozen occupied orbitals:'
  do iFre=1,nFre
    kOcc = iWork(ipPoint-1+iFre)
    jSym = Cho_iRange(kOcc,iOcc,nSym,.false.)
    jOcc = kOcc-iOcc(jSym)
    write(u6,'(1X,A,I5,A,I1,A,F15.8)') 'Occupied orbital',jOcc,' of symmetry ',jSym,' and energy ',-Work(ipEOcc-1+kOcc)
  end do
end if

! Free memory.
! ------------

call GetMem('Point','Free','Inte',ipPoint,lPoint)
call GetMem('Pivot','Free','Inte',ipPivot,lEOcc)
call GetMem('ScrOcc','Free','Real',ipEOcc,lEOcc)

end subroutine Freezer
