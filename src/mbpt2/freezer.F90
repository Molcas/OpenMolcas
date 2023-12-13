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

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: EAll(*)
integer(kind=iwp), intent(in) :: nFre, nSym, nFro(nSym), nOcc(nSym), nBas(nSym)
integer(kind=iwp), intent(out) :: nFro1(nSym)
logical(kind=iwp), intent(in) :: LocPrt
integer(kind=iwp) :: iCount, iFre, iOcc(8), iSym, jOcc, jSym, kAll, kOcc, lEOcc, NumFre
real(kind=wp) :: xMin
integer(kind=iwp), allocatable :: Pivot(:), Point(:)
real(kind=wp), allocatable :: EOcc(:)
integer(kind=iwp), external :: Cho_iRange
character(len=*), parameter :: SecNam = 'Freezer'

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
  nFro1(:) = 0
end if

! Set up array of active occupied orbital energies.
! -------------------------------------------------

iOcc(1) = 0
lEOcc = nOcc(1)
do iSym=2,nSym
  iOcc(iSym) = lEOcc
  lEOcc = lEOcc+nOcc(iSym)
end do

call mma_allocate(EOcc,lEOcc,label='ScrOcc')
call mma_allocate(Pivot,lEOcc,label='Pivot')
call mma_allocate(Point,nFre,label='Point')

iCount = 1
do iSym=1,nSym
  kAll = iCount+nFro(iSym)
  kOcc = iOcc(iSym)+1
  call dCopy_(nOcc(iSym),EAll(kAll),1,EOcc(kOcc),1)
  iCount = iCount+nBas(iSym)
end do

! Find pointers to lowest nFre occupied orbital energies.
! -------------------------------------------------------

xMin = -1.0e15_wp
NumFre = nFre
EOcc(:) = -EOcc(:) ! DiaMax finds MAX values
call CD_DiaMax(EOcc,lEOcc,Pivot,Point,NumFre,xMin)
if (NumFre /= nFre) then
  write(u6,*) SecNam,': an error occurred in CD_DiaMax!'
  write(u6,*) 'NumFre = ',NumFre,' != ',nFre,' = nFre'
  call SysAbendMsg(SecNam,'CD_DiaMax failure',' ')
end if

! Set up nFro1 array.
! -------------------

do iFre=1,nFre
  iSym = Cho_iRange(Point(iFre),iOcc,nSym,.false.)
  nFro1(iSym) = nFro1(iSym)+1
end do

! If requested, print.
! --------------------

if (LocPrt) then
  write(u6,'(/,3X,A,A,A)') 'Output from ',SecNam,':'
  write(u6,'(1X,A,I5,A)') 'The',nFre,' lowest occupied orbitals have been frozen.'
  write(u6,'(1X,A)') 'List of frozen occupied orbitals:'
  do iFre=1,nFre
    kOcc = Point(iFre)
    jSym = Cho_iRange(kOcc,iOcc,nSym,.false.)
    jOcc = kOcc-iOcc(jSym)
    write(u6,'(1X,A,I5,A,I1,A,F15.8)') 'Occupied orbital',jOcc,' of symmetry ',jSym,' and energy ',-EOcc(kOcc)
  end do
end if

! Free memory.
! ------------

call mma_deallocate(EOcc)
call mma_deallocate(Pivot)
call mma_deallocate(Point)

end subroutine Freezer
