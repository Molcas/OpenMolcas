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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine RPA_Frz(nFre,Prnt,nSym,EOcc,nFro,nOcc,nFro1)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nFre, nSym, nFro(nSym), nOcc(nSym)
logical(kind=iwp), intent(in) :: Prnt
real(kind=wp), intent(in) :: EOcc(*)
integer(kind=iwp), intent(out) :: nFro1(nSym)
integer(kind=iwp) :: l_Point, l_E, iCount, NumFre, iFre, i, iSym
real(kind=wp) :: xMin
integer(kind=iwp), allocatable :: Point(:), Pivot(:), iOcc(:)
real(kind=wp), allocatable :: E(:)
character(len=*), parameter :: SecNam = 'RPA_Frz'
integer(kind=iwp), external :: Cho_iRange

if ((nSym < 1) .or. (nSym > 8)) then
  write(u6,'(A,I6)') 'nSym=',nSym
  call RPA_Warn(3,SecNam//': illegal nSym')
else if (nSym == 1) then
  nFro1(1) = max(nFre,0)
  return
else
  nFro1(:) = 0
end if
if (nFre < 1) return

l_Point = nFre
l_E = nOcc(1)
do iSym=2,nSym
  l_E = l_E+nOcc(iSym)
end do
if (nFre > l_E) then
  call RPA_Warn(4,SecNam//': too many orbitals to freeze')
end if
call mma_allocate(Point,l_Point,label='ScrPnt')
call mma_allocate(iOcc,l_Point,label='iOcc')
call mma_allocate(E,l_Point,label='ScrOccE')
call mma_allocate(Pivot,l_Point,label='Pivot')

iCount = 0
do iSym=1,nSym
  iOcc(iSym) = iCount
  iCount = iCount+nOcc(iSym)
end do

iCount = 1
do iSym=1,nSym
  call dCopy_(nOcc(iSym),EOcc(iCount+nFro(iSym)),1,E(iOcc(iSym)+1),1)
  iCount = iCount+nFro(iSym)+nOcc(iSym)
end do

xMin = -1.0e15_wp
NumFre = nFre
E(:) = -E(:)
call CD_DiaMax(E,l_E,Pivot,Point,NumFre,xMin)
if (NumFre /= nFre) then
  write(u6,'(2(A,I12))') 'NumFre=',NumFre,'  nFre=',nFre
  call RPA_Warn(3,SecNam//': NumFre != nFre')
end if

do iFre=1,nFre
  iSym = Cho_iRange(Point(iFre),iOcc,nSym,.false.)
  nFro1(iSym) = nFro1(iSym)+1
end do

if (Prnt) then
  write(u6,'(/,3X,A,A,A)') 'Output from ',SecNam,':'
  write(u6,'(A,I5,A)') 'The',nFre,' lowest occupied orbitals have been frozen.'
  write(u6,'(A)') 'List of frozen occupied orbitals:'
  do iFre=1,nFre
    iCount = Point(iFre)
    iSym = Cho_iRange(iCount,iOcc,nSym,.false.)
    i = iCount-iOcc(iSym)
    write(u6,'(1X,A,I5,A,I1,A,F15.8)') 'Occupied orbital',i,' of symmetry ',iSym,' and energy ',-E(iCount)
  end do
  call xFlush(u6)
end if

call mma_deallocate(Point)
call mma_deallocate(iOcc)
call mma_deallocate(E)
call mma_deallocate(Pivot)

end subroutine RPA_Frz
