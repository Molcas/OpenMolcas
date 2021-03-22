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

implicit none
integer nFre
logical Prnt
integer nSym
real*8 EOcc(*)
integer nFro(nSym)
integer nOcc(nSym)
integer nFro1(nSym)
#include "WrkSpc.fh"

character*7 SecNam
parameter(SecNam='RPA_Frz')

integer Cho_iRange
external Cho_iRange

integer ip_Point, l_Point
integer ip_E, l_E
integer ip_Pivot, l_Pivot
integer ip_iOcc, l_iOcc
integer iCount
integer NumFre
integer iFre
integer i
integer iSym
integer ip1

real*8 xMin

integer j
integer iOcc
iOcc(j) = iWork(ip_iOcc-1+j)

if (nSym < 1 .or. nSym > 8) then
  write(6,'(A,I6)') 'nSym=',nSym
  call RPA_Warn(3,SecNam//': illegal nSym')
else if (nSym == 1) then
  nFro1(1) = max(nFre,0)
  return
else
  call iZero(nFro1,nSym)
end if
if (nFre < 1) return

l_Point = nFre
l_iOcc = nSym
l_E = nOcc(1)
do iSym=2,nSym
  l_E = l_E+nOcc(iSym)
end do
l_Pivot = l_E
if (nFre > l_E) then
  call RPA_Warn(4,SecNam//': too many orbitals to freeze')
end if
call GetMem('ScrPnt','Allo','Inte',ip_Point,l_Point)
call GetMem('iOcc','Allo','Inte',ip_iOcc,l_iOcc)
call GetMem('ScrOccE','Allo','Real',ip_E,l_E)
call GetMem('Pivot','Allo','Inte',ip_Pivot,l_Pivot)

iCount = 0
ip1 = ip_iOcc-1
do iSym=1,nSym
  iWork(ip1+iSym) = iCount
  iCount = iCount+nOcc(iSym)
end do

iCount = 1
do iSym=1,nSym
  call dCopy_(nOcc(iSym),EOcc(iCount+nFro(iSym)),1,Work(ip_E+iOcc(iSym)),1)
  iCount = iCount+nFro(iSym)+nOcc(iSym)
end do

xMin = -1.0d15
NumFre = nFre
call dScal_(l_E,-1.0d0,Work(ip_E),1)
call CD_DiaMax(Work(ip_E),l_E,iWork(ip_Pivot),iWork(ip_Point),NumFre,xMin)
if (NumFre /= nFre) then
  write(6,'(2(A,I12))') 'NumFre=',NumFre,'  nFre=',nFre
  call RPA_Warn(3,SecNam//': NumFre != nFre')
end if

do iFre=1,nFre
  iSym = Cho_iRange(iWork(ip_Point-1+iFre),iWork(ip_iOcc),nSym,.false.)
  nFro1(iSym) = nFro1(iSym)+1
end do

if (Prnt) then
  write(6,'(/,3X,A,A,A)') 'Output from ',SecNam,':'
  write(6,'(A,I5,A)') 'The',nFre,' lowest occupied orbitals have been frozen.'
  write(6,'(A)') 'List of frozen occupied orbitals:'
  do iFre=1,nFre
    iCount = iWork(ip_Point-1+iFre)
    iSym = Cho_iRange(iCount,iWork(ip_iOcc),nSym,.false.)
    i = iCount-iOcc(iSym)
    write(6,'(1X,A,I5,A,I1,A,F15.8)') 'Occupied orbital',i,' of symmetry ',iSym,' and energy ',-Work(ip_E-1+iCount)
  end do
  call xFlush(6)
end if

call GetMem('Pivot','Free','Inte',ip_Pivot,l_Pivot)
call GetMem('OccE','Free','Real',ip_E,l_E)
call GetMem('iOcc','Free','Inte',ip_iOcc,l_iOcc)
call GetMem('Point','Free','Inte',ip_Point,l_Point)

end subroutine RPA_Frz
