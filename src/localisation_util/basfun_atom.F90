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
! Copyright (C) Yannick Carissan                                       *
!***********************************************************************

subroutine BasFun_Atom(nBas_per_Atom,nBas_Start,BName,nBas,nAtoms,DoPrint)
! Author: Y. Carissan [put in separate subroutine by T.B. Pedersen]

use Definitions, only: iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: nAtoms, nBas
integer(kind=iwp), intent(out) :: nBas_per_Atom(nAtoms), nBas_Start(nAtoms)
character(len=LenIn8), intent(in) :: BName(nBas)
logical(kind=iwp), intent(in) :: DoPrint
integer(kind=iwp) :: iAt, iAt1, iBas, iCount, nBasAt
character(len=LenIn) :: Lbl, LblOld
character(len=80) :: Txt, Formt
character(len=*), parameter :: SecNam = 'BasFun_Atom'

! Counters.
! ---------

iAt = 1
nBasAt = 1
LblOld = BName(1)(1:LenIn)
do iBas=2,nBas
  Lbl = BName(iBas)(1:LenIn)
  if (Lbl /= LblOld) then
    nBas_per_Atom(iAt) = nBasAt
    iAt = iAt+1
    nBasAt = 0
    LblOld = Lbl
  end if
  nBasAt = nBasAt+1
end do
nBas_per_Atom(iAt) = nBasAt

if (iAt /= nAtoms) then ! centers without basis functions
  iAt1 = iAt+1
  do iAt=iAt1,nAtoms
    nBas_per_Atom(iAt) = 0
  end do
end if

! Offsets.
! --------

iCount = 0
do iAt=1,nAtoms
  nBas_Start(iAt) = iCount+1
  iCount = iCount+nBas_per_Atom(iAt)
end do
if (iCount /= nBas) then
  write(Txt,'(A,I9,A,I9)') 'iCount =',iCount,'  nBas =',nBas
  call SysAbendMsg(SecNam,'iCount /= nBas',Txt)
end if

! Print.
! ------

if (DoPrint) then
  write(Formt,'(3(a6,i3,a5))') '(/,a6,',nAtoms,'i5,/,','   a6,',nAtoms,'i5,/,','   a6,',nAtoms,'i5)'
  write(u6,Formt) 'Atom  ',(iAt,iAt=1,nAtoms),'Start ',nBas_Start(:),'nBas  ',nBas_per_Atom(:)
end if

end subroutine BasFun_Atom
