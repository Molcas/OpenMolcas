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

subroutine BasFun_Atom(nBas_per_Atom,nBas_Start,Name,nBas,nAtoms,DoPrint)
! Author: Y. Carissan [put in separate subroutine by T.B. Pedersen]

implicit none
#include "Molcas.fh"
integer nBas, nAtoms
integer nBas_per_Atom(nAtoms), nBas_Start(nAtoms)
character*(LENIN8) Name(nBas)
logical DoPrint
character*11 SecNam
parameter(SecNam='BasFun_Atom')
integer iAt, iAt1, nBasAt, iBas, iCount
character*(LENIN) Lbl, LblOld
character*80 Txt, Formt

! Counters.
! ---------

iAt = 1
nBasAt = 1
LblOld = Name(1)(1:LENIN)
do iBas=2,nBas
  Lbl = Name(iBas)(1:LENIN)
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
  write(6,Formt) 'Atom  ',(iAt,iAt=1,nAtoms),'Start ',(nBas_Start(iAt),iAt=1,nAtoms),'nBas  ',(nBas_per_Atom(iAt),iAt=1,nAtoms)
end if

end subroutine BasFun_Atom
