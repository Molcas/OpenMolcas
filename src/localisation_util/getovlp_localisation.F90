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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine GetOvlp_Localisation(S,Storage,nBas,nSym)
! Thomas Bondo Pedersen, Dec. 2005.
!
! Purpose: read the overlap matrix and return in S in lower triangular
! storage if Storage="Tri" else in full square storage.

implicit none
real*8 S(*)
character*3 Storage
integer nSym
integer nBas(nSym)
#include "WrkSpc.fh"
character*20 SecNam
parameter(SecNam='GetOvlp_Localisation')
logical Debug
parameter(Debug=.false.)
character*3 myStorage
character*8 Label
integer irc, iOpt, iComp, iSyLbl
integer iSym, l_Scr, ip_Scr, kTri, kSq, l_Tri

l_Tri = nBas(1)*(nBas(1)+1)/2
do iSym=2,nSym
  l_Tri = l_Tri+nBas(iSym)*(nBas(iSym)+1)/2
end do
l_Scr = l_Tri+4
call GetMem('OvlpScr','Allo','Real',ip_Scr,l_Scr)

irc = -1
iOpt = 2
iComp = 1
iSyLbl = 1
Label = 'Mltpl  0'
call RdOne(irc,iOpt,Label,iComp,Work(ip_Scr),iSyLbl)
if (irc /= 0) then
  write(6,*) SecNam,': RdOne returned ',irc
  write(6,*) 'Label = ',Label,'  iSyLbl = ',iSyLbl
  call SysAbendMsg(SecNam,'I/O error in RdOne',' ')
end if

myStorage = Storage
call UpCase(myStorage)
if (myStorage == 'TRI') then
  call dCopy_(l_Tri,Work(ip_Scr),1,S,1)
else
  kTri = ip_Scr
  kSq = 1
  do iSym=1,nSym
    call Tri2Rec(Work(kTri),S(kSq),nBas(iSym),Debug)
    kTri = kTri+nBas(iSym)*(nBas(iSym)+1)/2
    kSq = kSq+nBas(iSym)**2
  end do
end if

call GetMem('OvlpScr','Free','Real',ip_Scr,l_Scr)

end subroutine GetOvlp_Localisation
