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

subroutine Put_NucAttr()

use Constants, only: One
use Definitions, only: iwp, u6

implicit none
#include "WrkSpc.fh"
#include "embpcharg.fh"
character(len=8) :: Label
integer(kind=iwp) :: i, iComp, iOpt, ipAttr, ipXFdInt, irc, iSyLbl, nLT, nLT_, nSym, nBas(8)

call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)

nLT = nBas(1)*(nBas(1)+1)/2
do i=2,nSym
  nLT = nLT+nBas(i)*(nBas(i)+1)/2
end do
nLT_ = nLT
if (DoEMPC) nLT = 2*nLT
call Getmem('tempAtr','Allo','Real',ipAttr,nLT)
ipXFdInt = ipAttr+nLT_

irc = -1
iOpt = 6
iComp = 1
iSyLbl = 1
Label = 'Attract '
call RdOne(irc,iOpt,Label,iComp,Work(ipAttr),iSyLbl)
if (irc /= 0) then
  write(u6,*) 'Put_NucAttr: RdOne returned ',irc
  write(u6,*) 'Label = ',Label,'  iSyLbl = ',iSyLbl
  call SysAbendMsg('Put_NucAttr','I/O error in RdOne',' ')
end if

if (DoEMPC) then
  irc = -1
  iOpt = 2
  iComp = 1
  iSyLbl = 1
  Label = 'XFdInt  '
  call RdOne(irc,iOpt,Label,iComp,Work(ipXFdInt),iSyLbl)
  if (irc /= 0) then
    write(u6,*) 'Put_NucAttr: RdOne returned ',irc
    write(u6,*) 'Label = ',Label,'  iSyLbl = ',iSyLbl
    call SysAbendMsg('Put_NucAttr','I/O error in RdOne',' ')
  end if
  call daxpy_(nLT_,One,Work(ipXFdInt),1,Work(ipAttr),1)
end if

call Put_dArray('Nuc Potential',Work(ipAttr),nLT_)

#ifdef _DEBUGPRINT_
iAttr = ipAttr
do i=1,nSym
  call TriPrt('Attr Inte','',Work(iAttr),nBas(i))
  iAttr = iAttr+nBas(i)*(nBas(i)+1)/2
end do
#endif

call Getmem('tempAtr','Free','Real',ipAttr,nLT)

return

end subroutine Put_NucAttr
