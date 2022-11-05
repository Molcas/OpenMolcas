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

use Index_Functions, only: nTri_Elem
use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
#include "embpcharg.fh"
character(len=8) :: Label
integer(kind=iwp) :: i, iComp, iOpt, irc, iSyLbl, nC, nLT, nSym, nBas(8)
real(kind=wp), allocatable :: Attr(:,:)

call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)

nLT = 0
do i=1,nSym
  nLT = nLT+nTri_Elem(nBas(i))
end do
nC = 1
if (DoEMPC) nC = 2
call mma_allocate(Attr,nLT,nC,label='tempAtr')

irc = -1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
iComp = 1
iSyLbl = 1
Label = 'Attract '
call RdOne(irc,iOpt,Label,iComp,Attr(:,1),iSyLbl)
if (irc /= 0) then
  write(u6,*) 'Put_NucAttr: RdOne returned ',irc
  write(u6,*) 'Label = ',Label,'  iSyLbl = ',iSyLbl
  call SysAbendMsg('Put_NucAttr','I/O error in RdOne',' ')
end if

if (DoEMPC) then
  irc = -1
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  iComp = 1
  iSyLbl = 1
  Label = 'XFdInt  '
  call RdOne(irc,iOpt,Label,iComp,Attr(:,2),iSyLbl)
  if (irc /= 0) then
    write(u6,*) 'Put_NucAttr: RdOne returned ',irc
    write(u6,*) 'Label = ',Label,'  iSyLbl = ',iSyLbl
    call SysAbendMsg('Put_NucAttr','I/O error in RdOne',' ')
  end if
  Attr(:,1) = Attr(:,1)+Attr(:,2)
end if

call Put_dArray('Nuc Potential',Attr(:,1),nLT)

#ifdef _DEBUGPRINT_
iAttr = 1
do i=1,nSym
  call TriPrt('Attr Inte','',Attr(iAttr,1),nBas(i))
  iAttr = iAttr+nTri_Elem(nBas(i))
end do
#endif

call mma_deallocate(Attr)

return

end subroutine Put_NucAttr
