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

subroutine Integral_WrOut_Cho( &
#                             define _CALLING_
#                             include "int_wrout_interface.fh"
                             )
! calls the proper routines IndSft/PLF
! if IntOrd_jikl == .true. integral order within symblk: jikl
!                   else   integral order within symblk: ijkl

use Cholesky, only: IfcSew, nSym
use Definitions, only: wp, iwp, u6

implicit none
#include "int_wrout_interface.fh"
integer(kind=iwp) :: iAO(4), iAOst(4), iBas, iCmp(4), iShell(4), jBas, kBas, kOp(4), lBas
logical(kind=iwp) :: Shijij
character(len=*), parameter :: SecNam = 'Integral_WrOut_Cho'

#include "macros.fh"
unused_var(iSOSym)
unused_var(mSym)

iCmp(:) = iSD4(2,:)
iShell(:) = iSD4(11,:)
iAO(:) = iSD4(7,:)
iAOst(:) = iSD4(8,:)
iBas = iSD4(19,1)
jBas = iSD4(19,2)
kBas = iSD4(19,3)
lBas = iSD4(19,4)
Shijij = (iSD4(0,1) == iSD4(0,3)) .and. (iSD4(10,1) == iSD4(10,3)) .and. (iSD4(0,2) == iSD4(0,4)) .and. (iSD4(10,2) == iSD4(10,4))

! call sorting routine

kOp(:) = 0
if (IfcSew == 1) then
  if (nSym == 1) then
    call PLF_Cho(TInt,nTInt,AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),iAO,iAOst,iBas,jBas,kBas,lBas,kOp)
  else
    call IndSft_Cho(TInt,nTInt,iCmp,iShell,iBas,jBas,kBas,lBas,Shijij,iAO,iAOst,ijkl,SOInt,nSOint)
  end if
else if (IfcSew == 2) then
  if (nSym == 1) then
    call PLF_Cho_2(TInt,nTInt,AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),iAO,iAOst,iBas,jBas,kBas,lBas,kOp)
  else
    call IndSft_Cho_2(TInt,nTInt,iCmp,iShell,iBas,jBas,kBas,lBas,Shijij,iAO,iAOst,ijkl,SOInt,nSOint)
  end if
else if (IfcSew == 3) then
  if (nSym == 1) then
    call PLF_Cho_3(TInt,nTInt,AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),iAO,iAOst,iBas,jBas,kBas,lBas,kOp)
  else
    call IndSft_Cho_3(TInt,nTInt,iCmp,iShell,iBas,jBas,kBas,lBas,Shijij,iAO,iAOst,ijkl,SOInt,nSOint)
  end if
else
  write(u6,*)
  write(u6,*)
  write(u6,*) '!!!!!!!!!! IfcSew=',IfcSew,' !!!!!!!!!!'
  call Cho_Quit('IfcSew out of bounds in '//SecNam,105)
end if

end subroutine Integral_WrOut_Cho
