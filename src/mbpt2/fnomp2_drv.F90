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
! Copyright (C) 2008, Francesco Aquilante                              *
!***********************************************************************

subroutine FNOMP2_Drv(irc,EMP2,CMOI,EOcc,EVir)

use MBPT2_Global, only: nBas
use ChoMP2, only: ChoAlg, DoDens, DoFNO, DoMP2, vkept, XEMP2
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(out) :: EMP2
real(kind=wp), intent(inout) :: CMOI(*), EOcc(*), EVir(*)
integer(kind=iwp) :: ChoAlg_
logical(kind=iwp) :: DoDens_
real(kind=wp) :: Dum(2)
#include "corbinf.fh"

DoDens_ = DoDens
DoDens = .false.
ChoAlg_ = ChoAlg
ChoAlg = 2

call FNO_MP2(irc,nSym,nBas,nFro,nOcc,nExt,nDel,CMOI,EOcc,EVir,vkept,DoMP2,XEMP2)
if (irc /= 0) then
  write(u6,*) 'FNO_MP2 returned ',irc
  call SysAbendMsg('FNO_MP2','Non-zero return code from FNO_MP2',' ')
end if

ChoAlg = ChoAlg_
DoDens = DoDens_
DoFNO = .false.
call ChoMP2_Drv(irc,EMP2,CMOI,EOcc,EVir,Dum(1),Dum(2))
EMP2 = EMP2+XEMP2

return

end subroutine FNOMP2_Drv
