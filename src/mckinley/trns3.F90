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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine Trns3(Win,Wout,na,nb,nvec,nc,Temp)
!***********************************************************************
!                                                                      *
! Object: utility routine to transform a AO batch in case of redun-    *
!         dancy of type aA=bB or cC=dD.                                *
!                                                                      *
! Called from: TwoEl                                                   *
!                                                                      *
! Calling    : Trns2                                                   *
!              DGeTMO  (ESSL)                                          *
!              DCopy   (ESSL)                                          *
!                                                                      *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             June '90                                                 *
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: na, nb, nvec, nc
real(kind=wp) :: Win(na,nb), Wout(nb,na), Temp(na,nb)

!iRout = 71
!iPrint = nPrint(iRout)
!write(u6,*) ' In Trns1: na, nb, nVec, nc=',na,nb,nvec,nc
!call RecPrt(' Win',' ',Win,na,nb)
if (nc == 1) then
  call dcopy_(nvec,Win,1,Wout,1)
  return
end if
if ((na == 1) .or. (nb == 1)) then
  call Trns2(Win,Wout,nvec,nc)
else
  call DGeTMO(Win,na,na,nb,Wout,nb)
  !call RecPrt(' After first DGeTMO',' ',Wout,nb,na)
  call Trns2(Wout,Temp,nvec,nc)
  !call RecPrt(' After Trns2',' ',Temp,nb,na)
  call DGeTMO(Temp,nb,nb,na,Wout,na)
  !call RecPrt(' After second DGeTMO',' ',Wout,na,nb)
end if

return

end subroutine Trns3
