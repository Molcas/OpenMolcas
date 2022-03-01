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

subroutine unpackk_ic_1(i,vint,ndimv1,ndimv2,ndimv3,Vic,ndimvi)
! this routine vint(j,k,l) = <i,j|k,l>
! for given i from incore (nonreduced) expanded block Vic
!
! i      - value of pivot index (I)
! vint   - array of integrals (O)
! ndimv1 - first dimension of vint (norb(symj)) (I)
! ndimv2 - second dimension of vint (norb(symk)) (I)
! ndimv3 - third dimension of vint (norb(syml)) (I)
! Vic    - incore expanded block of integrals (I)
! ndimvi - first dimension of Vic norb(symi) (I)

#include "reorg.fh"
#include "SysDef.fh"
integer i, ndimv1, ndimv2, ndimv3, ndimvi
real*8 vint(1:ndimv1*ndimv2*ndimv3)
real*8 Vic(1:ndimvi,1:ndimv1*ndimv2*ndimv3)
! help variables
integer jkl

do jkl=1,ndimv1*ndimv2*ndimv3
  vint(jkl) = Vic(i,jkl)
end do

return

end subroutine unpackk_ic_1
