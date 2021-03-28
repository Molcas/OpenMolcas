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

subroutine Gen_Prim_Density_Matrix(nBas,nPrim,ip_D_p,nOcOb,oNum,oCof)

implicit real*8(a-h,o-z)
#include "WrkSpc.fh"
dimension oNum(nBas)
dimension oCof(nBas,nPrim)

call GetMem('D_p','ALLO','REAL',ip_D_p,nPrim*(nPrim+1)/2)
do K=1,nPrim
  do L=1,K
    Work(ip_D_p+k*(k-1)/2+l-1) = 0.0d0
  end do
end do
do K=1,nPrim
  do L=1,K
    do i=1,nOcOb
      Work(ip_D_p+k*(k-1)/2+l-1) = Work(ip_D_p+k*(k-1)/2+l-1)+oNum(I)*oCof(I,K)*oCof(I,L)
    end do
  end do
end do

return

end subroutine Gen_Prim_Density_Matrix
