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

subroutine Get_Prim_Density_Matrix(ip_D,nBas,ip_D_p,nPrim,TM)

use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ip_D, nBas, nPrim
integer(kind=iwp), intent(out) :: ip_D_p
real(kind=wp), intent(in) :: TM(nPrim,nBas)
#include "WrkSpc.fh"
integer(kind=iwp) :: i,j,k,l
real(kind=wp) :: Djl, TMij, TMkl, TmpDensity

do i=1,nBas
  do k=1,i-1
    Work(ip_D+i*(i-1)/2+k-1) = Work(ip_D+i*(i-1)/2+k-1)/Two
  end do
end do
call GetMem('D_p','ALLO','REAL',ip_D_p,nPrim*(nPrim+1)/2)
do i=1,nPrim
  do k=1,i
    TmpDensity = Zero
    do j=1,nBas
      do l=1,nBas
        TMij = TM(i,j)
        TMkl = TM(k,l)
        if (j < l) then
          Djl = Work(ip_D+l*(l-1)/2+j-1)
        else
          Djl = Work(ip_D+j*(j-1)/2+l-1)
        end if
        TmpDensity = TmpDensity+TMij*TMkl*Djl
      end do
    end do
    Work(ip_D_p+i*(i-1)/2+k-1) = TmpDensity
  end do
end do

return

end subroutine Get_Prim_Density_Matrix
