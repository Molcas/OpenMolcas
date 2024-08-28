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

!#define _DEBUGPRINT_
subroutine binte(k,alfa,beta,r0,a,ggrin,nz)

use Constants, only: Zero
use welcom, only: binom, fiint, kmax

implicit none
integer k, nz
real*8 grint(0:kmax,kmax), ggrin(nz,0:k,k/2+1,k/4+1), alfa(nz), a(nz)
real*8 beta, r0
integer iz, i, j, l, ind, k2, j2
real*8 tal

#ifdef _DEBUGPRINT_
call RecPrt(' In Binte: Alfa',' ',alfa,nz,1)
call RecPrt(' In Binte: A   ',' ',a,nz,1)
#endif
call dcopy_(nz*(k+1)*(k/2+1)*(k/4+1),[Zero],0,ggrin,1)
do iz=1,nz
  call dcopy_((kMax+1)*kMax,[Zero],0,grint,1)
  call rrint(k,alfa(iz),a(iz),beta,r0,grint,kmax)
  do i=0,k
    do j=0,i,2
      j2 = j/2+1
      ggrin(iz,i,j2,1) = Zero
      do l=j,i
        if (i-l == 0) then
          ggrin(iz,i,j2,1) = ggrin(iz,i,j2,1)+grint(l,j2)*binom(i-j,l-j)
        else
          ggrin(iz,i,j2,1) = ggrin(iz,i,j2,1)+grint(l,j2)*a(iz)**(i-l)*binom(i-j,l-j)
        end if
      end do
      ind = 1
      do k2=2,j2-1,2
        tal = fiint((j-k2)/2,k2/2)/fiint(j/2,0)
        ind = ind+1
        ggrin(iz,i,j2,ind) = ggrin(iz,i,j2,1)*tal
      end do
    end do
  end do
end do

#ifdef _DEBUGPRINT_
call RecPrt(' In Binte: Ggrin',' ',Ggrin,nz,(k+1)*(k/2+1)*(k/4+1))
#endif

return

end subroutine binte
