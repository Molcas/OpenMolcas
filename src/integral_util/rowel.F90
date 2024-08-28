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
! Copyright (C) 1992, Gunnar Karlstrom                                 *
!***********************************************************************

subroutine Rowel(nZeta,r0,Beta,K,alpha,P,a,gri,grin,jsum)
!***********************************************************************
! 1992                                                                 *
! Gunnar Karlstrom                                                     *
! Department of Theoretical Chemistry                                  *
! University of Lund                                                   *
! Lund                                                                 *
! Sweden                                                               *
!***********************************************************************

use Constants, only: One
use welcom, only: Fac, ipot3
use Definitions, only: wp

implicit none
integer nZeta, jSum, k
real*8 gri(nZeta*jsum), alpha(nZeta), a(nZeta), grin((k+1)*(k/2+1)*(k/4+1)*nZeta), P(nZeta,3)
real*8 r0, Beta
integer iZeta, i, iSum

call poti(k,ipot3)
!call IecPrt(' ipot3(0:k+1)',ipot3,k+2,1)
isum = ipot3(k+1)
call bino(k+6)
call fiin(k+1)
call tetin(k+1)
call ylmnor(k+1)

do iZeta=1,nZeta
  a(iZeta) = sqrt(P(iZeta,1)**2+P(iZeta,2)**2+P(iZeta,3)**2)
end do
!call RecPrt(' In Rowel: Distances',' ',a,nZeta,1)

fac(0) = One
do i=1,k+2
  fac(i) = fac(i-1)*real(i,kind=wp)
end do
call priwel(k,alpha,beta,r0,a,gri,nZeta,isum,grin)
!call RecPrt('Internal well integrals',' ',gri,nZeta,isum)

return

end subroutine Rowel
