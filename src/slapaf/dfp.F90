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

subroutine DFP(B,nDim,Bd,Delta,rGamma)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nDim
real(kind=wp) :: B(nDim,nDim), Bd(nDim), Delta(nDim), rGamma(nDim)
integer(kind=iwp) :: i, j
real(kind=wp) :: dBd, gd
real(kind=wp), parameter :: Thr = 1.0e-8_wp
real(kind=wp), external :: DDot_

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('DFP: B',' ',B,nDim,nDim)
!call RecPrt('DFP: Bd',' ',Bd,1,nDim)
call RecPrt('DFP: Gamma',' ',rGamma,1,nDim)
call RecPrt('DFP: Delta',' ',Delta,1,nDim)
#endif
call DGEMM_('N','N',nDim,1,nDim,One,B,nDim,Delta,nDim,Zero,Bd,nDim)
gd = DDot_(nDim,rGamma,1,Delta,1)
dBd = DDot_(nDim,Delta,1,Bd,1)
#ifdef _DEBUGPRINT_
call RecPrt('DFP: Bd',' ',Bd,1,nDim)
write(u6,*) 'gd=',gd
write(u6,*) 'dBd=',dBd
write(u6,*) 'Thr=',Thr
#endif
if (gd < Zero) then
  call MSP(B,rGamma,Delta,nDim)
else
  do i=1,nDim
    do j=1,nDim
      B(i,j) = B(i,j)+(rGamma(i)*rGamma(j))/gd-(Bd(i)*Bd(j))/dBd
    end do
  end do
end if

#ifdef _DEBUGPRINT_
call RecPrt('DFP: B',' ',B,nDim,nDim)
#endif

return

end subroutine DFP
