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

subroutine Hss_q()

use Slapaf_Info, only: Analytic_Hessian, Curvilinear, Degen, dqInt, iRef, lOld, nDimBC, Smmtrc
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: i, iAtom, ix, ixyz, nAtom, nQQ
real(kind=wp) :: rDum(1)
real(kind=wp), allocatable :: Hss_X(:), Degen2(:), Hss_Q_(:), KtB(:)

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
if (lOld) return

nQQ = size(dqInt,1)
nAtom = size(Degen,2)
!                                                                      *
!***********************************************************************
!                                                                      *
! Back-transform from cartesian to internals
!
! dQ/dx d^2E/dQ^2 dQ/dx + d^2Q/dx^2 dE/dQ = d^2E/dx^2

! Pickup d^2E/dx^2

call mma_allocate(Hss_x,nDimBC**2,Label='Hss_X')
call Get_dArray('Hss_X',Hss_x,nDimBC**2)
call mma_allocate(KtB,nDimBC*nQQ,Label='KtB')
call Get_dArray('KtB',KtB,nDimBC*nQQ)
#ifdef _DEBUGPRINT_
call RecPrt('Hss_x',' ',Hss_X,nDimBC,nDimBC)
#endif

call mma_allocate(Degen2,nDimBC,Label='Degen2')
i = 0
do ix=1,3*nAtom
  iAtom = (ix+2)/3
  ixyz = ix-(iAtom-1)*3
  if (Smmtrc(ixyz,iAtom)) then
    i = i+1
    Degen2(i) = Degen(ixyz,iAtom)
  end if
end do
#ifdef _DEBUGPRINT_
call RecPrt('Degen2',' ',Degen2,nDimBC,1)
#endif

if (Analytic_Hessian .and. Curvilinear) then

  ! Form u^(1/2) (Sum(i) d^2Q_i/dx^2 * dE/dQ_i) u^(1/2)
  !
  ! and form d^2E/dx^2 - d^2Q/dx^2 dE/dQ

  call dBuu(Degen2,nQQ,nDimBC,dqInt(:,iRef),Hss_X,.false.)
# ifdef _DEBUGPRINT_
  call RecPrt('H(X)-BtgQ',' ',Hss_X,nDimBC,nDimBC)
# endif
end if

call mma_allocate(Hss_Q_,nQQ**2,Label='Hss_Q_')
call Hess_Tra(Hss_X,nDimBC,Degen2,KtB,nQQ,Hss_Q_)

call Put_dArray('Hss_Q',Hss_Q_,nQQ**2)
call Put_dArray('Hss_upd',rDum,0)
#ifdef _DEBUGPRINT_
call RecPrt('Hss_Q: Hessian',' ',Hss_Q_,nQQ,nQQ)
#endif
call mma_deallocate(Hss_Q_)
call mma_deallocate(KtB)
call mma_deallocate(Degen2)
call mma_deallocate(Hss_X)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Hss_q
