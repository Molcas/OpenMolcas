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

subroutine W(R,iNQ,Weights,nNQ,nGrid,nRemoved)

use NQ_Structure, only: NQ_Data
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Three, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: iNQ, nNQ, nGrid
real(kind=wp), intent(inout) :: R(3,nGrid), Weights(nGrid)
integer(kind=iwp), intent(out) :: nRemoved
integer(kind=iwp) :: iGrid, jGrid, kNQ, lNQ
real(kind=wp) :: r_k, R_kl, r_l, rMU_kl, s, Sum_P_k, xdiff
real(kind=wp), allocatable :: P(:)
real(kind=wp), parameter :: Thrs = 1.0e-14_wp

!                                                                      *
!***********************************************************************
!                                                                      *
! iNQ is the index of the current atomic grid to which these grid points belong.

call mma_allocate(P,nNQ,label='P')

#ifdef _DEBUGPRINT_
write(u6,*) 'iNQ=',iNQ
write(u6,*) 'nGrid=',nGrid
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
jGrid = 0
nRemoved = 0
do iGrid=1,nGrid
# ifdef _DEBUGPRINT_
  write(u6,*) 'iGrid=',iGrid
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Becke's partitioning

  P(:) = One
  do kNQ=1,nNQ
    r_k = sqrt((R(1,iGrid)-NQ_Data(kNQ)%Coor(1))**2+(R(2,iGrid)-NQ_Data(kNQ)%Coor(2))**2+(R(3,iGrid)-NQ_Data(kNQ)%Coor(3))**2)
    do lNQ=1,kNQ-1
      r_l = sqrt((R(1,iGrid)-NQ_Data(lNQ)%Coor(1))**2+(R(2,iGrid)-NQ_Data(lNQ)%Coor(2))**2+(R(3,iGrid)-NQ_Data(lNQ)%Coor(3))**2)
      R_kl = sqrt((NQ_Data(kNQ)%Coor(1)-NQ_Data(lNQ)%Coor(1))**2+(NQ_Data(kNQ)%Coor(2)-NQ_Data(lNQ)%Coor(2))**2+ &
                  (NQ_Data(kNQ)%Coor(3)-NQ_Data(lNQ)%Coor(3))**2)
      rMU_kl = (r_k-r_l)/R_kl
      ! for abs(mu) > 0.986, s is 0 or 1 within ~1.0e-14
      if (abs(rMU_kl) <= 0.986_wp) then
        if (rMU_kl <= Half) then
          ! p(x) = 3/2*x-1/2*x**3
          ! p_i = p(...(p(mu)))
          xdiff = rMU_kl
          xdiff = (xdiff*Half)*(Three-xdiff**2) ! p_1
          xdiff = (xdiff*Half)*(Three-xdiff**2) ! p_2
          xdiff = (xdiff*Half)*(Three-xdiff**2) ! p_3
          s = Half*(One-xdiff)
        else
          ! q(x) = -3/2*x**2-1/2*x**3
          ! q_i = q(...(q(mu-1)))
          ! p_i = 1+q_i
          xdiff = rMU_kl-One
          xdiff = -(Three+xdiff)*(Half*xdiff**2) ! q_1
          xdiff = -(Three+xdiff)*(Half*xdiff**2) ! q_2
          xdiff = -(Three+xdiff)*(Half*xdiff**2) ! q_3
          s = -Half*xdiff
        end if
        P(kNQ) = P(kNQ)*s
        P(lNQ) = P(lNQ)*(One-s)
      else if (rMU_kl > Zero) then
        ! s = Zero
        P(kNQ) = Zero
      else   ! rMU_kl < Zero
        ! s = One
        P(lNQ) = Zero
      end if
    end do

  end do
  Sum_P_k = sum(P(:))
  Weights(iGrid) = Weights(iGrid)*P(iNQ)/Sum_P_k
  if (Weights(iGrid) >= Thrs) then
    jGrid = jGrid+1
    if (jGrid /= iGrid) then
      Weights(jGrid) = Weights(iGrid)
      R(1,jGrid) = R(1,iGrid)
      R(2,jGrid) = R(2,iGrid)
      R(3,jGrid) = R(3,iGrid)
    end if
  else
    nRemoved = nRemoved+1
  end if
# ifdef _DEBUGPRINT_
  write(u6,*) 'P_A,Z,Weights=',P(iNQ),Sum_P_k,Weights(jGrid)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end do
#ifdef _DEBUGPRINT_
write(u6,*) 'nRemoved=',nRemoved
#endif
call mma_deallocate(P)
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine W
