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

subroutine W(R,ilist_p,Weights,list_p,nlist_p,nGrid,nRemoved)

use NQ_Structure, only: NQ_Data
use Constants, only: Zero, One, Three, Half, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ilist_p, nlist_p, list_p(nlist_p), nGrid
real(kind=wp), intent(inout) :: R(3,nGrid), Weights(nGrid)
integer(kind=iwp), intent(out) :: nRemoved
integer(kind=iwp) :: iGrid, iNQ, jGrid, klist_p, kNQ, llist_p, lNQ
real(kind=wp) :: p1, p2, p3, P_i, P_k, r_k, R_kl, r_l, rMU_kl, s, Sum_P_k, xdiff
real(kind=wp), parameter :: Thrs = 1.0e-14_wp

!                                                                      *
!***********************************************************************
!                                                                      *
P_i = Zero ! dummy initialize

! iNQ is the index of the current atomic grid to which these grid
! points belong.

iNQ = list_p(ilist_p)
!write(u6,*) 'ilist_p=',ilist_p
!write(u6,*) 'nlist_p=',nlist_p
!write(u6,*) 'nGrid=',nGrid
!write(u6,*) 'iNQ=',iNQ
!                                                                      *
!***********************************************************************
!                                                                      *
jGrid = 0
nRemoved = 0
do iGrid=1,nGrid
  !write(u6,*) 'iGrid=',iGrid
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Becke's partitioning

  Sum_P_k = Zero
  do klist_p=1,nlist_p
    kNQ = list_p(klist_p)
    r_k = sqrt((R(1,iGrid)-NQ_Data(kNQ)%Coor(1))**2+(R(2,iGrid)-NQ_Data(kNQ)%Coor(2))**2+(R(3,iGrid)-NQ_Data(kNQ)%Coor(3))**2)
    P_k = One
    do llist_p=1,nlist_p
      lNQ = list_p(llist_p)

      if (kNQ /= lNQ) then

        r_l = sqrt((R(1,iGrid)-NQ_Data(lNQ)%Coor(1))**2+(R(2,iGrid)-NQ_Data(lNQ)%Coor(2))**2+(R(3,iGrid)-NQ_Data(lNQ)%Coor(3))**2)
        R_kl = sqrt((NQ_Data(kNQ)%Coor(1)-NQ_Data(lNQ)%Coor(1))**2+(NQ_Data(kNQ)%Coor(2)-NQ_Data(lNQ)%Coor(2))**2+ &
                    (NQ_Data(kNQ)%Coor(3)-NQ_Data(lNQ)%Coor(3))**2)
        rMU_kl = (r_k-r_l)/R_kl
        if (rMU_kl <= Half) then
          p1 = (rMU_kl*Half)*(Three-rMU_kl**2)
          p2 = (p1*Half)*(Three-p1**2)
          p3 = (p2*Half)*(Three-p2**2)
          s = Half*(One-p3)
        else
          xdiff = rMU_kl-One
          xdiff = (-OneHalf-Half*xdiff)*xdiff**2
          xdiff = (-OneHalf-Half*xdiff)*xdiff**2
          p3 = (OneHalf+Half*xdiff)*xdiff**2
          s = Half*p3
        end if
        P_k = P_k*s
      end if
    end do

    if (kNQ == iNQ) P_i = P_k
    Sum_P_k = Sum_P_k+P_k
  end do
  Weights(iGrid) = Weights(iGrid)*P_i/Sum_P_k
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
  !write(u6,*) 'P_A,Z,Weights=',P_i,Sum_P_k,Weights(jGrid)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end do
!write(u6,*) 'nRemoved=',nRemoved
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine W
