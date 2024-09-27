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

subroutine recover_exch_HAM_from_Stevens_ESO(n1,n2,S,HAM)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: cZero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n1, n2
complex(kind=wp), intent(in) :: S(n1-1,-(n1-1):n1-1,n2-1,-(n2-1):n2-1)
complex(kind=wp), intent(out) :: HAM(n1,n1,n2,n2)
integer(kind=iwp) :: k1, k2, l1, l2, q1, q2
complex(kind=wp) :: redME1, redME2
complex(kind=wp), allocatable :: O1(:,:), O2(:,:), OO(:,:,:,:), OW(:,:,:,:), W1(:,:), W2(:,:), WO(:,:,:,:), WW(:,:,:,:)

!-----------------------------------------------------------------------
! recover the original HAMILTONIAN using the S parameters
!=======================================================================
call mma_allocate(O1,n1,n1,'operator O1')
call mma_allocate(O2,n2,n2,'operator O2')
call mma_allocate(W1,n1,n1,'operator W1')
call mma_allocate(W2,n2,n2,'operator W2')
call mma_allocate(OO,n1,n1,n2,n2,'operator OO')
call mma_allocate(OW,n1,n1,n2,n2,'operator WO')
call mma_allocate(WO,n1,n1,n2,n2,'operator OW')
call mma_allocate(WW,n1,n1,n2,n2,'operator WW')
HAM(:,:,:,:) = cZero
do k1=1,n1-1
  do q1=0,k1
    do k2=1,n2-1
      do q2=0,k2
        ! generate the operator matrix K=ik, Q=iq, dimension=na
        call ESO(n1,k1,q1,O1,W1,redME1)
        call ESO(n2,k2,q2,O2,W2,redME2)
        ! generate coupled operators:
        do l1=1,n2
          do l2=1,n2
            OO(:,:,l1,l2) = O1(:,:)*O2(l1,l2)
            OW(:,:,l1,l2) = O1(:,:)*W2(l1,l2)
            WO(:,:,l1,l2) = W1(:,:)*O2(l1,l2)
            WW(:,:,l1,l2) = W1(:,:)*W2(l1,l2)
          end do
        end do
        ! compute the exchange Hamiltonian:
        if ((q1 == 0) .and. (q2 == 0)) then
          HAM(:,:,:,:) = HAM(:,:,:,:)+S(k1,0,k2,0)*OO(:,:,:,:)
        else if ((q1 == 0) .and. (q2 /= 0)) then
          HAM(:,:,:,:) = HAM(:,:,:,:)+S(k1,0,k2,q2)*OO(:,:,:,:)+S(k1,0,k2,-q2)*OW(:,:,:,:)
        else if ((q1 /= 0) .and. (q2 == 0)) then
          HAM(:,:,:,:) = HAM(:,:,:,:)+S(k1,q1,k2,0)*OO(:,:,:,:)+S(k1,-q1,k2,0)*WO(:,:,:,:)
        else if ((q1 /= 0) .and. (q2 /= 0)) then
          HAM(:,:,:,:) = HAM(:,:,:,:)+S(k1,q1,k2,q2)*OO(:,:,:,:)+S(k1,q1,k2,-q2)*OW(:,:,:,:)+S(k1,-q1,k2,q2)*WO(:,:,:,:)+ &
                         S(k1,-q1,k2,-q2)*WW(:,:,:,:)
        end if
      end do !q
    end do !k
  end do !q
end do !k

call mma_deallocate(O1)
call mma_deallocate(O2)
call mma_deallocate(W1)
call mma_deallocate(W2)
call mma_deallocate(OO)
call mma_deallocate(OW)
call mma_deallocate(WO)
call mma_deallocate(WW)

return

end subroutine recover_exch_HAM_from_Stevens_ESO
