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

subroutine JITO_Exchange_Int(MxR1,MxR2,imaxrank,n1,n2,JR,JI,HAM)
! this Subroutine calculates the anisotropic exchange interaction between
! two sites, of the one interacting pair on the basis of input ITO parameters

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: MxR1, MxR2, imaxrank(2), n1, n2
real(kind=wp), intent(in) :: JR(MxR1,-MxR1:MxR1,MxR2,-MxR2:MxR2), JI(MxR1,-MxR1:MxR1,MxR2,-MxR2:MxR2)
complex(kind=wp), intent(out) :: HAM(n1,n1,n2,n2)
integer(kind=iwp) :: ibuf, k1, k2, l1, l2, q1, q2
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: m1, m2
#endif
real(kind=wp) :: C01, C02, jpar
complex(kind=wp), allocatable :: J(:,:,:,:), O1(:,:), O2(:,:), OO(:,:,:,:), OW(:,:,:,:), W1(:,:), W2(:,:), WO(:,:,:,:), WW(:,:,:,:)
real(kind=wp), external :: dnrm2_

! ----  initial checks
if ((n1 <= 0) .or. (n2 <= 0)) return
HAM(:,:,:,:) = cZero
ibuf = MxR1*(2*MxR1+1)*MxR2*(2*MxR2+1)
if (ibuf == 0) return
jpar = dnrm2_(ibuf,JR,1)+dnrm2_(ibuf,JI,1)
if (jpar == Zero) return
! ---- end initial checks
call mma_allocate(J,[1,MxR1],[-MxR1,MxR1],[1,MxR2],[-MxR2,MxR2],label='J')
J(:,:,:,:) = cZero
do k1=1,MxR1,2
  do k2=1,MxR2,2
    J(k1,:,k2,:) = cmplx(JR(k1,:,k2,:),JI(k1,:,k2,:),kind=wp)
  end do
end do
!-----------------------------------------------------------------------
call mma_allocate(O1,n1,n1,'operator O1')
call mma_allocate(W1,n1,n1,'operator W1')
call mma_allocate(O2,n2,n2,'operator O2')
call mma_allocate(W2,n2,n2,'operator W2')
call mma_allocate(OO,n1,n1,n2,n2,'operator OO')
call mma_allocate(OW,n1,n1,n2,n2,'operator WO')
call mma_allocate(WO,n1,n1,n2,n2,'operator OW')
call mma_allocate(WW,n1,n1,n2,n2,'operator WW')
!-----------------------------------------------------------------------
do k1=1,imaxrank(1),2
  do q1=0,k1
    do k2=1,imaxrank(2),2
      do q2=0,k2
        ! generate the operator matrix K=ik, Q=iq, dimension=na
        call ITO(n1,k1,q1,C01,O1,W1)
        call ITO(n2,k2,q2,C02,O2,W2)
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
          HAM(:,:,:,:) = HAM(:,:,:,:)+J(k1,0,k2,0)*OO(:,:,:,:)
        else if ((q1 == 0) .and. (q2 /= 0)) then
          HAM(:,:,:,:) = HAM(:,:,:,:)+J(k1,0,k2,q2)*OO(:,:,:,:)+J(k1,0,k2,-q2)*OW(:,:,:,:)
        else if ((q1 /= 0) .and. (q2 == 0)) then
          HAM(:,:,:,:) = HAM(:,:,:,:)+J(k1,q1,k2,0)*OO(:,:,:,:)+J(k1,-q1,k2,0)*WO(:,:,:,:)
        else if ((q1 /= 0) .and. (q2 /= 0)) then
          HAM(:,:,:,:) = HAM(:,:,:,:)+J(k1,q1,k2,q2)*OO(:,:,:,:)+J(k1,q1,k2,-q2)*OW(:,:,:,:)+J(k1,-q1,k2,q2)*WO(:,:,:,:)+ &
                         J(k1,-q1,k2,-q2)*WW(:,:,:,:)
        end if
      end do
    end do
  end do
end do ! k1

#ifdef _DEBUGPRINT_
write(u6,'(A)') 'JITO_Exchange_Int: generated <m1,l1|HAM|m2,l2>'
do m1=1,n1
  do m2=1,n1
    do l1=1,n2
      do l2=1,n2
        write(u6,'(A,i2,A,i2,A,i2,A,i2,A,2ES22.14)') '<',m1,',',l1,'| HAM |',m2,',',l2,'> = ',HAM(m1,m2,l1,l2)
      end do
    end do
  end do
end do
#endif

call mma_deallocate(J)
call mma_deallocate(OO)
call mma_deallocate(OW)
call mma_deallocate(WO)
call mma_deallocate(WW)
call mma_deallocate(O1)
call mma_deallocate(W1)
call mma_deallocate(O2)
call mma_deallocate(W2)

return

end subroutine JITO_Exchange_Int
