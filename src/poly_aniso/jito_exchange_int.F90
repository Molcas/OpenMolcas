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

use Constants, only: Zero, cZero
use Definitions, only: wp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
#include "stdalloc.fh"
! input variables
integer, intent(in) :: imaxrank(2)
integer, intent(in) :: MxR1, MxR2
integer, intent(in) :: n1, n2
real(kind=8), intent(in) :: JR(MxR1,-MxR1:MxR1,MxR2,-MxR2:MxR2)
real(kind=8), intent(in) :: JI(MxR1,-MxR1:MxR1,MxR2,-MxR2:MxR2)
! output variables
complex(kind=8), intent(out) :: HAM(n1,n1,n2,n2)
! lcoal variables:
integer :: ibuf, k1, q1, k2, q2, m1, m2, l1, l2
real(kind=8) :: jpar, RR, RI
!real(kind=8) :: rK1, rK2, rQ1, rQ2, rM1, rM2, rJ1, rJ2, CGp1, CGp2, CGm1, CGm2, CG01, CG02
real(kind=8) :: C01, C02
complex(kind=8) :: J(MxR1,-MxR1:MxR1,MxR2,-MxR2:MxR2)
complex(kind=8), allocatable :: O1(:,:), O2(:,:)
complex(kind=8), allocatable :: W1(:,:), W2(:,:)
complex(kind=8), allocatable :: OO(:,:,:,:), WW(:,:,:,:), OW(:,:,:,:), WO(:,:,:,:)
real(kind=8) :: dnrm2_
external :: dnrm2_

! ----  initial checks
if ((n1 <= 0) .or. (n2 <= 0)) return
call zcopy_(n1*n1*n2*n2,[cZero],0,HAM,1)
ibuf = 0
ibuf = MxR1*(2*MxR1+1)*MxR2*(2*MxR2+1)
if (ibuf == 0) return
jpar = dnrm2_(ibuf,JR(1:MxR1,-MxR1:MxR1,1:MxR2,-MxR2:MxR2),1)+dnrm2_(ibuf,JI(1:MxR1,-MxR1:MxR1,1:MxR2,-MxR2:MxR2),1)
if (jpar == Zero) return
! ---- end initial checks
call zcopy_(ibuf,[cZero],0,J(1:MxR1,-MxR1:MxR1,1:MxR2,-MxR2:MxR2),1)
do k1=1,MxR1,2
  do k2=1,MxR2,2
    do q1=-k1,k1
      do q2=-k2,k2
        RR = JR(k1,q1,k2,q2)
        RI = JI(k1,q1,k2,q2)
        J(k1,q1,k2,q2) = cmplx(RR,RI,kind=wp)
      end do
    end do
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
        !generate coupled operators:
        call zcopy_(n1*n1*n2*n2,[cZero],0,OO,1)
        call zcopy_(n1*n1*n2*n2,[cZero],0,OW,1)
        call zcopy_(n1*n1*n2*n2,[cZero],0,WO,1)
        call zcopy_(n1*n1*n2*n2,[cZero],0,WW,1)
        do m1=1,n1
          do m2=1,n1
            do l1=1,n2
              do l2=1,n2
                OO(m1,m2,l1,l2) = O1(m1,m2)*O2(l1,l2)
                OW(m1,m2,l1,l2) = O1(m1,m2)*W2(l1,l2)
                WO(m1,m2,l1,l2) = W1(m1,m2)*O2(l1,l2)
                WW(m1,m2,l1,l2) = W1(m1,m2)*W2(l1,l2)
              end do
            end do
          end do
        end do !m1
        ! compute the exchange Hamiltonian:
        if ((q1 == 0) .and. (q2 == 0)) then
          call zaxpy_(n1*n1*n2*n2,J(k1,0,k2,0),OO,1,HAM,1)
        else if ((q1 == 0) .and. (q2 /= 0)) then
          call zaxpy_(n1*n1*n2*n2,J(k1,0,k2,q2),OO,1,HAM,1)
          call zaxpy_(n1*n1*n2*n2,J(k1,0,k2,-q2),OW,1,HAM,1)
        else if ((q1 /= 0) .and. (q2 == 0)) then
          call zaxpy_(n1*n1*n2*n2,J(k1,q1,k2,0),OO,1,HAM,1)
          call zaxpy_(n1*n1*n2*n2,J(k1,-q1,k2,0),WO,1,HAM,1)
        else if ((q1 /= 0) .and. (q2 /= 0)) then
          call zaxpy_(n1*n1*n2*n2,J(k1,q1,k2,q2),OO,1,HAM,1)
          call zaxpy_(n1*n1*n2*n2,J(k1,q1,k2,-q2),OW,1,HAM,1)
          call zaxpy_(n1*n1*n2*n2,J(k1,-q1,k2,q2),WO,1,HAM,1)
          call zaxpy_(n1*n1*n2*n2,J(k1,-q1,k2,-q2),WW,1,HAM,1)
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
