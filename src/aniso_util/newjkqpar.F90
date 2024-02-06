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

subroutine newjkqpar(n1,n2,H,J,B,S)

use Constants, only: Zero, One, cZero, cOne, Onei
use Definitions, only: wp, u6

implicit none
#include "stdalloc.fh"
integer, intent(in) :: n1, n2
complex(kind=8), intent(in) :: H(n1,n1,n2,n2)
complex(kind=8), intent(out) :: J(n1-1,-(n1-1):n1-1,n2-1,-(n2-1):n2-1)
complex(kind=8), intent(out) :: B(n1-1,-(n1-1):n1-1,n2-1,-(n2-1):n2-1)
complex(kind=8), intent(out) :: S(n1-1,-(n1-1):n1-1,n2-1,-(n2-1):n2-1)
! local variables:
integer :: k1, k2, q1, q2, m1, m2, m12, i1, i2, j1, j2
real(kind=8) :: cr1, cr2, C01, C02, r1, r2, F1, F2
complex(kind=8) :: cf1, cf2, c1, c2, c12, ci, cc1, cc2, trace_exch2
complex(kind=8), allocatable :: O1(:,:), W1(:,:)
complex(kind=8), allocatable :: O2(:,:), W2(:,:)
complex(kind=8), allocatable :: OO(:,:,:,:), WO(:,:,:,:)
complex(kind=8), allocatable :: OW(:,:,:,:), WW(:,:,:,:)
complex(kind=8), allocatable :: HAM(:,:,:,:)
real(kind=8) :: knm(12,0:12), dznrm2_
external :: dznrm2_, trace_exch2
logical :: dbg

!-------------------------------------------
if ((n1 < 1) .or. (n2 < 1)) return
!-------------------------------------------
dbg = .false.
call mma_allocate(O1,n1,n1,'operator O1')
call mma_allocate(W1,n1,n1,'operator W1')
call mma_allocate(O2,n2,n2,'operator O2')
call mma_allocate(W2,n2,n2,'operator W2')
call set_knm(knm)
!-------------------------------------------
J(:,:,:,:) = cZero
B(:,:,:,:) = cZero
S(:,:,:,:) = cZero
do k1=1,n1-1
  do q1=0,k1
    do k2=1,n2-1
      do q2=0,k2
        cr1 = Zero
        cr2 = Zero
        m1 = 0
        m2 = 0
        r1 = Zero
        r2 = Zero
        cf1 = cZero
        cf2 = cZero

        m1 = (-1)**q1
        m2 = (-1)**q2
        m12 = (-1)**(q1+q2)
        c1 = m1*cOne
        c2 = m2*cOne
        c12 = m12*cOne
        ci = -Onei

        ! generate the operators for site 1 and 2:
        call ITO(n1,k1,q1,C01,O1,W1)
        call ITO(n2,k2,q2,C02,O2,W2)

        ! compute the scaling factors for various Operators
        call coeff_redus_sub(n1,k1,cr1)
        call coeff_redus_sub(n2,k2,cr2)
        cc1 = One/(cr1*C01)*cOne
        cc2 = One/(cr2*C02)*cOne
        r1 = C01*C01*real(2*k1+1,kind=wp)/real(n1,kind=wp)
        r2 = C02*C02*real(2*k2+1,kind=wp)/real(n2,kind=wp)
        cf1 = c1*(r1*cOne)
        cf2 = c2*(r2*cOne)

        if (dbg) then
          ! use the old trace_exch function
          call mma_allocate(OO,n1,n1,n2,n2,'operator OO')
          call mma_allocate(WO,n1,n1,n2,n2,'operator WO')
          call mma_allocate(OW,n1,n1,n2,n2,'operator OW')
          call mma_allocate(WW,n1,n1,n2,n2,'operator WW')

          do i1=1,n1
            do j1=1,n1
              do i2=1,n2
                do j2=1,n2
                  OO(i1,j1,i2,j2) = O1(i1,j1)*O2(i2,j2)
                  OW(i1,j1,i2,j2) = O1(i1,j1)*W2(i2,j2)
                  WO(i1,j1,i2,j2) = W1(i1,j1)*O2(i2,j2)
                  WW(i1,j1,i2,j2) = W1(i1,j1)*W2(i2,j2)
                end do
              end do
            end do
          end do
          ! find the parameters:
          ! in the Naoya's operators
          J(k1,-q1,k2,-q2) = cf1*cf2*trace_exch2(n1,n2,H,O1,O2)
          J(k1,-q1,k2,q2) = cf1*cf2*trace_exch2(n1,n2,H,O1,W2)
          J(k1,q1,k2,-q2) = cf1*cf2*trace_exch2(n1,n2,H,W1,O2)
          J(k1,q1,k2,q2) = cf1*cf2*trace_exch2(n1,n2,H,W1,W2)

          call mma_deallocate(OO)
          call mma_deallocate(WO)
          call mma_deallocate(OW)
          call mma_deallocate(WW)
        else
          ! find the parameters:
          ! in the Naoya's operators
          J(k1,-q1,k2,-q2) = cf1*cf2*trace_exch2(n1,n2,H,O1,O2)
          J(k1,-q1,k2,q2) = cf1*cf2*trace_exch2(n1,n2,H,O1,W2)
          J(k1,q1,k2,-q2) = cf1*cf2*trace_exch2(n1,n2,H,W1,O2)
          J(k1,q1,k2,q2) = cf1*cf2*trace_exch2(n1,n2,H,W1,W2)
        end if

        ! in the Liviu operators
        B(k1,-q1,k2,-q2) = J(k1,-q1,k2,-q2)*cc1*cc2
        B(k1,-q1,k2,q2) = J(k1,-q1,k2,q2)*cc1*cc2
        B(k1,q1,k2,-q2) = J(k1,q1,k2,-q2)*cc1*cc2
        B(k1,q1,k2,q2) = J(k1,q1,k2,q2)*cc1*cc2

        ! generate real parameters for Extended Stevens Operators formalism:
        if ((q1 > 0) .and. (q2 > 0)) then
          ! BB = (Jmm + (-1)^q2 Jmp + (-1)^q1 Jpm + (-1)^(q1 + q2) Jpp);
          S(k1,q1,k2,q2) = B(k1,-q1,k2,-q2)+c2*B(k1,-q1,k2,q2)+c1*B(k1,q1,k2,-q2)+c12*B(k1,q1,k2,q2)
          ! BC = (Jmm - (-1)^q2 Jmp + (-1)^q1 Jpm - (-1)^(q1 + q2) Jpp) (-I);
          S(k1,q1,k2,-q2) = (B(k1,-q1,k2,-q2)-c2*B(k1,-q1,k2,q2)+c1*B(k1,q1,k2,-q2)-c12*B(k1,q1,k2,q2))*ci
          ! CB = (Jpp + (-1)^q2 Jpm - (-1)^q1 Jmp - (-1)^(q1 + q2) Jmm) (-I);
          S(k1,-q1,k2,q2) = (B(k1,-q1,k2,-q2)+c2*B(k1,-q1,k2,q2)-c1*B(k1,q1,k2,-q2)-c12*B(k1,q1,k2,q2))*ci
          ! CC = (-Jpp + (-1)^q2 Jpm + (-1)^q1 Jmp - (-1)^(q1 + q2) Jmm);
          S(k1,-q1,k2,-q2) = -B(k1,-q1,k2,-q2)+c2*B(k1,-q1,k2,q2)+c1*B(k1,q1,k2,-q2)-c12*B(k1,q1,k2,q2)
        else if ((q1 == 0) .and. (q2 > 0)) then
          S(k1,q1,k2,q2) = (B(k1,q1,k2,-q2)+c2*B(k1,q1,k2,q2))
          S(k1,q1,k2,-q2) = (B(k1,q1,k2,-q2)-c2*B(k1,q1,k2,q2))*ci
        else if ((q1 > 0) .and. (q2 == 0)) then
          S(k1,q1,k2,q2) = (B(k1,-q1,k2,q2)+c1*B(k1,q1,k2,q2))
          S(k1,-q1,k2,q2) = (B(k1,-q1,k2,q2)-c1*B(k1,q1,k2,q2))*ci
        else if ((q1 == 0) .and. (q2 == 0)) then
          S(k1,q1,k2,q2) = B(k1,q1,k2,q2)
        end if

      end do
    end do
  end do
end do

! scale the Stevens parameters to match the ESO operators:
do k1=1,n1-1
  do q1=-k1,k1
    do k2=1,n2-1
      do q2=-k2,k2
        F1 = Zero
        F2 = Zero
        if (k1 <= 12) F1 = knm(k1,abs(q1))
        if (k2 <= 12) F2 = knm(k2,abs(q2))
        S(k1,q1,k2,q2) = S(k1,q1,k2,q2)*(F1*F2*cOne)
      end do
    end do
  end do
end do

! in case verification is needed:
if (dbg) then
  call mma_allocate(HAM,n1,n1,n2,n2,'recovered HAM_S')
  write(u6,'(A)') 'Extracted exchange parameters: J(k1,q1,k2,q2):'
  write(u6,'(A)') 'using Naoya ITO:'
  do k1=1,n1-1
    do q1=-k1,k1
      do k2=1,n2-1
        do q2=-k2,k2
          !if (abs(S(k1,q1,k2,q2)) > 1.0e-20_wp) then
          write(u6,'(4(A,i2),A,2ES22.13)') 'J(',k1,',',q1,',',k2,',',q2,') = ',J(k1,q1,k2,q2)
          !end if
        end do
      end do
    end do
  end do
  call recover_exch_HAM_from_Naoya_ITO(n1,n2,J,HAM)
  write(u6,'(A,ES20.10)') 'recover from Naoya Jkqkq parameters'
  write(u6,'(A,ES24.14)') 'JKQP: total difference between HAM-H=',dznrm2_(n1*n1*n2*n2,HAM-H,1)

  write(u6,'(A)') 'Extracted exchange parameters: B(k1,q1,k2,q2):'
  write(u6,'(A)') 'using Liviu ITO:'
  do k1=1,n1-1
    do q1=-k1,k1
      do k2=1,n2-1
        do q2=-k2,k2
          !if (abs(S(k1,q1,k2,q2)) > 1.0e-20_wp) then
          write(u6,'(4(A,i2),A,2ES22.13)') 'B(',k1,',',q1,',',k2,',',q2,') = ',B(k1,q1,k2,q2)
          !end if
        end do
      end do
    end do
  end do
  call recover_exch_HAM_from_Liviu_ITO(n1,n2,B,HAM)
  write(u6,'(A,ES20.10)') 'recover from Liviu Bkqkq parameters'
  write(u6,'(A,ES24.14)') 'JKQP: total difference between HAM-H=',dznrm2_(n1*n1*n2*n2,HAM-H,1)

  write(u6,'(A)') 'Extracted exchange parameters: S(k1,q1,k2,q2):'
  write(u6,'(A)') 'using Stevens ESO:'
  do k1=1,n1-1
    do q1=-k1,k1
      do k2=1,n2-1
        do q2=-k2,k2
          !if (abs(S(k1,q1,k2,q2)) > 1.0e-20_wp) then
          write(u6,'(4(A,i2),A,2ES22.13)') 'S(',k1,',',q1,',',k2,',',q2,') = ',S(k1,q1,k2,q2)
          !end if
        end do
      end do
    end do
  end do
  call recover_exch_HAM_from_Stevens_ESO(n1,n2,S,HAM)
  write(u6,'(A,ES20.10)') 'recover from Stevens Skqkq parameters'
  write(u6,'(A,ES24.14)') 'JKQP: total difference between HAM-H=',dznrm2_(n1*n1*n2*n2,HAM-H,1)
  call mma_deallocate(HAM)
end if

!=======================================================================
call mma_deallocate(O1)
call mma_deallocate(O2)
call mma_deallocate(W1)
call mma_deallocate(W2)

return

end subroutine newjkqpar
