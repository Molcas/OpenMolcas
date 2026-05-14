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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, cZero, cOne, Onei
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: n1, n2
complex(kind=wp), intent(in) :: H(n1,n1,n2,n2)
complex(kind=wp), intent(out) :: J(n1-1,-(n1-1):n1-1,n2-1,-(n2-1):n2-1), B(n1-1,-(n1-1):n1-1,n2-1,-(n2-1):n2-1), &
                                 S(n1-1,-(n1-1):n1-1,n2-1,-(n2-1):n2-1)
integer(kind=iwp) :: k1, k2, m1, m12, m2, q1, q2
real(kind=wp) :: C01, C02, cr1, cr2, F1, F2, knm(12,0:12), r1, r2
complex(kind=wp) :: c1, c12, c2, cc1, cc2, cf1, cf2
complex(kind=wp), allocatable :: O1(:,:), O2(:,:), W1(:,:), W2(:,:)
complex(kind=wp), external :: trace_exch2
#ifdef _DEBUGPRINT_
complex(kind=wp), allocatable :: HAM(:,:,:,:)
real(kind=wp), external :: dznrm2_
#endif

!-------------------------------------------
if ((n1 < 1) .or. (n2 < 1)) return
!-------------------------------------------
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

        m1 = (-1)**q1
        m2 = (-1)**q2
        m12 = (-1)**(q1+q2)
        c1 = m1*cOne
        c2 = m2*cOne
        c12 = m12*cOne

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

        ! find the parameters:
        ! in the Naoya's operators
        J(k1,-q1,k2,-q2) = cf1*cf2*trace_exch2(n1,n2,H,O1,O2)
        J(k1,-q1,k2,q2) = cf1*cf2*trace_exch2(n1,n2,H,O1,W2)
        J(k1,q1,k2,-q2) = cf1*cf2*trace_exch2(n1,n2,H,W1,O2)
        J(k1,q1,k2,q2) = cf1*cf2*trace_exch2(n1,n2,H,W1,W2)

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
          S(k1,q1,k2,-q2) = -Onei*(B(k1,-q1,k2,-q2)-c2*B(k1,-q1,k2,q2)+c1*B(k1,q1,k2,-q2)-c12*B(k1,q1,k2,q2))
          ! CB = (Jpp + (-1)^q2 Jpm - (-1)^q1 Jmp - (-1)^(q1 + q2) Jmm) (-I);
          S(k1,-q1,k2,q2) = -Onei*(B(k1,-q1,k2,-q2)+c2*B(k1,-q1,k2,q2)-c1*B(k1,q1,k2,-q2)-c12*B(k1,q1,k2,q2))
          ! CC = (-Jpp + (-1)^q2 Jpm + (-1)^q1 Jmp - (-1)^(q1 + q2) Jmm);
          S(k1,-q1,k2,-q2) = -B(k1,-q1,k2,-q2)+c2*B(k1,-q1,k2,q2)+c1*B(k1,q1,k2,-q2)-c12*B(k1,q1,k2,q2)
        else if ((q1 == 0) .and. (q2 > 0)) then
          S(k1,q1,k2,q2) = (B(k1,q1,k2,-q2)+c2*B(k1,q1,k2,q2))
          S(k1,q1,k2,-q2) = -Onei*(B(k1,q1,k2,-q2)-c2*B(k1,q1,k2,q2))
        else if ((q1 > 0) .and. (q2 == 0)) then
          S(k1,q1,k2,q2) = (B(k1,-q1,k2,q2)+c1*B(k1,q1,k2,q2))
          S(k1,-q1,k2,q2) = -Onei*(B(k1,-q1,k2,q2)-c1*B(k1,q1,k2,q2))
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
#ifdef _DEBUGPRINT_
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
#endif

!=======================================================================
call mma_deallocate(O1)
call mma_deallocate(O2)
call mma_deallocate(W1)
call mma_deallocate(W2)

return

end subroutine newjkqpar
