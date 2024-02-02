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

implicit none
integer, parameter :: wp = kind(0.d0)
#include "stdalloc.fh"
integer, intent(in) :: n1, n2
complex(kind=8), intent(in) :: H(n1,n1,n2,n2)
complex(kind=8), intent(out) :: J((n1-1),-(n1-1):(n1-1),(n2-1),-(n2-1):(n2-1))
complex(kind=8), intent(out) :: B((n1-1),-(n1-1):(n1-1),(n2-1),-(n2-1):(n2-1))
complex(kind=8), intent(out) :: S((n1-1),-(n1-1):(n1-1),(n2-1),-(n2-1):(n2-1))
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
J(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)) = (0.0_wp,0.0_wp)
B(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)) = (0.0_wp,0.0_wp)
S(1:(n1-1),-(n1-1):(n1-1),1:(n2-1),-(n2-1):(n2-1)) = (0.0_wp,0.0_wp)
do k1=1,n1-1
  do q1=0,k1
    do k2=1,n2-1
      do q2=0,k2
        cr1 = 0.0_wp
        cr2 = 0.0_wp
        m1 = 0
        m2 = 0
        r1 = 0.0_wp
        r2 = 0.0_wp
        cf1 = (0.0_wp,0.0_wp)
        cf2 = (0.0_wp,0.0_wp)

        m1 = (-1)**q1
        m2 = (-1)**q2
        m12 = (-1)**(q1+q2)
        c1 = cmplx(m1,0.0_wp,wp)
        c2 = cmplx(m2,0.0_wp,wp)
        c12 = cmplx(m12,0.0_wp,wp)
        ci = (0.0_wp,-1.0_wp)

        ! generate the operators for site 1 and 2:
        call ITO(n1,k1,q1,C01,O1,W1)
        call ITO(n2,k2,q2,C02,O2,W2)

        ! compute the scaling factors for various Operators
        call coeff_redus_sub(n1,k1,cr1)
        call coeff_redus_sub(n2,k2,cr2)
        cc1 = cmplx(1.0_wp/(cr1*C01),0.0_wp,wp)
        cc2 = cmplx(1.0_wp/(cr2*C02),0.0_wp,wp)
        r1 = C01*C01*dble(2*k1+1)/dble(n1)
        r2 = C02*C02*dble(2*k2+1)/dble(n2)
        cf1 = c1*cmplx(r1,0.0_wp,wp)
        cf2 = c2*cmplx(r2,0.0_wp,wp)

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
        F1 = 0.0_wp
        F2 = 0.0_wp
        if (k1 <= 12) F1 = knm(k1,abs(q1))
        if (k2 <= 12) F2 = knm(k2,abs(q2))
        S(k1,q1,k2,q2) = S(k1,q1,k2,q2)*cmplx(F1*F2,0.0_wp,wp)
      end do
    end do
  end do
end do

! in case verification is needed:
if (dbg) then
  call mma_allocate(HAM,n1,n1,n2,n2,'recovered HAM_S')
  write(6,'(A)') 'Extracted exchange parameters: J(k1,q1,k2,q2):'
  write(6,'(A)') 'using Naoya ITO:'
  do k1=1,n1-1
    do q1=-k1,k1
      do k2=1,n2-1
        do q2=-k2,k2
          !if (abs(S(k1,q1,k2,q2)) > 1.0e-20_wp) then
          write(6,'(4(A,i2),A,2ES22.13)') 'J(',k1,',',q1,',',k2,',',q2,') = ',J(k1,q1,k2,q2)
          !end if
        end do
      end do
    end do
  end do
  call recover_exch_HAM_from_Naoya_ITO(n1,n2,J,HAM)
  write(6,'(A,ES20.10)') 'recover from Naoya Jkqkq parameters'
  write(6,'(A,ES24.14)') 'JKQP: total difference between HAM-H=',dznrm2_(n1*n1*n2*n2,HAM-H,1)

  write(6,'(A)') 'Extracted exchange parameters: B(k1,q1,k2,q2):'
  write(6,'(A)') 'using Liviu ITO:'
  do k1=1,n1-1
    do q1=-k1,k1
      do k2=1,n2-1
        do q2=-k2,k2
          !if (abs(S(k1,q1,k2,q2)) > 1.0e-20_wp) then
          write(6,'(4(A,i2),A,2ES22.13)') 'B(',k1,',',q1,',',k2,',',q2,') = ',B(k1,q1,k2,q2)
          !end if
        end do
      end do
    end do
  end do
  call recover_exch_HAM_from_Liviu_ITO(n1,n2,B,HAM)
  write(6,'(A,ES20.10)') 'recover from Liviu Bkqkq parameters'
  write(6,'(A,ES24.14)') 'JKQP: total difference between HAM-H=',dznrm2_(n1*n1*n2*n2,HAM-H,1)

  write(6,'(A)') 'Extracted exchange parameters: S(k1,q1,k2,q2):'
  write(6,'(A)') 'using Stevens ESO:'
  do k1=1,n1-1
    do q1=-k1,k1
      do k2=1,n2-1
        do q2=-k2,k2
          !if (abs(S(k1,q1,k2,q2)) > 1.0e-20_wp) then
          write(6,'(4(A,i2),A,2ES22.13)') 'S(',k1,',',q1,',',k2,',',q2,') = ',S(k1,q1,k2,q2)
          !end if
        end do
      end do
    end do
  end do
  call recover_exch_HAM_from_Stevens_ESO(n1,n2,S,HAM)
  write(6,'(A,ES20.10)') 'recover from Stevens Skqkq parameters'
  write(6,'(A,ES24.14)') 'JKQP: total difference between HAM-H=',dznrm2_(n1*n1*n2*n2,HAM-H,1)
  call mma_deallocate(HAM)
end if

!=======================================================================
call mma_deallocate(O1)
call mma_deallocate(O2)
call mma_deallocate(W1)
call mma_deallocate(W2)

return

end subroutine newjkqpar
