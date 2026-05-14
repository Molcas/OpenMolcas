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

subroutine JKQPar_Naoya(N1,N2,HEXCH,Jpar)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: cZero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: N1, N2
complex(kind=wp), intent(in) :: HEXCH(N1,N1,N2,N2)
complex(kind=wp), intent(out) :: Jpar(N1-1,(-N1+1):N1-1,N2-1,(-N2+1):N2-1)
integer(kind=iwp) :: is1, is2, j1, j2, js1, js2, k1, k2, ms1, ms2, ns1, ns2, q1, q2
real(kind=wp) :: FACT, OPER
complex(kind=wp) :: trace
complex(kind=wp), allocatable :: QMAT(:,:,:,:)
real(kind=wp), external :: WCG

call mma_allocate(QMAT,N1,N1,N2,N2,label='QMAT')

! we need to project now the HEXCH: in products of ITOs
!  HEXCH = SUM(rank1,proj1,rank2,proj2) = { B(rank1,proj1,rank2,proj2)* O1(rank1,proj1) * O2(rank2,proj2) }
! Naoya definition
! eq.40 in Doi:10.1103/PhysRevB.91.174438
Jpar(:,:,:,:) = cZero
J1 = N1-1 ! i.e. double of the dimension of the spin on site 1
J2 = N2-1 ! i.e. double of the dimension of the spin on site 2
do k1=0,J1
  do q1=-k1,k1
    do k2=0,J2
      do q2=-k2,k2
        ! the rank of individual spins must be even
        if (mod(k1,2) /= 1) cycle
        ! the rank of individual spins must be even
        if (mod(k2,2) /= 1) cycle
        ! If the total rank is odd, Then it is a local ZFS contribution; ==> to be Done later
        ! compute the qmat:
        ! projections q and q' are with opposite sign:
        !  -is1 and  -is2
        do is1=1,N1
          ms1 = 2*is1-N1-1 ! spin projection on site 1
          do js1=1,N1
            ns1 = 2*js1-N1-1 ! spin projection on site 1
            do is2=1,N2
              ms2 = 2*is2-N2-1 ! spin projection on site 2
              do js2=1,N2
                ns2 = 2*js2-N2-1 ! spin projection on site 2
                if (WCG(J1,J1,2*k1,0,J1,J1) == 0) cycle
                if (WCG(J2,J2,2*k2,0,J2,J2) == 0) cycle

                QMAT(is1,js1,is2,js2) = WCG(J1,ns1,2*k1,-2*q1,J1,ms1)*WCG(J2,ns2,2*k2,-2*q2,J2,ms2)/WCG(J1,J1,2*k1,0,J1,J1)/ &
                                        WCG(J2,J2,2*k2,0,J2,J2)

#               ifdef _DEBUGPRINT_
                if (abs(QMAT(is1,js1,is2,js2)) > 1.0e-12_wp) &
                  write(u6,'(8(A,i3),A,2F20.14)') ' QMAT(',k1,',',q1,',',k2,',',q2,'|||',is1,',',js1,',',is2,',',js2,') = ', &
                                                  QMAT(is1,js1,is2,js2)
#               endif

              end do
            end do
          end do
        end do
        ! compute the trace Tr[QMAT * HEXCH]
        trace = cZero
        do is1=1,N1
          do js1=1,N1
            do is2=1,N2
              do js2=1,N2
                trace = trace+QMAT(is1,js1,is2,js2)*HEXCH(js1,is1,js2,is2)
              end do
            end do
          end do
        end do

#       ifdef _DEBUGPRINT_
        if (abs(trace) > 1.0e-12_wp) write(u6,'(4(A,i3),A,2F20.14)') 'trace(',k1,',',q1,',',k2,',',q2,') = ',trace
#       endif

        OPER = WCG(J1,J1,2*k1,0,J1,J1)*WCG(J1,J1,2*k1,0,J1,J1)*WCG(J2,J2,2*k2,0,J2,J2)*WCG(J2,J2,2*k2,0,J2,J2)
        FACT = real(((-1)**(q1+q2))*((2*k1+1)*(2*k2+1)),kind=wp)/real(N1*N2,kind=wp)

        Jpar(k1,q1,k2,q2) = FACT*OPER*trace

#       ifdef _DEBUGPRINT_
        if (abs(Jpar(k1,q1,k2,q2)) > 0.5e-13_wp) &
          write(u6,'(4(A,i3),A,2F20.14)') '    J(',k1,',',q1,',',k2,',',q2,') = ',Jpar(k1,q1,k2,q2)
#       endif

      end do
    end do
  end do
end do

call mma_deallocate(QMAT)

return

end subroutine JKQPar_Naoya
