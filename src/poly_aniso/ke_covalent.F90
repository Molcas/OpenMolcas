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

subroutine KE_Covalent(N,lant,t,u,OPT,HCOV)
! this function computes the covalent CF Hamiltonian ofr a given Lanthanide

use jcoeff, only: dE, init_Jx, Jx
use Constants, only: Zero, cOne

implicit none
integer N, OPT, lant
real(kind=8) :: t, u
real(kind=8) :: WCG ! Clebsch_Gordan Coefficients
complex(kind=8) :: HCOV(N,N)
! local variables
integer i, j, JLn, ms1, ns1
integer :: iJ, iLS
integer :: iK, ika
real(kind=8) :: HCOV1(N,N)
external WCG

call init_Jx()
HCOV1(:,:) = Zero
JLn = N-1

do i=1,N
  ms1 = -(N-1)+2*(i-1)
  do j=1,N
    ns1 = -(N-1)+2*(j-1)

    if (OPT == 1) then !FULL calculation

      do iLS=1,4
        do iJ=1,17
          do iK=0,6,2
            do ika=-4,4,4
              if (WCG(JLn,JLn,2*iK,0,JLn,JLn) == Zero) cycle
              HCOV1(i,j) = HCOV1(i,j)+t*t/(u+dE(lant,iLS,iJ))*Jx(lant,iLS,iJ,iK,ika,0,0)*WCG(JLn,ns1,2*iK,2*ika,JLn,ms1)/ &
                           WCG(JLn,JLn,2*iK,0,JLn,JLn)
            end do
          end do
        end do
      end do

    else if (OPT == 2) then ! 1/U approximation

      do iLS=1,4
        do iJ=1,17
          do iK=0,6,2
            do ika=-4,4,4
              if (WCG(JLn,JLn,2*iK,0,JLn,JLn) == Zero) cycle
              HCOV1(i,j) = HCOV1(i,j)+t*t/u*Jx(lant,iLS,iJ,iK,ika,0,0)*WCG(JLn,ns1,2*iK,2*ika,JLn,ms1)/WCG(JLn,JLn,2*iK,0,JLn,JLn)
            end do
          end do
        end do
      end do
    end if
    ! kind=8, complex double precision
    HCOV(i,j) = HCOV1(i,j)*cOne
  end do !j
end do !i

return

end subroutine KE_Covalent
