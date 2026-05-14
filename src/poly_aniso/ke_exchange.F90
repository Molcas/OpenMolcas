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

subroutine KE_exchange(N1,N2,lant,t,u,OPT,HEXC)
! this function computes the exchange+covalent contributions to Hamiltonian of a given Lanthanide

use jcoeff, only: dE, init_Jx, Jx
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cOne, auTocm, auToeV
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N1, N2, lant, OPT
real(kind=wp), intent(in) :: t, u
complex(kind=wp), intent(out) :: HEXC(N1,N1,N2,N2)
integer(kind=iwp) :: i1, i2, iJ, iK, ika, iLS, iP, iph, j1, j2, JLn, ms1, ms2, ns1, ns2, SR
real(kind=wp) :: Jfinal(0:7,-5:5,0:1,-1:1)
real(kind=wp), allocatable :: HEXC1(:,:,:,:)
real(kind=wp), parameter :: conv = 1.0e3_wp*auToeV/auTocm ! cm-1 to meV?
real(kind=wp), external :: WCG ! Clebsch-Gordan Coefficients

call init_Jx()
call mma_allocate(HEXC1,N1,N1,N2,N2,label='HEXC1')
HEXC1(:,:,:,:) = Zero
JLn = N1-1  !nexch(iLn)-1
SR = N2-1  !nexch(iRad)-1
do i1=1,N1
  ms1 = -(N1-1)+2*(i1-1)
  do j1=1,N1
    ns1 = -(N1-1)+2*(j1-1)
    do i2=1,N2
      ms2 = -(N2-1)+2*(i2-1)
      do j2=1,N2
        ns2 = -(N2-1)+2*(j2-1)

        if (OPT == 1) then ! FULL calculation

          do iLS=1,4
            do iJ=1,17
              do iK=0,7
                do ika=-5,5
                  do iP=0,1
                    do iph=-1,1
                      if (WCG(JLn,JLn,2*iK,0,JLn,JLn)*WCG(SR,SR,2*iP,0,SR,SR) == Zero) cycle
                      if (WCG(JLn,ns1,2*iK,2*ika,JLn,ms1)*WCG(SR,ns2,2*iP,2*iph,SR,ms2) == Zero) cycle

                      HEXC1(i1,j1,i2,j2) = HEXC1(i1,j1,i2,j2)+t*t/(u+dE(lant,iLS,iJ))*Jx(lant,iLS,iJ,iK,ika,iP,iph)* &
                                           WCG(JLn,ns1,2*iK,2*ika,JLn,ms1)*WCG(SR,ns2,2*iP,2*iph,SR,ms2)/ &
                                           WCG(JLn,JLn,2*iK,0,JLn,JLn)/WCG(SR,SR,2*iP,0,SR,SR)
                    end do
                  end do
                end do
              end do
            end do
          end do

          if ((i1 == 1) .and. (j1 == 1) .and. (i2 == 1) .and. (j2 == 1)) then
            Jfinal = Zero
            write(u6,'(A)') 'Jfinal(iK,ika,iP,iph):'
            do iK=0,7
              do ika=-5,5
                do iP=0,1
                  do iph=-1,1

                    do iLS=1,4
                      do iJ=1,17

                        Jfinal(iK,ika,iP,iph) = Jfinal(iK,ika,iP,iph)+t*t/(u+dE(lant,iLS,iJ))*Jx(lant,iLS,iJ,iK,ika,iP,iph)

                      end do
                    end do
                    if (abs(Jfinal(iK,ika,iP,iph)) > 1.0e-13_wp) &
                      write(u6,'(4i4,4x,2ES24.14)') iK,ika,iP,iph,Jfinal(iK,ika,iP,iph),Jfinal(iK,ika,iP,iph)*conv
                  end do
                end do
              end do
            end do
          end if

        else if (OPT == 2) then ! 1/U approximation

          do iLS=1,4
            do iJ=1,17
              do iK=0,7
                do ika=-5,5
                  do iP=0,1
                    do iph=-1,1
                      if (WCG(JLn,JLn,2*iK,0,JLn,JLn)*WCG(SR,SR,2*iP,0,SR,SR) == Zero) cycle
                      if (WCG(JLn,ns1,2*iK,2*ika,JLn,ms1)*WCG(SR,ns2,2*iP,2*iph,SR,ms2) == Zero) cycle

                      HEXC1(i1,j1,i2,j2) = HEXC1(i1,j1,i2,j2)+t*t/u*Jx(lant,iLS,iJ,iK,ika,iP,iph)*WCG(JLn,ns1,2*iK,2*ika,JLn,ms1)* &
                                           WCG(SR,ns2,2*iP,2*iph,SR,ms2)/WCG(JLn,JLn,2*iK,0,JLn,JLn)/WCG(SR,SR,2*iP,0,SR,SR)
                    end do
                  end do
                end do
              end do
            end do
          end do

        end if !( OPT )

        HEXC(i1,j1,i2,j2) = HEXC1(i1,j1,i2,j2)*cOne

      end do !j2
    end do !i2
  end do !j1
end do !i1

call mma_deallocate(HEXC1)

return

end subroutine KE_exchange
