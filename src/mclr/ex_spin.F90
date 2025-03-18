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

subroutine Ex_spin(rD,Fock,Temp1,ntemp,Temp2)

use Constants, only: Zero, Half
use MCLR_Data, only: nDens2, ipCM, nNA
use input_mclr, only: nSym, nAsh, nIsh, nBas

implicit none
integer nTemp
real*8 rD(*), Fock(*), Temp1(nTemp), Temp2(nDens2)
integer jS, kS, llB, lB, jjB, jB
real*8 rDens

Temp2(:) = Zero
do jS=1,nSym
  call FZero(Fock(ipCm(js)),nbas(js)**2)
  do kS=1,nsym

    ! To be debugged!

    if (nBas(js)*nash(ks) > 0) then

      do llB=1,nAsh(kS)
        lB = nIsh(kS)+llB
        do jjB=1,nAsh(kS)
          jB = nIsh(kS)+jjB

          call Exch(jS,kS,jS,kS,lB,jB,Temp2,Temp2)
          rDens = -Half*rD(nna*(jjB-1)+llB)
          call DaXpY_(nBas(jS)**2,rDens,Temp1,1,Fock(ipCM(jS)),1)

        end do
      end do

    end if
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Ex_spin
