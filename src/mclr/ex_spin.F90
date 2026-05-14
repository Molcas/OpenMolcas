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

use MCLR_Data, only: ipCM, nDens, nNA
use input_mclr, only: nAsh, nBas, nIsh, nSym
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: rD(*)
real(kind=wp), intent(_OUT_) :: Fock(*)
integer(kind=iwp), intent(in) :: nTemp
real(kind=wp), intent(out) :: Temp1(nTemp), Temp2(nDens)
integer(kind=iwp) :: jB, jjB, jS, kS, lB, llB
real(kind=wp) :: rDens

Temp2(:) = Zero
do jS=1,nSym
  Fock(ipCm(js):ipCm(js)+nbas(js)**2-1) = Zero
  do kS=1,nsym

    ! To be debugged!

    if (nBas(js)*nash(ks) > 0) then

      do llB=1,nAsh(kS)
        lB = nIsh(kS)+llB
        do jjB=1,nAsh(kS)
          jB = nIsh(kS)+jjB

          call Exch(jS,kS,jS,kS,lB,jB,Temp1,Temp2)
          rDens = -Half*rD(nna*(jjB-1)+llB)
          Fock(ipCM(jS):ipCM(jS)+nBas(jS)**2-1) = Fock(ipCM(jS):ipCM(jS)+nBas(jS)**2-1)+rDens*Temp1(1:nBas(jS)**2)

        end do
      end do

    end if
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Ex_spin
