!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************

subroutine creq2(q,G2,idSym,Temp,Scr,n2)
! Constructs the Q matrix

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use MCLR_Data, only: ipMatBA, nA, nDens
use input_mclr, only: nAsh, nIsh, nOrb, nSym
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idSym, n2
real(kind=wp), intent(out) :: Q(nDens), Temp(n2), Scr(n2)
real(kind=wp), intent(in) :: G2(*)
integer(kind=iwp) :: iAsh, iij, ijS, ikl, ipG, ipI, ipQ, ipS, iS, jAsh, jS, kAsh, kS, lAsh, lS

!                                                                      *
!***********************************************************************
!                                                                      *
!  Q = (pj|kl)d
!   pi         ijkl

Q(:) = Zero

do iS=1,nSym
  ipS = Mul(is,idsym)
  if (norb(ips) /= 0) then
    do jS=1,nsym
      ijS = Mul(is,js)
      do kS=1,nSym
        ls = Mul(ijs,ks)

        do kAsh=1,nAsh(ks)
          do lAsh=1,nAsh(ls)
            ikl = iTri(lAsh+nA(lS),kAsh+nA(kS))

            call Coul(ipS,jS,kS,lS,nIsh(kS)+kAsh,nIsh(lS)+lAsh,Temp,Scr)

            do iAsh=1,nAsh(is)
              ipQ = ipMatba(ipS,iS)+nOrb(ipS)*(iAsh-1)
              do jAsh=1,nAsh(jS)
                iij = iTri(iAsh+nA(iS),jAsh+nA(jS))
                ipG = iTri(iij,ikl)
                ipI = (nIsh(jS)+jAsh-1)*nOrb(ipS)+1

                Q(ipQ:ipQ+nOrb(ipS)-1) = Q(ipQ:ipQ+nOrb(ipS)-1)+G2(ipG)*Temp(ipI:ipI+nOrb(ipS)-1)

              end do
            end do
          end do
        end do

      end do
    end do
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine creq2
