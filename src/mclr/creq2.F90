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
use MCLR_Data, only: nDens, ipMatBA, nA
use input_mclr, only: nSym, nAsh, nIsh, nOrb
use Constants, only: Zero

implicit none
integer idSym, n2
real*8 Q(nDens), G2(*), Temp(n2), Scr(n2)
integer iS, jS, kS, lS, ijS, ipS, kAsh, lAsh, ikl, iAsh, jAsh, ipQ, iij, ipG, ipI

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
