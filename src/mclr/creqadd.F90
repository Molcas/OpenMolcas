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

subroutine creqadd(q,G2,idSym,MO,Scr,n2)
! Constructs the Q matrix

use MCLR_Data, only: nDens2, nNA, ipMat, nA
use input_mclr, only: nSym, nAsh, nIsh, nOrb

implicit none
integer idSym, n2
real*8 Q(nDens2), G2(*), MO(n2), Scr(n2)
integer iS, jS, kS, lS, ipS, ijS, iAsh, jAsh, kAsh, lAsh, kAA, lAA, iAA, ikl, ipQ, ipM, iij, ipG
real*8 P_ijkl
! Statement function
integer i, j, itri
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

!                                                                      *
!***********************************************************************
!                                                                      *
! Q = (pj|kl)P
!  pi         ijkl
!                                                                      *
!***********************************************************************
!                                                                      *
do iS=1,nSym
  ipS = ieor(is-1,idsym-1)+1
  if (norb(ips) /= 0) then
    do jS=1,nsym
      ijS = ieor(is-1,js-1)+1
      do kS=1,nSym
        ls = ieor(ijs-1,ieor(ks-1,idsym-1))+1
        !                                                              *
        !***************************************************************
        !                                                              *
        do kAsh=1,nAsh(kS)
          kAA = kAsh+nIsh(kS)
          do lAsh=1,nAsh(lS)
            lAA = lAsh+nIsh(lS)
            ikl = nna*(lAsh+nA(lS)-1)+kAsh+nA(kS)

            ! Pick up (pj|kl)

            call Coul(ipS,jS,kS,lS,kAA,lAA,MO,Scr)

            do iAsh=1,nAsh(iS)
              iAA = iAsh+nIsh(iS)
              ipQ = ipMat(ipS,iS)+nOrb(ipS)*(iAA-1)
              ipM = 1+nIsh(jS)*nOrb(ipS)
              do jAsh=1,nAsh(jS)
                iij = nna*(jAsh+nA(jS)-1)+iAsh+nA(iS)
                ipG = itri(iij,ikl)
                P_ijkl = G2(ipG)

                call DaXpY_(nOrb(ipS),P_ijkl,MO(ipM),1,Q(ipQ),1)
                ipM = ipM+nOrb(ipS)

              end do
            end do

          end do
        end do
        !                                                              *
        !***************************************************************
        !                                                              *
      end do  ! kS
    end do    ! jS
  end if
end do        ! iS
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine creqadd
