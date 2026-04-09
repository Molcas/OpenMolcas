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

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use MCLR_Data, only: ipMat, nA, nDens, nNA
use input_mclr, only: nAsh, nIsh, nOrb, nSym
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: Q(nDens)
real(kind=wp), intent(in) :: G2(*)
integer(kind=iwp), intent(in) :: idSym, n2
real(kind=wp), intent(out) :: MO(n2), Scr(n2)
integer(kind=iwp) :: iAA, iAsh, iij, ijS, ikl, ipG, ipM, ipQ, ipS, iS, jAsh, jS, kAA, kAsh, kS, lAA, lAsh, lS
real(kind=wp) :: P_ijkl

!                                                                      *
!***********************************************************************
!                                                                      *
! Q = (pj|kl)P
!  pi         ijkl
!                                                                      *
!***********************************************************************
!                                                                      *
do iS=1,nSym
  ipS = Mul(is,idsym)
  if (norb(ips) /= 0) then
    do jS=1,nsym
      ijS = Mul(is,js)
      do kS=1,nSym
        ls = Mul(ijs,Mul(ks,idsym))
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
                ipG = iTri(iij,ikl)
                P_ijkl = G2(ipG)

                Q(ipQ:ipQ+nOrb(ipS)-1) = Q(ipQ:ipQ+nOrb(ipS)-1)+P_ijkl*MO(ipM:ipM+nOrb(ipS)-1)
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
