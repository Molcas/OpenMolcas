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

subroutine creqadd_sp(q,G2,idsym,Temp,Scr,n2)
! Constructs the Q matrix

use Symmetry_Info, only: Mul
use MCLR_Data, only: ipMat, nA, nDens, nNA
use input_mclr, only: nAsh, nIsh, nOrb, nSym
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: Q(nDens)
real(kind=wp), intent(in) :: G2(nna,nna,nna,nna)
integer(kind=iwp), intent(in) :: idSym, n2
real(kind=wp), intent(out) :: Temp(n2), Scr(n2)
integer(kind=iwp) :: iAA, iAsh, ijS, ipM, ipQ, ipS, iS, jAA, jAsh, jS, kAA, kAsh, kS, lAA, lAsh, lS
real(kind=wp) :: RD

!                                                                      *
!***********************************************************************
!                                                                      *
! Q = (pj|kl)d
!  pi         ijkl

do iS=1,nSym
  ipS = Mul(is,idsym)
  if (norb(ips) /= 0) then
    do jS=1,nsym
      ijS = Mul(is,js)
      do kS=1,nSym
        ls = Mul(ijs,Mul(ks,idsym))

        do kAsh=1,nAsh(ks)
          kAA = kAsh+nIsh(ks)
          do lAsh=1,nAsh(ls)
            lAA = lAsh+nIsh(ls)

            call Coul(ipS,jS,kS,lS,kAA,lAA,Temp,Scr)

            do iAsh=1,nAsh(is)
              iAA = iAsh+nIsh(is)
              ipQ = ipMat(ips,is)+(iAA-1)*nOrb(ips)
              do jAsh=1,nAsh(js)
                jAA = jAsh+nIsh(js)
                ipM = (jAA-1)*nOrb(ipS)+1

                rd = G2(iAsh+na(is),jAsh+na(js),kAsh+na(ks),lash+na(ls))
                Q(ipQ:ipQ+nOrb(ipS)-1) = Q(ipQ:ipQ+nOrb(ipS)-1)+rd*Temp(ipM:ipM+nOrb(ipS)-1)
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

end subroutine creqadd_sp
