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

subroutine creq_sp(q,rint,G2,idsym)
! Constructs the Q matrix

use Symmetry_Info, only: Mul
use MCLR_Data, only: ipMat, ipMO, nA, nDens, nNA
use input_mclr, only: nAsh, nIsh, nOrb, nSym
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: Q(nDens)
real(kind=wp), intent(in) :: rint(*), G2(nna,nna,nna,nna)
integer(kind=iwp), intent(in) :: idSym
integer(kind=iwp) :: iAsh, ijS, ipi, ipQ, ipS, iS, jAsh, jS, kAsh, kS, lAsh, lS
real(kind=wp) :: rd

! Q = (pj|kl)d
!  pi         ijkl

!Q(:) = Zero
do iS=1,nSym
  ipS = Mul(is,idsym)
  if (norb(ips) /= 0) then

    do jS=1,nsym
      ijS = Mul(is,js)
      do kS=1,nSym
        ls = Mul(ijs,ks)
        do iAsh=1,nAsh(is)
          do jAsh=1,nAsh(js)
            do kAsh=1,nAsh(ks)
              do lAsh=1,nAsh(ls)
                ipQ = ipMat(ips,is)+norb(ips)*(nish(is)+iAsh-1)
                ipi = ipMO(js,ks,ls)+(norb(ips)*(jAsh-1+nAsh(js)*(kAsh-1+nAsh(ks)*(lAsh-1))))
                rd = G2(na(is)+iash,na(js)+jash,na(ks)+kash,na(ls)+lash)
                Q(ipQ:ipQ+norb(is)-1) = Q(ipQ:ipQ+norb(is)-1)+rd*rint(ipI:ipI+norb(is)-1)
              end do
            end do
          end do
        end do
      end do
    end do
  end if
end do

end subroutine creq_sp
