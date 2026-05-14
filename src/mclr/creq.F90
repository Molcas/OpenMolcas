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

subroutine creq(q,rint,G2,idsym)
! Constructs the Q matrix

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use MCLR_Data, only: ipMatBA, ipMO, nA, nDens
use input_mclr, only: nAsh, nOrb, nSym
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: Q(nDens)
real(kind=wp), intent(in) :: rint(*), G2(*)
integer(kind=iwp), intent(in) :: idSym
integer(kind=iwp) :: iAsh, iij, ijS, ikl, ipG, ipi, ipQ, ipS, iS, jAsh, jS, kAsh, kS, lAsh, lS

! Q = (pj|kl)d
!  pi         ijkl

Q(:) = Zero
do iS=1,nSym
  ipS = Mul(is,idsym)
  if (norb(ips) /= 0) then

    do jS=1,nsym
      ijS = Mul(is,js)
      do kS=1,nSym
        ls = Mul(ijs,ks)
        do iAsh=1,nAsh(is)
          do jAsh=1,nAsh(js)
            iij = iTri(iAsh+nA(is),jAsh+nA(jS))
            do kAsh=1,nAsh(ks)
              do lAsh=1,nAsh(ls)
                ikl = iTri(lAsh+nA(lS),kAsh+nA(kS))
                ipQ = ipMatba(ips,is)+norb(ips)*(iAsh-1)
                ipG = iTri(iij,ikl)
                ipi = ipMO(js,ks,ls)+(norb(ips)*(jAsh-1+nAsh(js)*(kAsh-1+nAsh(ks)*(lAsh-1))))
                Q(ipQ:ipQ+norb(ips)-1) = Q(ipQ:ipQ+norb(ips)-1)+G2(ipG)*rint(ipI:ipI+norb(ips)-1)
              end do
            end do
          end do
        end do
      end do
    end do
  end if
end do

end subroutine creq
