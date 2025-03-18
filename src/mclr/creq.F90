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

use Constants, only: Zero
use MCLR_Data, only: nDens2, ipMatBA, ipMO, nA
use input_mclr, only: nSym, nAsh, nOrb

implicit none
integer idSym
real*8 Q(nDens2), rint(*), G2(*)
integer iS, jS, kS, lS, ijS, iAsh, jAsh, kAsh, lAsh, iij, ikl, ipS, ipQ, ipG, ipi
integer i, j, itri
! Statement function
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

! Q = (pj|kl)d
!  pi         ijkl

Q(:) = Zero
do iS=1,nSym
  ipS = ieor(is-1,idsym-1)+1
  if (norb(ips) /= 0) then

    do jS=1,nsym
      ijS = ieor(is-1,js-1)+1
      do kS=1,nSym
        ls = ieor(ijs-1,ks-1)+1
        do iAsh=1,nAsh(is)
          do jAsh=1,nAsh(js)
            iij = itri(iAsh+nA(is),jAsh+nA(jS))
            do kAsh=1,nAsh(ks)
              do lAsh=1,nAsh(ls)
                ikl = itri(lAsh+nA(lS),kAsh+nA(kS))
                ipQ = ipMatba(ips,is)+norb(ips)*(iAsh-1)
                ipG = itri(iij,ikl)
                ipi = ipMO(js,ks,ls)+(norb(ips)*(jAsh-1+nAsh(js)*(kAsh-1+nAsh(ks)*(lAsh-1))))
                call daxpy_(norb(ips),G2(ipG),rint(ipI),1,Q(ipQ),1)
              end do
            end do
          end do
        end do
      end do
    end do
  end if
end do

end subroutine creq
