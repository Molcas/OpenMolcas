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
use MCLR_Data, only: nDens2, nNA, ipMat, ipMO, nA
use input_mclr, only: nSym, nAsh, nIsh, nOrb

implicit none
integer idSym
real*8 Q(nDens2), rint(*), G2(nna,nna,nna,nna)
integer iS, jS, kS, lS, ipS, ijS, iAsh, jAsh, kAsh, lAsh, ipQ, ipi
real*8 rd

! Q = (pj|kl)d
!  pi         ijkl

!call dcopy_(ndens2,[Zero],0,Q,1)
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
                call daxpy_(norb(ips),rd,rint(ipI),1,Q(ipQ),1)
              end do
            end do
          end do
        end do
      end do
    end do
  end if
end do

end subroutine creq_sp
