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
use MCLR_Data, only: nDens2, nNA, ipMat, nA
use input_mclr, only: nSym, nAsh, nIsh, nOrb

implicit none
integer idSym, n2
real*8 Q(nDens2), G2(nna,nna,nna,nna), Temp(n2), Scr(n2)
integer iS, jS, kS, lS, ipS, ijS, iAsh, jAsh, kAsh, lAsh, kAA, lAA, iAA, jAA, ipQ, ipM
real*8 RD

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
                call daxpy_(nOrb(ipS),rd,Temp(ipM),1,Q(ipQ),1)
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
