!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine FockGen_td(d_0,rDens1,rdens2,fock,idsym)
!***********************************************************************
!                                                                      *
!   Constructs active fockmatrix and Q matrix                          *
!                                                                      *
!   Input: rkappa : Rotation matrix                                    *
!          idsym  : symmetry of perturbation                           *
!                                                                      *
!                                                                      *
!   Output:MO     :MO integrals                                        *
!          Fock   :Fock matrix (one index transformed integrals)       *
!          MOtilde:MO (one index transformed integrals)                *
!                                                                      *
!***********************************************************************

use Arrays, only: FIMO
use MCLR_Data, only: nDens2, nNA, ipMat, ipCM, nA
use input_mclr, only: nSym, nAsh, nIsh, nBas
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Half

implicit none
real*8 d_0
integer idSym
real*8 Fock(nDens2), rdens2(*), rDens1(nna,nna)
real*8, allocatable :: MO(:), Scr(:), TQ(:)
integer n1, iS, n2, ipS, kS, jS, iA, iAA, jA, jAA, ipF, ipM, kA, ip1, ip2, ip3
real*8 rd, rd1, rd2

!                                                                      *
!***********************************************************************
!                                                                      *
Fock(:) = Zero

n1 = 0
do iS=1,nSym
  n1 = max(n1,nBas(iS))
end do
n2 = n1**2
call mma_allocate(MO,n2,Label='MO')
call mma_allocate(Scr,n2,Label='Scr')

do ipS=1,nSym
  do kS=1,nSym
    do iS=1,nSym
      jS = ieor(ieor(ipS-1,kS-1),iS-1)+1
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Coulomb term: F  = 2(pk|ji)d
      !                kp           ij
      !                                                                *
      !*****************************************************************
      !                                                                *
      if ((ieor(ipS-1,kS-1)+1 == iDsym) .and. (nBas(ipS)*nIsh(kS) > 0)) then
        do iA=1,nAsh(iS)
          iAA = iA+nIsh(iS)
          do jA=1,nAsh(jS)
            jAA = jA+nIsh(jS)

            call Coul(ipS,kS,iS,jS,iAA,jAA,MO,Scr)

            rD = rDens1(iA+nA(iS),jA+nA(jS))*Two
            call DaXpY_(nBas(ipS)*nIsh(kS),rd,MO,1,Fock(ipMat(ipS,Ks)),1)

          end do
        end do
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Exchange term   F = -(pk|ji)d     (i>j)  OBS KAN VARA FEL
      !                  pl          kj

      if ((ieor(ips-1,iS-1)+1 == iDsym) .and. (nBas(ipS) > 0)) then
        do iA=1,nIsh(iS)
          ipF = ipMat(ipS,iS)+nBas(ipS)*(iA-1)
          do jA=1,nAsh(jS)
            jAA = jA+nIsh(js)

            call Coul(ipS,kS,iS,jS,iA,jAA,MO,Scr)

            ipM = 1+nIsh(kS)*nBas(ipS)
            do kA=1,nAsh(ks)

              ! Two different densities for the exchange term.

              rd1 = rDens1(jA+nA(jS),kA+nA(ks))*Two
              rd2 = rDens1(kA+nA(kS),jA+nA(js))*Two
              call DaXpY_(nBas(ipS),-rd1*Half,MO(ipM),1,Fock(ipF),1)
              call DaXpY_(nBas(ipS),rd2*Half,MO(ipM),1,Fock(ipMat(is,ips)+iA-1),nbas(is))
              ipM = ipM+nBas(ipS)
            end do
          end do
        end do
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! KAN GAA FEL
!
! $F_{pb}=F_{pa}^ID_{ba}$

do iS=1,nSym
  if (nBas(iS) > 0) then
    jS = ieor(is-1,iDSym-1)+1
    do iA=1,nAsh(is)
      do jA=1,nAsh(js)
        rd2 = rDens1(iA+nA(iS),jA+nA(js))
        rd1 = rDens1(jA+nA(jS),iA+nA(is))
        ip1 = nBas(iS)*(nIsh(is)+iA-1)+ipCM(is)
        ip2 = nBas(iS)*(nIsh(js)+jA-1)+ipmat(is,js)
        ip3 = nIsh(js)+jA-1+ipmat(js,is)
        call DaxPy_(nBas(iS),Rd1,FIMO(ip1),1,Fock(ip2),1)
        call DaxPy_(nBAs(iS),-Rd2,FIMO(ip1),1,Fock(ip3),nbas(js))
      end do
    end do
  end if
end do
! QB is calc here

call CreQADD(Fock,rdens2,idsym,MO,Scr,n2)

call mma_allocate(TQ,ndens2,Label='TQ')
TQ(:) = Zero

! QA here
call CreQADD2(TQ,rdens2,idsym,MO,Scr,n2)

call mma_deallocate(Scr)
call mma_deallocate(MO)

do iS=1,nsym
  jS = ieor(is-1,idsym-1)+1
  if (nBas(iS)*nBas(jS) > 0) &
    call DGeSub(Fock(ipMat(is,js)),nbas(is),'N',TQ(ipMat(js,is)),nbas(js),'T',Fock(ipmat(is,js)),nbas(is),nbas(is),nbas(js))
end do

if (idSym == 1) call AddGrad2(Fock,d_0)

call DScal_(nDens2,Two,Fock,1)

call mma_deallocate(TQ)
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine FockGen_td
