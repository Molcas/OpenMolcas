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

use Symmetry_Info, only: Mul
use MCLR_Data, only: FIMO, ipCM, ipMat, nA, nDens, nNA
use input_mclr, only: nAsh, nBas, nIsh, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: d_0, rDens1(nna,nna), rdens2(*)
real(kind=wp), intent(out) :: Fock(nDens)
integer(kind=iwp), intent(in) :: idSym
integer(kind=iwp) :: iA, iAA, ip1, ip2, ip3, ipF, ipM, ipS, iS, jA, jAA, jS, kA, kS, n1, n2
real(kind=wp) :: rd, rd1, rd2
real(kind=wp), allocatable :: MO(:), Scr(:), TQ(:)

!                                                                      *
!***********************************************************************
!                                                                      *
Fock(:) = Zero

n1 = max(0,maxval(nBas(1:nSym)))
n2 = n1**2
call mma_allocate(MO,n2,Label='MO')
call mma_allocate(Scr,n2,Label='Scr')

do ipS=1,nSym
  do kS=1,nSym
    do iS=1,nSym
      jS = Mul(Mul(ipS,kS),iS)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Coulomb term: F  = 2(pk|ji)d
      !                kp           ij
      !                                                                *
      !*****************************************************************
      !                                                                *
      if ((Mul(ipS,kS) == iDsym) .and. (nBas(ipS)*nIsh(kS) > 0)) then
        do iA=1,nAsh(iS)
          iAA = iA+nIsh(iS)
          do jA=1,nAsh(jS)
            jAA = jA+nIsh(jS)

            call Coul(ipS,kS,iS,jS,iAA,jAA,MO,Scr)

            rD = rDens1(iA+nA(iS),jA+nA(jS))*Two
            Fock(ipMat(ipS,Ks):ipMat(ipS,Ks)+nBas(ipS)*nIsh(kS)-1) = Fock(ipMat(ipS,Ks):ipMat(ipS,Ks)+nBas(ipS)*nIsh(kS)-1)+ &
                                                                     rd*MO(1:nBas(ipS)*nIsh(kS))

          end do
        end do
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Exchange term   F = -(pk|ji)d     (i>j)  OBS KAN VARA FEL
      !                  pl          kj

      if ((Mul(ips,iS) == iDsym) .and. (nBas(ipS) > 0)) then
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
              Fock(ipF:ipF+nBas(ipS)-1) = Fock(ipF:ipF+nBas(ipS)-1)-rd1*Half*MO(ipM:ipM+nBas(ipS)-1)
              Fock(ipMat(is,ips)+iA-1:ipMat(is,ips)+iA-1+nBas(ipS)*nBas(iS)-1:nBas(iS)) = &
                Fock(ipMat(is,ips)+iA-1:ipMat(is,ips)+iA-1+nBas(ipS)*nBas(iS)-1:nBas(iS))+rd2*Half*MO(ipM:ipM+nBas(ipS)-1)
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
    jS = Mul(is,iDSym)
    do iA=1,nAsh(is)
      do jA=1,nAsh(js)
        rd2 = rDens1(iA+nA(iS),jA+nA(js))
        rd1 = rDens1(jA+nA(jS),iA+nA(is))
        ip1 = nBas(iS)*(nIsh(is)+iA-1)+ipCM(is)
        ip2 = nBas(iS)*(nIsh(js)+jA-1)+ipmat(is,js)
        ip3 = nIsh(js)+jA-1+ipmat(js,is)
        Fock(ip2:ip2+nBas(iS)-1) = Fock(ip2:ip2+nBas(iS)-1)+Rd1*FIMO(ip1:ip1+nBas(iS)-1)
        Fock(ip3:ip3+nBas(iS)*nBas(jS)-1:nBas(jS)) = Fock(ip3:ip3+nBas(iS)*nBas(jS)-1:nBas(jS))-Rd2*FIMO(ip1:ip1+nBas(iS)-1)
      end do
    end do
  end if
end do
! QB is calc here

call CreQADD(Fock,rdens2,idsym,MO,Scr,n2)

call mma_allocate(TQ,nDens,Label='TQ')
TQ(:) = Zero

! QA here
call CreQADD2(TQ,rdens2,idsym,MO,Scr,n2)

call mma_deallocate(Scr)
call mma_deallocate(MO)

call mma_allocate(Scr,nDens,Label='Scr')
Scr(:) = Fock(:)

do iS=1,nsym
  jS = Mul(is,idsym)
  if (nBas(iS)*nBas(jS) > 0) &
    call DGeSub(Fock(ipMat(is,js)),nbas(is),'N',TQ(ipMat(js,is)),nbas(js),'T',Scr(ipmat(is,js)),nbas(is),nbas(is),nbas(js))
end do

if (idSym == 1) call AddGrad2(Scr,d_0)

Fock(:) = Two*Scr(:)

call mma_deallocate(TQ)
call mma_deallocate(Scr)
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine FockGen_td
