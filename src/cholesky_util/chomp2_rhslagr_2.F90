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

subroutine ChoMP2_RHSlagr_2(EOcc,EVir,EFro,EDel,Xaibj,LnPQRSprod,LiPQRSprod,iBatch,jBatch,nOccLeftI,nOccLeftJ,nOrbLeftI,nOrbLeftJ, &
                            nFroLeftI,nFroLeftJ)
! This will calculate the righthandside of the mp2lagrangian.

use ChoMP2, only: iFirstS, LnBatOrb, LnPQprod, LiPQprod

implicit real*8(a-h,o-z)
real*8 EOcc(*), EVir(*), EFro(*), EDel(*), Xaibj(LnPQRSprod)
integer LiPQRSprod(8)
integer nOccLeftI(8), nOccLeftJ(8)
integer nOrbLeftI(8), nOrbLeftJ(8)
integer nFroLeftI(8), nFroLeftJ(8)
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"
! Statement functions
MulD2h(i,j) = ieor(i-1,j-1)+1
iTri(i,j) = max(i,j)*(max(i,j)-3)/2+i+j
iDensVir(i,j,k) = ip_Density(k)+nFro(k)+nOcc(k)+j-1+(nOrb(k)+nDel(k))*(i+nFro(k)+nOcc(k)-1)
iDensOcc(i,j,k) = ip_Density(k)+j-1+(nOrb(k)+nDel(k))*(i-1)
iMp2Lagr(i,j,k) = ip_Mp2Lagr(k)+j-1+(nOcc(k)+nFro(k))*(i-1)
iDiaA(i,j,k) = ip_DiaA(k)+j-1+(nOcc(k)+nFro(k))*(i-1)

do iSymIP=1,nSym
  iSymAQ = iSymIP
  do iSymI=1,nSym
    iSymP = MulD2h(iSymI,iSymIP)
    LiOrb = min(LnBatOrb(iSymI,iBatch),nOccLeftI(iSymI))
    do Li=1,LiOrb
      if (nFroLeftI(iSymI) > 0) then
        iI = Li
      else
        iI = iFirstS(iSymI,iBatch)+Li-1
      end if
      Lii = LiPQprod(iSymI,iSymI,iBatch)+(nFro(iSymI)+nOcc(iSymI)+nVir(iSymI)+nDel(iSymI))*(Li-1)+iI

      do iSymA=1,nSym
        iSymAI = MulD2h(iSymI,iSymA)
        iSymQ = MulD2h(iSymA,iSymAQ)
        iSymIQ = MulD2h(iSymI,iSymQ)
        LaOrb = LnBatOrb(iSymA,jBatch)-nOccLeftJ(iSymA)
        do La=1,LaOrb
          if (nOccLeftJ(iSymA) > 0) then
            iA = La
          else
            iA = iFirstS(iSymA,jBatch)-nFro(iSymA)-nOcc(iSymA)+La-1
          end if
          ! Calculating the diagonal of A: 3(ia|ia) - (ii|aa)
          ! (Only needed once for each iA hence iP = 1)

          if ((iSymI == iSymA) .and. (iSymIP == 1)) then
            Lii = LiPQprod(iSymI,iSymI,iBatch)+(nFro(iSymI)+nOcc(iSymI)+nVir(iSymI)+nDel(iSymI))*(Li-1)+iI
            Laa = LiPQprod(iSymA,iSymA,jBatch)+(nFro(iSymA)+nOcc(iSymA)+nVir(iSymA)+nDel(iSymA))*(La+nOccLeftJ(iSymA)-1)+iA+ &
                  nFro(iSymA)+nOcc(iSymA)
            Lai = LiPQprod(iSymA,iSymI,iBatch)+(nFro(iSymA)+nOcc(iSymA)+nVir(iSymA)+nDel(iSymA))*(Li-1)+iA+nFro(iSymA)+nOcc(iSymA)
            Lia = LiPQprod(iSymI,iSymA,jBatch)+(nFro(iSymA)+nOcc(iSymA)+nVir(iSymA)+nDel(iSymA))*(La+nOccLeftJ(iSymA)-1)+iI
            if (iBatch == jBatch) then
              ip_iiaa = LiPQRSprod(iSymAI)+iTri(Laa,Lii)
              ip_aiai = LiPQRSprod(iSymAI)+iTri(Lia,Lai)
            else
              ip_iiaa = LiPQRSprod(iSymAI)+LnPQprod(iSymAI,iBatch)*(Laa-1)+Lii
              ip_aiai = LiPQRSprod(iSymAI)+LnPQprod(iSymAI,iBatch)*(Lia-1)+Lai
            end if
            DiagInt = 3.0d0*Xaibj(ip_aiai)-Xaibj(ip_iiaa)
            if (iI <= nFro(iSymI)) then
              E_Occ = EFro(iFro(iSymI)+iI)
            else
              E_Occ = EOcc(iOcc(iSymI)+iI-nFro(iSymI))
            end if
            if (iA > nVir(iSymA)) then
              E_Vir = EDel(iDel(iSymA)+iA-nVir(iSymA))
            else
              E_Vir = EVir(iVir(iSymA)+iA)
            end if

            DiagEn = E_Vir-E_Occ
            Work(iDiaA(iA,iI,iSymI)) = Work(iDiaA(iA,iI,iSymI))+1.0d0/(DiagInt+DiagEn)
            !-----------------------------------------------------------
            !write(6,*) 'IA',iI,iA
            !write(6,*) 'Symm',iSymI,iSymA
            !write(6,*) 'iiaa',Xaibj(ip_iiaa)
            !write(6,*) 'iaia',Xaibj(ip_aiai)
            !write(6,*) 'DiagInt',DiagInt
            !write(6,*) 'DiagEn',DiagEn
            !-----------------------------------------------------------
          end if
          do iP=1,nOrb(iSymP)+nDel(iSymP)
            Lip = LiPQprod(iSymP,iSymI,iBatch)+(nFro(iSymP)+nOcc(iSymP)+nVir(iSymP)+nDel(iSymP))*(Li-1)+iP
            Lap = LiPQProd(iSymP,iSymA,jBatch)+(nFro(iSymP)+nOcc(iSymP)+nVir(iSymP)+nDel(iSymP))*(La+nOccLeftJ(iSymA)-1)+iP
            if ((iP <= nFro(iSymP)+nOcc(iSymP)) .and. (iSymI == iSymP)) then
              do iK=1,nFro(iSymQ)+nOcc(iSymQ)
                iJ = iP
                Lak = LiPQProd(iSymQ,iSymA,jBatch)+(nFro(iSymQ)+nOcc(iSymQ)+nVir(iSymQ)+nDel(iSymQ))*(La+nOccLeftJ(iSymA)-1)+iK
                Lik = LiPQProd(iSymQ,iSymI,iBatch)+(nFro(iSymQ)+nOcc(iSymQ)+nVir(iSymQ)+nDel(iSymQ))*(Li-1)+iK

                if (iBatch == jBatch) then
                  ip_ijak = LiPQRSprod(iSymIP)+iTri(Lip,Lak)
                  ip_ikaj = LiPQRSprod(iSymIQ)+iTri(Lik,Lap)
                else
                  ip_ijak = LiPQRSprod(iSymIP)+LnPQprod(iSymIP,iBatch)*(Lak-1)+Lip
                  ip_ikaj = LiPQRSprod(iSymIQ)+LnPQprod(iSymIQ,iBatch)*(Lap-1)+Lik
                end if
                A = 2.0d0*Xaibj(ip_ijak)-Xaibj(ip_ikaj)
                Work(iMp2Lagr(iA,iK,iSymQ)) = Work(iMp2Lagr(iA,iK,iSymQ))+Work(iDensOcc(iJ,iI,iSymI))*A
                !----- Debug Comments -----------------------
                !if (iP <= nFro(iSymP)+nOcc(iSymP)) then
                !  write(6,*) 'AIJK',iA,iI,iJ,iK
                !  write(6,*) 'aijk',Xaibj(ip_ijak)
                !  write(6,*) 'akji',Xaibj(ip_ikaj)
                !  write(6,*) 'A',A
                !  write(6,*) 'Density',Work(iDensOcc(iJ,iI,iSymI))
                !end if
                !--------------------------------------------
              end do
            else if ((iP > nFro(iSymP)+nOcc(iSymP)) .and. (iSymA == iSymQ)) then
              do iB=1,nVir(iSymQ)+nDel(iSymQ)
                iC = iP-nFro(iSymP)-nOcc(iSymP)
                Lab = LiPQProd(iSymQ,iSymA,jBatch)+(nFro(iSymQ)+nOcc(iSymQ)+nVir(iSymQ)+nDel(iSymQ))*(La+nOccLeftJ(iSymA)-1)+iB+ &
                      nFro(iSymQ)+nOcc(iSymQ)
                Lib = LiPQProd(iSymQ,iSymI,iBatch)+(nFro(iSymQ)+nOcc(iSymQ)+nVir(iSymQ)+nDel(iSymQ))*(Li-1)+iB+nFro(iSymQ)+ &
                      nOcc(iSymQ)
                if (iBatch == jBatch) then
                  ip_icab = LiPQRSprod(iSymIP)+iTri(Lip,Lab)
                  ip_ibac = LiPQRSprod(iSymIQ)+iTri(Lib,Lap)
                else
                  ip_icab = LiPQRSprod(iSymIP)+LnPQprod(iSymIP,iBatch)*(Lab-1)+Lip
                  ip_ibac = LiPQRSprod(iSymIQ)+LnPQprod(iSymIQ,iBatch)*(Lap-1)+Lib
                end if
                A = 2.0d0*Xaibj(ip_icab)-Xaibj(ip_ibac)
                Work(iMp2Lagr(iC,iI,iSymI)) = Work(iMp2Lagr(iC,iI,iSymI))+Work(iDensVir(iB,iA,iSymA))*A
                !------ Debug Comments ---------------------------------
                !write(6,*) 'AIBC',iA,iI,iB,iC
                !write(6,*) 'Symm',iSymA,iSymI,iSymQ,iSymP
                !write(6,*) 'icab',Xaibj(ip_icab)
                !write(6,*) 'ibac',Xaibj(ip_ibac)
                !write(6,*) 'adress',iMp2Lagr(iC,iI,iSymI)-ip_Mp2Lagr(1)
                !write(6,*) 'DensAdress',iDensVir(iB,iA,iSymA)-ip_Density(1)
                !write(6,*) 'Dens',Work(iDensVir(iB,iA,iSymA))
                !write(6,*) 'A',A
                !write(6,*) 'Bidrag',A*Work(iDensVir(iB,iA,iSymA))
                !-------------------------------------------------------
              end do
            end if
          end do
        end do
      end do
    end do
  end do
end do

! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(nOrbLeftI)
  call Unused_integer_array(nOrbLeftJ)
  call Unused_integer_array(nFroLeftJ)
end if

end subroutine ChoMP2_RHSlagr_2
