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

subroutine ChoMP2_rhslagr_1(EOcc,EVir,EFro,EDel,Xaibj,LnPQRSprod,LiPQRSprod,iBatch,jBatch,nOccLeftI,nOccLeftJ,nOrbLeftI,nOrbLeftJ, &
                            nFroLeftI,nFroLeftJ)
! This will calculate the righthandside of the mp2lagrangian.

use ChoMP2, only: iFirstS, LiPQprod, LnBatOrb, LnPQprod
use Constants, only: Two, Four, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: LnPQRSprod, LiPQRSprod(8), iBatch, jBatch, nOccLeftI(8), nOccLeftJ(8), nOrbLeftI(8), nOrbLeftJ(8), &
                    nFroLeftI(8), nFroLeftJ(8)
real(kind=wp) :: EOcc(*), EVir(*), EFro(*), EDel(*), Xaibj(LnPQRSprod)
#include "cholesky.fh"
#include "chomp2.fh"
#include "chomp2_cfg.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iA, iB, iC, iCFroz, iI, iJ, iK, iKFroz, iP, ip_iajb, ip_iapb, ip_ibja, ip_ipjb, iSymA, iSymAI, iSymAJ, iSymB, &
                     iSymBI, iSymBJ, iSymI, iSymJ, iSymP, LA, Lai, Laj, LaOrb, Lb, Lbi, Lbj, LbOrb, Lbp, Li, LiOrb, Lj, LjOrb, Lpi
real(kind=wp) :: Dnom1, Dnom2, T1, X2
! Statement functions
integer(kind=iwp) :: MulD2h, iTri, iDensVir, iWDensVir, iDensVirFro, iWDensVirFro, iDensFroVir, iWDensOcc, iWDensOccFro, iDensOcc, &
                     iDensFroOcc, iDensOccFro, iMp2Lagr, iWDensVactOcc, i, j, k
MulD2h(i,j) = ieor(i-1,j-1)+1
iTri(i,j) = max(i,j)*(max(i,j)-3)/2+i+j
iDensVir(i,j,k) = ip_Density(k)+nFro(k)+nOcc(k)+j-1+(nOrb(k)+nDel(k))*(i+nFro(k)+nOcc(k)-1)
iWDensVir(i,j,k) = ip_WDensity(k)+nFro(k)+nOcc(k)+j-1+(nOrb(k)+nDel(k))*(i+nFro(k)+nOcc(k)-1)
iDensVirFro(i,j,k) = ip_Density(k)+nFro(k)+nOcc(k)+nVir(k)+j-1+(nOrb(k)+nDel(k))*(i+nFro(k)+nOcc(k)-1)
iWDensVirFro(i,j,k) = ip_WDensity(k)+nFro(k)+nOcc(k)+nVir(k)+j-1+(nOrb(k)+nDel(k))*(i+nFro(k)+nOcc(k)-1)
iDensFroVir(i,j,k) = ip_Density(k)+nFro(k)+nOcc(k)+j-1+(nOrb(k)+nDel(k))*(i+nFro(k)+nOcc(k)+nVir(k)-1)
iWDensOcc(i,j,k) = ip_WDensity(k)+j+nFro(k)-1+(nOrb(k)+nDel(k))*(i+nFro(k)-1)
iWDensOccFro(i,j,k) = ip_WDensity(k)+j-1+(nOrb(k)+nDel(k))*(i+nFro(k)-1)
iDensOcc(i,j,k) = ip_Density(k)+nFro(k)+j-1+(nOrb(k)+nDel(k))*(i+nFro(k)-1)
iDensFroOcc(i,j,k) = ip_Density(k)+j-1+nFro(k)+(nOrb(k)+nDel(k))*(i-1)
iDensOccFro(i,j,k) = ip_Density(k)+j-1+(nOrb(k)+nDel(k))*(i+nFro(k)-1)
iMp2Lagr(i,j,k) = ip_Mp2Lagr(k)+j-1+(nOcc(k)+nFro(k))*(i-1)
iWDensVactOcc(i,j,k) = ip_WDensity(k)+j-1+(nOrb(k)+nDel(k))*(i+nFro(k)+nOcc(k)-1)

do iSymBJ=1,nSym
  iSymAI = iSymBJ
  !***************************************************************
  ! Common code for Pab, Wab, Lagr(3) and PaB
  !***************************************************************
  do iSymJ=1,nSym
    iSymB = MulD2h(iSymJ,iSymBJ)
    LjOrb = min(LnBatOrb(iSymJ,jBatch)-nFroLeftJ(iSymJ),nOccLeftJ(iSymJ)-nFroLeftJ(iSymJ))
    do Lj=1,LjOrb
      if (nFroLeftJ(iSymJ) > 0) then
        iJ = Lj
      else
        iJ = iFirstS(iSymJ,jBatch)+Lj-nFro(iSymJ)-1
      end if
      do iB=1,nVir(iSymB)
        Lbj = LiPQprod(iSymB,iSymJ,jBatch)+(nFro(iSymB)+nOcc(iSymB)+nVir(iSymB)+nDel(iSymB))*(Lj+nFroLeftJ(iSymJ)-1)+iB+ &
              nFro(iSymB)+nOcc(iSymB)
        do iSymI=1,nSym
          iSymA = MulD2h(iSymI,iSymAI)
          iSymAJ = MulD2h(iSymA,iSymJ)
          iSymBI = iSymAJ
          iSymP = iSymA
          LiOrb = min(LnBatOrb(iSymI,iBatch)-nFroLeftI(iSymI),nOccleftI(iSymI)-nFroLeftI(iSymI))
          do Li=1,LiOrb
            if (nFroLeftI(iSymI) > 0) then
              iI = Li
            else
              iI = iFirstS(iSymI,iBatch)-nFro(iSymI)+Li-1
            end if
            Lbi = LiPQprod(iSymB,iSymI,iBatch)+(nFro(iSymB)+nOcc(iSymB)+nVir(iSymB)+nDel(iSymB))*(Li+nFroLeftI(iSymI)-1)+iB+ &
                  nFro(iSymB)+nOcc(iSymB)
            do iA=1,nVir(iSymA)
              Lai = LiPQprod(iSymA,iSymI,iBatch)+(nOcc(iSymA)+nVir(iSymA)+nFro(iSymA)+nDel(iSymA))*(Li+nFroLeftI(iSymI)-1)+iA+ &
                    nFro(iSymA)+nOcc(iSymA)
              Laj = LiPQprod(iSymA,iSymJ,jBatch)+(nOcc(iSymA)+nVir(iSymA)+nFro(iSymA)+nDel(iSymA))*(Lj+nFroLeftJ(iSymJ)-1)+iA+ &
                    nFro(iSymA)+nOcc(iSymA)

              ! Construct the T_iajb (common for all terms)

              Dnom1 = EOcc(iOcc(iSymI)+iI)-EVir(iVir(iSymA)+iA)+EOcc(iOcc(iSymJ)+iJ)-EVir(iVir(iSymB)+iB)
              if (iBatch == jBatch) then
                ip_iajb = LiPQRSprod(iSymAI)+iTri(Lai,Lbj)
                ip_ibja = LiPQRSprod(iSymBI)+iTri(Lbi,Laj)
              else
                ip_iajb = LiPQRSprod(iSymAI)+LnPQprod(iSymAI,iBatch)*(Lbj-1)+Lai
                ip_ibja = LiPQRSprod(iSymBI)+LnPQprod(iSymBI,iBatch)*(Laj-1)+Lbi
              end if
              T1 = Four*Xaibj(ip_iajb)
              T1 = (T1-Two*Xaibj(ip_ibja))/Dnom1
              ! Here we calculate the MP2-energy
              EMP2_Dens = EMP2_Dens+Half*T1*Xaibj(ip_iajb)
              do iP=1,nOrb(iSymP)+nDel(iSymP)
                Lpi = LiPQprod(iSymP,iSymI,iBatch)+(nOcc(iSymP)+nVir(iSymP)+nFro(iSymP)+nDel(iSymP))*(Li+nFroLeftI(iSymI)-1)+iP
                !Lpj = LiPQprod(iSymP,iSymJ,jBatch)+(nOcc(iSymP)+nVir(iSymP)+nFro(iSymP)+nDel(iSymP))*(Lj+nFroLeftJ(iSymJ)-1)+iP
                if (iBatch == jBatch) then
                  ip_ipjb = LiPQRSprod(iSymBJ)+iTri(Lpi,Lbj)
                  !ip_jpib = 0
                else
                  ip_ipjb = LiPQRSprod(iSymBJ)+LnPQprod(iSymBJ,iBatch)*(Lbj-1)+Lpi
                  !ip_jpib = LiPQRSprod(iSymBI)+LnPQprod(iSymBI,iBatch)*(Lpj-1)+Lbi
                end if
                X2 = Xaibj(ip_ipjb)
                !**************************************************
                ! Calculate P_Ca frozen virtual - active virtual
                !**************************************************
                if (iP > nFro(iSymP)+nOcc(iSymP)+nVir(iSymP)) then
                  iCFroz = iP-nFro(iSymP)-nOcc(iSymP)-nVir(iSymP)
                  Dnom2 = EVir(iVir(iSymA)+iA)-EDel(iDel(iSymP)+iCFroz)
                  Work(iDensVirFro(iA,iCFroz,iSymP)) = Work(iDensVirFro(iA,iCFroz,iSymP))+T1*X2/Dnom2
                  Work(iDensFroVir(iCFroz,iA,iSymA)) = Work(iDensFroVir(iCFroz,iA,iSymA))+T1*X2/Dnom2
                  Work(iWDensVirFro(iA,iCFroz,iSymP)) = Work(iWDensVirFro(iA,iCFroz,iSymP))-Two*T1*X2
                  !*** The 2 ON THE ROW ABOVE IS QUESTIONABLE *****'''' ***
                  !write(u6,*) 'index Dens',iDensFroVir(iCFroz,iA,iSymA)-ip_Density(iSymA)
                  !write(u6,*) 'Value',Work(iDensFroVir(iCFroz,iA,iSymA))+Two*T1*X2/Dnom2
                  !write(u6,*) 'xicjb',X2
                  !write(u6,*) 'EDiffbc',Dnom2
                  !write(u6,*) 'E_a',EVir(iVir(iSymA)+iA)
                  !write(u6,*) 'E:B',EDel(iDel(iSymP)+iCFroz)
                  !write(u6,*) 'Tij',T1
                  !write(u6,*) 'iB, iC',iCFroz+nOcc(iSymP)-nFro(iSymP),iB
                  !write(u6,*) 'iI, iJ',iI,iJ
                  !************************************************
                  ! Calculate P_ab active vir - active vir
                  !************************************************
                else if (iP > nFro(iSymP)+nOcc(iSymP)) then
                  iC = iP-nFro(iSymP)-nOcc(iSymP)
                  Dnom2 = EOcc(iOcc(iSymI)+iI)-EVir(iVir(iSymP)+iC)+EOcc(iOcc(iSymJ)+iJ)-EVir(iVir(iSymB)+iB)
                  Work(iDensVir(iA,iC,iSymP)) = Work(iDensVir(iA,iC,iSymP))+T1*X2/Dnom2
                  Work(iWDensVir(iA,iC,iSymP)) = Work(iWDensVir(iA,iC,iSymP))-Two*T1*X2
                  !************************************************
                  ! Calculate Lagr(3)_ai  active virtual - all occupied
                  !************************************************
                else
                  iK = iP
                  Work(iMp2Lagr(iA,iK,iSymP)) = Work(iMp2Lagr(iA,iK,iSymP))-T1*X2
                  Work(iWDensVactOcc(iA,iK,iSymP)) = Work(iWDensVactOcc(iA,iK,iSymP))-Four*T1*X2
                end if

                !---------- Debugging comments --------------------------
                !if ((iP > nFro(iSymP)+nOcc(iSymP)) .and. (ip <= nFro(iSymP)+nOcc(iSymP)+nVir(iSymP))) then
                !  write(u6,*) 'AIBJC',iA,iI,iB,iJ,iC
                !  write(u6,*) 'Symm',iSymA,iSymI,iSymB,iSymJ,iSymP
                !  write(u6,*) 'Dnom1',Dnom1
                !  write(u6,*) 'Dnom2',Dnom2
                !  write(u6,*) 'iajb',Xaibj(ip_iajb)
                !  write(u6,*) 'ibja',Xaibj(ip_ibja)
                !  write(u6,*) 'icjb',Xaibj(ip_ipjb)
                !  write(u6,*) 'Bidrag T',T1*X2/dnom2
                !  write(u6,*) 'adress',iDensVir(iA,iC,iSymP)-ip_Density(1)
                !end if
                !--------------------------------------------------------
                !---------- Debugging comments --------------------------
                !if (iP <= nFro(iSymP)+nOcc(iSymP)) then
                !  write(u6,*) 'AIBJK',iA,iI,iB,iJ,iK
                !  write(u6,*) 'Symm',iSymA,iSymI,iSymB,iSymJ,iSymP
                !  write(u6,*) 'Dnom2',Dnom2
                !  write(u6,*) 'iajb',Xaibj(ip_iajb)
                !  write(u6,*) 'ibja',Xaibj(ip_ibja)
                !  write(u6,*) 'ibjk',Xaibj(ip_ipjb)
                !  if (iBatch /= jBatch) then
                !    write(u6,*) 'ikjb' Xaibj(ip_jpib)
                !    write(u6,*) 'Bidrag U',U1*U2
                !  end if
                !  write(u6,*) 'Bidrag T',T1*T2
                !  write(u6,*) 'adress',iMp2Lagr(iA,iK,iSymP)-ip_Mp2Lagr(1)
                !end if
                !--------------------------------------------------------
              end do
            end do
          end do
        end do
      end do
    end do
  end do

  !*********************************************************************
  !  Common code for Pij, Wij(I), Lagr(4) and PiJ
  !*********************************************************************

  ! To avoid compiler warnings
  !iKfroz = 0
  !iK = 0

  do iSymA=1,nSym
    iSymI = MulD2h(iSymA,iSymAI)
    LaOrb = min(LnBatOrb(iSymA,iBatch)-nOccLeftI(iSymA),nOrbLeftI(iSymA)-nOccLeftI(iSymA))
    do La=1,LaOrb
      if (nOccLeftI(iSymA) > 0) then
        iA = La
      else
        iA = iFirstS(iSymA,iBatch)-nFro(iSymA)-nOcc(iSymA)+La-1
      end if
      do iI=1,nOcc(iSymI)
        Lai = LiPQprod(iSymI,iSymA,iBatch)+(nOcc(iSymI)+nVir(iSymI)+nFro(iSymI)+nDel(iSymI))*(La+nOccLeftI(iSymA)-1)+iI+nFro(iSymI)
        do iSymJ=1,nSym
          iSymB = MulD2h(iSymJ,iSymBJ)
          iSymAJ = MulD2h(iSymA,iSymJ)
          iSymBI = iSymAJ
          iSymP = iSymJ
          LbOrb = min(LnBatOrb(iSymB,jBatch)-nOccLeftJ(iSymB),nOrbLeftJ(iSymB)-nOccLeftJ(iSymB))
          do Lb=1,LbOrb
            if (nOccLeftJ(iSymB) > 0) then
              iB = Lb
            else
              iB = iFirstS(iSymB,jBatch)-nFro(iSymB)-nOcc(iSymB)+Lb-1
            end if
            Lbi = LiPQprod(iSymI,iSymB,jBatch)+(nOcc(iSymI)+nVir(iSymI)+nFro(iSymI)+nDel(iSymI))*(Lb+nOccLeftJ(iSymB)-1)+iI+ &
                  nFro(iSymI)
            do iJ=1,nOcc(iSymJ)
              Lbj = LiPQprod(iSymJ,iSymB,jBatch)+(nOcc(iSymJ)+nVir(iSymJ)+nFro(iSymJ)+nDel(iSymJ))*(Lb+nOccLeftJ(iSymB)-1)+iJ+ &
                    nFro(iSymJ)
              Laj = LiPQprod(iSymJ,iSymA,iBatch)+(nOcc(iSymJ)+nVir(iSymJ)+nFro(iSymJ)+nDel(iSymJ))*(La+nOccLeftI(iSymA)-1)+iJ+ &
                    nFro(iSymJ)
              Dnom1 = EOcc(iOcc(iSymI)+iI)-EVir(iVir(iSymA)+iA)+EOcc(iOcc(iSymJ)+iJ)-EVir(iVir(iSymB)+iB)
              if (iBatch == jBatch) then
                ip_iajb = LiPQRSprod(iSymAI)+iTri(Lai,Lbj)
                ip_ibja = LiPQRSprod(iSymBI)+iTri(Lbi,Laj)
              else
                ip_iajb = LiPQRSprod(iSymAI)+LnPQprod(iSymAI,iBatch)*(Lbj-1)+Lai
                ip_ibja = LiPQRSprod(iSymBI)+LnPQprod(iSymBI,iBatch)*(Lbi-1)+Laj
              end if
              T1 = Four*Xaibj(ip_iajb)
              T1 = (T1-Two*Xaibj(ip_ibja))/Dnom1
              do iP=1,nOrb(iSymP)+nDel(iSymP)

                Lbp = LiPQprod(iSymP,iSymB,jBatch)+(nOcc(iSymP)+nVir(iSymP)+nFro(iSymP)+nDel(iSymP))*(Lb+nOccLeftJ(iSymB)-1)+iP
                if (iBatch == jBatch) then
                  ip_iapb = LiPQRSprod(iSymAI)+iTri(Lai,Lbp)
                else
                  ip_iapb = LiPQRSprod(iSymAI)+LnPQprod(iSymAI,iBatch)*(Lbp-1)+Lai
                end if

                X2 = Xaibj(ip_iapb)
                !**************************************************
                ! Calculate P_Pj frozen occupied - active occupied
                !**************************************************
                if (ip <= nFro(iSymP)) then
                  iKFroz = iP
                  Dnom2 = EOcc(iOcc(iSymI)+iI)-EFro(iFro(iSymP)+iKFroz)
                  Work(iDensFroOcc(iKFroz,iJ,iSymJ)) = Work(iDensFroOcc(iKFroz,iJ,iSymJ))+T1*X2/Dnom2
                  Work(iDensOccFro(iJ,iKFroz,iSymJ)) = Work(iDensOccFro(iJ,iKFroz,iSymJ))+T1*X2/Dnom2
                  Work(iWDensOccFro(iJ,iKFroz,iSymP)) = Work(iWDensOccFro(iJ,iKFroz,iSymP))-Two*T1*X2
                  !************************************************
                  ! Calculate P_ij active occ - active occ
                  !************************************************
                else if (iP <= nFro(iSymP)+nOcc(iSymP)) then
                  iK = iP-nFro(iSymP)
                  Dnom2 = EOcc(iOcc(iSymI)+iI)-EVir(iVir(iSymA)+iA)+EOcc(iOcc(iSymP)+iK)-EVir(iVir(iSymB)+iB)
                  Work(iDensOcc(iJ,iK,iSymP)) = Work(iDensOcc(iJ,iK,iSymP))-T1*X2/Dnom2
                  Work(iWDensOcc(iJ,iK,iSymP)) = Work(iWDensOcc(iJ,iK,iSymP))-Two*T1*X2
                  !************************************************
                  ! Calculate Lagr(4)_ai all virtual - active occupied
                  !************************************************
                else
                  iC = iP-nFro(iSymP)-nOcc(iSymP)
                  Work(iMp2Lagr(iC,iJ+nFro(iSymJ),iSymJ)) = Work(iMp2Lagr(iC,iJ+nFro(iSymJ),iSymJ))+T1*X2
                end if

                !---------- Debugging comments -------------------------
                !if ((iP > nFro(iSymP)) .and. (iP <= nFro(iSymP)+nOcc(iSymP))) then
                !  write(u6,*) 'AIBJK',iA,iI,iB,iJ,iK
                !  write(u6,*) 'Symm',iSymA,iSymI,iSymB,iSymJ,iSymP
                !  write(u6,*) 'Dnom1',Dnom1
                !  write(u6,*) 'Dnom2',Dnom2
                !  write(u6,*) 'iajb',Xaibj(ip_iajb)
                !  write(u6,*) 'ibja',Xaibj(ip_ibja)
                !  write(u6,*) 'iakb',Xaibj(ip_iapb)
                !  write(u6,*) 'Bidrag',T1*T2
                !  write(u6,*) 'adress',iDensOcc(iJ,iK,iSymP)-ip_Density(1)
                !end if
                !-------------------------------------------------------
                !---------- Debugging comments -------------------------
                !write(u6,*) 'AIBJC',iA,iI,iB,iJ,iC
                !write(u6,*) 'Symm',iSymA,iSymI,iSymB,iSymJ,iSymC
                !write(u6,*) 'Dnom1',Dnom1
                !write(u6,*) 'aibj',Xaibj(ip_iajb)
                !write(u6,*) 'biaj',Xaibj(ip_ibja)
                !write(u6,*) 'iacb',Xaibj(ip_iacb)
                !if (iBatch /= jBatch) then
                !   write(u6,*) 'iacb',Xaibj(ip_ibca)
                !   write(u6,*) 'Bidrag U',U1*U2
                !end if
                !write(u6,*) 'Bidrag T',T1*T2
                !write(u6,*) 'adress',iMp2Lagr(iC,iJ,iSymJ)-ip_Mp2Lagr(1)
                !-------------------------------------------------------

              end do
            end do
          end do
        end do
      end do
    end do
  end do
end do

return

end subroutine ChoMP2_rhslagr_1
