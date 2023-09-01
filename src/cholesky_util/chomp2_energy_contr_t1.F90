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
! Copyright (C) 2004,2005, Thomas Bondo Pedersen                       *
!               Francesco Aquilante                                    *
!***********************************************************************

subroutine ChoMP2_Energy_Contr_T1(EMP2,EOcc,EVir,Xaibj,LnT2am,LiT2am,iBatch,jBatch)
!
! Thomas Bondo Pedersen, Dec. 2004 / Feb. 2005.
!
! Purpose: compute (MINUS the) energy contribution from a
!          batch of (ai|bj) integrals.
!
! Modified by F. Aquilante to add contributions from T1 amplitudes
!                          determined from Thouless formula

use Symmetry_Info, only: Mul
use Index_Functions, only: iTri
use Cholesky, only: nSym
use ChoMP2, only: ChoAlg, EOSMP2, iFirstS, iMatab, iOcc, iOffT1, iVir, LiMatij, LiT1am, LnOcc, LnT1am, nMatab, nOcc, nVir, T1amp, &
                  Wref
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LnT2am, LiT2am(8), iBatch, jBatch
real(kind=wp), intent(inout) :: EMP2, Xaibj(LnT2am)
real(kind=wp), intent(in) :: EOcc(*), EVir(*)
integer(kind=iwp) :: a, abij, aibj, b, baij, biaj, i, ia, ij, iSym1, iSym2, iSyma, iSymab, iSymai, iSymaj, iSymb, iSymbi, iSymbj, &
                     iSymi, iSymij, iSymj, j, jb, Lai, Laj, Lbi, Lbj, Li, Lj
real(kind=wp) :: Dnom, Eaibj, EMP2_sav, EOSMP2_sav, Taibj, Waibj, WREF_sav, Xtmp

call Cho_GAdGOp(Xaibj,LnT2am,'+')

if (iBatch == jBatch) then

  if (ChoAlg == 2) then ! M(ab,ij)=(ai|bj) with i<=j.
    do iSymj=1,nSym
      iSymi = iSymj
      do Lj=1,LnOcc(iSymj,jBatch)
        j = iFirstS(iSymj,jBatch)+Lj-1
        do Li=1,LnOcc(iSymi,iBatch)
          i = iFirstS(iSymi,iBatch)+Li-1
          ij = LiMatij(iSymi,iSymj,iBatch)+iTri(Li,Lj)
          do iSymb=1,nSym
            iSyma = iSymb
            if (iSymj == iSymb) then
              do b=1,nVir(iSymb)
                abij = LiT2am(1)+nMatab(1)*(ij-1)+iMatab(iSyma,iSymb)+nVir(iSyma)*(b-1)
                baij = LiT2am(1)+nMatab(1)*(ij-1)+iMatab(iSymb,iSyma)-nVir(iSymb)+b
                do a=1,nVir(iSyma)
                  abij = abij+1
                  baij = baij+nVir(iSymb)
                  Dnom = EVir(iVir(iSyma)+a)-EOcc(iOcc(iSymi)+i)+EVir(iVir(iSymb)+b)-EOcc(iOcc(iSymj)+j)
                  ia = iOffT1(iSymi)+nOcc(iSymi)*(a-1)+i
                  jb = iOffT1(iSymj)+nOcc(iSymj)*(b-1)+j
                  Xtmp = Xaibj(abij)
                  Taibj = Xtmp/Dnom-T1amp(ia)*T1amp(jb)
                  Waibj = Two*Xtmp
                  EOSMP2 = EOSMP2+Taibj*Waibj
                  Waibj = Waibj-Xaibj(baij)
                  Eaibj = Taibj*Waibj
                  EMP2 = EMP2+Eaibj
                  WREF = WREF+Eaibj/Dnom
                end do
              end do
            else
              do b=1,nVir(iSymb)
                abij = LiT2am(1)+nMatab(1)*(ij-1)+iMatab(iSyma,iSymb)+nVir(iSyma)*(b-1)
                baij = LiT2am(1)+nMatab(1)*(ij-1)+iMatab(iSymb,iSyma)-nVir(iSymb)+b
                do a=1,nVir(iSyma)
                  abij = abij+1
                  baij = baij+nVir(iSymb)
                  Dnom = EVir(iVir(iSyma)+a)-EOcc(iOcc(iSymi)+i)+EVir(iVir(iSymb)+b)-EOcc(iOcc(iSymj)+j)
                  Xtmp = Xaibj(abij)
                  Taibj = Xtmp/Dnom
                  Waibj = Two*Xtmp
                  EOSMP2 = EOSMP2+Taibj*Waibj
                  Waibj = Waibj-Xaibj(baij)
                  Eaibj = Taibj*Waibj
                  EMP2 = EMP2+Eaibj
                  WREF = WREF+Eaibj/Dnom
                end do
              end do
            end if
          end do
        end do
      end do
    end do
    do iSymij=2,nSym
      iSymab = iSymij
      do iSym2=1,nSym
        iSym1 = Mul(iSym2,iSymij)
        iSymj = max(iSym1,iSym2)
        iSymi = min(iSym1,iSym2)
        do Lj=1,LnOcc(iSymj,jBatch)
          j = iFirstS(iSymj,jBatch)+Lj-1
          do Li=1,LnOcc(iSymi,iBatch)
            i = iFirstS(iSymi,iBatch)+Li-1
            ij = LiMatij(iSymi,iSymj,iBatch)+LnOcc(iSymi,iBatch)*(Lj-1)+Li
            do iSymb=1,nSym
              iSyma = Mul(iSymb,iSymab)
              if ((iSymi == iSyma) .and. (iSymj == iSymb)) then
                do b=1,nVir(iSymb)
                  abij = LiT2am(iSymij)+nMatab(iSymab)*(ij-1)+iMatab(iSyma,iSymb)+nVir(iSyma)*(b-1)
                  baij = LiT2am(iSymij)+nMatab(iSymab)*(ij-1)+iMatab(iSymb,iSyma)-nVir(iSymb)+b
                  do a=1,nVir(iSyma)
                    abij = abij+1
                    baij = baij+nVir(iSymb)
                    Dnom = EVir(iVir(iSyma)+a)-EOcc(iOcc(iSymi)+i)+EVir(iVir(iSymb)+b)-EOcc(iOcc(iSymj)+j)
                    ia = iOffT1(iSymi)+nOcc(iSymi)*(a-1)+i
                    jb = iOffT1(iSymj)+nOcc(iSymj)*(b-1)+j
                    Xtmp = Xaibj(abij)
                    Taibj = Xtmp/Dnom-T1amp(ia)*T1amp(jb)
                    Waibj = Two*Xtmp
                    EOSMP2 = EOSMP2+Taibj*Waibj
                    Waibj = Waibj-Xaibj(baij)
                    Eaibj = Taibj*Waibj
                    EMP2 = EMP2+Eaibj
                    WREF = WREF+Eaibj/Dnom
                  end do
                end do
              else
                do b=1,nVir(iSymb)
                  abij = LiT2am(iSymij)+nMatab(iSymab)*(ij-1)+iMatab(iSyma,iSymb)+nVir(iSyma)*(b-1)
                  baij = LiT2am(iSymij)+nMatab(iSymab)*(ij-1)+iMatab(iSymb,iSyma)-nVir(iSymb)+b
                  do a=1,nVir(iSyma)
                    abij = abij+1
                    baij = baij+nVir(iSymb)
                    Dnom = EVir(iVir(iSyma)+a)-EOcc(iOcc(iSymi)+i)+EVir(iVir(iSymb)+b)-EOcc(iOcc(iSymj)+j)
                    Xtmp = Xaibj(abij)
                    Taibj = Xtmp/Dnom
                    Waibj = Two*Xtmp
                    EOSMP2 = EOSMP2+Taibj*Waibj
                    Waibj = Waibj-Xaibj(baij)
                    Eaibj = Taibj*Waibj
                    EMP2 = EMP2+Eaibj
                    WREF = WREF+Eaibj/Dnom
                  end do
                end do
              end if
            end do
          end do
        end do
      end do
    end do
  else ! triangular storage, (ai|bj) with ai<=bj.
    do iSymbj=1,nSym
      iSymai = iSymbj
      do iSymj=1,nSym
        iSymb = Mul(iSymj,iSymbj)
        do Lj=1,LnOcc(iSymj,jBatch)
          j = iFirstS(iSymj,jBatch)+Lj-1
          do b=1,nVir(iSymb)
            Lbj = LiT1am(iSymb,iSymj,jBatch)+nVir(iSymb)*(Lj-1)+b
            do iSymi=1,nSym
              iSyma = Mul(iSymi,iSymai)
              iSymaj = Mul(iSyma,iSymj)
              iSymbi = iSymaj
              if ((iSymi == iSyma) .and. (iSymj == iSymb)) then
                do Li=1,LnOcc(iSymi,iBatch)
                  i = iFirstS(iSymi,iBatch)+Li-1
                  Lbi = LiT1am(iSymb,iSymi,iBatch)+nVir(iSymb)*(Li-1)+b
                  do a=1,nVir(iSyma)
                    Lai = LiT1am(iSyma,iSymi,iBatch)+nVir(iSyma)*(Li-1)+a
                    Laj = LiT1am(iSyma,iSymj,jBatch)+nVir(iSyma)*(Lj-1)+a
                    aibj = LiT2am(iSymai)+iTri(Lai,Lbj)
                    biaj = LiT2am(iSymbi)+iTri(Lbi,Laj)
                    Dnom = EVir(iVir(iSyma)+a)-EOcc(iOcc(iSymi)+i)+EVir(iVir(iSymb)+b)-EOcc(iOcc(iSymj)+j)
                    ia = iOffT1(iSymi)+nOcc(iSymi)*(a-1)+i
                    jb = iOffT1(iSymj)+nOcc(iSymj)*(b-1)+j
                    Xtmp = Xaibj(aibj)
                    Taibj = Xtmp/Dnom-T1amp(ia)*T1amp(jb)
                    Waibj = Two*Xtmp
                    EOSMP2 = EOSMP2+Taibj*Waibj
                    Waibj = Waibj-Xaibj(biaj)
                    Eaibj = Taibj*Waibj
                    EMP2 = EMP2+Eaibj
                    WREF = WREF+Eaibj/Dnom
                  end do
                end do
              else
                do Li=1,LnOcc(iSymi,iBatch)
                  i = iFirstS(iSymi,iBatch)+Li-1
                  Lbi = LiT1am(iSymb,iSymi,iBatch)+nVir(iSymb)*(Li-1)+b
                  do a=1,nVir(iSyma)
                    Lai = LiT1am(iSyma,iSymi,iBatch)+nVir(iSyma)*(Li-1)+a
                    Laj = LiT1am(iSyma,iSymj,jBatch)+nVir(iSyma)*(Lj-1)+a
                    aibj = LiT2am(iSymai)+iTri(Lai,Lbj)
                    biaj = LiT2am(iSymbi)+iTri(Lbi,Laj)
                    Dnom = EVir(iVir(iSyma)+a)-EOcc(iOcc(iSymi)+i)+EVir(iVir(iSymb)+b)-EOcc(iOcc(iSymj)+j)
                    Xtmp = Xaibj(aibj)
                    Taibj = Xtmp/Dnom
                    Waibj = Two*Xtmp
                    EOSMP2 = EOSMP2+Taibj*Waibj
                    Waibj = Waibj-Xaibj(biaj)
                    Eaibj = Taibj*Waibj
                    EMP2 = EMP2+Eaibj
                    WREF = WREF+Eaibj/Dnom
                  end do
                end do
              end if
            end do
          end do
        end do
      end do
    end do
  end if

else ! rectangular storage (ai|bj) with ai<bj.

  EMP2_sav = EMP2
  EMP2 = Zero
  WREF_sav = WREF
  WREF = Zero
  EOSMP2_sav = EOSMP2
  EOSMP2 = Zero
  do iSymbj=1,nSym
    iSymai = iSymbj
    do iSymj=1,nSym
      iSymb = Mul(iSymj,iSymbj)
      do Lj=1,LnOcc(iSymj,jBatch)
        j = iFirstS(iSymj,jBatch)+Lj-1
        do b=1,nVir(iSymb)
          Lbj = LiT1am(iSymb,iSymj,jBatch)+nVir(iSymb)*(Lj-1)+b
          do iSymi=1,nSym
            iSyma = Mul(iSymi,iSymai)
            iSymaj = Mul(iSyma,iSymj)
            iSymbi = iSymaj
            if ((iSymi == iSyma) .and. (iSymj == iSymb)) then
              do Li=1,LnOcc(iSymi,iBatch)
                i = iFirstS(iSymi,iBatch)+Li-1
                Lbi = LiT1am(iSymb,iSymi,iBatch)+nVir(iSymb)*(Li-1)+b
                do a=1,nVir(iSyma)
                  Lai = LiT1am(iSyma,iSymi,iBatch)+nVir(iSyma)*(Li-1)+a
                  Laj = LiT1am(iSyma,iSymj,jBatch)+nVir(iSyma)*(Lj-1)+a
                  aibj = LiT2am(iSymai)+LnT1am(iSymai,iBatch)*(Lbj-1)+Lai
                  biaj = LiT2am(iSymbi)+LnT1am(iSymbi,iBatch)*(Laj-1)+Lbi
                  Dnom = EVir(iVir(iSyma)+a)-EOcc(iOcc(iSymi)+i)+EVir(iVir(iSymb)+b)-EOcc(iOcc(iSymj)+j)
                  ia = iOffT1(iSymi)+nOcc(iSymi)*(a-1)+i
                  jb = iOffT1(iSymj)+nOcc(iSymj)*(b-1)+j
                  Xtmp = Xaibj(aibj)
                  Taibj = Xtmp/Dnom-T1amp(ia)*T1amp(jb)
                  Waibj = Two*Xtmp
                  EOSMP2 = EOSMP2+Taibj*Waibj
                  Waibj = Waibj-Xaibj(biaj)
                  Eaibj = Taibj*Waibj
                  EMP2 = EMP2+Eaibj
                  WREF = WREF+Eaibj/Dnom
                end do
              end do
            else
              do Li=1,LnOcc(iSymi,iBatch)
                i = iFirstS(iSymi,iBatch)+Li-1
                Lbi = LiT1am(iSymb,iSymi,iBatch)+nVir(iSymb)*(Li-1)+b
                do a=1,nVir(iSyma)
                  Lai = LiT1am(iSyma,iSymi,iBatch)+nVir(iSyma)*(Li-1)+a
                  Laj = LiT1am(iSyma,iSymj,jBatch)+nVir(iSyma)*(Lj-1)+a
                  aibj = LiT2am(iSymai)+LnT1am(iSymai,iBatch)*(Lbj-1)+Lai
                  biaj = LiT2am(iSymbi)+LnT1am(iSymbi,iBatch)*(Laj-1)+Lbi
                  Dnom = EVir(iVir(iSyma)+a)-EOcc(iOcc(iSymi)+i)+EVir(iVir(iSymb)+b)-EOcc(iOcc(iSymj)+j)
                  Xtmp = Xaibj(aibj)
                  Taibj = Xtmp/Dnom
                  Waibj = Two*Xtmp
                  EOSMP2 = EOSMP2+Taibj*Waibj
                  Waibj = Waibj-Xaibj(biaj)
                  Eaibj = Taibj*Waibj
                  EMP2 = EMP2+Eaibj
                  WREF = WREF+Eaibj/Dnom
                end do
              end do
            end if
          end do
        end do
      end do
    end do
  end do
  EMP2 = EMP2_sav+Two*EMP2
  EOSMP2 = EOSMP2_sav+Two*EOSMP2
  WREF = WREF_sav+Two*WREF

end if

EOSMP2 = Half*EOSMP2

end subroutine ChoMP2_Energy_Contr_T1
