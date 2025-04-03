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

subroutine UnC(ArrayIn,ArrayOut,dsym,Sign)
! Uncompresses the PCG vector to a orbital rotation matrix
!
! The redundant rotations are set to zero

use Symmetry_Info, only: Mul
use MCLR_Data, only: nDensC, nDens2, ipMat, nB
use input_mclr, only: nSym, nIsh, nRS1, nRS2, nRS3, nOrb, TimeDep
use Constants, only: Zero, One

implicit none
real*8 ArrayIn(nDensC), ArrayOut(nDens2)
integer dsym
real*8 sign
integer IndexC, iSym, jSym, jBas, jT, iBas, iT, Index1, Index2
real*8 Fact

indexC = 0
Fact = One
if (dsym < 0) Fact = -Fact
dsym = abs(dsym)
ArrayOut(:) = Zero
do iSym=1,nSym
  do jSym=1,nSym
    if (Mul(iSym,jSym) == dSym) then
      do jBas=1,nB(jSym)
        if (jBas <= nIsh(jsym)) then
          jT = 0
        else if (jBas <= nIsh(jsym)+nRs1(jsym)) then
          jT = 1
        else if (jBas <= nIsh(jsym)+nRs2(jsym)) then
          jT = 2
        else if (jBas <= nIsh(jsym)+nRs3(jsym)) then
          jT = 3
        else
          jT = 4
        end if
        do iBas=1,nOrb(iSym)
          if (iBas <= nIsh(isym)) then
            iT = 0
          else if (iBas <= nIsh(isym)+nRs1(isym)) then
            iT = 1
          else if (iBas <= nIsh(isym)+nRs2(isym)) then
            iT = 2
          else if (iBas <= nIsh(isym)+nRs3(isym)) then
            iT = 3
          else
            iT = 4
          end if
          if (Timedep) then
            if (iT /= jT) then
              indexC = indexc+1
              Index1 = ipMat(iSym,jSym)+(jBas-1)*nOrb(iSym)+iBas-1
              ArrayOut(Index1) = Fact*ArrayIn(indexC)
            end if
          else
            if (iT > jT) then
              indexC = indexc+1
              Index1 = ipMat(iSym,jSym)+(jBas-1)*nOrb(iSym)+iBas-1
              Index2 = ipMat(jSym,iSym)+(iBas-1)*nOrb(jSym)+jBas-1
              ArrayOut(Index1) = Fact*ArrayIn(indexC)
              ArrayOut(Index2) = Sign*Fact*ArrayIn(indexC)
            end if
          end if
        end do
      end do
    end if
  end do
end do

end subroutine UnC
