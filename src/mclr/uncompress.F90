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

subroutine UnCompress(ArrayIn,ArrayOut,idsym)
! Uncompresses the PCG vector to a orbital rotation matrix
!
! The redundant rotations are set to zero

use Symmetry_Info, only: Mul
use MCLR_Data, only: ipMat, nB, nDens, nDensC
use input_mclr, only: nBas, nIsh, nOrb, nRS1, nRS2, nRS3, nSym, TimeDep
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: ArrayIn(nDensC)
real(kind=wp), intent(out) :: ArrayOut(nDens)
integer(kind=iwp), intent(in) :: idsym
integer(kind=iwp) :: Bas(8), dsym, iBas, Index1, Index2, IndexC, iSym, iT, jBas, jSym, jT
real(kind=wp) :: Fact

indexC = 0
Fact = One
if (idsym < 0) Fact = -Fact
dsym = abs(idsym)
ArrayOut(:) = Zero
if (TimeDep) then
  Bas(1:nSym) = nBas(1:nSym)
else
  Bas(1:nSym) = nB(1:nSym)
end if
do iSym=1,nSym
  do jSym=1,nSym
    if (Mul(iSym,jSym) == dSym) then
      do jBas=1,Bas(jSym)
        if (jBas <= nIsh(jsym)) then
          jT = 0
        else if (jBas <= nIsh(jsym)+nRs1(jsym)) then
          jT = 1
        else if (jBas <= nIsh(jsym)+nRs1(jsym)+nRs2(jsym)) then
          jT = 2
        else if (jBas <= nIsh(jsym)+nRs1(jsym)+nRs2(jsym)+nRs3(jsym)) then
          jT = 3
        else
          jT = 4
        end if
        do iBas=1,nOrb(iSym)
          if (iBas <= nIsh(isym)) then
            iT = 0
          else if (iBas <= nIsh(isym)+nRs1(isym)) then
            iT = 1
          else if (iBas <= nIsh(isym)+nRs1(isym)+nRs2(isym)) then
            iT = 2
          else if (iBas <= nIsh(isym)+nRs1(isym)+nRs2(isym)+nRs3(isym)) then
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
              ArrayOut(Index2) = -Fact*ArrayIn(indexC)
            end if
          end if
        end do
      end do
    end if
  end do
end do

end subroutine UnCompress
