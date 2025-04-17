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
! Copyright (C) 1997, Anders Bernhardsson                              *
!***********************************************************************

subroutine UnCompress2(ArrayIn,ArrayOut,dsym)
!************************************************************
!
! Change from vector to kappa   notation
!                            ip
!
! where i is occupied and p is general
!
! used by the preconditioner.
!
!************************************************************

use Symmetry_Info, only: Mul
use MCLR_Data, only: ipMat, nB, nDens, nDensC
use input_mclr, only: nIsh, nOrb, nRS1, nRS2, nRS3, nSym, TimeDep
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: ArrayIn(nDensC)
real(kind=wp), intent(out) :: ArrayOut(nDens)
integer(kind=iwp), intent(inout) :: dsym
integer(kind=iwp) :: i1, iBas, ij, Index1, Index2, IndexC, iSym, iT, j1, jBas, ji, jSym, jT
real(kind=wp) :: Fact

indexC = 0
Fact = One
if (dsym < 0) Fact = -Fact
dsym = abs(dsym)
ArrayOut(:) = Zero
jT = 0 ! dummy initialize
i1 = 0 ! dummy initialize
j1 = 0 ! dummy initialize
do iSym=1,nSym
  do jSym=1,nSym
    if (Mul(iSym,jSym) == dSym) then
      do jBas=1,nB(jSym)
        if (jBas <= nIsh(jsym)) then
          jT = 0
          i1 = nIsh(isym)
          j1 = nIsh(jsym)
        else if (jBas <= nIsh(jsym)+nRs1(jsym)) then
          jT = 1
          i1 = nRs1(isym)
          j1 = nRs1(jsym)
        else if (jBas <= nIsh(jsym)+nRs1(jsym)+nRs2(jsym)) then
          jT = 2
          i1 = nRs2(isym)
          j1 = nRs2(jsym)
        else if (jBas <= nIsh(jsym)+nRs1(jsym)+nRs2(jsym)+nRs3(jsym)) then
          jT = 3
          i1 = nRs3(isym)
          j1 = nRs3(jsym)
        end if
        do iBas=1,nOrb(iSym)
          if (iBas <= nIsh(isym)) then
            iT = 0
            ij = 0
            ji = j1
          else if (iBas <= nIsh(isym)+nRs1(isym)) then
            iT = 1
            ij = 0
            ji = j1
            if (it > jt) then
              ij = i1
              ji = 0
            end if
          else if (iBas <= nIsh(isym)+nRs1(isym)+nRs2(isym)) then
            iT = 2
            ij = 0
            ji = j1
            if (it > jt) then
              ij = i1
              ji = 0
            end if
          else if (iBas <= nIsh(isym)+nRs1(isym)+nRs2(isym)+nRs3(isym)) then
            iT = 3
            ij = 0
            ji = j1
            if (it > jt) then
              ij = i1
              ji = 0
            end if
          else
            iT = 4
            ij = i1
            ji = 0
          end if
          if (Timedep) then
            if (iT /= jT) then
              indexC = indexc+1
              Index1 = ipMat(iSym,jSym)+(jBas-1)*norb(iSym)+iBas-1
              ArrayOut(Index1) = Fact*ArrayIn(indexC)
            end if
          else
            if (iT > jT) then
              indexC = indexc+1
              Index1 = ipMat(iSym,jSym)+(jBas-1)*norb(iSym)+iBas-ij-1
              Index2 = ipMat(jSym,iSym)+(iBas-1)*norb(jSym)+jBas-ji-1
              ArrayOut(Index1) = Fact*ArrayIn(indexC)
              if (iBas <= nB(isym)) ArrayOut(Index2) = Fact*ArrayIn(indexC)
            end if
          end if
        end do
      end do
    end if
  end do
end do

end subroutine UnCompress2
