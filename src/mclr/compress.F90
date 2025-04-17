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

subroutine Compress(ArrayIn,ArrayOut,dsym)
! Compresses the orbital rotation matrix to
! the vector that is used in the PCG routines
! the indexes are ordered to fit the preconditioner
!
! The redundant rotations are set to zero

use Symmetry_Info, only: Mul
use MCLR_Data, only: ipMat, nDens, nDensC
use input_mclr, only: nIsh, nOrb, nRs1, nRs2, nRs3, nSym, TimeDep
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: ArrayIn(nDens)
real(kind=wp), intent(out) :: ArrayOut(nDensC)
integer(kind=iwp), intent(in) :: dsym
integer(kind=iwp) :: iBas, index1, indexC, isym, iT, jBas, jsym, jT

indexC = 0
ArrayOut(:) = Zero
do iSym=1,nSym
  do jSym=1,nSym
    if (Mul(iSym,jSym) == abs(dSym)) then
      do jBas=1,nOrb(jSym)
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
          if (TimeDep) then
            if (iT /= jT) then !
              indexC = indexc+1 !
              Index1 = ipMat(iSym,jSym)+(jBas-1)*nOrb(iSym)+iBas-1 !
              ArrayOut(IndexC) = ArrayIn(index1) !
            end if !
          else
            if (iT > jT) then
              indexC = indexc+1
              Index1 = ipMat(iSym,jSym)+(jBas-1)*nOrb(iSym)+iBas-1
              ArrayOut(IndexC) = ArrayIn(index1)
            end if
          end if
        end do
      end do
    end if
  end do
end do
if (indexc /= nDensC) call SysAbendMsg('compress','indexc /= nDensC',' ')

end subroutine Compress
