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

module RI_glob

use Data_Structures, only: Alloc1DiArray_Type, DSBA_Type
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: iMP2prpt, iAdrCVec(8,8,2), iOff_Ymnij(8,5), iOffA(4,0:7), iOpt, ip_Chunk = 0, iRsv, klS, Lu_A(0:7), &
                     Lu_Q(0:7), LuAVector(2), LuBVector(2), LuCVector(8,2), MxChVInShl, nAdens, nAuxVe, nAvec, nChOrb(0:7,5), &
                     nChV(0:7), nIJ1(8,8,2), nIJR(8,8,2), nJdens, nKdens, nKvec, nScreen, nSO, nSkal_Valence, nTask, NumAuxVec(8), &
                     nYmnij(8,5)
real(kind=wp) :: dmpK = Zero, tavec(2), tbvec(2)
logical(kind=iwp) :: DoCholExch, Timings_default
type(Alloc1DiArray_Type) :: Ymnij(5)
type(DSBA_Type), target :: CMOi(5), DMLT(5)
integer(kind=iwp), allocatable :: iBDsh(:), iMap(:), iShij(:,:), iSSOff(:,:,:), nBasSh(:,:), ShlSO(:), SO2Ind(:), SOShl(:), &
                                  TskList(:)
real(kind=wp), allocatable :: A(:), AMP2(:,:), BMP2(:,:), Chunk(:)
real(kind=wp), allocatable, target :: BklK(:), CijK(:), CilK(:), Yij(:,:,:)

public :: A, AMP2, BklK, BMP2, Chunk, CijK, CilK, CMOi, DMLT, dmpK, DoCholExch, iAdrCVec, iBDsh, iMap, iMP2prpt, iOff_Ymnij, &
          iOffA, iOpt, ip_Chunk, iRsv, iShij, iSSOff, klS, Lu_A, Lu_Q, LuAVector, LuBVector, LuCVector, MxChVInShl, nAdens, &
          nAuxVe, nAvec, nBasSh, nChOrb, nChV, nIJ1, nIJR, nJdens, nKdens, nKvec, nScreen, nSkal_Valence, nSO, nTask, NumAuxVec, &
          nYmnij, ShlSO, SO2Ind, SOShl, tavec, tbvec, Timings_default, TskList, Yij, Ymnij

end module RI_glob
