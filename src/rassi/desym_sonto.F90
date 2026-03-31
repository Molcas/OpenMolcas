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
! Copyright (C) 2021, Rulin Feng                                       *
!***********************************************************************

!****************************************************
!          Desymmetrize a given array(matrix)
!****************************************************
! This routine is made to expand a given symmetry-adapted array
! into a C1 symmetry. SYMLAB is a bit flag,
! e.g., for symmetry 3, symlab=4=2^(3-1). A is input array with
! size SIZA, B is output array with size of NBST**2.
! Not tested for general cases, so for use of SO-NTOs only.
!
!                                               -RF 8/24,2021
subroutine DESYM_SONTO(A,SIZA,B,SYMLAB)

use Symmetry_Info, only: MUL, nIrrep
use rassi_data, only: NBASF, NBST
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: SIZA, SYMLAB
real(kind=wp) :: A(SIZA), B(NBST**2)
integer(kind=iwp) :: I, IJ, IOF, ISY1, ISY12_MA, ISY12_MA_BI, ISY2, ITD, J, JI, NB1, NB1_F, NB1_I, NB2, NB2_F, NB2_I
real(kind=wp) :: ME, TDM
real(kind=wp), allocatable :: SCR(:)

! Initialize
B(:) = Zero

call mma_allocate(SCR,SIZA,Label='SCR')
SCR(:) = Zero

if (SYMLAB == 1) then
  ! Diagonal symmetry blocks.
  ! Dont need to do anything, just leave it be
  call DCOPY_(SIZA,A(:),1,SCR,1)
else
  ! Non-diagonal symmetry blocks
  ! note that only half of the total matrix has been stored
  ITD = 0
  IOF = 0
  do ISY1=1,nIrrep
    NB1 = NBASF(ISY1)
    do ISY2=1,nIrrep
      NB2 = NBASF(ISY2)
      ISY12_ma = MUL(ISY1,ISY2)
      ISY12_ma_bi = 2**(ISY12_ma-1)
      if (ISY12_ma_bi == SYMLAB) then
        if (ISY1 > ISY2) then
          do J=1,NB2
            do I=1,NB1
              ITD = ITD+1
              TDM = A(ITD)
              IJ = IOF+J+NB2*(I-1)
              SCR(IJ) = TDM
            end do
          end do
          IOF = IOF+NB1*NB2
        end if
      end if
    end do
  end do
end if
! Expand into C1

ITD = 0
NB1_i = 0
NB1_f = 0
do ISY1=1,nIrrep
  NB1 = NBASF(ISY1)
  NB1_f = NB1_i+NB1
  NB2_i = 0
  NB2_f = 0
  do ISY2=1,nIrrep
    ISY12_ma = MUL(ISY1,ISY2)
    ISY12_ma_bi = 2**(ISY12_ma-1)
    NB2 = NBASF(ISY2)
    NB2_f = NB2_i+NB2
    if (ISY12_ma_bi == SYMLAB) then
      do J=NB2_i+1,NB2_f
        do I=NB1_i+1,NB1_f
          if (I <= J) then
            ITD = ITD+1
            ME = SCR(ITD)
            IJ = I+NBST*(J-1)
            JI = J+NBST*(I-1)
            B(IJ) = ME
            B(JI) = ME
          end if
        end do
      end do
    end if
    NB2_i = NB2_i+NB2
  end do
  NB1_i = NB1_i+NB1
end do
call mma_deallocate(SCR)

end subroutine DESYM_SONTO
