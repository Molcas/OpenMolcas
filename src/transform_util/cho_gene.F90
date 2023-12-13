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
! Copyright (C) 2005, Giovanni Ghigo                                   *
!***********************************************************************
!  Cho_GenE
!
!> @brief
!>   Routine for the generation of the ``A,B`` block of exchange integrals (in symmetries \p iSymA, \p iSymB)
!>   for occupied MO \p iI, \p iJ in symmetries \p iSymI, \p iSymJ
!> @author Giovanni Ghigo
!>
!> @details
!> The routine generates the ``A,B`` block of integrals gathering 9
!> sub-blocks. These are combination of inactive, active, and
!> secondary ``A,B`` MO.
!>
!> @param[in] iSymI,iSymJ,iSymA,iSymB Symmetry block of the two-electrons integrals
!> @param[in] iI,iJ
!> @param[in] NumV                    Number of Cholesky vectors to transform in the current batch
!> @param[in,out] AddEx               Array of the ``A,B`` integrals block
!> @param[in] LenEx                   Length of the ``A,B`` integrals block
!***********************************************************************

subroutine Cho_GenE(iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV,AddEx,LenEx)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           January-June 2005                                          *
!----------------------------------------------------------------------*
! This is the routine that really generates the A,B block of exchange  *
! integrals (in symmetries iSymA,iSymB) for occupied MO iI,iJ in       *
! symmetries iSymI,iSymJ. The 3 x 3 A,B block is built gathering 9     *
! sub-blocks. These are combination of inactive, active, and secondary *
! A,B MO                                                               *
!***********************************************************************

use Cho_Tra, only: DoTCVA, nAsh, nIsh, nSsh, SubBlocks
use MkSubs, only: MkExSB11, MkExSB12, MkExSB13, MkExSB21, MkExSB22, MkExSB23, MkExSB31, MkExSB32, MkExSB33
use Data_Structures, only: Alloc1DArray_Type
use stdalloc, only: mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV, LenEx
real(kind=wp), intent(inout) :: AddEx(LenEx)
integer(kind=iwp) :: iA, iAddExSB, iAddSBi, iB, iLenA, iLenAa, iLenB, iLenBb, iSB_A, iSB_B, LenA(3), LenB(3)
type(Alloc1DArray_Type) :: AddSB(3,3)

!                                                                      *
!***********************************************************************
!                                                                      *
! GENERATION of SubBlocks
if (SubBlocks(1,1)) call MkExSB11(AddSB(1,1)%A,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
if (SubBlocks(1,2)) call MkExSB12(AddSB(1,2)%A,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
if (SubBlocks(1,3)) call MkExSB13(AddSB(1,3)%A,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
if (SubBlocks(2,1)) call MkExSB21(AddSB(2,1)%A,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV,AddSB(1,2)%A)
if (SubBlocks(2,2)) call MkExSB22(AddSB(2,2)%A,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
if (SubBlocks(2,3)) call MkExSB23(AddSB(2,3)%A,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
if (SubBlocks(3,1)) call MkExSB31(AddSB(3,1)%A,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV,AddSB(1,3)%A)
if (SubBlocks(3,2)) call MkExSB32(AddSB(3,2)%A,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV,AddSB(2,3)%A)
if (SubBlocks(3,3)) call MkExSB33(AddSB(3,3)%A,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
! END GENERATION of SubBlocks

! GATHERING of SubBlocks
iAddExSB = 1
if (DoTCVA) then

  LenB(1) = nIsh(iSymB)
  LenB(2) = nAsh(iSymB)
  LenB(3) = nSsh(iSymB)
  LenA(1) = nIsh(iSymA)
  LenA(2) = nAsh(iSymA)
  LenA(3) = nSsh(iSymA)
  if (iSymA /= iSymB) then

    do iSB_A=1,3
      do iA=1,LenA(iSB_A)
        do iSB_B=1,3
          if (LenB(iSB_B) == 0) cycle ! SB(1,iSB_A)

          iAddSBi = 1+LenB(iSB_B)*(iA-1)
          call daxpy_(LenB(iSB_B),One,AddSB(iSB_A,iSB_B)%A(iAddSBi),1,AddEx(iAddExSB),1)
          iAddExSB = iAddExSB+LenB(iSB_B)
        end do
      end do  ! iA
    end do ! iSB_A

  else

    do iSB_B=1,3
      do iB=1,LenB(iSB_B)
        do iSB_A=1,3
          if (LenA(iSB_A) == 0) cycle ! SB(1,iSB_B)

          iAddSBi = 1+LenA(iSB_A)*(iB-1)
          call daxpy_(LenA(iSB_A),One,AddSB(iSB_B,iSB_A)%A(iAddSBi),1,AddEx(iAddExSB),1)
          iAddExSB = iAddExSB+LenA(iSB_A)
        end do  ! iSB_A

      end do  ! iB
    end do ! iSB_B

  end if

else

  if (iSymA /= iSymB) then
    iLenBb = nSsh(iSymB)
    if (iLenBb > 0) then ! SB(3,3)
      iLenA = nSsh(iSymA)
      do iA=1,iLenA
        iAddSBi = 1+iLenBb*(iA-1)
        call daxpy_(iLenBb,One,AddSB(3,3)%A(iAddSBi),1,AddEx(iAddExSB),1)
        iAddExSB = iAddExSB+iLenBb
      end do  ! iA
    end if
  else
    iLenAa = nSsh(iSymA)
    if (iLenAa > 0) then ! SB(3,3)
      iLenB = nSsh(iSymB)
      do iB=1,iLenB
        iAddSBi = 1+iLenAa*(iB-1)
        call daxpy_(iLenAa,One,AddSB(3,3)%A(iAddSBi),1,AddEx(iAddExSB),1)
        iAddExSB = iAddExSB+iLenAa
      end do  ! iB
    end if
  end if

end if
! END GATHERING of SubBlocks

do iSB_A=1,3
  do iSB_B=1,3
    if (allocated(AddSB(iSB_A,iSB_B)%A)) call mma_deallocate(AddSB(iSB_A,iSB_B)%A)
  end do
end do

return

end subroutine Cho_GenE
