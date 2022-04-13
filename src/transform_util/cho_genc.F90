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
!  Cho_GenC
!
!> @brief
!>   Routine for the generation of the ``A,B`` block of Coulomb integrals (in symmetries \p iSymA, \p iSymB)
!>   for occupied MO \p iI, \p iJ in symmetries \p iSymI, \p iSymJ.
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
!> @param[in,out] AddCou              Array of the ``A,B`` integrals block
!> @param[in] LenCou                  Length of the ``A,B`` integrals block
!> @param[in] LenEx
!***********************************************************************

subroutine Cho_GenC(iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV,AddCou,LenCou,LenEx)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           July 2005                                                  *
!----------------------------------------------------------------------*
! This is the routine that really generates the A,B block of coulomb   *
! integrals (in symmetries iSymA,iSymB) for occupied MO iI,iJ in       *
! symmetries iSymI,iSymJ. The 3 x 3 A,B block is built gathering 9     *
! sub-blocks. These are combination of inactive, active, and secondary *
! A,B MO                                                               *
! OBS !!!!!  By now, it works only for iSymA == iSymB  !!!             *
!***********************************************************************

use Cho_Tra, only: IfTest, nAsh, nIsh, nOrb, nSsh, SubBlocks
use MkSubs, only: MkCouSB11, MkCouSB21, MkCouSB22, MkCouSB31, MkCouSB32, MkCouSB33
use Data_Structures, only: Alloc1DArray_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB, iI, iJ, numV, LenCou, LenEx
real(kind=wp), intent(inout) :: AddCou(LenCou)
integer(kind=iwp) :: iAddCouSB, iAddSBi, iB, iSB_A, iSB_B, LenA(3), LenB(3), LenSB, nOrbA
type(Alloc1DArray_Type) :: AddSB(3,3)
real(kind=wp), allocatable :: AddSq(:)

!                                                                      *
!***********************************************************************
!                                                                      *
!-----------------------------------------------------------------------
!IfTest = .true.
!-----------------------------------------------------------------------

LenA(1) = nIsh(iSymA)
LenA(2) = nAsh(iSymA)
LenA(3) = nSsh(iSymA)
LenB(1) = nIsh(iSymB)
LenB(2) = nAsh(iSymB)
LenB(3) = nSsh(iSymB)

! GENERATION of SubBlocks
!-----------------------------------------------------------------------
if (SubBlocks(1,1)) call MkCouSB11(AddSB(1,1)%A,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!-----------------------------------------------------------------------
if (IfTest .and. allocated(AddSB(1,1)%A)) then
  write(u6,*) '       SB_11 :',nIsh(iSymA),' x',nIsh(iSymB)
  write(u6,'(8F10.6)') AddSB(1,1)%A(:)
  call XFlush(u6)
end if
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
if (SubBlocks(1,2)) then
  if (iSymA /= iSymB) then
    ! excluded
    !call MkCouSB12(AddSB(1,2)%A,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
  else
    LenSB = nIsh(iSymA)*nAsh(iSymB)
    call mma_Allocate(AddSB(1,2)%A,LenSB,Label='AddSB12')
    AddSB(1,2)%A(:) = Zero
  end if
end if
!-----------------------------------------------------------------------
if (IfTest .and. allocated(AddSB(1,2)%A)) then
  write(u6,*) '       SB_12 :',nIsh(iSymA),' x',nAsh(iSymB)
  write(u6,'(8F10.6)') AddSB(1,2)%A(:)
  call XFlush(u6)
end if
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
if (SubBlocks(1,3)) then
  if (iSymA /= iSymB) then
    ! excluded
    !call MkCouSB13(AddSB(1,3)%A,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
  else
    LenSB = nIsh(iSymA)*nSsh(iSymB)
    call mma_allocate(AddSB(1,3)%A,LenSB,Label='AddSB13')
    AddSB(1,3)%A(:) = Zero
  end if
end if
!-----------------------------------------------------------------------
if (IfTest .and. allocated(AddSB(1,3)%A)) then
  write(u6,*) '       SB_13 :',nIsh(iSymA),' x',nSsh(iSymB)
  write(u6,'(8F10.6)') AddSB(1,3)%A(:)
  call XFlush(u6)
end if
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
if (SubBlocks(2,1)) call MkCouSB21(AddSB(2,1)%A,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!-----------------------------------------------------------------------
if (IfTest .and. allocated(AddSB(2,1)%A)) then
  write(u6,*) '       SB_21 :',nAsh(iSymA),' x',nIsh(iSymB)
  write(u6,'(8F10.6)') AddSB(2,1)%A(:)
  call XFlush(u6)
end if
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
if (SubBlocks(2,2)) call MkCouSB22(AddSB(2,2)%A,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!-----------------------------------------------------------------------
if (IfTest .and. allocated(AddSB(2,2)%A)) then
  write(u6,*) '       SB_22 :',nAsh(iSymA),' x',nAsh(iSymB)
  write(u6,'(8F10.6)') AddSB(2,2)%A(:)
  call XFlush(u6)
end if
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
if (SubBlocks(2,3)) then
  if (iSymA /= iSymB) then
    ! excluded
    !call MkCouSB23(AddSB(2,3)%A,iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
  else
    LenSB = nAsh(iSymA)*nSsh(iSymB)
    call mma_allocate(AddSB(2,3)%A,LenSB,Label='AddSB23')
    AddSB(2,3)%A(:) = Zero
  end if
end if
!-----------------------------------------------------------------------
if (IfTest .and. allocated(AddSB(2,3)%A)) then
  write(u6,*) '       SB_23 :',nAsh(iSymA),' x',nSsh(iSymB)
  write(u6,'(8F10.6)') AddSB(2,3)%A(:)
  call XFlush(u6)
end if
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
if (SubBlocks(3,1)) call MkCouSB31(AddSB(3,1)%A,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!-----------------------------------------------------------------------
if (IfTest .and. allocated(AddSB(3,1)%A)) then
  write(u6,*) '       SB_31 :',nSsh(iSymA),' x',nIsh(iSymB)
  write(u6,'(8F10.6)') AddSB(3,1)%A(:)
  call XFlush(u6)
end if
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
if (SubBlocks(3,2)) call MkCouSB32(AddSB(3,2)%A,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!-----------------------------------------------------------------------
if (IfTest .and. allocated(AddSB(3,2)%A)) then
  write(u6,*) '       SB_32 :',nSsh(iSymA),' x',nAsh(iSymB)
  write(u6,'(8F10.6)') AddSB(3,2)%A(:)
  call XFlush(u6)
end if
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
if (SubBlocks(3,3)) call MkCouSB33(AddSB(3,3)%A,iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV)
!-----------------------------------------------------------------------
if (IfTest .and. allocated(AddSB(3,3)%A)) then
  write(u6,*) '       SB_33 :',nSsh(iSymA),' x',nSsh(iSymB)
  write(u6,'(8F10.6)') AddSB(3,3)%A(:)
  call XFlush(u6)
end if
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
if (IfTest) then
  write(u6,*) '     END GENERATION of SubBlocks'
  call XFlush(u6)
end if
!-----------------------------------------------------------------------
! END GENERATION of SubBlocks

! GATHERING of SubBlocks
call mma_allocate(AddSq,LenEx,Label='AddSq')

iAddCouSB = 1

do iSB_B=1,3
  do iB=1,LenB(iSB_B)
    do iSB_A=1,3

      if (LenA(iSB_A) == 0) cycle

      ! SB(iSB_A,iSB_B)
      iAddSBi = 1+LenA(iSB_A)*(iB-1)
      call dCopy_(LenA(iSB_A),AddSB(iSB_B,iSB_A)%A(iAddSBi),1,AddSq(iAddCouSB),1)
      iAddCouSB = iAddCouSB+LenA(iSB_A)

    end do ! iSB_B
  end do  ! iB
end do ! iSB_B
nOrbA = nOrb(iSymA)
!-----------------------------------------------------------------------
if (IfTest) then
  write(u6,*)
  write(u6,*) '        The Square Gatered matrix'
  call PrintSquareMat(nOrbA,AddSq)
  call XFlush(u6)
end if
!-----------------------------------------------------------------------

call Local_Triang(nOrbA,AddSq)
call daxpy_(LenCou,One,AddSq,1,AddCou,1)
call mma_deallocate(AddSq)
!-----------------------------------------------------------------------
if (IfTest) then
  write(u6,*)
  write(u6,*) '        The Triangular Integrals matrix'
  call PrintTriangMat(nOrbA,AddCou)
  call XFlush(u6)
end if
!-----------------------------------------------------------------------
! END GATHERING of SubBlocks

do iSB_A=1,3
  do iSB_B=1,3
    if (allocated(AddSB(iSB_A,iSB_B)%A)) call mma_deallocate(AddSB(iSB_A,iSB_B)%A)
  end do
end do

!-----------------------------------------------------------------------
!IfTest = .False.
!-----------------------------------------------------------------------

return

end subroutine Cho_GenC
