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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

!> The pieces shared by the loops over coupled pairs of excitation cases
!> (SIGMA_CASPT2, SIGDER, TRDNS2O): the outer loop runs over the blocks
!> (iCase1,iSym1), the inner one over the blocks (iCase2,iSym2) with
!> iCase2 > iCase1 that couple to it (Cpl_Code in EQSOLV).
module Case_Pairs

use Definitions, only: wp, iwp

implicit none
private

public :: Cpl_Block, Cpl_D_Arrays, Cpl_GSum, Cpl_N1, Cpl_Read_Vec, Cpl_Spec1

contains

!> .true. if the block (iCase,iSym) of the first-order vector holds
!> anything (independent functions and a non-empty nAS x nIS array);
!> nAS, nIS (out): its dimensions.
function Cpl_Block(iCase,iSym,nAS,nIS) result(Active)

  use caspt2_module, only: nASup, nInDep, nISup

  logical(kind=iwp) :: Active
  integer(kind=iwp), intent(in) :: iCase, iSym
  integer(kind=iwp), intent(out) :: nAS, nIS

  nAS = 0
  nIS = 0
  Active = .false.
  if (nInDep(iSym,iCase) == 0) return
  nAS = nASup(iSym,iCase)
  nIS = nISup(iSym,iCase)
  Active = (nAS*nIS /= 0)

end function Cpl_Block

!> Size of the one-electron part (SPEC1A/C/D) of the block (iCase,iSym)
!> with nIS non-active superindices: cases A (1) and C (4), and D (5) in
!> the totally symmetric irrep; 0 otherwise.
function Cpl_N1(iCase,iSym,nIS) result(n1)

  use Orbital_Spaces, only: nAsh
  use caspt2_module, only: nIsh, nSsh

  integer(kind=iwp) :: n1
  integer(kind=iwp), intent(in) :: iCase, iSym, nIS

  n1 = 0
  if (iCase == 1) then
    n1 = nAsh(iSym)*nIsh(iSym)
  else if (iCase == 4) then
    n1 = nAsh(iSym)*nSsh(iSym)
  else if ((iCase == 5) .and. (iSym == 1)) then
    n1 = nIS
  end if

end function Cpl_N1

!> The one-electron part of block (iCase,iSym), the case dispatch of
!> SPEC1A/SPEC1C/SPEC1D with the factor 1/nActEl. iMltOp, X(nX), Y(nY):
!> as in those routines; nothing is done for the other blocks.
subroutine Cpl_Spec1(iMltOp,iCase,iSym,X,nX,Y,nY)

  use CI_Lib_Info, only: CI_Lib_Space
  use Constants, only: One

  integer(kind=iwp), intent(in) :: iMltOp, iCase, iSym, nX, nY
  real(kind=wp), intent(inout) :: X(nX), Y(nY)
  integer(kind=iwp), parameter :: iState = 1
  real(kind=wp) :: Fact

  Fact = One/real(max(1,CI_Lib_Space(iState)%nActEl),kind=wp)
  if (iCase == 1) then
    call SPEC1A(iMltOp,Fact,iSym,X,nX,Y,nY)
  else if (iCase == 4) then
    call SPEC1C(iMltOp,Fact,iSym,X,nX,Y,nY)
  else if ((iCase == 5) .and. (iSym == 1)) then
    call SPEC1D(iMltOp,Fact,X,nX,Y,nY)
  end if

end subroutine Cpl_Spec1

!> The operands of the transposed coupling loops (iMltOp = 1) for block
!> (iCase,iSym): D2 (nAS*nIS) from the RHS array lg_D2, which is freed, and
!> its one-electron part D1 (nD1 elements, Cpl_N1; a dummy of one element
!> if nD1 = 0). D1 and D2 are allocated here, the caller frees them.
subroutine Cpl_D_Arrays(iCase,iSym,nAS,nIS,lg_D2,D2,D1,nD1)

  use stdalloc, only: mma_allocate
  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: iCase, iSym, nAS, nIS
  integer(kind=iwp), intent(inout) :: lg_D2
  real(kind=wp), allocatable, intent(out) :: D2(:), D1(:)
  integer(kind=iwp), intent(out) :: nD1

  call mma_allocate(D2,nAS*nIS,Label='D2')
  call RHS_GET(nAS,nIS,lg_D2,D2)
  call RHS_FREE(lg_D2)

  nD1 = Cpl_N1(iCase,iSym,nIS)
  if (nD1 > 0) then
    call mma_allocate(D1,nD1,Label='D1')
    D1(:) = Zero
    call Cpl_Spec1(1,iCase,iSym,D2,size(D2),D1,size(D1))
  else
    call mma_allocate(D1,1,Label='D1')
  end if

end subroutine Cpl_D_Arrays

!> Read block (iCase,iSym) of vector iVec into the RHS array lg (nAS x
!> nIS, allocated by the caller); with Add_Lambda the gradient form
!> T + lambda: scaled by Scal, plus the same block of vector jVec.
subroutine Cpl_Read_Vec(nAS,nIS,lg,iCase,iSym,iVec,jVec,Add_Lambda,Scal)

  use Constants, only: One

  integer(kind=iwp), intent(in) :: nAS, nIS, lg, iCase, iSym, iVec, jVec
  logical(kind=iwp), intent(in) :: Add_Lambda
  real(kind=wp), intent(in) :: Scal
  integer(kind=iwp) :: lg_V

  call RHS_READ(nAS,nIS,lg,iCase,iSym,iVec)
  if (.not. Add_Lambda) return
  if (Scal /= One) call RHS_SCAL(nAS,nIS,lg,Scal)
  call RHS_ALLO(nAS,nIS,lg_V)
  call RHS_READ(nAS,nIS,lg_V,iCase,iSym,jVec)
  call RHS_DAXPY(nAS,nIS,One,lg_V,lg)
  call RHS_FREE(lg_V)

end subroutine Cpl_Read_Vec

!> Sum the replicated array A(n) over the processes, in messages of at
!> most 2**27 elements.
subroutine Cpl_GSum(A,n)

  integer(kind=iwp), intent(in) :: n
  real(kind=wp), intent(inout) :: A(n)
  integer(kind=iwp), parameter :: Max_Mesg_Size = 2**27
  integer(kind=iwp) :: iSta, nBlk

  do iSta=1,n,Max_Mesg_Size
    nBlk = min(Max_Mesg_Size,n-iSta+1)
    call GADGOP(A(iSta:),nBlk,'+')
  end do

end subroutine Cpl_GSum

end module Case_Pairs
