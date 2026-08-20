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
! Copyright (C) 2026, Yoshio Nishimoto                                 *
!***********************************************************************
!
! Expansion of Cholesky vectors from the reduced set into square symmetry blocks.
!
! Usage:
!   call CHO_SQ_SETUP(USE_SQ)            ! decides NSUB and allocates LSQ
!   ... call CHO_RED2SQ(...) per sub-batch, read LSQ/IPLSQ ...
!   call CHO_SQ_CLOSE()

module Cho_Square

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: NSUB = 0        ! vectors expanded to square form at a time
integer(kind=iwp) :: NSQ1 = 0        ! size of the square blocks of a single vector
integer(kind=iwp) :: IPLSQ(8) = 0    ! offset of each symmetry block within LSQ
integer(kind=iwp) :: JSYMPRV = -1    ! vector symmetry of the previous expansion
integer(kind=iwp) :: JREDPRV = -1    ! reduced set of the previous expansion
integer(kind=iwp) :: NJPRV = -1      ! number of vectors in the previous expansion
real(kind=wp), allocatable :: LSQ(:) ! the expanded vectors, LSQ(a,j,b)

public :: Cho_Sq_Setup, Cho_Sq_Close
public :: Cho_Red2Sq
public :: LSQ, IPLSQ, NSUB

contains

!-----------------------------------------------------------------------

subroutine Cho_Sq_Setup(USE_SQ)

  use Symmetry_Info, only: Mul
  use ChoCASPT2, only: MXCHARR, NCHSPC
  use caspt2_module, only: nBas, nSym
  use stdalloc, only: mma_allocate, mma_MaxDBLE

  logical(kind=iwp), intent(inout) :: USE_SQ

  integer(kind=iwp) :: ISYMA, ISYMB, JSYM, MAXSUB, MXAVAIL, NPB
  integer(kind=iwp), parameter :: MAXSUB_DEF = 32 ! cap on NSUB

  MAXSUB = MAXSUB_DEF

  ! A fresh set of blocks: the first expansion has to zero them.
  NSUB = 0
  JSYMPRV = -1
  JREDPRV = -1
  NJPRV = -1
  if (.not. USE_SQ) return

  ! Largest square-block set over all vector symmetries, for one vector.
  NSQ1 = 0
  do JSYM=1,NSYM
    NPB = 0
    do ISYMA=1,NSYM
      ISYMB = Mul(ISYMA,JSYM)
      NPB = NPB+NBAS(ISYMA)*NBAS(ISYMB)
    end do
    NSQ1 = max(NSQ1,NPB)
  end do

  ! Spend at most half of what is left, and never more vectors than a batch holds anyway
  call mma_MaxDBLE(MXAVAIL)
  NSUB = min(NCHSPC/MXCHARR,max(1,MXAVAIL/(2*NSQ1)),MAXSUB)
  call mma_allocate(LSQ,NSUB*NSQ1,LABEL='LSQ')

end subroutine Cho_Sq_Setup

!-----------------------------------------------------------------------

subroutine Cho_Sq_Close()

  use stdalloc, only: mma_deallocate

  if (allocated(LSQ)) call mma_deallocate(LSQ)
  NSUB = 0

end subroutine Cho_Sq_Close

!-----------------------------------------------------------------------

subroutine Cho_Red2Sq(SCR,LSCR,JSUB,NJ,JSYM,JREDC,JVGLB)

  use Symmetry_Info, only: Mul
  use Cholesky, only: nBas, nDimRS, nSym
  use Constants, only: Zero

  integer(kind=iwp), intent(in) :: LSCR, JSUB, NJ, JSYM, JREDC, JVGLB
  real(kind=wp), intent(in) :: SCR(LSCR)

  integer(kind=iwp) :: IOFF, IREDCL, IRC, ISKIP(8), ISYMA, ISYMB, NRS, NTOT
  logical(kind=iwp) :: do_clear

  ! Expand NJ Cholesky vectors from the reduced set storage into square symmetry blocks,
  ! laid out as LSQ(a,j,b) with the vector index in the middle.

  IOFF = 1
  do ISYMA=1,NSYM
    ISYMB = Mul(ISYMA,JSYM)
    IPLSQ(ISYMA) = IOFF
    IOFF = IOFF+NBAS(ISYMA)*NJ*NBAS(ISYMB)
  end do
  NTOT = IOFF-1

  do_clear = (JSYM /= JSYMPRV) .or. (JREDC /= JREDPRV) .or. (NJ /= NJPRV)
  JSYMPRV = JSYM
  JREDPRV = JREDC
  NJPRV = NJ
  if (do_clear) LSQ(1:NTOT) = Zero

  ! CHO_REORDR reads its input from the start, so hand it the sub-batch.
  NRS = nDimRS(JSYM,JREDC)
  ISKIP(:) = 1
  IREDCL = JREDC
  call CHO_REORDR(IRC,SCR(1+(JSUB-1)*NRS),NJ*NRS,1,JVGLB,NJ,NJ,JSYM,IREDCL,3,IPLSQ,LSQ,ISKIP)

end subroutine Cho_Red2Sq

end module Cho_Square
