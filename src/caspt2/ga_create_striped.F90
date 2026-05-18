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
!SVC: DGA can't distribute stripes for some small dimensions so we have
!to explicitly use the irregular versions. These is a wrapper routine to
!create horizontal (H) and vertical (V) stripes. The dimensions are
!divided evenly, with the remainder spread over the leading stripes.

#include "compiler_features.h"
#ifdef _MOLCAS_MPP_

subroutine GA_CREATE_STRIPED(ORI,NROW,NCOL,LABEL,LG_M)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
character, intent(in) :: ORI
integer(kind=iwp), intent(in) :: NROW, NCOL
character(len=*), intent(in) :: LABEL
integer(kind=iwp), intent(in) :: LG_M
integer(kind=iwp) :: I, IOFF, NBASE, NBLOCK1, NBLOCK2, NDIM, NPROCS, NREST
logical(kind=iwp) :: BSTAT
integer(kind=iwp), allocatable :: MAP1(:), MAP2(:)
#include "global.fh"
#include "mafdecls.fh"

NPROCS = GA_NNODES()

NBLOCK1 = 1
call MMA_ALLOCATE(MAP1,NBLOCK1,Label='MAP1')
MAP1(1) = 1

NDIM = 0
if (ORI == 'H') then
  NDIM = NROW
else if (ORI == 'V') then
  NDIM = NCOL
end if
NBLOCK2 = min(NDIM,NPROCS)
NBASE = NDIM/NPROCS
NREST = mod(NDIM,NPROCS)
call MMA_ALLOCATE(MAP2,NBLOCK2,Label='MAP2')
IOFF = 1
do I=1,NBLOCK2
  MAP2(I) = IOFF
  if (I <= NREST) then
    IOFF = IOFF+NBASE+1
  else
    IOFF = IOFF+NBASE
  end if
end do

BSTAT = .false.
if (ORI == 'H') then
  BSTAT = GA_CREATE_IRREG(MT_DBL,NROW,NCOL,LABEL,MAP2,NBLOCK2,MAP1,NBLOCK1,LG_M)
else if (ORI == 'V') then
  BSTAT = GA_CREATE_IRREG(MT_DBL,NROW,NCOL,LABEL,MAP1,NBLOCK1,MAP2,NBLOCK2,LG_M)
end if

call MMA_DEALLOCATE(MAP1)
call MMA_DEALLOCATE(MAP2)

if (.not. bStat) then
  write(u6,*) 'GA_CREATE_HS: could not create array, abort'
  call AbEnd()
end if

end subroutine GA_CREATE_STRIPED

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(GA_CREATE_STRIPED)

#endif
