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
integer(kind=iwp) :: I, IOFF, NBASE, NBLOCK, NDIM, NPROCS, NREST
logical(kind=iwp) :: BSTAT
integer(kind=iwp), allocatable :: MAP(:)
#include "global.fh"
#include "mafdecls.fh"

NPROCS = GA_NNODES()

select case (ORI)
  case ('H')
    NDIM = NROW
  case ('V')
    NDIM = NCOL
  case default
    write(u6,*) 'GA_CREATE_STRIPED: unknown orientation: ',ORI
    call AbEnd()
    NDIM = 0
    BSTAT = .false.
end select
NBLOCK = min(NDIM,NPROCS)
NBASE = NDIM/NPROCS
NREST = mod(NDIM,NPROCS)
call MMA_ALLOCATE(MAP,NBLOCK,Label='MAP')
IOFF = 1
do I=1,NBLOCK
  MAP(I) = IOFF
  if (I <= NREST) then
    IOFF = IOFF+NBASE+1
  else
    IOFF = IOFF+NBASE
  end if
end do

select case (ORI)
  case ('H')
    BSTAT = GA_CREATE_IRREG(MT_DBL,NROW,NCOL,LABEL,MAP,NBLOCK,[1],1,LG_M)
  case ('V')
    BSTAT = GA_CREATE_IRREG(MT_DBL,NROW,NCOL,LABEL,[1],1,MAP,NBLOCK,LG_M)
end select

call MMA_DEALLOCATE(MAP)

if (.not. bStat) then
  write(u6,*) 'GA_CREATE_STRIPED: could not create array, abort'
  call AbEnd()
end if

end subroutine GA_CREATE_STRIPED

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(GA_CREATE_STRIPED)

#endif
