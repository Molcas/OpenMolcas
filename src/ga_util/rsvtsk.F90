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
!***********************************************************************
! Integer Function RsvTsk(igaTsk,iTskLs,nTsk,iStart)                   *
!  -> reserve a task for my node and mark it on the global task list   *
!     as reserved.                                                     *
!     igaTsk:   global array handle to global task list                *
!     iTskLs:   my private task list (my favourite sequence of tasks)  *
!     nTsk:     # of tasks                                             *
!     iStart:   starting value of index of private task list           *
!***********************************************************************

function RsvTsk(igaTsk,iTskLs,nTsk,mTsk,iStart,iS,iE)

#if defined (_MOLCAS_MPP_) && !defined (_GA_)
use stdalloc, only: mma_allocate, mma_deallocate
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: RsvTsk
integer(kind=iwp) :: igaTsk, nTsk, iTskLs(nTsk,2), mTsk, iStart, iS, iE
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: iCnt, iTsk
logical(kind=iwp) :: Reserved
#include "global.fh"

#ifdef _GA_
if (iStart > mTsk) then
  iTsk = 0
  iCnt = nTsk
else
  do iCnt=iStart,mTsk
    iTsk = iTskLs(iCnt,1)
    ! try to reserve iTsk on global task list...
    Reserved = ga_read_inc(igaTsk,iTsk,1,1) /= 0
    if (Reserved) then
      iE = iE-1
      iTskLs(iE,2) = iTsk
    else
      iS = iS+1
      iTskLs(iS,2) = iTsk
      Go To 100
    end if
  end do
  iTsk = 0
100 continue

end if
RsvTsk = iTsk
iStart = iCnt

#else
integer(kind=iwp), allocatable :: TSKR(:)

if (iStart > mTsk) then
  iTsk = 0
  iCnt = nTsk
else
  call mma_allocate(TSKR,mTsk,Label='TSKR')
  call ga_readb_inc(igaTsk,mTsk,iTskLs(iStart,1),TSKR)
  do iCnt=iStart,mTsk
    iTsk = iTskLs(iCnt,1)
    Reserved = TSKR(1+iCnt-iStart) /= 0
    if (Reserved) then
      iE = iE-1
      iTskLs(iE,2) = iTsk
    else
      iS = iS+1
      iTskLs(iS,2) = iTsk
      Go To 100
    end if
  end do
  iTsk = 0
100 continue
  call mma_deallocate(TSKR)
end if
RsvTsk = iTsk
iStart = iCnt

#endif
#else

#include "macros.fh"
unused_var(igaTsk)
unused_var(iTskLs)
unused_var(mTsk)
unused_var(iS)
unused_var(iE)

RsvTsk = 0
iStart = nTsk

#endif

end function RsvTsk
