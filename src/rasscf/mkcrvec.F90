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
! Copyright (C) 2019, Per Ake Malmqvist                                *
!***********************************************************************

subroutine MKCRVEC(CMO_0,CRVEC)

use OneDat, only: sNoNuc, sNoOri
use rasscf_global, only: ITCORE
use general_data, only: NTOT, NTOT2, NTOT1, NBAS, NFRO, NISH
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: u6

implicit none
real*8 CRVEC(NTOT), CMO_0(NTOT2)
character(len=8) LABEL
real*8, allocatable :: STRI(:), SAO(:,:)
integer IRC, IOPT, ICOMP, ISYMLBL, NB, NFI
#include "warnings.h"

! Note: Nbas,etc are in included common. So is ITCORE.
! CMO_0 is the starting CMO vectors. Active orbital nr ITCORE in
! symmetry 1 is a specially prepared core orbital, which will be used
! to define a projector on any CI vector. This projector is that part
! of the state vector where the core orbital would be doubly occupied.
! This orbital is computed as a covariant vector CRVEC, to allow the
! projector to be invariant to the orbital basis in each interation.

call mma_allocate(STRI,NTOT1+4,Label='STRI')
IRC = 0
IOPT = ibset(ibset(0,sNoOri),sNoNuc)
LABEL = 'Mltpl  0'
ICOMP = 1
ISYMLBL = 1
call RdOne(IRC,IOPT,LABEL,ICOMP,STRI,ISYMLBL)
if (iRc /= 0) then
  write(u6,*) ' MKCRVEC could not read overlaps from ONEINT.'
  write(u6,*) ' Something is wrong with that file, or possibly'
  write(u6,*) ' with the program. Please check.'
  call Quit(_RC_IO_ERROR_READ_)
end if
NB = NBAS(1)
NFI = NFRO(1)+NISH(1)
call mma_allocate(SAO,NB,NB,Label='SAO')
call Square(STRI,SAO,1,NB,NB)
call mma_deallocate(STRI)
call DGEMV_('N',NB,NB,One,SAO,NB,CMO_0(NB*(NFI+ITCORE-1)+1),1,Zero,CRVEC,1)
call mma_deallocate(SAO)

! Test:
!write(u6,*) 'MKCRVEC test: Overlaps all orbs/core :'
!do it=1,nb
!  write(u6,'(1x,i5,f16.8)') it,ddot_(NB,CMO_0(NB*(IT-1)+1),1,CRVEC,1)
!end do

end subroutine MKCRVEC
