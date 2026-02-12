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

subroutine MKPROJ(CRVEC,CMO,TUVX)

use rasscf_global, only: CORESHIFT
use general_data, only: NCRVEC, NTOT2, NASH, NBAS
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
real*8 CRVEC(NCRVEC), CMO(NTOT2)
real*8 TUVX(*)
real*8, allocatable :: CS_TMP(:)
real*8 CS_TU, CS_TUV, COREPROJ
integer NA, NB, IPOS, ITU, IT, IU, IVX, IV, IXMX, IX
real*8, external :: DDot_

NA = NASH(1)
NB = NBAS(1)

call mma_allocate(CS_TMP,NB,Label='CS_TMP')
do IT=1,NA
  CS_TMP(IT) = DDOT_(NB,CMO((IT-1)*NB+1),1,CRVEC,1)
end do

IPOS = 0
ITU = 0
do IT=1,NA
  do IU=1,IT
    CS_TU = CS_TMP(IT)*CS_TMP(IU)
    ITU = ITU+1
    IVX = 0
    do IV=1,IT
      CS_TUV = CS_TU*CS_TMP(IV)
      IXMX = IV
      if (IT == IV) IXMX = IU
      do IX=1,IXMX
        IVX = IVX+1
        IPOS = IPOS+1
        COREPROJ = CS_TUV*CS_TMP(IX)
        TUVX(IPOS) = TUVX(IPOS)+CORESHIFT*COREPROJ
      end do
    end do
  end do
end do

call mma_deallocate(CS_TMP)

end subroutine MKPROJ
