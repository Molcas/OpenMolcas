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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

#include "compiler_features.h"
#ifdef _MOLCAS_MPP_

subroutine DF3_DEPSA_MPP(NG3,NASHT,DF3,DEPSA,lg_S,idxG3)

use Symmetry_Info, only: Mul
use SUPERINDEX, only: KTUV
use Para_Info, only: nProcs
use caspt2_module, only: IASYM, NTUVES
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, byte

implicit none
integer(kind=iwp), intent(in) :: NG3, NASHT, lg_S
real(kind=wp), intent(in) :: DF3(NG3)
real(kind=wp), intent(inout) :: DEPSA(NASHT,NASHT)
integer(kind=byte), intent(in) :: idxG3(6,NG3)
integer(kind=iwp) :: iG3, IHI, ILO, iRank, iST, iSU, ISUP1, ISUP2, iSV, iSX, iSY, isym, iSZ, iT, ituvs, iU, iV, IW, iX, ixyzs, iY, &
                     iZ, JHI, JLO, JSUP1, JSUP2, NCOL, NROW, NSEQ
real(kind=wp) :: F3VAL
real(kind=wp), allocatable :: WRK(:)

!! do depsa
isym = 1
do iRank=0,NPROCS-1
  call GA_Distribution(lg_S,iRank,ILO,IHI,JLO,JHI)
  NROW = iHi-iLo+1
  NCOL = jHi-jLo+1
  call mma_allocate(WRK,NROW*NCOL,Label='WRK')
  call GA_Get(lg_S,iLo,iHi,jLo,jHi,WRK,NROW)
  do iG3=1,NG3
    iT = idxG3(1,iG3)
    iU = idxG3(2,iG3)
    iV = idxG3(3,iG3)
    iX = idxG3(4,iG3)
    iY = idxG3(5,iG3)
    iZ = idxG3(6,iG3)
    iST = IASYM(iT)
    iSU = IASYM(iU)
    iSV = IASYM(iV)
    iSX = IASYM(iX)
    iSY = IASYM(iY)
    iSZ = IASYM(iZ)
    ituvs = Mul(IST,Mul(ISU,ISV))
    ixyzs = Mul(ISX,Mul(ISY,ISZ))
    if (ituvs /= ixyzs) cycle

    F3VAL = DF3(iG3)

    ISUP1 = KTUV(iV,iU,iT)-nTUVES(iSYM)
    if ((ISUP1 >= iLo) .and. (ISUP1 <= iHi)) then
      do iW=1,nAshT
        JSUP1 = KTUV(iX,iW,iZ)-nTUVES(iSYM)
        NSEQ = 1+ISUP1-iLo+NROW*(JSUP1-jLo)
        DEPSA(iW,iY) = DEPSA(iW,iY)-F3VAL*WRK(NSEQ)
      end do
    end if
    JSUP2 = KTUV(iX,iY,iZ)-nTUVES(iSYM)
    if ((JSUP2 >= iLo) .and. (JSUP2 <= iHi)) then
      do iW=1,nAshT
        ISUP2 = KTUV(iV,iW,iT)-nTUVES(iSYM)
        NSEQ = 1+JSUP2-iLo+NROW*(ISUP2-jLo)
        DEPSA(iW,iU) = DEPSA(iW,iU)-F3VAL*WRK(NSEQ)
      end do
    end if
  end do
  call mma_deallocate(WRK)
end do

return

end subroutine DF3_DEPSA_MPP

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(DF3_DEPSA_MPP)

#endif
