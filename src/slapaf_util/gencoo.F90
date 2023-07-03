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

subroutine GenCoo(Cart,nsAtom,Coor,mTtAtm,Vctrs,nDim,jAnr,iTabAI)

use Symmetry_Info, only: nIrrep, iOper
use Slapaf_Info, only: Degen, Smmtrc, ANr

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 Cart(3,nsAtom), Coor(3,mTtAtm), Vctrs(3*mTtAtm,nDim), r(3)
integer jAnr(mTtAtm), iTabAI(2,mTtAtm)
logical New

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt('GenCoo: Cart',' ',Cart,3,nsAtom)
call RecPrt('GenCoo: Degen',' ',Degen,3,nsAtom)
#endif

! Loop over list of symmetry unique centers

iSt = 1
iDim = 0
do iAtom=1,nsAtom
  Fact = One/sqrt(Degen(1,iAtom))
  iEnd = iSt
  jDim = iDim
  call dcopy_(3,Cart(1,iAtom),1,Coor(1,iSt),1)
  iTabAI(1,iEnd) = iAtom
  iTabAI(2,iEnd) = iOper(0)
  jAnr(iEnd) = Anr(iAtom)
  do ix=1,3
    if (Smmtrc(ix,iAtom)) then
      jDim = jDim+1
      call dcopy_(3*mTtAtm,[Zero],0,Vctrs(1,jDim),1)
      Vctrs((iEnd-1)*3+ix,jDim) = Fact
    end if
  end do

  ! Loop over the operators of the point group

  iElem = 1
  do ig=1,nIrrep-1
    r(1) = One
    if (iand(iOper(ig),1) /= 0) r(1) = -One
    r(2) = One
    if (iand(iOper(ig),2) /= 0) r(2) = -One
    r(3) = One
    if (iand(iOper(ig),4) /= 0) r(3) = -One
    x = r(1)*Cart(1,iAtom)
    y = r(2)*Cart(2,iAtom)
    z = r(3)*Cart(3,iAtom)

    New = .true.
    do iGo=iSt,iEnd
      if (New .and. (x == Coor(1,iGo)) .and. (y == Coor(2,iGo)) .and. (z == Coor(3,iGo))) New = .false.
    end do
    if (New) then
      iElem = iElem+1
      iEnd = iEnd+1
      Coor(1,iEnd) = x
      Coor(2,iEnd) = y
      Coor(3,iEnd) = z
      iTabAI(1,iEnd) = iAtom
      iTabAI(2,iEnd) = iOper(ig)
      jAnr(iEnd) = Anr(iAtom)
      jDim = iDim
      do ix=1,3
        if (Smmtrc(ix,iAtom)) then
          jDim = jDim+1
          Vctrs((iEnd-1)*3+ix,jDim) = r(ix)*Fact
        end if
      end do
    end if
  end do      ! End loop over operators

  do ix=1,3
    if (Smmtrc(ix,iAtom)) then
      iDim = iDim+1
    end if
  end do
  iSt = iEnd+1
end do         ! End loop over centers

#ifdef _DEBUGPRINT_
call RecPrt(' In GenCoo: Coor',' ',Coor,3,mTtAtm)
call RecPrt(' In GenCoo: Vctrs',' ',Vctrs,3*mTtAtm,nDim)
write(6,*)
write(6,*) ' iTabAI'
write(6,*)
do iAtom=1,mTtAtm
  write(6,*) iTabAI(1,iAtom),iTabAI(2,iAtom)
end do
#endif

return

end subroutine GenCoo
