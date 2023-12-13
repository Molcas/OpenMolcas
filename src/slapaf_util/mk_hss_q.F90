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

subroutine Mk_Hss_Q()

use Slapaf_Info, only: BMx, BSet, Coor, Cx, Delta, DipM, HSet, iter, lNmHss, mTROld, nDimBC, NmIter, dqInt, mRowH, qInt
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iNeg, mInt, nsAtom

! Compute the Hessian in internal coordinates.

!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('Mk_Hss_Q: DipM',' ',DipM,size(DipM,1),size(DipM,2))
#endif
if ((lNmHss .or. allocated(mRowH)) .and. (iter == NmIter)) then
  mInt = nDimBC-mTROld
  nsAtom = size(Coor,2)
  call Put_dArray('Unique Coordinates',Cx,3*nsAtom)
  call Put_Coord_New(Cx,nsAtom)
  if (allocated(mRowH)) then
    if (BSet .and. HSet) call Hss_Q()
    call RowHessian(NmIter,mInt,Delta/2.5_wp)
  else
    call FormNumHess(iter,mInt,Delta,nsAtom,iNeg,DipM)
  end if

  Coor(:,:) = Cx(:,:,1)
  call Get_dArray('BMxOld',BMx,size(Coor)*size(qInt,1))
  qInt(:,Iter) = qInt(:,1)
  dqInt(:,Iter) = dqInt(:,1)
else if (BSet .and. HSet) then
  call Hss_Q()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Mk_Hss_Q
