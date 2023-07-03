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

subroutine Get_dDipM(dDipM,DipM,nDoF,nInter)
!***********************************************************************
!                                                                      *
!     Objective: to compute the dipole moment derivative in Cartesians *
!                                                                      *
!***********************************************************************

use Slapaf_Info, only: Coor
use Slapaf_Parameters, only: mTROld

implicit real*8(a-h,o-z)
#include "stdalloc.fh"
real*8 dDipM(3,nDoF), DipM(3)
logical Found
real*8, allocatable :: Tmp2(:), BOld(:), TROld(:)

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
nX = 3*size(Coor,2)

call mma_allocate(Tmp2,nX**2,Label='Tmp2')
call mma_allocate(BOld,nX*nInter,Label='BOld')
call Qpg_dArray('BMxOld',Found,nBMx)
if (Found .and. (nBMx == nX*nInter)) then
  call Get_dArray('BMxOld',BOld,nX*nInter)
else
  call Get_dArray('BMtrx',BOld,nX*nInter)
end if
if (mTROld > 0) then
  call mma_allocate(TROld,nX*mTROld,Label='TROld')
  call Qpg_dArray('TROld',Found,nTR)
  if (Found .and. (nTR == nX*mTROld)) then
    call Get_dArray('TROld',TROld,nX*mTROld)
  else
    call Get_dArray('TR',TROld,nX*mTROld)
  end if
else
  call mma_allocate(TROld,1,Label='TROld')
end if

call Get_dDipM_(nX,BOld,TROld,nDoF,nInter,Tmp2,dDipM,mTROld,size(Coor,2),DipM)

call mma_deallocate(TROld)
call mma_deallocate(BOld)
call mma_deallocate(Tmp2)

return

end subroutine Get_dDipM

subroutine Get_dDipM_(nX,BMtrx,TRVec,nDoF,nInter,Tmp2,dDipM,mTR,nAtom,DipM)

use Slapaf_Info, only: Cx, Degen, Smmtrc

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 TRVec(nX,mTR), BMtrx(nX,nInter), Tmp2(nX**2), dDipM(3,nInter+mTR), DipM(3)
real*8 CM(3)
parameter(thr=1.0D-12)

#ifdef _DEBUGPRINT_
call RecPrt('TRVec',' ',TRVec,nX,mTR)
call RecPrt('BMtrx',' ',BMtrx,nX,nInter)
call RecPrt('dDipM',' ',dDipM,3,nInter+mTR)
call RecPrt('DipM',' ',DipM,3,1)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Analysis of the translational and rotational modes.
! Compute the dipole moment derivative with respect to the
! translational and rotational modes.

do i=1,3
  CM(i) = Zero
  rNorm = Zero
  do iAtom=1,nAtom
    rNorm = rNorm+Degen(i,iAtom)
    if (Smmtrc(i,iAtom)) then
      CM(i) = CM(i)+Degen(i,iAtom)*Cx(i,iAtom,1)
    end if
  end do
  CM(i) = CM(i)/rNorm
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over the rotations and translations

iTR = nInter+mTR
do iX=mTR,1,-1
  Tx = Zero
  Ty = Zero
  Tz = Zero
  Rx = Zero
  Ry = Zero
  Rz = Zero
  do i=1,nAtom
    Tx = Tx+TRVec((i-1)*3+1,iX)*Degen(1,i)
    Ty = Ty+TRVec((i-1)*3+2,iX)*Degen(2,i)
    Tz = Tz+TRVec((i-1)*3+3,iX)*Degen(3,i)
    Rx = Rx+(TRVec((i-1)*3+2,iX)*(Cx(3,i,1)-CM(3))-TRVec((i-1)*3+3,iX)*(Cx(2,i,1)-CM(2)))*Degen(1,i)
    Ry = Ry+(TRVec((i-1)*3+3,iX)*(Cx(1,i,1)-CM(1))-TRVec((i-1)*3+1,iX)*(Cx(3,i,1)-CM(3)))*Degen(2,i)
    Rz = Rz+(TRVec((i-1)*3+1,iX)*(Cx(2,i,1)-CM(2))-TRVec((i-1)*3+2,iX)*(Cx(1,i,1)-CM(1)))*Degen(3,i)
  end do
#ifdef _DEBUGPRINT_
  write(6,*) 'Tx,Ty,Tz=',Tx,Ty,Tz
  write(6,*) 'Rx,Ry,Rz=',Rx,Ry,Rz
#endif

  if ((Rx**2+Ry**2+Rz**2 < thr) .and. (Tx**2+Ty**2+Tz**2 > thr)) then

  ! Translation, dipole moment invariant to translation

#   ifdef _DEBUGPRINT_
    write(6,*) 'Translation'
#   endif
    call dcopy_(3,[Zero],0,dDipM(1,iTR),1)

  else if ((Tx**2+Ty**2+Tz**2 < thr) .and. (Rx**2+Ry**2+Rz**2 > thr)) then
    rNorm = (Rx**2+Ry**2+Rz**2)

    ! Rotation, dipole moment variant to rotation

    if (rNorm > thr) then

      ! General axis

      dDipM(1,iTR) = (DipM(2)*Rz-DipM(3)*Ry)/rNorm
      dDipM(2,iTR) = (DipM(3)*Rx-DipM(1)*Rz)/rNorm
      dDipM(3,iTR) = (DipM(1)*Ry-DipM(2)*Rx)/rNorm
    else
      call WarningMessage(2,' GF: too small rNorm!')
      call Abend()

    end if
  end if
  iTR = iTR-1
end do
#ifdef _DEBUGPRINT_
call RecPrt('dDipM(Original)',' ',dDipM,3,nInter+mTR)
#endif
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Now backtransform to cartesian coordinates.
!
! dmu/dx = dq/dx  dmu/dq

do ix=1,3

  jj = 0
  jjj = 0
  do jAtom=1,nAtom
    do jx=1,3
      jjj = jjj+1
      if (Smmtrc(jx,jAtom)) then
        jj = jj+1
        ij = (jj-1)*3+ix

        tmp_ij = 0.0d0
        do k=1,nInter
          tmp_ij = tmp_ij+dDipM(ix,k)*BMtrx(jjj,k)
        end do
        do k=1,mTR
          tmp_ij = tmp_ij+dDipM(ix,nInter+k)*TRVec(jjj,k)
        end do
        Tmp2(ij) = tmp_ij

      end if
    end do
  end do

end do
call dcopy_(3*nDoF,Tmp2,1,dDipM,1)
#ifdef _DEBUGPRINT_
call RecPrt('dDipM(cartesian)',' ',dDipM,3,nDoF)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Get_dDipM_
