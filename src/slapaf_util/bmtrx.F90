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

!#define _DEBUGPRINT_
subroutine BMtrx(nsAtom,Coor,nIter,mTtAtm,nWndw)

use Index_Functions, only: nTri_Elem
use Slapaf_Info, only: BMx, BSet, Cx, Curvilinear, HSet, iRef, KtB, Lbl, lOld, MaxItr, nDimBC, nLambda, Numerical, qInt, &
                       Redundant, Shift, Smmtrc, User_Def
use Slapaf_procedures, only: Box, Hidden
use UnixInfo, only: SuperName
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nsAtom, nIter, mTtAtm, nWndw
real(kind=wp), intent(in) :: Coor(3,nsAtom)
integer(kind=iwp) :: i, iAtom, ix, ixyz, Lu, mTR, mTtAtm_, nBonds, nHidden, nMax, nQQ
integer(kind=iwp), allocatable :: AN(:), TabA(:,:,:), TabAI(:,:), TabB(:,:)
real(kind=wp), allocatable :: Coor2(:,:), EVal(:), Hss_X(:), Scr2(:), TR(:), TRNew(:), TROld(:), Vec(:,:)

!                                                                      *
!***********************************************************************
!                                                                      *
nQQ = 0

Lu = u6
!                                                                      *
!***********************************************************************
!                                                                      *
! iRef: point at the reference geometry, used for non-redundant
!       internal coordinate to define the K-matrix and to generate
!       the raw model Hessian and TR vectors.

if (Numerical) then
  iRef = 1              ! Numerical Hessian Computation
else if (iRef == 0) then
  if (.not. BSet) then
    iRef = nIter-1      ! Compute cartesian Structure
  else
    iRef = nIter        ! Normal Computation
  end if
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' Actual structure from iteration',iRef
write(u6,*) ' Last structure from iteration',nIter
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Find the translational and rotational eigenvectors for the
! current structure.

call mma_allocate(TR,18*nsAtom,Label='TR')
TR(:) = Zero

call TRPGen(nDimBC,nsAtom,Cx(1,1,iRef),mTR,.false.,TR)

call mma_allocate(TRnew,3*nsAtom*mTR,Label='TRNew')
TRNew(:) = Zero
i = 0
do ix=1,3*nsAtom
  iAtom = (ix+2)/3
  ixyz = ix-(iAtom-1)*3
  if (Smmtrc(ixyz,iAtom)) then
    i = i+1
    call dcopy_(mTR,TR(i),-nDimBC,TRNew(ix),3*nsAtom)
  end if
end do
call Put_dArray('TR',TRnew,3*nsAtom*mTR)
call mma_deallocate(TRnew)

#ifdef _DEBUGPRINT_
call RecPrt('TR',' ',TR,nDimBC,mTR)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(TabAI,2,mTtAtm,Label='TabAI')
call mma_allocate(Vec,3*mTtAtm,nDimBC,Label='Vec')
call mma_allocate(AN,mTtAtm,Label='AN')
call mma_allocate(Coor2,3,mTtAtm,Label='Coor2')

! Generate Grand atoms list

call GenCoo(Cx(:,:,iRef),nsAtom,Coor2,mTtAtm,Vec,nDimBC,AN,TabAI)

! Are there some hidden frozen atoms ?

call Hidden(Coor2,AN,nHidden)

! Generate bond list

mTtAtm_ = mTtAtm+nHidden
call Box(Coor2,mTtAtm_,AN,TabB,TabA,nBonds,nMax)
!                                                                      *
!***********************************************************************
!                                                                      *
! First compute the approximate Hessian in cartesians
!
! OBSERVE if the analytic Hessian is available it will be used
! rather than the model Hessian.
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the raw Cartesian Hessian

call mma_allocate(EVal,nTri_Elem(3*mTtAtm),Label='EVal')
call mma_Allocate(Hss_X,(3*mTtAtm)**2,Label='Hss_X')
call mma_allocate(Scr2,(3*mTtAtm)**2,Label='Scr2')

if (HSet .or. (.not. (Curvilinear .or. User_Def))) call LNM(Coor2,mTtAtm,EVal,Hss_X,Scr2,Vec,nsAtom,nDimBC,AN,nIter,TabB,TabA, &
                                                            nBonds,nMax,nHidden)

call mma_deallocate(Scr2)
call mma_deallocate(Coor2)
!                                                                      *
!***********************************************************************
!                                                                      *
! The internal coordinates are of either two sets.
!
!          i)  user supplied internal coordinates
!
!         or
!
!         ii)  automatic internal coordinates
!                a) cartesian coordinates (lnm)
!                b) non-redundant internal coordinates (nrc)
!                   1) Conventional
!                   2) Curvature weighted   (HWRS)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*) 'User_Def   :',User_Def
write(u6,*) 'Curvilinear:',Curvilinear
#endif
if (User_Def) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call BMtrx_User_Defined(nsAtom,Coor,nDimBC,nIter,mTR,nQQ)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if (Curvilinear) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (Redundant) then
    call WarningMessage(2,' Bmtrx: Redundant option not implemented yet.')
    call Abend()
  end if

  ! Re-generate the bonds if there were hidden atoms

  if (nHidden /= 0) call Box(Coor2,mTtAtm,AN,TabB,TabA,nBonds,nMax)
  call BMtrx_Internal(nsAtom,nDimBC,nIter,mTtAtm,iRef,mTR,TR,TabAI,TabA,TabB,nBonds,nMax,iRef,nQQ,nWndw)

  ! Set the Labels for internal coordinates.

  do i=1,nQQ
    write(Lbl(i),'(A,I3.3,A)') 'nrc',i,'  '
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else   ! Cartesian coordinates
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call BMtrx_Cartesian(nsAtom,nDimBC,nIter,mTtAtm,mTR,TR,EVal,Hss_X,nQQ,nWndw)

  ! Set the Labels for cartesian normal modes.

  do i=1,nQQ
    write(Lbl(i),'(A,I3.3,A)') 'lnm',i,'  '
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
if (BSet .and. HSet .and. (.not. lOld)) then
  call Put_dArray('Hss_X',Hss_X,nDimBC**2)
  call Put_dArray('KtB',KtB,nDimBC*nQQ)
  call mma_deallocate(KtB)
end if
call mma_deallocate(Hss_X)
call mma_deallocate(EVal)
call mma_deallocate(TabA)
call mma_deallocate(TabB)
call mma_deallocate(AN)
call mma_deallocate(Vec)
call mma_deallocate(TabAI)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the shift vector in the basis.

if (Bset) then
  call mma_allocate(Shift,nQQ,MaxItr,Label='Shift')
  Shift(:,:) = Zero
  call ShfANM(nQQ,nIter,qInt,Shift)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Store B matrices to be used to transform the numerical
! Hessian to the new basis as the optimization proceeds.

if (((nIter == 1) .and. BSet) .and. (SuperName /= 'numerical_gradient')) then

  call Put_dArray('BMxOld',BMx,3*nsAtom*nQQ)

  if (mTR /= 0) then
    call mma_allocate(TROld,3*nsAtom*mTR,Label='TROld')
    TROld(:) = Zero
#   ifdef _DEBUGPRINT_
    call RecPrt('TRVec',' ',TR,3*nsAtom,mTR)
#   endif
    i = 0
    do ix=1,3*nsAtom
      iAtom = (ix+2)/3
      ixyz = ix-(iAtom-1)*3
      if (Smmtrc(ixyz,iAtom)) then
        i = i+1
        call dcopy_(mTR,TR(i),-nDimBC,TROld(ix),3*nsAtom)
      end if
    end do
    call Put_dArray('TROld',TROld,3*nsAtom*mTR)
    call mma_deallocate(TROld)
  end if
end if

! Print the B-matrix

#ifdef _DEBUGPRINT_
call RecPrt(' The BMtrx',' ',BMx,3*nsAtom,nQQ)
#endif
call mma_deallocate(TR)
!                                                                      *
!***********************************************************************
!                                                                      *
! Print out the values of the internal coordinates

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' Internal coordinates'
write(u6,*)
write(u6,'(1X,A,2X,F10.4)') (Lbl(iInter),qInt(iInter,nIter),iInter=1,nQQ)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Too many constraints?

if (nLambda > nQQ) then
  call WarningMessage(2,'Error in RlxCtl')
  write(Lu,*)
  write(Lu,*) '********************************************'
  write(Lu,*) ' ERROR: nLambda > nQQ'
  write(Lu,*) ' nLambda=',nLambda
  write(Lu,*) ' nQQ=',nQQ
  write(Lu,*) ' There are more constraints than coordinates'
  write(Lu,*) '********************************************'
  call Quit_OnUserError()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine BMtrx
