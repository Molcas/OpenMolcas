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

subroutine LNM(Cart,mTtAtm,Hess,Scrt1,Scrt2,Vctrs,nsAtom,nDim,iAnr,nIter,iTabBonds,iTabAtoms,nBonds,nMax,nHidden)

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: nIrrep
use Slapaf_Info, only: Analytic_Hessian, Degen, iOptH, Smmtrc
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: mTtAtm, nsAtom, nDim, nHidden, iANr(mTtAtm+nHidden), nIter, nBonds, iTabBonds(3,nBonds), nMax, &
                                 iTabAtoms(2,0:nMax,mTtAtm+nHidden)
real(kind=wp), intent(in) :: Cart(3,mTtAtm+nHidden), Vctrs(3*mTtAtm,nDim)
real(kind=wp), intent(out) :: Hess(nTri_Elem(3*mTtAtm)), Scrt1((3*mTtAtm)**2), Scrt2((3*mTtAtm)**2)
integer(kind=iwp) :: i, iAtom, ii, ij, ijTri, IterHess, ix, ixyz, j, jAtom, ji, jj, jxyz, Len3, Length, nRP
real(kind=wp) :: eigen, Tmp
logical(kind=iwp) :: Found, RunOld
real(kind=wp), allocatable :: HTanVec(:), TanVec(:)

!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('In LNM: Cart',' ',Cart,3,mTtAtm)
if (nHidden /= 0) call RecPrt('In LNM: Cart(hidden atoms)',' ',Cart(1,mTtAtm+1),3,nHidden)
call RecPrt('In LNM: Vctrs',' ',Vctrs,3*mTtAtm,nDim)
write(u6,*) 'iAnr=',iANr
write(u6,*) 'Analytic Hessian=',Analytic_Hessian
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! CARTESIAN HESSIAN is generated approximately here or the analytic    *
! is read from file.                                                   *
!                                                                      *
! Retrieve the analytic Hessian or compute the Hessian in cartesians
! from the Hessian model function.

!                                                                      *
!***********************************************************************
!                                                                      *
if (Analytic_Hessian) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! The Hessian matrix is read in the basis of the symmetric displacements.

  Len3 = ndim
  Len3 = nTri_Elem(Len3)

  call qpg_dArray('Analytic Hessian',Found,Length)
  if (Found) then
    call Get_dArray_chk('Analytic Hessian',Hess,Len3)
    RunOld = .false.
  else
    call NameRun('RUNOLD')
    call qpg_dArray('Analytic Hessian',Found,Length)
    if (.not. Found) then
      call WarningMessage(2,' Error in LNM: Analytic Hessian not found')
      call Abend()
    end if
    call Get_dArray_chk('Analytic Hessian',Hess,Len3)
    call NameRun('#Pop')
    RunOld = .true.
  end if
# ifdef _DEBUGPRINT_
  call TriPrt('LNM: Analytic Hessian',' ',Hess,nDim)
  call RecPrt('LNM: Degen',' ',Degen,size(Degen,1),size(Degen,2))
# endif

  ! Modify matrix with degeneracy factors and square the matrix.

  Scrt1(1:nDim**2) = Zero
  ii = 0
  do i=1,3*nsAtom
    iAtom = (i+2)/3
    ixyz = i-(iAtom-1)*3
    if (Smmtrc(ixyz,iAtom)) then
      ii = ii+1
      jj = 0
      do j=1,i
        jAtom = (j+2)/3
        jxyz = j-(jAtom-1)*3
        if (Smmtrc(jxyz,jAtom)) then
          jj = jj+1
          ijTri = nTri_Elem(ii-1)+jj
          ij = (jj-1)*ndim+ii
          ji = (ii-1)*ndim+jj
          Tmp = Hess(ijTri)*sqrt(degen(ixyz,iAtom)*degen(jxyz,jAtom))
          Scrt1(ij) = Tmp
          Scrt1(ji) = Tmp
        end if
      end do
    end if
  end do
# ifdef _DEBUGPRINT_
  call TriPrt('Hessian(anal.)',' ',Hess,nDim)
  write(u6,*) 'nDim,mTtAtm,nsAtom=',nDim,mTtAtm,nsAtom
  call RecPrt('Hessian(anal.)',' ',Scrt1,nDim,nDim)
# endif

  ! If the analytic Hessian corresponds to the current iteration,
  ! disable Hessian updating

  if (RunOld) call NameRun('RUNOLD')
  call Get_iScalar('HessIter',IterHess)
  if (RunOld) call NameRun('#Pop')
  if (IterHess == nIter) iOptH = ibset(0,3)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else   ! Use the Hessian Model Function
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call ddV(Cart,mTtAtm,Hess,iANr,iTabBonds,iTabAtoms,nBonds,nMax,nHidden)
# ifdef _DEBUGPRINT_
  call TriPrt(' The Model Hessian','(12f9.5)',Hess,3*mTtAtm)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Square the matrix

  do i=1,3*mTtAtm
    do j=1,i
      ijTri = nTri_Elem(i-1)+j
      ij = (j-1)*(3*mTtAtm)+i
      ji = (i-1)*(3*mTtAtm)+j
      Scrt1(ij) = Hess(ijTri)
      Scrt1(ji) = Hess(ijTri)
    end do
  end do
  !call RecPrt('Scrt1',' ',Scrt1,3*mTtAtm,3*mTtAtm)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
# ifdef _DEBUGPRINT_
  call RecPrt(' Scrt1',' ',Scrt1,3*mTtAtm,3*mTtAtm)
# endif
  if (nIrrep /= 1) then

    ! Now project out the total symmetric part of the Hessian

    call DGEMM_('N','N',3*mTtAtm,nDim,3*mTtAtm,One,Scrt1,3*mTtAtm,Vctrs,3*mTtAtm,Zero,Scrt2,3*mTtAtm)
#   ifdef _DEBUGPRINT_
    call RecPrt(' Scrt2',' ',Scrt2,3*mTtAtm,nDim)
#   endif
    call DGEMM_('T','N',nDim,nDim,3*mTtAtm,One,Vctrs,3*mTtAtm,Scrt2,3*mTtAtm,Zero,Scrt1,nDim)
#   ifdef _DEBUGPRINT_
    call RecPrt(' The Symmetrized Hessian',' ',Scrt1,nDim,nDim)
#   endif
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Now project hessian using reaction path vector if TS search

  call qpg_darray('TanVec',Found,nRP)
  if (Found) then
    if (nRP /= 3*nsAtom) then
      call WarningMessage(2,' Error in LNM: nRP /= 3*nsAtom')
      write(u6,*) 'nRP,3*nsAtom=',nRP,nsAtom
      call Abend()
    end if

    call mma_allocate(TanVec,nRP,Label='TanVec')
    call mma_allocate(HTanVec,nRP**2,Label='HTanVec')
    TanVec(:) = Zero
    HTanVec(:) = Zero

    call Get_dArray('TanVec',TanVec,nRP)

    !call RecPrt('TanVec',' ',TanVec,nRP,1)
    i = 0
    do ix=1,nRP
      iAtom = (ix+2)/3
      ixyz = ix-(iAtom-1)*3
      if (Smmtrc(ixyz,iAtom)) then
        i = i+1
        TanVec(i) = TanVec(ix)
      end if
    end do
    nRP = nDim
    !call RecPrt('TanVec',' ',TanVec,nRP,1)

    ! Compute H|i>

    call dGeMV_('N',nRP,nRP,One,scrt1,nRP,TanVec,1,Zero,HTanVec,1)

    ! Compute <i|H|i>

    eigen = Zero
    do i=1,nRP
      eigen = eigen+HTanVec(i)*TanVec(i)
    end do
    !write(u6,*) 'Eigen=',eigen
    if (eigen >= Zero) then

      ! Form H - H|i><i| - |i><i|H

      do i=1,nRP
        scrt1((i-1)*nRP+1:i*nRP) = scrt1((i-1)*nRP+1:i*nRP)-HTanVec(i)*TanVec(:)-TanVec(i)*HTanVec(:)
      end do
    end if

    call mma_deallocate(HTanVec)
    call mma_deallocate(TanVec)
  end if
  !call RecPrt('Scrt1(final)',' ',Scrt1,nDim,nDim)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
do j=1,nDim
  do i=1,j
    ij = (j-1)*nDim+i
    ijTri = nTri_Elem(i-1)+j
    Hess(ijTri) = Scrt1(ij)
  end do
end do
!call RecPrt('Scrt1(final)',' ',Scrt1,nDim,nDim)
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine LNM
