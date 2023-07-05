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

subroutine BMtrx_Cartesian(nsAtom,nDimBC,nIter,mTtAtm,mTR,TRVec,EVal,Hss_x,nQQ,nWndw)

use Slapaf_Info, only: AtomLbl, BMx, Cx, Degen, dqInt, dqInt_Aux, Gx, Gx0, KtB, NAC, qInt, Smmtrc
use Slapaf_Parameters, only: BSet, HSet, lOld, MaxItr, PrQ, Redundant
use Kriging_Mod, only: nSet
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nsAtom, nDimBC, nIter, mTtAtm, mTR, nQQ, nWndw
real(kind=wp) :: TRVec(nDimBC,mTR), Eval(3*mTtAtm*(3*mTtAtm+1)/2), Hss_x((3*mTtAtm)**2)
integer(kind=iwp) :: i, i_Dim, iAtom, iInd, iInter, ij, ijTri, ik, ipFrom, iTR, ix, ixyz, j, jAtom, ji, jx, jxyz, k, kAtom, kx, kxyz
real(kind=wp) :: Hii, Omega, Temp
integer(kind=iwp), allocatable :: Ind(:)
real(kind=wp), allocatable :: Degen2(:), EVec(:), Hi(:,:), iHi(:)
real(kind=wp), external :: DDot_

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('BMtrx_Cartesian: TRVec',' ',TRVec,nDimBC,mTR)
call RecPrt('BMtrx_Cartesian: Degen',' ',Degen,3,size(Degen,2))
call RecPrt('BMtrx_Cartesian: Hss_x',' ',Hss_x,3*mTtAtm,3*mTtAtm)
#endif

! Recompute the B matrix once each macroIteration, this is
! not done if a numerical Hessian is computed.
!                                                                      *
!***********************************************************************
!                                                                      *
!     R E D U N D A N T  C A R T E S I A N  C O O R D S

if (Redundant) then
  nQQ = nDimBC
  if (allocated(qInt)) then
    if (size(qInt,1) /= nQQ) then
      call mma_deallocate(qInt)
      call mma_deallocate(dqInt)
    end if
  end if
  if (.not. allocated(qInt)) then
    call mma_allocate(qInt,nQQ,MaxItr,Label='qInt')
    call mma_allocate(dqInt,nQQ,MaxItr,Label='dqInt')
    qInt(:,:) = Zero
    dqInt(:,:) = Zero
  end if
  if (allocated(dqInt_Aux)) then
    if (size(dqInt_Aux,1) /= nQQ) call mma_deallocate(dqInt_Aux)
  end if
  if ((.not. allocated(dqInt_Aux)) .and. (nSet > 1)) then
    call mma_allocate(dqInt_Aux,nQQ,MaxItr,nSet-1,Label='dqInt_Aux')
    dqInt_Aux(:,:,:) = Zero
  end if
  call mma_allocate(EVec,nDimBC**2,Label='EVec')
  EVec(:) = Zero
  call dcopy_(nDimBC,[One],0,EVec,nDimBC+1)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Move over the eigenvectors putting to BMx

  call mma_allocate(BMx,3*nsAtom,nQQ,Label='BMx')
  BMX(:,:) = Zero
  ipFrom = 1
  call BPut(EVec(ipFrom),nDimBC,BMx,3*nsAtom,Smmtrc,nQQ,Degen)
# ifdef _DEBUGPRINT_
  call RecPrt('In Bmtrx: B',' ',BMx,3*nsAtom,nQQ)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (PrQ .and. (nsAtom <= 5)) call List2('Cartesian Redundant',AtomLbl,BMx,nsAtom,nQQ,Smmtrc)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Project the model Hessian with respect to rotations and
  ! translations. The eigenvalues are shifted to large positive
  ! eigenvalues to effectively remove any displacements in the
  ! rotational and translations directions and to make sure that
  ! the matrix is not singular.

  call mma_allocate(Ind,nDimBC,Label='Ind')
  iInd = 0
  do i=1,3*nsAtom
    iAtom = (i+2)/3
    ixyz = i-(iAtom-1)*3
    if (Smmtrc(ixyz,iAtom)) then
      iInd = iInd+1
      Ind(iInd) = i
    end if
  end do

  ! Compute H|i>

  call mma_allocate(Hi,nDimBC,mTR,Label='Hi')
  call mma_allocate(iHi,mTR,Label='iHi')
  Hi(:,:) = Zero

  do j=1,mTR
    do i=1,nDimBC
      Temp = Zero
      do k=1,nDimBC
        kx = Ind(k)
        iAtom = (kx+2)/3
        ixyz = kx-(iAtom-1)*3
        ik = (i-1)*nDimBC+k
        Temp = Temp+Hss_X(ik)*sqrt(Degen(ixyz,iAtom))*TRVec(k,j)
      end do
      Hi(i,j) = Temp
    end do
  end do
  !call RecPrt('Hi',' ',Hi,nDimBC,mTR)
  do iTR=1,mTR
    iHi(iTR) = DDot_(nDimBC,TRVec(1,iTR),1,Hi(:,iTR),1)
  end do
  !call RecPrt('iHi',' ',iHi,mTR,1)

  do i=1,nDimBC
    ix = Ind(i)
    iAtom = (ix+2)/3
    ixyz = ix-(iAtom-1)*3
    do j=1,i
      jx = Ind(j)
      jAtom = (jx+2)/3
      jxyz = jx-(jAtom-1)*3
      ij = (j-1)*nDimBC+i
      ji = (i-1)*nDimBC+j
      Temp = Half*(Hss_x(ij)+Hss_x(ji))
#     ifndef UNIT_MM

      ! Here we shift the eigenvectors corresponding to
      ! translations and rotations to a large positive values.

      do iTR=1,mTR
        Omega = 1.0e5_wp
        Hii = iHi(iTR)
        Temp = Temp+sqrt(Degen(ixyz,iAtom))* &
               (-TRVec(i,iTR)*Hi(j,iTR)-Hi(i,iTR)*TRVec(j,iTR)+TRVec(i,iTR)*(Omega+Hii)*TRVec(j,iTR))*sqrt(Degen(jxyz,jAtom))
      end do
#     endif

      Hss_X(ij) = Temp
      Hss_X(ji) = Temp
    end do
  end do
  call mma_deallocate(iHi)
  call mma_deallocate(Hi)
  call mma_deallocate(Ind)

  ! Clean up the gradient wrt translational and rotational component.
  !
  ! |g> = |g> - Sum(TR) |i><i|g>

  if (BSet) then

    !call RecPrt('Gx',' ',Gx(:,:,nIter),1,3*nsAtom)
    !call RecPrt('TRVec',' ',TRVec,nDimBC,mTR)

    do iTR=1,mTR

      ! <i|g>

      Temp = Zero
      iInd = 0
      do iAtom=1,nsAtom
        do j=1,3
          if (Smmtrc(j,iAtom)) then
            iInd = iInd+1
            Temp = Temp+Degen(j,iAtom)*Gx(j,iAtom,nIter)*TRVec(iInd,iTR)
          end if
        end do
      end do

      iInd = 0
      do iAtom=1,nsAtom
        do j=1,3
          if (Smmtrc(j,iAtom)) then
            iInd = iInd+1
            Gx(j,iAtom,nIter) = Gx(j,iAtom,nIter)-TRVec(iInd,iTR)*Temp
          end if
        end do
      end do

    end do

    !call RecPrt('Gx',' ',Gx(:,:,nIter),1,3*nsAtom)

  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !   N O N - R E D U N D A N T  C A R T E S I A N  C O O R D S

  nQQ = nDimBC-mTR
  if (allocated(qInt)) then
    if (size(qInt,1) /= nQQ) then
      call mma_deallocate(qInt)
      call mma_deallocate(dqInt)
    end if
  end if
  if (.not. allocated(qInt)) then
    call mma_allocate(qInt,nQQ,MaxItr,Label='qInt')
    call mma_allocate(dqInt,nQQ,MaxItr,Label='dqInt')
    qInt(:,:) = Zero
    dqInt(:,:) = Zero
  end if
  if (allocated(dqInt_Aux)) then
    if (size(dqInt_Aux,1) /= nQQ) call mma_deallocate(dqInt_Aux)
  end if
  if ((.not. allocated(dqInt_Aux)) .and. (nSet > 1)) then
    call mma_allocate(dqInt_Aux,nQQ,MaxItr,nSet-1,Label='dqInt_Aux')
    dqInt_Aux(:,:,:) = Zero
  end if

  ! Project the model Hessian with respect to rotations and
  ! translations. The eigenvalues are shifted to negative eigenvalues.

  call mma_allocate(Ind,nDimBC,Label='Ind')
  iInd = 0
  do i=1,3*nsAtom
    iAtom = (i+2)/3
    ixyz = i-(iAtom-1)*3
    if (Smmtrc(ixyz,iAtom)) then
      iInd = iInd+1
      Ind(iInd) = i
    end if
  end do

  ! Compute H|i>

  call mma_allocate(Hi,nDimBC,mTR,Label='Hi')
  call mma_allocate(iHi,mTR,Label='iHi')
  Hi(:,:) = Zero

  do j=1,mTR
    do i=1,nDimBC
      Temp = Zero
      do k=1,nDimBC
        kx = Ind(k)
        kAtom = (kx+2)/3
        kxyz = kx-(kAtom-1)*3
        ik = (i-1)*nDimBC+k
        Temp = Temp+Hss_X(ik)*sqrt(Degen(kxyz,kAtom))*TRVec(k,j)
      end do
      Hi(i,j) = Temp
    end do
  end do
  do iTR=1,mTR
    iHi(iTR) = DDot_(nDimBC,TRVec(1,iTR),1,Hi(:,iTR),1)
  end do
# ifdef _DEBUGPRINT_
  call RecPrt('Hi',' ',Hi,nDimBC,mTR)
  call RecPrt('iHi',' ',iHi,mTR,1)
# endif

  do i=1,nDimBC
    ix = Ind(i)
    iAtom = (ix+2)/3
    ixyz = ix-(iAtom-1)*3
    do j=1,i
      jx = Ind(j)
      jAtom = (jx+2)/3
      jxyz = jx-(jAtom-1)*3
      ijTri = i*(i-1)/2+j
      ij = (j-1)*nDimBC+i
      ji = (i-1)*nDimBC+j
      EVal(ijTri) = Half*(Hss_X(ij)+Hss_X(ji))

      ! Here we shift the eigenvectors corresponding to translations
      ! and rotations down to negative faked eigenvalues.

      do iTR=1,mTR
        Omega = -real(iTR,kind=wp)
        Hii = iHi(iTR)
        Eval(ijTri) = Eval(ijTri)+sqrt(Degen(ixyz,iAtom))* &
                      (-TRVec(i,iTR)*Hi(j,iTR)-Hi(i,iTR)*TRVec(j,iTR)+TRVec(i,iTR)*(Omega+Hii)*TRVec(j,iTR))*sqrt(Degen(jxyz,jAtom))
      end do

      Hss_X(ij) = EVal(ijTri)
      Hss_X(ji) = EVal(ijTri)
    end do
  end do
  call mma_deallocate(iHi)
  call mma_deallocate(Hi)
  call mma_deallocate(Ind)
# ifdef _DEBUGPRINT_
  call TriPrt(' The Projected Model Hessian','(5G20.10)',EVal,nDimBC)
  call RecPrt(' The Projected Model Hessian','(5G20.10)',Hss_x,nDimBC,nDimBC)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute the eigenvectors for the Cartesian Hessian

  call mma_allocate(EVec,(3*mTtAtm)**2,Label='EVec')
  call Hess_Vec(mTtAtm,EVal,EVec,nsAtom,nDimBC)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Move over the eigenvectors putting to BMx

  call mma_allocate(BMx,3*nsAtom,3*nsAtom,Label='BMx')
  BMx(:,:) = Zero
  ipFrom = 1+mTR*nDimBC
  call BPut(EVec(ipFrom),nDimBC,BMx,3*nsAtom,Smmtrc,nQQ,Degen)
# ifdef _DEBUGPRINT_
  call RecPrt('In Bmtrx: B',' ',BMx,3*nsAtom,nQQ)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (PrQ .and. (nsAtom <= 5)) call List2('Cartesian Approximate Normal Modes',AtomLbl,BMx,nsAtom,nQQ,Smmtrc)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
if (HSet .and. (.not. lOld)) then
  call mma_allocate(KtB,nDimBC,nQQ,Label='KtB')

  call mma_allocate(Degen2,nDimBC,Label='Degen2')
  i = 0
  do ix=1,3*nsAtom
    iAtom = (ix+2)/3
    ixyz = ix-(iAtom-1)*3
    if (Smmtrc(ixyz,iAtom)) then
      i = i+1
      Degen2(i) = Degen(ixyz,iAtom)
    end if
  end do

  call dcopy_(nDimBC*nQQ,EVec(ipFrom),1,KtB,1)
  do iInter=1,nQQ
    do i_Dim=1,nDimBC
      KtB(i_Dim,iInter) = KtB(i_Dim,iInter)/sqrt(Degen2(i_Dim))
    end do
  end do
  call mma_deallocate(Degen2)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(EVec)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the value and gradient vectors in the new basis.

call ValANM(nsAtom,nQQ,nIter,BMx,Degen,qInt,Cx,'Values',nWndw)
if (BSet) then
  call ValANM(nsAtom,nQQ,nIter,BMx,Degen,dqInt,Gx,'Gradients',nWndw)
  if (nSet > 1) then
    call ValANM(nsAtom,nQQ,nIter,BMx,Degen,dqInt_Aux(:,:,1),Gx0,'Gradients',nWndw)
  end if
  if (nSet > 2) then
    call ValANM(nsAtom,nQQ,nIter,BMx,Degen,dqInt_Aux(:,:,2),NAC,'Gradients',nWndw)
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine BMtrx_Cartesian
