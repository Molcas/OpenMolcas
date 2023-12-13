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
subroutine rotder(nmass,xmass,currxyz,ref123,trans,rotang,rotvec,rotmat,norder,dRVdXYZ,d2RVdXYZ2)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nmass, norder
real(kind=wp), intent(in) :: xMass(nMass), CurrXYZ(3,nMass), Ref123(3,nMass)
real(kind=wp), intent(out) :: trans(3), RotAng, RotMat(3,3), dRVdXYZ(3,3*nMass)
real(kind=wp), intent(inout) :: RotVec(3), d2RVdXYZ2(3,3*nMass,3*nMass)
integer(kind=iwp) :: i, imass, ip, ip1, ip2, ipk, ipk1, ipk2, iter, j, j1, k, k1, k2
real(kind=wp) :: d2RVdA2(3,3,3), d3RVdA3(3,3,3,3), d4RVdA4(3,3,3,3,3), dRVdA(3,3), MOI(3,3), MOIInv(3,3), RotErr, rSum, &
                 Sayvetz2(3), SmallRot(3), SMat(3,3), sval(3), TotMass, Trace, umat(3,3), vmat(3,3), wTmp(100) !, &
                 !det, detinv, G(3,3)
real(kind=wp), allocatable :: Curr123(:,:), dAdXYZ(:,:,:), tmp(:,:)

#ifdef _DEBUGPRINT_
call RecPrt('Ref123',' ',Ref123,3,nMass)
#endif
call mma_allocate(Curr123,3,nMass,label='Curr123')
call mma_allocate(tmp,3,3*nMass,label='tmp')
call mma_allocate(dAdXYZ,3,3,nMass,label='dAdXYZ')

! Compute the Center-of-Mass coordinates.
! These are the same as the translation vector.
Trans(:) = Zero
TotMass = Zero
do imass=1,nmass
  TotMass = TotMass+xMass(imass)
  Trans(:) = Trans(:)+xMass(imass)*CurrXYZ(:,imass)
end do
Trans(:) = Trans(:)/TotMass
! A lengthy piece of code determines the exact orientation
! of the internal frame. This is defined by Sayvetz conditions
! using a reference conformation. The procedure is iterative.
! At least an approximate rotation vector is assumed to be
! known when entering the Subroutine.
! Given rotation vector (''BigOmega'' in formulas),
! compute scalar rotation angle, and rotation matrix.
#ifdef _DEBUGPRINT_
call RecPrt('RotVec(00)',' ',RotVec,1,3)
#endif
RotAng = sqrt(RotVec(1)**2+RotVec(2)**2+RotVec(3)**2)
call mkRotMat(RotVec,RotMat)
#ifdef _DEBUGPRINT_
call RecPrt('RotVec(0)',' ',RotVec,1,3)
#endif
iter = 0
do
  iter = iter+1
  if (iter >= 100) then
    !call Quit_OnConvError()
    call WarningMessage(1,'Warning: Convergence problem in the internal frame')
    exit
  end if
  ! compute cartesian position vectors in the co-moving frame:
  ! (Note: Inverse rotation matrix = transpose)
  do imass=1,nmass
    do i=1,3
      rSum = Zero
      do j=1,3
        rSum = rSum+RotMat(j,i)*(CurrXYZ(j,imass)-Trans(j))
      end do
      Curr123(i,imass) = rSum
    end do
  end do
  ! Compute the Moments-of-Inertia matrix. The asymmetrical
  ! one, used for Sayvetz conditions.
  tmp(:,1:nmass) = Ref123(:,:)
# ifdef _DEBUGPRINT_
  call RecPrt('Curr123',' ',Curr123,3,nMass)
  call RecPrt('tmp',' ',tmp,3,nMass)
# endif
  do imass=1,nMass
    tmp(:,imass) = xMass(iMass)*tmp(:,imass)
  end do
  call DGEMM_('N','T',3,3,nmass,One,tmp,3,Curr123,3,Zero,SMat,3)
  Trace = (SMat(1,1)+SMat(2,2)+SMat(3,3))
  MOI(:,:) = -SMat(:,:)
  do i=1,3
    MOI(i,i) = MOI(i,i)+Trace
  end do
# ifdef _DEBUGPRINT_
  call RecPrt('MOI',' ',MOI,size(MOI,1),size(MOI,2))
# endif
  ! Invert the MOI-matrix.
  !     G(1,1) = MOI(2,2)*MOI(3,3)-MOI(2,3)*MOI(3,2)
  !     G(2,1) = MOI(3,2)*MOI(1,3)-MOI(3,3)*MOI(1,2)
  !     G(3,1) = MOI(1,2)*MOI(2,3)-MOI(1,3)*MOI(2,2)
  !     G(1,2) = MOI(2,3)*MOI(3,1)-MOI(2,1)*MOI(3,3)
  !     G(2,2) = MOI(3,3)*MOI(1,1)-MOI(3,1)*MOI(1,3)
  !     G(3,2) = MOI(1,3)*MOI(2,1)-MOI(1,1)*MOI(2,3)
  !     G(1,3) = MOI(2,1)*MOI(3,2)-MOI(2,2)*MOI(3,1)
  !     G(2,3) = MOI(3,1)*MOI(1,2)-MOI(3,2)*MOI(1,1)
  !     G(3,3) = MOI(1,1)*MOI(2,2)-MOI(1,2)*MOI(2,1)
  !     Det = MOI(1,1)*G(1,1)+MOI(2,1)*G(2,1)+MOI(3,1)*G(3,1)
  !     DetInv = One/Det
  !     MOIInv(1,1) = G(1,1)*DetInv
  !     MOIInv(2,1) = G(1,2)*DetInv
  !     MOIInv(3,1) = G(1,3)*DetInv
  !     MOIInv(1,2) = G(2,1)*DetInv
  !     MOIInv(2,2) = G(2,2)*DetInv
  !     MOIInv(3,2) = G(2,3)*DetInv
  !     MOIInv(1,3) = G(3,1)*DetInv
  !     MOIInv(2,3) = G(3,2)*DetInv
  !     MOIInv(3,3) = G(3,3)*DetInv
  ! (use the Moore-Penrose pseudoinverse, to deal with linear systems)
  call dgesvd_('A','A',3,3,MOI,3,sval,umat,3,vmat,3,wTmp,100,i)
  do i=1,3
    if (abs(sval(i)) > 1.0e-12_wp) then
      umat(:,i) = umat(:,i)/sval(i)
    else
      umat(:,i) = Zero
    end if
  end do
  call dgemm_('T','T',3,3,3,One,vmat,3,umat,3,Zero,MOIInv,3)
# ifdef _DEBUGPRINT_
  call RecPrt('vmat',' ',vmat,size(vmat,1),size(vmat,2))
  call RecPrt('umat',' ',umat,size(umat,1),size(umat,2))
  call RecPrt('MOI',' ',MOIInv,size(MOIInv,1),size(MOIInv,2))
  call RecPrt('xMass',' ',xMass,1,nMass)
  call RecPrt('Ref123',' ',Ref123,3,nMass)
  call RecPrt('Curr123',' ',Curr123,3,nMass)
# endif
  ! Determine small rotation that zeroes the Sayvetz rotational conditions
  !write(u6,*) ' Determine small rotation. First compute Sayvetz2:'
  Sayvetz2(:) = Zero
  do i=1,3
    j = 1+mod(i,3)
    k = 1+mod(j,3)
    do imass=1,nMass
      Sayvetz2(i) = Sayvetz2(i)+xMass(imass)*(Ref123(j,imass)*Curr123(k,imass)-Ref123(k,imass)*Curr123(j,imass))
    end do
  end do
# ifdef _DEBUGPRINT_
  call RecPrt('Sayvetz2',' ',Sayvetz2,1,3)
# endif
  RotErr = Zero
  do i=1,3
    rSum = Zero
    do j=1,3
      rSum = rSum+MOIInv(j,i)*Sayvetz2(j)
    end do
    SmallRot(i) = rSum
    RotErr = RotErr+rSum**2
  end do
# ifdef _DEBUGPRINT_
  call RecPrt('SmallRot',' ',SmallRot,1,3)
# endif
  RotErr = sqrt(RotErr)
  ! Limiting step size:
  if (RotErr > One) then
    SmallRot(1) = SmallRot(1)/RotErr
    SmallRot(2) = SmallRot(2)/RotErr
    SmallRot(3) = SmallRot(3)/RotErr
    RotErr = One
  end if
  ! Apply this rotation to the rotation matrix:
  call updRotMat(SmallRot,RotMat)
# ifdef _DEBUGPRINT_
  call RecPrt('RotMat(i)',' ',RotMat,3,3)
  call RecPrt('SmallRot',' ',SmallRot,1,3)
# endif
  if (RotErr <= 1.0e-12_wp) exit
end do
#ifdef _DEBUGPRINT_
call RecPrt('RotMat(Final)',' ',RotMat,3,3)
#endif
! Now RotMat is converged. Recompute RotVec, RotAng.
call Mat2Vec(RotMat,RotVec,RotAng)
#ifdef _DEBUGPRINT_
call RecPrt('RotVec(Recomputed)',' ',RotVec,1,3)
#endif
! Derivatives w.r.t A(i), of global rotation parameters RotVec
! Note: A(i) is the dual vector of the antisymmetric part of the
! asymmetrically defined MOI matrix, with frozen orientation.
! RotVec is the rotation parameter array necessary to make MOI
! symmetrical.
call rotder4(norder,SMat,RotVec,dRVdA,d2RVdA2,d3RVdA3,d4RVdA4)
! Note that A(i) is exactly linear in atom coordinates.
! Definition of A(i) is as follows:
!    A(1) = sum(xMass(imass)*(Ref123(2,imass)*Curr123(3,imass)-Ref123(3,imass)*Curr123(2,imass)) )
!    A(2) = sum(xMass(imass)*(Ref123(3,imass)*Curr123(1,imass)-Ref123(1,imass)*Curr123(3,imass)) )
!    A(3) = sum(xMass(imass)*(Ref123(1,imass)*Curr123(2,imass)-Ref123(2,imass)*Curr123(1,imass)) )
! First, compute derivatives without accounting for COM motion:
do ip=1,nmass
  dAdXYZ(1,:,ip) = xmass(ip)*(RotMat(:,3)*Ref123(2,ip)-RotMat(:,2)*Ref123(3,ip))
  dAdXYZ(2,:,ip) = xmass(ip)*(RotMat(:,1)*Ref123(3,ip)-RotMat(:,3)*Ref123(1,ip))
  dAdXYZ(3,:,ip) = xmass(ip)*(RotMat(:,2)*Ref123(1,ip)-RotMat(:,1)*Ref123(2,ip))
end do
! Correction for COM motion: (Theoretically, this makes no
! difference to the results, but why not...)
do i=1,3
  do k=1,3
    rSum = Zero
    do ip=1,nmass
      rSum = rSum+dAdXYZ(i,k,ip)
    end do
    dAdXYZ(i,k,:) = dAdXYZ(i,k,:)-(xMass(:)/TotMass)*rSum
  end do
end do
! Finally,transform dXdA, etc, into derivatives w.r.t atom coordinates:
!    dRVdXYZ(i,k,ip) = sum dRVdA(i,j)*dAdXYZ(j,k,ip)
!    d2RVdXYZ2(i,k1,ip1,k2,ip2) = sum d2RVdA2(i,j1,j2)*dAdXYZ(j1,k1,ip1)*dAdXYZ(j2,k2,ip2)
! and so on. Note: A is linear in nuclear coordinates.
if (nOrder >= 1) then
  do i=1,3
    do ip=1,nMass
      do k=1,3
        ipk = k+3*(ip-1)
        dRVdXYZ(i,ipk) = sum(dRVdA(i,:)*dAdXYZ(:,k,ip))
      end do
    end do
  end do
end if
if (nOrder >= 2) then
  do i=1,3
    do j1=1,3
      do ip2=1,nMass
        do k2=1,3
          ipk2 = k2+3*(ip2-1)
          tmp(j1,ipk2) = sum(d2RVdA2(i,j1,:)*dAdXYZ(:,k2,ip2))
        end do
      end do
    end do
    do ipk2=1,3*nMass
      do ip1=1,nMass
        do k1=1,3
          ipk1 = k1+3*(ip1-1)
          d2RVdXYZ2(i,ipk1,ipk2) = sum(tmp(:,ipk2)*dAdXYZ(:,k1,ip1))
        end do
      end do
    end do
  end do
end if
call mma_deallocate(Curr123)
call mma_deallocate(tmp)
call mma_deallocate(dAdXYZ)

return

end subroutine rotder
