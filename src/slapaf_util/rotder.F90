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

implicit none
integer nmass, norder
real*8 xMass(nMass), CurrXYZ(3,nMass), Ref123(3,nMass)
real*8 trans(3), RotAng, RotVec(3), RotMat(3,3)
real*8 dRVdXYZ(3,3*nMass)
real*8 d2RVdXYZ2(3,3*nMass,3*nMass)
! Local variables:
integer i, imass, ip, ip1, ip2, ipk, ipk1, ipk2
integer iter, j, j1, j2, k, k1, k2
real*8 dRVdA(3,3), d2RVdA2(3,3,3), d3RVdA3(3,3,3,3)
real*8 d4RVdA4(3,3,3,3,3)
real*8 TotMass, MOI(3,3), MOIInv(3,3)
real*8 Sum, Trace, SMat(3,3)
!real*8 G(3,3),det,detinv
real*8 Sayvetz2(3), RotErr, SmallRot(3)
real*8 umat(3,3), vmat(3,3), sval(3), wTmp(100)
real*8, dimension(:,:), allocatable :: Curr123, tmp
real*8, dimension(:,:,:), allocatable :: dAdXYZ

#ifdef _DEBUGPRINT_
call RecPrt('Ref123',' ',Ref123,3,nMass)
#endif
call mma_allocate(Curr123,3,nMass,label='Curr123')
call mma_allocate(tmp,3,3*nMass,label='tmp')
call mma_allocate(dAdXYZ,3,3,nMass,label='dAdXYZ')

! Compute the Center-of-Mass coordinates.
! These are the same as the translation vector.
Trans(1) = 0.0d0
Trans(2) = 0.0d0
Trans(3) = 0.0d0
TotMass = 0.0d0
do imass=1,nmass
  TotMass = TotMass+xMass(imass)
  Trans(1) = Trans(1)+xMass(imass)*CurrXYZ(1,imass)
  Trans(2) = Trans(2)+xMass(imass)*CurrXYZ(2,imass)
  Trans(3) = Trans(3)+xMass(imass)*CurrXYZ(3,imass)
end do
Trans(1) = Trans(1)/TotMass
Trans(2) = Trans(2)/TotMass
Trans(3) = Trans(3)/TotMass
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
10 continue
iter = iter+1
if (iter >= 100) then
  !call Quit_OnConvError()
  call WarningMessage(1,'Warning: Convergence problem in the internal frame')
  goto 11
end if
! compute cartesian position vectors in the co-moving frame:
! (Note: Inverse rotation matrix = transpose)
do imass=1,nmass
  do i=1,3
    sum = 0.0d0
    do j=1,3
      sum = sum+RotMat(j,i)*(CurrXYZ(j,imass)-Trans(j))
    end do
    Curr123(i,imass) = sum
  end do
end do
! Compute the Moments-of-Inertia matrix. The asymmetrical
! one, used for Sayvetz conditions.
call dcopy_(3*nmass,Ref123,1,tmp,1)
#ifdef _DEBUGPRINT_
call RecPrt('Curr123',' ',Curr123,3,nMass)
call RecPrt('tmp',' ',tmp,3,nMass)
#endif
do imass=1,nMass
  do i=1,3
    tmp(i,imass) = xMass(iMass)*tmp(i,imass)
  end do
end do
call DGEMM_('N','T',3,3,nmass,1.0d0,tmp,3,Curr123,3,0.0d0,SMat,3)
Trace = (SMat(1,1)+SMat(2,2)+SMat(3,3))
do i=1,3
  do j=1,3
    MOI(i,j) = -SMat(i,j)
  end do
  MOI(i,i) = Trace-SMat(i,i)
end do
#ifdef _DEBUGPRINT_
call RecPrt('MOI',' ',MOI,size(MOI,1),size(MOI,2))
#endif
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
!     DetInv = 1.0D0/Det
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
  if (abs(sval(i)) > 1.0d-12) then
    call dscal_(3,1.0d0/sval(i),umat(1,i),1)
  else
    call dcopy_(3,[0.0d0],0,umat(1,i),1)
  end if
end do
call dgemm_('T','T',3,3,3,1.0d0,vmat,3,umat,3,0.0d0,MOIInv,3)
#ifdef _DEBUGPRINT_
call RecPrt('vmat',' ',vmat,size(vmat,1),size(vmat,2))
call RecPrt('umat',' ',umat,size(umat,1),size(umat,2))
call RecPrt('MOI',' ',MOIInv,size(MOIInv,1),size(MOIInv,2))
call RecPrt('xMass',' ',xMass,1,nMass)
call RecPrt('Ref123',' ',Ref123,3,nMass)
call RecPrt('Curr123',' ',Curr123,3,nMass)
#endif
! Determine small rotation that zeroes the Sayvetz rotational conditions
!write(6,*) ' Determine small rotation. First compute Sayvetz2:'
do i=1,3
  j = 1+mod(i,3)
  k = 1+mod(j,3)
  Sayvetz2(i) = 0.0d0
  do imass=1,nMass
    Sayvetz2(i) = Sayvetz2(i)+xMass(imass)*(Ref123(j,imass)*Curr123(k,imass)-Ref123(k,imass)*Curr123(j,imass))
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('Sayvetz2',' ',Sayvetz2,1,3)
#endif
RotErr = 0.0d0
do i=1,3
  sum = 0.0d0
  do j=1,3
    sum = sum+MOIInv(j,i)*Sayvetz2(j)
  end do
  SmallRot(i) = +sum
  RotErr = RotErr+sum**2
end do
#ifdef _DEBUGPRINT_
call RecPrt('SmallRot',' ',SmallRot,1,3)
#endif
RotErr = sqrt(RotErr)
! Limiting step size:
if (RotErr > 1.0d0) then
  SmallRot(1) = SmallRot(1)/RotErr
  SmallRot(2) = SmallRot(2)/RotErr
  SmallRot(3) = SmallRot(3)/RotErr
  RotErr = 1.0d0
end if
! Apply this rotation to the rotation matrix:
call updRotMat(SmallRot,RotMat)
#ifdef _DEBUGPRINT_
call RecPrt('RotMat(i)',' ',RotMat,3,3)
call RecPrt('SmallRot',' ',SmallRot,1,3)
#endif
if (RotErr > 1.0D-12) goto 10
11 continue
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
  do k=1,3
    dAdXYZ(1,k,ip) = xmass(ip)*(RotMat(k,3)*Ref123(2,ip)-RotMat(k,2)*Ref123(3,ip))
    dAdXYZ(2,k,ip) = xmass(ip)*(RotMat(k,1)*Ref123(3,ip)-RotMat(k,3)*Ref123(1,ip))
    dAdXYZ(3,k,ip) = xmass(ip)*(RotMat(k,2)*Ref123(1,ip)-RotMat(k,1)*Ref123(2,ip))
  end do
end do
! Correction for COM motion: (Theoretically, this makes no
! difference to the results, but why not...)
do i=1,3
  do k=1,3
    sum = 0.0d0
    do ip=1,nmass
      sum = sum+dAdXYZ(i,k,ip)
    end do
    do ip=1,nmass
      dAdXYZ(i,k,ip) = dAdXYZ(i,k,ip)-(xMass(ip)/TotMass)*sum
    end do
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
        sum = 0.0d0
        do j=1,3
          sum = sum+dRVdA(i,j)*dAdXYZ(j,k,ip)
        end do
        dRVdXYZ(i,ipk) = sum
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
          sum = 0.0d0
          do j2=1,3
            sum = sum+d2RVdA2(i,j1,j2)*dAdXYZ(j2,k2,ip2)
          end do
          tmp(j1,ipk2) = sum
        end do
      end do
    end do
    do ipk2=1,3*nMass
      do ip1=1,nMass
        do k1=1,3
          ipk1 = k1+3*(ip1-1)
          sum = 0.0d0
          do j1=1,3
            sum = sum+tmp(j1,ipk2)*dAdXYZ(j1,k1,ip1)
          end do
          d2RVdXYZ2(i,ipk1,ipk2) = sum
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
