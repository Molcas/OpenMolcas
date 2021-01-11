************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine rotder(nmass,xmass,currxyz,ref123,trans,rotang,
     &   rotvec,rotmat,norder,dRVdXYZ,d2RVdXYZ2)
      implicit none
#include "stdalloc.fh"
      Integer nmass,norder
      real*8 xMass(nMass),CurrXYZ(3,nMass),Ref123(3,nMass)
      real*8 trans(3),RotAng,RotVec(3), RotMat(3,3)
      real*8 dRVdXYZ(3,3*nMass)
      real*8 d2RVdXYZ2(3,3*nMass,3*nMass)
* Local variables:
      Integer i,imass,ip,ip1,ip2,ipk,ipk1,ipk2
      Integer iter,j,j1,j2,k,k1,k2
      real*8 dRVdA(3,3),d2RVdA2(3,3,3),d3RVdA3(3,3,3,3)
      real*8 d4RVdA4(3,3,3,3,3)
      real*8 TotMass,MOI(3,3),MOIInv(3,3)
      real*8 Sum,Trace,SMat(3,3)
c     real*8 G(3,3),det,detinv
      real*8 Sayvetz2(3),RotErr,SmallRot(3)
      real*8 umat(3,3),vmat(3,3),sval(3),wTmp(100)
      real*8, dimension(:,:), allocatable :: Curr123,tmp
      real*8, dimension(:,:,:), allocatable :: dAdXYZ
*
      Call mma_allocate(Curr123,3,nMass,label='Curr123')
      Call mma_allocate(tmp,3,3*nMass,label='tmp')
      Call mma_allocate(dAdXYZ,3,3,nMass,label='dAdXYZ')

* Compute the Center-of-Mass coordinates.
* These are the same as the translation vector.
      Trans(1)=0.0D0
      Trans(2)=0.0D0
      Trans(3)=0.0D0
      TotMass=0.0D0
      do imass=1,nmass
        TotMass=TotMass+xMass(imass)
        Trans(1)=Trans(1)+xMass(imass)*CurrXYZ(1,imass)
        Trans(2)=Trans(2)+xMass(imass)*CurrXYZ(2,imass)
        Trans(3)=Trans(3)+xMass(imass)*CurrXYZ(3,imass)
      end do
      Trans(1)=Trans(1)/TotMass
      Trans(2)=Trans(2)/TotMass
      Trans(3)=Trans(3)/TotMass
* A lengthy piece of code determines the exact orientation
* of the internal frame. This is defined by Sayvetz conditions
* using a reference conformation. The procedure is iterative.
* At least an approximate rotation vector is assumed to be
* known when entering the Subroutine.
* Given rotation vector (''BigOmega'' in formulas),
* compute scalar rotation angle, and rotation matrix.
      RotAng=sqrt(RotVec(1)**2+RotVec(2)**2+RotVec(3)**2)
      Call mkRotMat(RotVec,RotMat)
      iter=0
  10  Continue
      iter=iter+1
      if(iter.ge.100) then
c       Call Quit_OnConvError()
        Call WarningMessage(1,
     &       'Warning: Convergence problem in the internal frame')
        goto 11
      end if
* compute cartesian position vectors in the co-moving frame:
* (Note: Inverse rotation matrix = transpose)
      do imass=1,nmass
       do i=1,3
        sum=0.0D0
        do j=1,3
         sum=sum+RotMat(j,i)*(CurrXYZ(j,imass)-Trans(j))
        end do
        Curr123(i,imass)=sum
       end do
      end do
* Compute the Moments-of-Inertia matrix. The asymmetrical
* one, used for Sayvetz conditions.
      call dcopy_(3*nmass,Ref123,1,tmp,1)
      Do imass=1,nMass
       Do i=1,3
        tmp(i,imass)=xMass(iMass)*tmp(i,imass)
       End Do
      End Do
      Call DGEMM_('N','T',3,3,nmass,1.0D0,tmp,3,Curr123,3,0.0D0,SMat,3)
      Trace=(SMat(1,1)+SMat(2,2)+SMat(3,3))
      do i=1,3
       do j=1,3
        MOI(i,j)=-SMat(i,j)
       end do
       MOI(i,i)=Trace-SMat(i,i)
      end do
* Invert the MOI-matrix.
c     G(1,1)=MOI(2,2)*MOI(3,3)-MOI(2,3)*MOI(3,2)
c     G(2,1)=MOI(3,2)*MOI(1,3)-MOI(3,3)*MOI(1,2)
c     G(3,1)=MOI(1,2)*MOI(2,3)-MOI(1,3)*MOI(2,2)
c     G(1,2)=MOI(2,3)*MOI(3,1)-MOI(2,1)*MOI(3,3)
c     G(2,2)=MOI(3,3)*MOI(1,1)-MOI(3,1)*MOI(1,3)
c     G(3,2)=MOI(1,3)*MOI(2,1)-MOI(1,1)*MOI(2,3)
c     G(1,3)=MOI(2,1)*MOI(3,2)-MOI(2,2)*MOI(3,1)
c     G(2,3)=MOI(3,1)*MOI(1,2)-MOI(3,2)*MOI(1,1)
c     G(3,3)=MOI(1,1)*MOI(2,2)-MOI(1,2)*MOI(2,1)
c     Det=MOI(1,1)*G(1,1)+MOI(2,1)*G(2,1)+MOI(3,1)*G(3,1)
c     DetInv=1.0D0/Det
c     MOIInv(1,1)=G(1,1)*DetInv
c     MOIInv(2,1)=G(1,2)*DetInv
c     MOIInv(3,1)=G(1,3)*DetInv
c     MOIInv(1,2)=G(2,1)*DetInv
c     MOIInv(2,2)=G(2,2)*DetInv
c     MOIInv(3,2)=G(2,3)*DetInv
c     MOIInv(1,3)=G(3,1)*DetInv
c     MOIInv(2,3)=G(3,2)*DetInv
c     MOIInv(3,3)=G(3,3)*DetInv
* (use the Moore-Penrose pseudoinverse, to deal with linear systems)
      call dgesvd_('A','A',3,3,MOI,3,sval,umat,3,vmat,3,wTmp,100,i)
      Do i=1,3
        If (Abs(sval(i)).gt.1.0d-12) Then
          call dscal_(3,1.0d0/sval(i),umat(1,i),1)
        Else
          call dcopy_(3,[0.0d0],0,umat(1,i),1)
        End If
      End Do
      Call dgemm_('T','T',3,3,3,1.0d0,vmat,3,umat,3,0.0d0,MOIInv,3)
* Determine small rotation that zeroes the Sayvetz rotational conditions
*      write(*,*)' Determine small rotation. First compute Sayvetz2:'
      do i=1,3
       j=1+mod(i,3)
       k=1+mod(j,3)
       sum=0.0D0
       do imass=1,nMass
        sum=sum+xMass(imass)*(Ref123(j,imass)*Curr123(k,imass)-
     &                       Ref123(k,imass)*Curr123(j,imass))
       end do
       Sayvetz2(i)=sum
      end do
      RotErr=0.0D0
      do i=1,3
       sum=0.0D0
       do j=1,3
        sum=sum+MOIInv(j,i)*Sayvetz2(j)
       end do
       SmallRot(i)=+sum
       RotErr=RotErr+sum**2
      end do
      RotErr=sqrt(RotErr)
* Limiting step size:
      If(RotErr.gt.1.0D0) then
        SmallRot(1)=SmallRot(1)/RotErr
        SmallRot(2)=SmallRot(2)/RotErr
        SmallRot(3)=SmallRot(3)/RotErr
        RotErr=1.0D0
      End If
* Apply this rotation to the rotation matrix:
      Call updRotMat(SmallRot,RotMat)
      if(RotErr.gt.1.0D-12) goto 10
  11  continue
* Now RotMat is converged. Recompute RotVec, RotAng.
      Call Mat2Vec(RotMat,RotVec,RotAng)
* Derivatives w.r.t A(i), of global rotation parameters RotVec
* Note: A(i) is the dual vector of the antisymmetric part of the
* asymmetrically defined MOI matrix, with frozen orientation.
* RotVec is the rotation parameter array necessary to make MOI
* symmetrical.
      call rotder4(norder,SMat,RotVec,dRVdA,d2RVdA2,d3RVdA3,d4RVdA4)
* Note that A(i) is exactly linear in atom coordinates.
* Definition of A(i) is as follows:
*    A(1)=sum( xMass(imass)*(Ref123(2,imass)*Curr123(3,imass)-
*                            Ref123(3,imass)*Curr123(2,imass)) )
*    A(2)=sum( xMass(imass)*(Ref123(3,imass)*Curr123(1,imass)-
*                            Ref123(1,imass)*Curr123(3,imass)) )
*    A(3)=sum( xMass(imass)*(Ref123(1,imass)*Curr123(2,imass)-
*                            Ref123(2,imass)*Curr123(1,imass)) )
* First, compute derivatives without accounting for COM motion:
      do ip=1,nmass
       do k=1,3
        dAdXYZ(1,k,ip)=xmass(ip)*
     &          (RotMat(k,3)*Ref123(2,ip)-RotMat(k,2)*Ref123(3,ip))
        dAdXYZ(2,k,ip)=xmass(ip)*
     &          (RotMat(k,1)*Ref123(3,ip)-RotMat(k,3)*Ref123(1,ip))
        dAdXYZ(3,k,ip)=xmass(ip)*
     &          (RotMat(k,2)*Ref123(1,ip)-RotMat(k,1)*Ref123(2,ip))
       end do
      end do
* Correction for COM motion: (Theoretically, this makes no
* difference to the results, but why not...)
      do i=1,3
       do k=1,3
        sum=0.0d0
        do ip=1,nmass
         sum=sum+dAdXYZ(i,k,ip)
        end do
        do ip=1,nmass
         dAdXYZ(i,k,ip)=dAdXYZ(i,k,ip)-(xMass(ip)/TotMass)*sum
        end do
       end do
      end do
* Finally,transform dXdA, etc, into derivatives w.r.t atom coordinates:
*  dRVdXYZ(i,k,ip)= sum dRVdA(i,j)*dAdXYZ(j,k,ip)
*  d2RVdXYZ2(i,k1,ip1,k2,ip2)=
*    sum d2RVdA2(i,j1,j2)*dAdXYZ(j1,k1,ip1)*dAdXYZ(j2,k2,ip2)
*  and so on.. Note: A is linear in nuclear coordinates.
      If (nOrder.ge.1) Then
      do i=1,3
       do ip=1,nMass
        do k=1,3
         ipk=k+3*(ip-1)
         sum=0.0D0
         do j=1,3
          sum=sum+dRVdA(i,j)*dAdXYZ(j,k,ip)
         end do
         dRVdXYZ(i,ipk)=sum
        end do
       end do
      end do
      End If
      If (nOrder.ge.2) Then
      do i=1,3
       do j1=1,3
        do ip2=1,nMass
         do k2=1,3
          ipk2=k2+3*(ip2-1)
          sum=0.0D0
          do j2=1,3
           sum=sum+d2RVdA2(i,j1,j2)*dAdXYZ(j2,k2,ip2)
          end do
          tmp(j1,ipk2)=sum
         end do
        end do
       end do
       do ipk2=1,3*nMass
        do ip1=1,nMass
         do k1=1,3
          ipk1=k1+3*(ip1-1)
          sum=0.0D0
          do j1=1,3
           sum=sum+tmp(j1,ipk2)*dAdXYZ(j1,k1,ip1)
          end do
          d2RVdXYZ2(i,ipk1,ipk2)=sum
         end do
        end do
       end do
      end do
      End If
      Call mma_deallocate(Curr123)
      Call mma_deallocate(tmp)
      Call mma_deallocate(dAdXYZ)
      return
      end
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine rotder4(norder,S,X,dXdA,d2XdA2,d3XdA3,d4XdA4)
      implicit none
* Call arguments:
      integer norder
      real*8 S(3,3)
      real*8 X(3),dXdA(3,3),d2XdA2(3,3,3),d3XdA3(3,3,3,3)
      real*8 d4XdA4(3,3,3,3,3)
* Local variables:
      integer i,i1,i2,i3,ia,ib,ic,id,j,k,l,m,n
      real*8 c0,c1,c2,Q,cs,sn,sum,XN
      real*8 dc0dQ,d2c0dQ2,d3c0dQ3,d4c0dQ4
      real*8 dc1dQ,d2c1dQ2,d3c1dQ3,d4c1dQ4
      real*8 dc2dQ,d2c2dQ2,d3c2dQ3,d4c2dQ4
      real*8 dc0dX(3),d2c0dX2(3,3),d3c0dX3(3,3,3),d4c0dX4(3,3,3,3)
      real*8 dc1dX(3),d2c1dX2(3,3),d3c1dX3(3,3,3),d4c1dX4(3,3,3,3)
      real*8 dc2dX(3),d2c2dX2(3,3),d3c2dX3(3,3,3),d4c2dX4(3,3,3,3)
      real*8 U(3,3),dUdX(3,3,3),d2UdX2(3,3,3,3),d3UdX3(3,3,3,3,3)
      real*8 d4UdX4(3,3,3,3,3,3)
      real*8 S0(3,3)
      real*8 dPdX(3,3,3),d2PdX2(3,3,3,3),d3PdX3(3,3,3,3,3)
      real*8 d4PdX4(3,3,3,3,3,3)
      real*8 dAdX(3,3),d2AdX2(3,3,3),d3AdX3(3,3,3,3)
      real*8 d4AdX4(3,3,3,3,3)
c     real*8 T(3,3),det,detinv
      real*8 tmp1(3),tmp2(3,3),tmp3(3,3,3),tmp4(3,3,3,3)
      real*8 tmp5(3,3,3,3,3),tmp5A(3,3,3,3,3)
      real*8 umat(3,3),vmat(3,3),sval(3),wTmp(100)
*
      d2c0dQ2=0.0D0 ! Dummy initialize
      d2c1dQ2=0.0D0 ! Dummy initialize
      d2c2dQ2=0.0D0 ! Dummy initialize
      d3c0dQ3=0.0D0 ! Dummy initialize
      d3c1dQ3=0.0D0 ! Dummy initialize
      d3c2dQ3=0.0D0 ! Dummy initialize
      d4c0dQ4=0.0D0 ! Dummy initialize
      d4c1dQ4=0.0D0 ! Dummy initialize
      d4c2dQ4=0.0D0 ! Dummy initialize
*
      Q=X(1)**2+X(2)**2+X(3)**2
      XN=sqrt(Q)
*                                                                      *
************************************************************************
*                                                                      *
      if (Q.lt.0.5D0) then
* Use a Taylor expansion:
        c0=1.D0-(Q/2.D0)
     &    *(1.D0-(Q/12.D0)*(1.D0-(Q/30.D0)*(1.D0-(Q/56.D0))))
        c1=(1.D0-(Q/6.D0)
     &    *(1.D0-(Q/20.D0)*(1.D0-(Q/42.D0)*(1.D0-(Q/72.D0)))))
        c2=(1.D0-(Q/12.D0)
     &    *(1.D0-(Q/30.D0)*(1.D0-(Q/56.D0)*(1.D0-(Q/90.D0)))))/2.D0
        dc0dQ=-(1.D0-(Q/6.D0)
     &       *(1.D0-(Q/20.D0)*(1.D0-(Q/42.D0)*(1.D0-(Q/72.D0)))))/2.D0
        dc1dQ=-(1.D0-(Q/10.D0)*(1.D0-(Q/28.D0)*(1.D0-(Q/54.D0))))/6.D0
        dc2dQ=-(1.D0-(Q/30.D0)
     &       *(2.D0-(Q/56.D0)*(3.D0-(Q/90.D0)*4.D0)))/24.D0
       if (norder.le.1) goto 101
        d2c0dQ2=(1.D0-(Q/10.D0)*(1.D0-(Q/28.D0)*(1.D0-(Q/54.D0))))/12.D0
        d2c1dQ2=(1.D0-(Q/14.D0)*(1.D0-(Q/36.D0)*(1.D0-(Q/66.D0))))/60.D0
        d2c2dQ2= (1.D0-(Q/56.D0)
     &         *(3.D0-(Q/90.D0)*(6.D0-(Q/132.D0)*10.D0)))/360.D0
       if (norder.le.2) goto 101
        d3c0dQ3=-(1.D0-(Q/14.D0)
     &         *(1.D0-(Q/36.D0)*(1.D0-(Q/66.D0))))/120.D0
        d3c1dQ3=-(1.D0-(Q/18.D0)
     &         *(1.D0-(Q/44.D0)*(1.D0-(Q/78.D0))))/840.D0
        d3c2dQ3=-(1.D0-(Q/90.D0)*(4.D0-(Q/132.D0)*10.D0))/6720.D0
       if (norder.le.3) goto 101
        d4c0dQ4=(1.D0-(Q/18.D0)
     &         *(1.D0-(Q/44.D0)*(1.D0-(Q/78.D0))))/1680.D0
        d4c1dQ4=(1.D0-(Q/22.D0)
     &         *(1.D0-(Q/52.D0)*(1.D0-(Q/90.D0))))/15120.D0
        d4c2dQ4= (1.D0-(Q/132.D0)*(5.D0-(Q/182.D0)*15.D0))/151200.D0
 101  continue
      else
       cs=cos(XN)
       sn=sin(XN)
       c0=cs
       c1=sn/XN
       c2=(1.0D0-cs)/(XN**2)
       dc0dQ=-sn/(2.D0*XN)
       dc1dQ=-(sn-XN*cs)/(2.D0*XN**3)
       dc2dQ=(-2.D0+2.D0*cs+XN*sn)/(2.D0*XN**4)
       if (norder.le.1) goto 102
       d2c0dQ2=(sn-XN*cs)/(4.D0*XN**3)
       d2c1dQ2=-(3.D0*XN*cs-(3.D0-Q)*sn)/(4.D0*XN**5)
       d2c2dQ2=(8.D0-(8.D0-XN**2)*cs-5.D0*XN*sn)/(4.D0*XN**6)
       if (norder.le.2) goto 102
       d3c0dQ3=(3.D0*XN*cs-(3.D0-Q)*sn)/(8.D0*XN**5)
       d3c1dQ3=-((15.D0-6.D0*Q)*sn-XN*(15.D0-Q)*cs)/(8.D0*XN**7)
       d3c2dQ3=(-48.D0+(48.D0-9.D0*XN**2)*cs+XN*(33.D0-XN**2)*sn)
     &        /(8.D0*XN**8)
       if (norder.le.3) goto 102
       d4c0dQ4=((15.D0-6.D0*Q)*sn-XN*(15.D0-Q)*cs)/(16.D0*XN**7)
       d4c1dQ4=((105.D0-Q*(45.D0-Q))*sn-XN*(105.D0-10.D0*Q)*cs)
     &        /(16.D0*XN**9)
       d4c2dQ4=(384.D0-(384.D0-87.D0*XN**2+XN**4)*cs
     &        -XN*(279.D0-14.D0*XN**2)*sn)/(16.D0*XN**10)
 102  continue
      end if
      do i=1,3
       dc0dX(i)=2.D0*dc0dQ*X(i)
       dc1dX(i)=2.D0*dc1dQ*X(i)
       dc2dX(i)=2.D0*dc2dQ*X(i)
      end do
      if (norder.le.1) goto 103
      do j=1,3
       do i=1,3
        d2c0dX2(i,j)=4.D0*d2c0dQ2*X(i)*X(j)
        d2c1dX2(i,j)=4.D0*d2c1dQ2*X(i)*X(j)
        d2c2dX2(i,j)=4.D0*d2c2dQ2*X(i)*X(j)
       end do
      end do
      do i=1,3
       d2c0dX2(i,i)=d2c0dX2(i,i)+2.D0*dc0dQ
       d2c1dX2(i,i)=d2c1dX2(i,i)+2.D0*dc1dQ
       d2c2dX2(i,i)=d2c2dX2(i,i)+2.D0*dc2dQ
      end do
      if (norder.le.2) goto 103
      do k=1,3
       do j=1,3
        do i=1,3
         d3c0dX3(i,j,k)=8.D0*d3c0dQ3*X(i)*X(j)*X(k)
         d3c1dX3(i,j,k)=8.D0*d3c1dQ3*X(i)*X(j)*X(k)
         d3c2dX3(i,j,k)=8.D0*d3c2dQ3*X(i)*X(j)*X(k)
        end do
       end do
      end do
      do i=1,3
       do j=1,3
        d3c0dX3(i,i,j)=d3c0dX3(i,i,j)+4.D0*d2c0dQ2*X(j)
        d3c0dX3(i,j,i)=d3c0dX3(i,j,i)+4.D0*d2c0dQ2*X(j)
        d3c0dX3(j,i,i)=d3c0dX3(j,i,i)+4.D0*d2c0dQ2*X(j)
        d3c1dX3(i,i,j)=d3c1dX3(i,i,j)+4.D0*d2c1dQ2*X(j)
        d3c1dX3(i,j,i)=d3c1dX3(i,j,i)+4.D0*d2c1dQ2*X(j)
        d3c1dX3(j,i,i)=d3c1dX3(j,i,i)+4.D0*d2c1dQ2*X(j)
        d3c2dX3(i,i,j)=d3c2dX3(i,i,j)+4.D0*d2c2dQ2*X(j)
        d3c2dX3(i,j,i)=d3c2dX3(i,j,i)+4.D0*d2c2dQ2*X(j)
        d3c2dX3(j,i,i)=d3c2dX3(j,i,i)+4.D0*d2c2dQ2*X(j)
       end do
      end do
      if (norder.le.3) goto 103
      do l=1,3
       do k=1,3
        do j=1,3
         do i=1,3
          d4c0dX4(i,j,k,l)=16.D0*d4c0dQ4*X(i)*X(j)*X(k)*X(l)
          d4c1dX4(i,j,k,l)=16.D0*d4c1dQ4*X(i)*X(j)*X(k)*X(l)
          d4c2dX4(i,j,k,l)=16.D0*d4c2dQ4*X(i)*X(j)*X(k)*X(l)
         end do
        end do
       end do
      end do
      do k=1,3
       do j=1,3
        do i=1,3
         d4c0dX4(i,i,j,k)=d4c0dX4(i,i,j,k)+8.D0*d3c0dQ3*X(j)*X(k)
         d4c0dX4(i,j,i,k)=d4c0dX4(i,j,i,k)+8.D0*d3c0dQ3*X(j)*X(k)
         d4c0dX4(i,j,k,i)=d4c0dX4(i,j,k,i)+8.D0*d3c0dQ3*X(j)*X(k)
         d4c0dX4(j,i,i,k)=d4c0dX4(j,i,i,k)+8.D0*d3c0dQ3*X(j)*X(k)
         d4c0dX4(j,i,k,i)=d4c0dX4(j,i,k,i)+8.D0*d3c0dQ3*X(j)*X(k)
         d4c0dX4(j,k,i,i)=d4c0dX4(j,k,i,i)+8.D0*d3c0dQ3*X(j)*X(k)
         d4c1dX4(i,i,j,k)=d4c1dX4(i,i,j,k)+8.D0*d3c1dQ3*X(j)*X(k)
         d4c1dX4(i,j,i,k)=d4c1dX4(i,j,i,k)+8.D0*d3c1dQ3*X(j)*X(k)
         d4c1dX4(i,j,k,i)=d4c1dX4(i,j,k,i)+8.D0*d3c1dQ3*X(j)*X(k)
         d4c1dX4(j,i,i,k)=d4c1dX4(j,i,i,k)+8.D0*d3c1dQ3*X(j)*X(k)
         d4c1dX4(j,i,k,i)=d4c1dX4(j,i,k,i)+8.D0*d3c1dQ3*X(j)*X(k)
         d4c1dX4(j,k,i,i)=d4c1dX4(j,k,i,i)+8.D0*d3c1dQ3*X(j)*X(k)
         d4c2dX4(i,i,j,k)=d4c2dX4(i,i,j,k)+8.D0*d3c2dQ3*X(j)*X(k)
         d4c2dX4(i,j,i,k)=d4c2dX4(i,j,i,k)+8.D0*d3c2dQ3*X(j)*X(k)
         d4c2dX4(i,j,k,i)=d4c2dX4(i,j,k,i)+8.D0*d3c2dQ3*X(j)*X(k)
         d4c2dX4(j,i,i,k)=d4c2dX4(j,i,i,k)+8.D0*d3c2dQ3*X(j)*X(k)
         d4c2dX4(j,i,k,i)=d4c2dX4(j,i,k,i)+8.D0*d3c2dQ3*X(j)*X(k)
         d4c2dX4(j,k,i,i)=d4c2dX4(j,k,i,i)+8.D0*d3c2dQ3*X(j)*X(k)
        end do
       end do
      end do
      do j=1,3
       do i=1,3
        d4c0dX4(i,i,j,j)=d4c0dX4(i,i,j,j)+4.D0*d2c0dQ2
        d4c0dX4(i,j,i,j)=d4c0dX4(i,j,i,j)+4.D0*d2c0dQ2
        d4c0dX4(i,j,j,i)=d4c0dX4(i,j,j,i)+4.D0*d2c0dQ2
        d4c1dX4(i,i,j,j)=d4c1dX4(i,i,j,j)+4.D0*d2c1dQ2
        d4c1dX4(i,j,i,j)=d4c1dX4(i,j,i,j)+4.D0*d2c1dQ2
        d4c1dX4(i,j,j,i)=d4c1dX4(i,j,j,i)+4.D0*d2c1dQ2
        d4c2dX4(i,i,j,j)=d4c2dX4(i,i,j,j)+4.D0*d2c2dQ2
        d4c2dX4(i,j,i,j)=d4c2dX4(i,j,i,j)+4.D0*d2c2dQ2
        d4c2dX4(i,j,j,i)=d4c2dX4(i,j,j,i)+4.D0*d2c2dQ2
       end do
      end do
 103  continue
* The unitary matrix U=exp(XMat), and its derivatives.
* Here, X is the matrix with elements XMat(i,j)=eps(i,k,j)*X(k)
* First term, cos(X)*delta(i,j):
      do i=1,3
       do j=1,3
        U(i,j)=0.0d0
       end do
       U(i,i)=c0
      end do
      do k=1,3
      do i=1,3
       do j=1,3
        dUdX(i,j,k)=0.0d0
       end do
       dUdX(i,i,k)=dc0dX(k)
      end do
      end do
      if (norder.le.1) goto 201

      do k=1,3
      do l=1,3
      do i=1,3
       do j=1,3
        d2UdX2(i,j,k,l)=0.0d0
       end do
       d2UdX2(i,i,k,l)=d2c0dX2(k,l)
      end do
      end do
      end do
      if (norder.le.2) goto 201
      do k=1,3
      do l=1,3
      do m=1,3
      do i=1,3
       do j=1,3
        d3UdX3(i,j,k,l,m)=0.0d0
       end do
       d3UdX3(i,i,k,l,m)=d3c0dX3(k,l,m)
      end do
      end do
      end do
      end do
      if (norder.le.3) goto 201
      do k=1,3
      do l=1,3
      do m=1,3
      do n=1,3
      do i=1,3
       do j=1,3
        d4UdX4(i,j,k,l,m,n)=0.0d0
       end do
       d4UdX4(i,i,k,l,m,n)=d4c0dX4(k,l,m,n)
      end do
      end do
      end do
      end do
      end do
 201  continue
* Second term, (sin(X)/X)*eps(i,k,j)*X(k):
      do i1=1,3
       i2=1+mod(i1,3)
       i3=1+mod(i2,3)
       U(i1,i3)=U(i1,i3)+c1*X(i2)
       U(i3,i1)=U(i3,i1)-c1*X(i2)
       do k=1,3
        dUdX(i1,i3,k)=dUdX(i1,i3,k)+dc1dX(k)*X(i2)
        dUdX(i3,i1,k)=dUdX(i3,i1,k)-dc1dX(k)*X(i2)
       end do
       dUdX(i1,i3,i2)=dUdX(i1,i3,i2)+c1
       dUdX(i3,i1,i2)=dUdX(i3,i1,i2)-c1
       if (norder.le.1) goto 202
       do k=1,3
        do l=1,3
         d2UdX2(i1,i3,k,l)=d2UdX2(i1,i3,k,l)+d2c1dX2(k,l)*X(i2)
         d2UdX2(i3,i1,k,l)=d2UdX2(i3,i1,k,l)-d2c1dX2(k,l)*X(i2)
        end do
       end do
       do k=1,3
        d2UdX2(i1,i3,i2,k)=d2UdX2(i1,i3,i2,k)+dc1dX(k)
        d2UdX2(i3,i1,i2,k)=d2UdX2(i3,i1,i2,k)-dc1dX(k)
        d2UdX2(i1,i3,k,i2)=d2UdX2(i1,i3,k,i2)+dc1dX(k)
        d2UdX2(i3,i1,k,i2)=d2UdX2(i3,i1,k,i2)-dc1dX(k)
       end do
       if (norder.le.2) goto 202
       do k=1,3
        do l=1,3
         do m=1,3
          d3UdX3(i1,i3,k,l,m)=d3UdX3(i1,i3,k,l,m)+d3c1dX3(k,l,m)*X(i2)
          d3UdX3(i3,i1,k,l,m)=d3UdX3(i3,i1,k,l,m)-d3c1dX3(k,l,m)*X(i2)
         end do
        end do
       end do
       do k=1,3
        do l=1,3
         d3UdX3(i1,i3,i2,k,l)=d3UdX3(i1,i3,i2,k,l)+d2c1dX2(k,l)
         d3UdX3(i3,i1,i2,k,l)=d3UdX3(i3,i1,i2,k,l)-d2c1dX2(k,l)
         d3UdX3(i1,i3,k,i2,l)=d3UdX3(i1,i3,k,i2,l)+d2c1dX2(k,l)
         d3UdX3(i3,i1,k,i2,l)=d3UdX3(i3,i1,k,i2,l)-d2c1dX2(k,l)
         d3UdX3(i1,i3,k,l,i2)=d3UdX3(i1,i3,k,l,i2)+d2c1dX2(k,l)
         d3UdX3(i3,i1,k,l,i2)=d3UdX3(i3,i1,k,l,i2)-d2c1dX2(k,l)
        end do
       end do
       if (norder.le.3) goto 202
       do k=1,3
        do l=1,3
         do m=1,3
          do n=1,3
           d4UdX4(i1,i3,k,l,m,n)=
     &                   d4UdX4(i1,i3,k,l,m,n)+d4c1dX4(k,l,m,n)*X(i2)
           d4UdX4(i3,i1,k,l,m,n)=
     &                   d4UdX4(i3,i1,k,l,m,n)-d4c1dX4(k,l,m,n)*X(i2)
          end do
         end do
        end do
       end do
       do k=1,3
        do l=1,3
         do m=1,3
          d4UdX4(i1,i3,i2,k,l,m)=d4UdX4(i1,i3,i2,k,l,m)+d3c1dX3(k,l,m)
          d4UdX4(i3,i1,i2,k,l,m)=d4UdX4(i3,i1,i2,k,l,m)-d3c1dX3(k,l,m)
          d4UdX4(i1,i3,k,i2,l,m)=d4UdX4(i1,i3,k,i2,l,m)+d3c1dX3(k,l,m)
          d4UdX4(i3,i1,k,i2,l,m)=d4UdX4(i3,i1,k,i2,l,m)-d3c1dX3(k,l,m)
          d4UdX4(i1,i3,k,l,i2,m)=d4UdX4(i1,i3,k,l,i2,m)+d3c1dX3(k,l,m)
          d4UdX4(i3,i1,k,l,i2,m)=d4UdX4(i3,i1,k,l,i2,m)-d3c1dX3(k,l,m)
          d4UdX4(i1,i3,k,l,m,i2)=d4UdX4(i1,i3,k,l,m,i2)+d3c1dX3(k,l,m)
          d4UdX4(i3,i1,k,l,m,i2)=d4UdX4(i3,i1,k,l,m,i2)-d3c1dX3(k,l,m)
         end do
        end do
       end do
 202   continue
      end do
* Third term,((1-cos(X))/X**2)*X(i)*X(j):
      do i=1,3
       do j=1,3
        U(i,j)=U(i,j)+c2*X(i)*X(j)
        do k=1,3
         dUdX(i,j,k)=dUdX(i,j,k)+dc2dX(k)*X(i)*X(j)
        end do
        dUdX(i,j,i)=dUdX(i,j,i)+c2*X(j)
        dUdX(i,j,j)=dUdX(i,j,j)+c2*X(i)
       end do
      end do
      if (norder.le.1) goto 203
      do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3
          d2UdX2(i,j,k,l)=d2UdX2(i,j,k,l)+d2c2dX2(k,l)*X(i)*X(j)
         end do
         d2UdX2(i,j,i,k)=d2UdX2(i,j,i,k)+dc2dX(k)*X(j)
         d2UdX2(i,j,j,k)=d2UdX2(i,j,j,k)+dc2dX(k)*X(i)
         d2UdX2(i,j,k,i)=d2UdX2(i,j,k,i)+dc2dX(k)*X(j)
         d2UdX2(i,j,k,j)=d2UdX2(i,j,k,j)+dc2dX(k)*X(i)
        end do
        d2UdX2(i,j,i,j)=d2UdX2(i,j,i,j)+c2
        d2UdX2(i,j,j,i)=d2UdX2(i,j,j,i)+c2
       end do
      end do
      if (norder.le.2) goto 203
      do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3
          do m=1,3
           d3UdX3(i,j,k,l,m)=d3UdX3(i,j,k,l,m)+d3c2dX3(k,l,m)*X(i)*X(j)
          end do
          d3UdX3(i,j,i,k,l)=d3UdX3(i,j,i,k,l)+d2c2dX2(k,l)*X(j)
          d3UdX3(i,j,j,k,l)=d3UdX3(i,j,j,k,l)+d2c2dX2(k,l)*X(i)
          d3UdX3(i,j,k,i,l)=d3UdX3(i,j,k,i,l)+d2c2dX2(k,l)*X(j)
          d3UdX3(i,j,k,j,l)=d3UdX3(i,j,k,j,l)+d2c2dX2(k,l)*X(i)
          d3UdX3(i,j,k,l,i)=d3UdX3(i,j,k,l,i)+d2c2dX2(k,l)*X(j)
          d3UdX3(i,j,k,l,j)=d3UdX3(i,j,k,l,j)+d2c2dX2(k,l)*X(i)
         end do
         d3UdX3(i,j,i,j,k)=d3UdX3(i,j,i,j,k)+dc2dX(k)
         d3UdX3(i,j,j,i,k)=d3UdX3(i,j,j,i,k)+dc2dX(k)
         d3UdX3(i,j,i,k,j)=d3UdX3(i,j,i,k,j)+dc2dX(k)
         d3UdX3(i,j,j,k,i)=d3UdX3(i,j,j,k,i)+dc2dX(k)
         d3UdX3(i,j,k,i,j)=d3UdX3(i,j,k,i,j)+dc2dX(k)
         d3UdX3(i,j,k,j,i)=d3UdX3(i,j,k,j,i)+dc2dX(k)
        end do
       end do
      end do
      if (norder.le.3) goto 203
      do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3
          do m=1,3
           do n=1,3
            d4UdX4(i,j,k,l,m,n)=d4UdX4(i,j,k,l,m,n)+
     &                          d4c2dX4(k,l,m,n)*X(i)*X(j)
           end do
           d4UdX4(i,j,i,k,l,m)=d4UdX4(i,j,i,k,l,m)+d3c2dX3(k,l,m)*X(j)
           d4UdX4(i,j,j,k,l,m)=d4UdX4(i,j,j,k,l,m)+d3c2dX3(k,l,m)*X(i)
           d4UdX4(i,j,k,i,l,m)=d4UdX4(i,j,k,i,l,m)+d3c2dX3(k,l,m)*X(j)
           d4UdX4(i,j,k,j,l,m)=d4UdX4(i,j,k,j,l,m)+d3c2dX3(k,l,m)*X(i)
           d4UdX4(i,j,k,l,i,m)=d4UdX4(i,j,k,l,i,m)+d3c2dX3(k,l,m)*X(j)
           d4UdX4(i,j,k,l,j,m)=d4UdX4(i,j,k,l,j,m)+d3c2dX3(k,l,m)*X(i)
           d4UdX4(i,j,k,l,m,i)=d4UdX4(i,j,k,l,m,i)+d3c2dX3(k,l,m)*X(j)
           d4UdX4(i,j,k,l,m,j)=d4UdX4(i,j,k,l,m,j)+d3c2dX3(k,l,m)*X(i)
          end do
          d4UdX4(i,j,i,j,k,l)=d4UdX4(i,j,i,j,k,l)+d2c2dX2(k,l)
          d4UdX4(i,j,i,k,j,l)=d4UdX4(i,j,i,k,j,l)+d2c2dX2(k,l)
          d4UdX4(i,j,i,k,l,j)=d4UdX4(i,j,i,k,l,j)+d2c2dX2(k,l)
          d4UdX4(i,j,k,i,j,l)=d4UdX4(i,j,k,i,j,l)+d2c2dX2(k,l)
          d4UdX4(i,j,k,i,l,j)=d4UdX4(i,j,k,i,l,j)+d2c2dX2(k,l)
          d4UdX4(i,j,k,l,i,j)=d4UdX4(i,j,k,l,i,j)+d2c2dX2(k,l)
          d4UdX4(i,j,j,i,k,l)=d4UdX4(i,j,j,i,k,l)+d2c2dX2(k,l)
          d4UdX4(i,j,j,k,i,l)=d4UdX4(i,j,j,k,i,l)+d2c2dX2(k,l)
          d4UdX4(i,j,j,k,l,i)=d4UdX4(i,j,j,k,l,i)+d2c2dX2(k,l)
          d4UdX4(i,j,k,j,i,l)=d4UdX4(i,j,k,j,i,l)+d2c2dX2(k,l)
          d4UdX4(i,j,k,j,l,i)=d4UdX4(i,j,k,j,l,i)+d2c2dX2(k,l)
          d4UdX4(i,j,k,l,j,i)=d4UdX4(i,j,k,l,j,i)+d2c2dX2(k,l)
         end do
        end do
       end do
      end do
 203  continue
* The matrix S (which should be symmetrical) is the product
* S=S0*U=S0*exp(X) with X given as input arguments to this Subroutine.
* Need to compute S0=S*U(transpose):
      do i=1,3
       do j=1,3
        sum=0.0D0
        do k=1,3
         sum=sum+S(i,k)*U(j,k)
        end do
        S0(i,j)=sum
       end do
      end do
* The matrix product P=S0*U, where now U is used as an expansion in
* the variables X, and its derivatives:
      do i=1,3
       do j=1,3
*       sum=0.0D0
*       do i1=1,3
*        sum=sum+S0(i,i1)*U(i1,j)
*       end do
*       P(i,j)=sum
        do k=1,3
         sum=0.0D0
         do i1=1,3
          sum=sum+S0(i,i1)*dUdX(i1,j,k)
         end do
         dPdX(i,j,k)=sum
        end do
       end do
      end do
      if (norder.le.1) goto 301
      do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3
          sum=0.0D0
          do i1=1,3
           sum=sum+S0(i,i1)*d2UdX2(i1,j,k,l)
          end do
          d2PdX2(i,j,k,l)=sum
         end do
        end do
       end do
      end do
      if (norder.le.2) goto 301
      do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3
          do m=1,3
           sum=0.0D0
           do i1=1,3
            sum=sum+S0(i,i1)*d3UdX3(i1,j,k,l,m)
           end do
           d3PdX3(i,j,k,l,m)=sum
          end do
         end do
        end do
       end do
      end do
      if (norder.le.3) goto 301
      do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3
          do m=1,3
           do n=1,3
            sum=0.0D0
            do i1=1,3
             sum=sum+S0(i,i1)*d4UdX4(i1,j,k,l,m,n)
            end do
            d4PdX4(i,j,k,l,m,n)=sum
           end do
          end do
         end do
        end do
       end do
      end do
 301  continue
* The vector A is the dual of the antisymmetric part of P:
*       A(1)=P(3,2)-P(2,3)
*       A(2)=P(1,3)-P(3,1)
*       A(3)=P(2,1)-P(1,2)
        do i=1,3
         dAdX(1,i)=dPdX(3,2,i)-dPdX(2,3,i)
         dAdX(2,i)=dPdX(1,3,i)-dPdX(3,1,i)
         dAdX(3,i)=dPdX(2,1,i)-dPdX(1,2,i)
         if (norder.le.1) goto 304
         do j=1,3
          d2AdX2(1,i,j)=d2PdX2(3,2,i,j)-d2PdX2(2,3,i,j)
          d2AdX2(2,i,j)=d2PdX2(1,3,i,j)-d2PdX2(3,1,i,j)
          d2AdX2(3,i,j)=d2PdX2(2,1,i,j)-d2PdX2(1,2,i,j)
          if (norder.le.2) goto 303
          do k=1,3
           d3AdX3(1,i,j,k)=d3PdX3(3,2,i,j,k)-d3PdX3(2,3,i,j,k)
           d3AdX3(2,i,j,k)=d3PdX3(1,3,i,j,k)-d3PdX3(3,1,i,j,k)
           d3AdX3(3,i,j,k)=d3PdX3(2,1,i,j,k)-d3PdX3(1,2,i,j,k)
           if (norder.le.3) goto 302
           do l=1,3
            d4AdX4(1,i,j,k,l)=d4PdX4(3,2,i,j,k,l)-d4PdX4(2,3,i,j,k,l)
            d4AdX4(2,i,j,k,l)=d4PdX4(1,3,i,j,k,l)-d4PdX4(3,1,i,j,k,l)
            d4AdX4(3,i,j,k,l)=d4PdX4(2,1,i,j,k,l)-d4PdX4(1,2,i,j,k,l)
           end do
 302       continue
          end do
 303      continue
         end do
 304     continue
        end do
* Finally, we have obtained derivatives of A w.r.t X.
* So now we can obtain derivatives of the inverse mapping.
* This will require the inverse T of the matrix dAdX:
c     T(1,1)=dAdX(2,2)*dAdX(3,3)-dAdX(3,2)*dAdX(2,3)
c     T(2,1)=dAdX(3,2)*dAdX(1,3)-dAdX(1,2)*dAdX(3,3)
c     T(3,1)=dAdX(1,2)*dAdX(2,3)-dAdX(2,2)*dAdX(1,3)
c     T(1,2)=dAdX(2,3)*dAdX(3,1)-dAdX(3,3)*dAdX(2,1)
c     T(2,2)=dAdX(3,3)*dAdX(1,1)-dAdX(1,3)*dAdX(3,1)
c     T(3,2)=dAdX(1,3)*dAdX(2,1)-dAdX(2,3)*dAdX(1,1)
c     T(1,3)=dAdX(2,1)*dAdX(3,2)-dAdX(3,1)*dAdX(2,2)
c     T(2,3)=dAdX(3,1)*dAdX(1,2)-dAdX(1,1)*dAdX(3,2)
c     T(3,3)=dAdX(1,1)*dAdX(2,2)-dAdX(2,1)*dAdX(1,2)
c     det=dAdX(1,1)*T(1,1)+dAdX(2,1)*T(2,1)+dAdX(3,1)*T(3,1)
c     detInv=1.0D0/det
*  First derivatives dXdA(j,i):
c     do i=1,3
c      do ia=1,3
c       dXdA(ia,i)=DetInv*T(i,ia)
c      end do
c     end do
* Use the Moore-Penrose pseudoinverse instead
      call dgesvd_('A','A',3,3,dAdX,3,sval,umat,3,vmat,3,wTmp,100,i)
      Do i=1,3
        If (Abs(sval(i)).gt.1.0d-12) Then
          call dscal_(3,1.0d0/sval(i),umat(1,i),1)
        Else
          call dcopy_(3,[0.0d0],0,umat(1,i),1)
        End If
      End Do
      Call dgemm_('T','T',3,3,3,1.0d0,vmat,3,umat,3,0.0d0,dXdA,3)
*  Second derivatives d2XdA(ic,j,k)
      if (norder.le.1) goto 401
      do i=1,3
       do k=1,3
        do ia=1,3
         sum=0.0D0
         do ib=1,3
          sum=sum+d2AdX2(i,ia,ib)*dXdA(ib,k)
         end do
         tmp1(ia)=sum
        end do
        do j=1,3
         sum=0.0D0
         do ia=1,3
          sum=sum+tmp1(ia)*dXdA(ia,j)
         end do
         tmp3(i,j,k)=sum
        end do
       end do
      end do
      do ic=1,3
       do j=1,3
        do k=1,3
         sum=0.0D0
         do i=1,3
          sum=sum+dXdA(ic,i)*tmp3(i,j,k)
         end do
         d2XdA2(ic,j,k)=-sum
        end do
       end do
      end do
*  Third derivatives d3XdA3(id,j,k,l)
      if (norder.le.2) goto 401
      do i=1,3
       do l=1,3
        do ia=1,3
         do ib=1,3
          sum=0.0D0
          do ic=1,3
           sum=sum+d3AdX3(i,ia,ib,ic)*dXdA(ic,l)
          end do
          tmp1(ib)=sum
         end do
         do k=1,3
          sum=0.0D0
          do ib=1,3
           sum=sum+tmp1(ib)*dXdA(ib,k)
          end do
          tmp2(ia,k)=sum
         end do
        end do
        do j=1,3
         do k=1,3
          sum=0.0D0
          do ia=1,3
           sum=sum+tmp2(ia,k)*dXdA(ia,j)
          end do
          tmp4(i,j,k,l)=sum
         end do
        end do
       end do
      end do
      do id=1,3
       do j=1,3
        do k=1,3
         do l=1,3
          sum=0.0D0
          do i=1,3
           sum=sum+dXdA(id,i)*tmp4(i,j,k,l)
          end do
          d3XdA3(id,j,k,l)=-sum
         end do
        end do
       end do
      end do
      do i=1,3
       do j=1,3
        do ib=1,3
         sum=0.0D0
         do ia=1,3
          sum=sum+d2AdX2(i,ia,ib)*dXdA(ia,j)
         end do
         tmp1(ib)=sum
        end do
        do k=1,3
         do l=1,3
          sum=0.0D0
          do ib=1,3
           sum=sum+tmp1(ib)*d2XdA2(ib,k,l)
          end do
          tmp4(i,j,k,l)=sum
         end do
        end do
       end do
      end do
      do ic=1,3
       do j=1,3
        do k=1,3
         do l=1,3
          sum=d3XdA3(ic,j,k,l)
          do i=1,3
           sum=sum-DXDA(ic,i)*
     &             (tmp4(i,j,k,l)+tmp4(i,k,l,j)+tmp4(i,l,j,k))
          end do
          d3XdA3(ic,j,k,l)=sum
         end do
        end do
       end do
      end do
*  Fourth derivatives d4XdA4(id,j,k,l,m)
      if (norder.le.3) goto 401
      do i=1,3
       do ia=1,3
        do ib=1,3
         do ic=1,3
          do m=1,3
           sum=0.0d0
           do id=1,3
            sum=sum+d4AdX4(i,ia,ib,ic,id)*dXdA(id,m)
           end do
           tmp2(ic,m)=sum
          end do
         end do
         do m=1,3
          do l=1,3
           sum=0.0d0
           do ic=1,3
            sum=sum+tmp2(ic,m)*dXdA(ic,l)
           end do
           tmp3(ib,l,m)=sum
          end do
         end do
        end do
        do m=1,3
         do l=1,3
          do k=1,3
           sum=0.0d0
           do ib=1,3
            sum=sum+tmp3(ib,l,m)*dXdA(ib,k)
           end do
           tmp4(ia,k,l,m)=sum
          end do
         end do
        end do
       end do
       do m=1,3
        do l=1,3
         do k=1,3
          do j=1,3
           sum=0.0d0
           do ia=1,3
            sum=sum+tmp4(ia,k,l,m)*dXdA(ia,j)
           end do
           tmp5(i,j,k,l,m)=sum
          end do
         end do
        end do
       end do
       do l=1,3
        do m=1,3
         do ia=1,3
          do ib=1,3
           sum=0.0d0
           do ic=1,3
            sum=sum+d3AdX3(i,ia,ib,ic)*d2XdA2(ic,l,m)
           end do
           tmp1(ib)=sum
          end do
          do k=1,3
           sum=0.0d0
           do ib=1,3
            sum=sum+tmp1(ib)*dXdA(ib,k)
           end do
           tmp2(ia,k)=sum
          end do
         end do
         do k=1,3
          do j=1,3
           sum=0.0d0
           do ia=1,3
            sum=sum+tmp2(ia,k)*dXdA(ia,j)
           end do
           tmp5A(i,j,k,l,m)=sum
          end do
         end do
        end do
       end do
       do j=1,3
        do k=1,3
         do l=1,3
          do m=1,3
           tmp5(i,j,k,l,m)=tmp5(i,j,k,l,m)+tmp5A(i,j,k,l,m)+
     &          tmp5A(i,j,l,k,m)+tmp5A(i,k,l,j,m)+tmp5A(i,j,m,k,l)+
     &          tmp5A(i,k,m,j,l)+tmp5A(i,l,m,j,k)
          end do
         end do
        end do
       end do
       do j=1,3
        do ib=1,3
         sum=0.0d0
         do ia=1,3
          sum=sum+d2AdX2(i,ia,ib)*dXdA(ia,j)
         end do
         tmp2(j,ib)=sum
        end do
        do k=1,3
         do l=1,3
          do m=1,3
           sum=0.0d0
           do ib=1,3
            sum=sum+tmp2(j,ib)*d3XdA3(ib,k,l,m)
           end do
           tmp5A(i,j,k,l,m)=sum
          end do
         end do
        end do
       end do
       do j=1,3
        do m=1,3
         do l=1,3
          do k=1,3
       tmp5(i,j,k,l,m)=tmp5(i,j,k,l,m)+tmp5A(i,j,k,l,m)+
     &       tmp5A(i,k,l,m,j)+tmp5A(i,l,m,j,k)+tmp5A(i,m,j,k,l)
          end do
         end do
        end do
       end do
       do j=1,3
        do m=1,3
         do ib=1,3
          sum=0.0D0
          do ia=1,3
           sum=sum+d2AdX2(i,ia,ib)*d2XdA2(ia,j,m)
          end do
          tmp1(ib)=sum
         end do
         do l=1,3
          do k=1,3
           sum=0.0D0
           do ib=1,3
            sum=sum+tmp1(ib)*d2XdA2(ib,k,l)
           end do
           tmp5A(i,j,k,l,m)=sum
          end do
         end do
        end do
       end do
       do m=1,3
        do j=1,3
         do l=1,3
          do k=1,3
           tmp5(i,j,k,l,m)=tmp5(i,j,k,l,m)+tmp5A(i,j,k,l,m)+
     &                  tmp5A(i,j,k,m,l)+tmp5A(i,j,m,l,k)
          end do
         end do
        end do
       end do
      end do
      do id=1,3
       do j=1,3
        do k=1,3
         do l=1,3
          do m=1,3
           sum=0.0D0
           do i=1,3
            sum=sum+dXdA(id,i)*tmp5(i,j,k,l,m)
           end do
           d4XdA4(id,j,k,l,m)=-sum
          end do
         end do
        end do
       end do
      end do
 401  continue
      return
      end
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine mat2vec(RotMat,RotVec,RotAng)
      implicit real*8 (a-h,o-z)
      parameter (Pi=3.14159265358979323846D0)
      dimension RotMat(3,3),RotVec(3)
      trace=RotMat(1,1)+RotMat(2,2)+RotMat(3,3)
      x1=0.5D0*(RotMat(3,2)-RotMat(2,3))
      x2=0.5D0*(RotMat(1,3)-RotMat(3,1))
      x3=0.5D0*(RotMat(2,1)-RotMat(1,2))
      sr=sqrt(x1**2+x2**2+x3**2)
      cr=0.5D0*(trace-1.0d0)
      if(0.05D0*cr.gt.sr) then
        tn=sr/cr
        tn2=tn**2
* cos(theta) is positive, C is theta/sin(theta)
        C=(45045.0D0-tn2*(15015.0D0-tn2*(9009.0D0-tn2*(6435.0D0-
     &   tn2*(5005.0D0-tn2*(4095.0D0-tn2*3465.0D0))))))/(45045.0D0*cr)
        RotVec(1)=x1*C
        RotVec(2)=x2*C
        RotVec(3)=x3*C
        RotAng   =sr*C
      else if(-0.05D0*abs(cr).gt.sr.and.sr.gt.0.0D0) then
        tn=-sr/cr
        tn2=tn**2
* cos(theta) is negative, sin(theta) is positive (always),
* Pimth is (Pi-theta)
        Pimth=tn*(45045.0D0-tn2*(15015.0D0-tn2*(9009.0D0-tn2*(6435.0D0-
     &   tn2*(5005.0D0-tn2*(4095.0D0-tn2*3465.0D0))))))/45045.0D0
        RotAng   =Pi-Pimth
        RotVec(1)=(x1/sr)*RotAng
        RotVec(2)=(x2/sr)*RotAng
        RotVec(3)=(x3/sr)*RotAng
      else if (sr.ne.0.0D0) then
        RotAng=atan2(sr,cr)
        RotVec(1)=RotAng*(x1/sr)
        RotVec(2)=RotAng*(x2/sr)
        RotVec(3)=RotAng*(x3/sr)
      else
        RotVec(1)=0.0D0
        RotVec(2)=0.0D0
        RotVec(3)=0.0D0
        RotAng   =0.0D0
      end if
      return
      end
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine mkRotMat(rotvec,rotmat)
      implicit real*8 (a-h,o-z)
      real*8 RotMat(3,3)
      real*8 RotVec(3)
*
C     Call RecPrt('mkRotMat: RotVec',' ',RotVec,1,3)
*
      Q=RotVec(1)**2+RotVec(2)**2+RotVec(3)**2
      XN=sqrt(Q)
      if (Q.lt.0.01D0) then
        c0=1.d0-(Q/2.d0)
     &    *(1.d0-(Q/12.d0)*(1.d0-(Q/30.d0)*(1.d0-(Q/56.d0))))
        c1=1.d0-(Q/6.d0)
     &    *(1.d0-(Q/20.d0)*(1.d0-(Q/42.d0)*(1.d0-(Q/72.d0))))
        c2=(1.d0-(Q/12.d0)
     &    *(1.d0-(Q/30.d0)*(1.d0-(Q/56.d0)*(1.d0-(Q/90.d0)))))/2.d0
      else
       cs=cos(XN)
       sn=sin(XN)
       c0=cs
       c1=sn/XN
       c2=(1.0D0-cs)/(XN**2)
      end if
      RotMat(1,1)=c0
      RotMat(2,2)=c0
      RotMat(3,3)=c0
      RotMat(3,2)= c1*RotVec(1)
      RotMat(1,3)= c1*RotVec(2)
      RotMat(2,1)= c1*RotVec(3)
      RotMat(2,3)=-c1*RotVec(1)
      RotMat(3,1)=-c1*RotVec(2)
      RotMat(1,2)=-c1*RotVec(3)
      do i=1,3
       do j=1,3
        RotMat(i,j)=RotMat(i,j)+c2*RotVec(i)*RotVec(j)
       end do
      end do
      do i=1,3
       do j=1,3
        sum=0.0d0
        if (i.eq.j) sum=-1.0D0
        do k=1,3
         sum=sum+RotMat(i,k)*RotMat(j,k)
        end do
        if(abs(sum).gt.1.0D-10) then
         Call WarningMessage(2,'Error in RotDer')
         write(6,*)' MKROTMAT: ON check sum error=',sum
         Call Abend()
        end if
       end do
      end do
      return
      end
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine updRotMat(SmallRot,RotMat)
      implicit real*8 (a-h,o-z)
      dimension SmallRot(3),RotMat(3,3)
      dimension tmp(3,3),SmallMat(3,3)
      Call mkRotMat(SmallRot,SmallMat)
      do i=1,3
       do j=1,3
        sum=0.0D0
        do k=1,3
         sum=sum+RotMat(i,k)*SmallMat(k,j)
        end do
        tmp(i,j)=sum
       end do
      end do
      do i=1,3
       do j=1,3
        RotMat(i,j)=tmp(i,j)
       end do
      end do
* Check for orthonormality:
      do i=1,3
       do j=1,3
        sum=0.0d0
        if (i.eq.j) sum=-1.0D0
        do k=1,3
         sum=sum+RotMat(i,k)*RotMat(j,k)
        end do
        if(abs(sum).gt.1.0D-10) then
         write(6,*)' UPDROTMAT ON check sum error:',sum
        end if
       end do
      end do
      return
      end
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine moveit(nmass,currxyz,ref123,trans,rotvec)
      implicit real*8 (a-h,o-z)
      real*8 U(3,3)
      real*8 CurrXYZ(3,nMass)
      real*8 Ref123(3,nMass)
      real*8 RotVec(3)
      real*8 Trans(3)
      Call mkRotMat(RotVec,U)
      do imass=1,nmass
       do i=1,3
        sum=Trans(i)
        do j=1,3
         sum=sum+U(i,j)*Ref123(j,imass)
        end do
        CurrXYZ(i,imass)=sum
       end do
      end do
      return
      end
