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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************
      Subroutine SumVeeNew(SV,A,GD,I1,I2,G,V1,V2,Update)
      use stdalloc, only : mma_allocate, mma_deallocate
      use rasscf_global, only: lRoots, NAC
      Implicit None

#include "warnings.h"

      Real*8 SV,A,V1,V2
      INTEGER I1,I2
      Real*8,DIMENSION(LRoots*(LRoots+1)/2,NAC,NAC)::GD
      Real*8,DIMENSION(NAC,NAC,NAC,NAC)::G
      Real*8,DIMENSION(:,:),Allocatable::D11,D22
      Real*8,DIMENSION(:,:,:),Allocatable::D1J,D2J
      Logical Update

      INTEGER t,u,v,x,i11,i22,i12
      INTEGER J,I1J,I2J
      IF(Update) THEN
       CALL mma_allocate(D1J,lRoots,NAC,NAC)
       CALL mma_allocate(D2J,lRoots,NAC,NAC)
!      calculating
       DO J=1,I2-1                           !(J<I2<I1)
        I1J=(I1-1)*I1/2+J
        I2J=(I2-1)*I2/2+J
        Do t=1,NAC
        Do u=1,NAC
         D2J(J,t,u)= cos(A)*GD(I2J,t,u)+sin(A)*GD(I1J,t,u)
         D1J(J,t,u)=-sin(A)*GD(I2J,t,u)+cos(A)*GD(I1J,t,u)
        End Do
        End Do
       END DO
       J=I2                                  !(J=I2<I1)
      i11=(I1+1)*I1/2
      i22=(I2+1)*I2/2
      i12=(I1-1)*I1/2+I2
        Do t=1,NAC
        Do u=1,NAC
         D2J(i2,t,u)=GD(i11,t,u)*sin(A)**2+GD(i22,t,u)*cos(A)**2        &
     &+cos(A)*sin(A)*(GD(i12,u,t)+GD(i12,t,u))
         D1J(i2,t,u)=cos(A)*sin(A)*(GD(i11,t,u)-GD(i22,t,u))            &
     &+GD(i12,t,u)*cos(A)**2-GD(i12,u,t)*sin(A)**2
        End Do
        End Do
       DO J=I2+1,I1-1                        !(I2<J<I1)
        I1J=(I1-1)*I1/2+J
        I2J=(J-1)*J/2+I2
        Do t=1,NAC
        Do u=1,NAC
         D2J(J,t,u)= cos(A)*GD(I2J,u,t)+sin(A)*GD(I1J,t,u)
         D1J(J,t,u)=-sin(A)*GD(I2J,u,t)+cos(A)*GD(I1J,t,u)
        End Do
        End Do
       END DO
       J=I1                                  !(I2<J=I1)
      i11=(I1+1)*I1/2
      i22=(I2+1)*I2/2
      i12=(I1-1)*I1/2+I2
        Do t=1,NAC
        Do u=1,NAC
         D1J(i1,t,u)=GD(i11,t,u)*cos(A)**2+GD(i22,t,u)*sin(A)**2        &
     &   -cos(A)*sin(A)*(GD(i12,t,u)+GD(i12,u,t))
        End Do
        End Do

       DO J=I1+1,lRoots                      !(I2<I1<J)
        I1J=(J-1)*J/2+I1
        I2J=(J-1)*J/2+I2
        Do t=1,NAC
        Do u=1,NAC
         D2J(J,t,u)= cos(A)*GD(I2J,u,t)+sin(A)*GD(I1J,u,t)
         D1J(J,t,u)=-sin(A)*GD(I2J,u,t)+cos(A)*GD(I1J,u,t)
        End Do
        End Do
       END DO
!      updating
       DO J=1,I2-1                           !(J<I2<I1)
        I1J=(I1-1)*I1/2+J
        I2J=(I2-1)*I2/2+J
        Do t=1,NAC
        Do u=1,NAC
         GD(I2J,t,u)=D2J(J,t,u)
         GD(I1J,t,u)=D1J(J,t,u)
        End Do
        End Do
       END DO
       J=I2                                  !(J=I2<I1)
!      i11=(I1+1)*I1/2
      i22=(I2+1)*I2/2
      i12=(I1-1)*I1/2+I2
        Do t=1,NAC
        Do u=1,NAC
         GD(I22,t,u)=D2J(I2,t,u)
         GD(I12,t,u)=D1J(I2,t,u)
        End Do
        End Do
       DO J=I2+1,I1-1                        !(I2<J<I1)
        I1J=(I1-1)*I1/2+J
        I2J=(J-1)*J/2+I2
        Do t=1,NAC
        Do u=1,NAC
         GD(I2J,t,u)=D2J(J,u,t)
         GD(I1J,t,u)=D1J(J,t,u)
        End Do
        End Do
       END DO
       J=I1                                  !(I2<J=I1)
       i11=(I1+1)*I1/2
        Do t=1,NAC
        Do u=1,NAC
         GD(i11,t,u)=D1J(I1,t,u)
        End Do
        End Do

       DO J=I1+1,lRoots                      !(I2<I1<J)
        I1J=(J-1)*J/2+I1
        I2J=(J-1)*J/2+I2
        Do t=1,NAC
        Do u=1,NAC
         GD(I2J,t,u)=D2J(J,u,t)
         GD(I1J,t,u)=D1J(J,u,t)
        End Do
        End Do
       END DO
       Call mma_deallocate(D1J)
       Call mma_deallocate(D2J)
      ELSE
       CALL mma_allocate(D11,NAC,NAC)
       CALL mma_allocate(D22,NAC,NAC)
       i11=(I1+1)*I1/2
       i22=(I2+1)*I2/2
       i12=(I1-1)*I1/2+I2
       V1=0.0d0
       V2=V1
       DO t=1,NAC
        Do u=1,NAC
         D11(t,u)=GD(i11,t,u)*cos(A)**2+GD(i22,t,u)*sin(A)**2           &
     &   -cos(A)*sin(A)*(GD(i12,t,u)+GD(i12,u,t))
         D22(t,u)=GD(i11,t,u)*sin(A)**2+GD(i22,t,u)*cos(A)**2           &
     &   +cos(A)*sin(A)*(GD(i12,u,t)+GD(i12,t,u))
        End Do
       END DO
       DO t=1,NAC
       DO u=1,NAC
        Do v=1,NAC
        Do x=1,NAC
         V1=V1+D11(t,u)*D11(v,x)*G(t,u,v,x)
         V2=V2+D22(t,u)*D22(v,x)*G(t,u,v,x)
        End Do
        End Do
       END DO
       END DO
       V1=V1/2.0d0
       V2=V2/2.0d0
       SV=V1+V2
       Call mma_deallocate(D11)
       Call mma_deallocate(D22)
      END IF
      End Subroutine SumVeeNew
