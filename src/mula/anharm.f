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
! Copyright (C) 1996, Niclas Forsberg                                  *
!***********************************************************************
!!-----------------------------------------------------------------------!
!!
      Subroutine Anharm(eigenVec,harmfreq,D3,D4,Gprime,                 &
     &       Gdbleprime,x,max_term,nOsc,C,Temp,V3,T3,V4,T4)
!!
!!  Purpose:
!!    Calculate the anharmonicity constants.
!!    This routine assumes that the curvilinear coordinates are such
!!    that they coincide with dimensionless normal coordinates for
!!    small displacements.
!!
!!  Input:
!!    eigenVec   : Real*8 two dimensional array - eigenvectors.
!!    harmfreq   : Real*8 array - harmonic frequencies.
!!    D3         : Real*8 three dimensional array - cubic force
!!                 constants.
!!    D4         : Real*8 four dimensional array - quartic force
!!                 constants.
!!    Gprime     : Real*8 three dimensional array - first
!!                 derivatives of the inverse mass tensor.
!!    Gdbleprime : Real*8 four dimensional array - second
!!                 derivatives of the inverse mass tensor.
!!    max_term   : Integer - highest power of term in polynomial fit.
!!
!!  Output:
!!    x          : Real*8 two dimensional array - anharmonicity
!!                 constants.
!!
!!  Written by:
!!    Niclas Forsberg,
!!    Dept. of Theoretical Chemistry, Lund University, 1996.
!!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
#include "dims.fh"
      Real*8 eigenVec (nosc,nosc)
      Real*8 harmfreq (nosc)
      Real*8 D3(ngdim,ngdim,ngdim)
      Real*8 D4 (ngdim,ngdim,ngdim,ngdim)
      Real*8 Gprime (ngdim,ngdim,ngdim)
      Real*8 Gdbleprime(ngdim,ngdim,ngdim,ngdim)
      Real*8 x (nosc,nosc)
      Real*8  C(nOsc,nOsc)
      Real*8   V3(nOsc,nOsc,nOsc)
      Real*8   V4(nOsc,nOsc,nOsc,nOsc)
      Real*8   T3(nOsc,nOsc,nOsc)
      Real*8   T4(nOsc,nOsc,nOsc,nOsc)
      Real*8   Temp(nOsc,nOsc)

!!
      NumInt = nOsc
!!
!!---- Calculate the eigenvector matrix, C, in dimensionless normal coordinates.
      call dcopy_(NumInt**2,[0.0d0],0,Temp,1)
      Do i = 1,NumInt
      Temp(i,i) = 1.0d0/sqrt(harmfreq(i))
      End Do
      Call DGEMM_('n','n',                                              &
     &            NumInt,NumInt,NumInt,                                 &
     &            1.0d0,eigenVec,NumInt,                                &
     &            Temp,NumInt,                                          &
     &            0.0d0,C,NumInt)
!!
!!---- Transform cubic force constants to dimensionless normal coordinates
      Do i = 1,NumInt
      Do j = 1,NumInt
      Do k = 1,NumInt
      sum1 = 0.0d0
      sum2 = 0.0d0
      Do i1 = 1,NumInt
      Do j1 = 1,NumInt
      Do k1 = 1,NumInt
      coef = C(i1,i)*C(j1,j)*C(k1,k)
      sum1 = sum1+coef*D3(i1,j1,k1)
      sum2 = sum2+coef*Gprime(i1,j1,k1)
      End Do
      End Do
      End Do
      V3(i,j,k) = sum1
      T3(i,j,k) = sum2
      End Do
      End Do
      End Do
!!
!!---- Transform quartic force constants to dimensionless normal coordinates
      Do i = 1,NumInt
      Do j = 1,NumInt
      Do k = 1,NumInt
      Do l = 1,NumInt
      sum1 = 0.0d0
      sum2 = 0.0d0
      Do i1 = 1,NumInt
      Do j1 = 1,NumInt
      Do k1 = 1,NumInt
      Do l1 = 1,NumInt
      coef = C(i1,i)*C(j1,j)*C(k1,k)*C(l1,l)
      sum1 = sum1+coef*D4(i1,j1,k1,l1)
      sum2 = sum2+coef*Gdbleprime(i1,j1,k1,l1)
      End Do
      End Do
      End Do
      End Do
      V4(i,j,k,l) = sum1
      T4(i,j,k,l) = sum2
      End Do
      End Do
      End Do
      End Do
!!
!!
!!---- Calculate diagonal anharmonicity constants.
      Do i = 1,NumInt
      x(i,i) = V4(i,i,i,i)/16.0d0
      tmp = 1.0d0/(48.0d0*harmfreq(i))
      x(i,i) = x(i,i)-tmp*(V3(i,i,i)*                                   &
     &    (5.0d0*V3(i,i,i)+6.0d0*T3(i,i,i))+9.0d0*T3(i,i,i)**2)
      Do j = 1,NumInt
      If ( j.ne.i ) Then
      tmp =-1.0d0/(16.0d0*harmfreq(j)*                                  &
     &          (4.0d0*harmfreq(i)**2-harmfreq(j)**2))
      x(i,i) = x(i,i)+                                                  &
     &          tmp*(8.0d0*harmfreq(i)**2-3.0d0*harmfreq(j)**2)*        &
     &                     (V3(i,i,j)**2+T3(i,i,j)**2)
      x(i,i) = x(i,i)+                                                  &
     &          tmp*(8.0d0*harmfreq(i)**2-harmfreq(j)**2)*              &
     &         (2.0d0*V3(i,i,j)*T3(i,i,j))
      x(i,i) = x(i,i)-                                                  &
     &          tmp*8.0d0*harmfreq(i)*harmfreq(j)*T3(i,j,i)*            &
     &                     (V3(i,i,j)-T3(i,i,j))
      x(i,i) = x(i,i)-tmp*4.0d0*harmfreq(j)**2*T3(i,j,i)**2
      End If
      End Do
      End Do
!!
!!---- Calculate off-diagonal anharmonicity constants.
      Do i = 1,NumInt
      Do j = 1,NumInt
      If ( i.ne.j ) Then
      x(i,j) = V4(i,i,j,j)/4.0d0
      x(i,j) = x(i,j)+T4(i,i,j,j)/2.0d0
      x(i,j) = x(i,j)-                                                  &
     &                     (V3(i,i,i)+T3(i,i,i))*                       &
     &                  (V3(i,j,j)+T3(j,j,i))/(4.0d0*harmfreq(i))
      x(i,j) = x(i,j)-                                                  &
     &                     (V3(i,i,j)+T3(i,i,j))*                       &
     &                  (V3(j,j,j)+T3(j,j,j))/(4.0d0*harmfreq(j))
      Do k = 1,NumInt
      If (( k.ne.i ).and.( k.ne.j )) Then
      x(i,j) = x(i,j)-                                                  &
     &                           (V3(i,i,k)+T3(i,i,k))*                 &
     &                   (V3(j,j,k)+T3(j,j,k))/(2.0d0*harmfreq(k))
      End If
      End Do
      tmp = 0.5d0*harmfreq(i)/                                          &
     &          (4.0d0*harmfreq(i)**2-harmfreq(j)**2)
      x(i,j) = x(i,j)-tmp*((V3(i,i,j)-T3(i,i,j))**2+                    &
     &                     2.0d0*(V3(j,j,i)-T3(j,j,i))*T3(i,j,j)+       &
     &                    4.0d0*T3(i,j,i)**2)
      tmp = 0.5d0*harmfreq(j)/                                          &
     &                          (4.0d0*harmfreq(j)**2-harmfreq(i)**2)
      x(i,j) = x(i,j)-tmp*((V3(j,j,i)-T3(j,j,i))**2+                    &
     &                     2.0d0*(V3(i,i,j)-T3(i,i,j))*                 &
     &         T3(i,j,i)+4.0d0*T3(i,j,j)**2)
      Do k = 1,NumInt
      If (( k.ne.i ).and.( k.ne.j )) Then
      delta = 1.0d0
      sum =  harmfreq(i)+harmfreq(j)+harmfreq(k)
      delta = delta*sum
      sum =  harmfreq(i)-harmfreq(j)-harmfreq(k)
      delta = delta*sum
      sum = -harmfreq(i)+harmfreq(j)-harmfreq(k)
      delta = delta*sum
      sum = -harmfreq(i)-harmfreq(j)+harmfreq(k)
      delta = delta*sum
      deltaInv = 1.0d0/delta
      x(i,j) = x(i,j)+deltaInv*0.5d0*harmfreq(k)*                       &
     &                           (harmfreq(i)**2+harmfreq(j)**2-        &
     &                           harmfreq(k)**2)*                       &
     &                           (V3(i,j,k)**2+T3(i,j,k)**2+            &
     &                           T3(j,k,i)**2+T3(k,i,j)**2)
      x(i,j) = x(i,j)+deltaInv*2.0d0*                                   &
     &                      (harmfreq(i)*harmfreq(j)*harmfreq(k))*      &
     &                      (V3(i,j,k)*T3(i,j,k)-T3(j,k,i)*T3(k,i,j))
      x(i,j) = x(i,j)-deltaInv*harmfreq(i)*                             &
     &                      (harmfreq(k)**2+harmfreq(j)**2-             &
     &                       harmfreq(i)**2)*                           &
     &                       (V3(i,j,k)*T3(k,i,j)-T3(i,j,k)*T3(j,k,i))
      x(i,j) = x(i,j)-deltaInv*harmfreq(j)*                             &
     &                           (harmfreq(k)**2+harmfreq(i)**2-        &
     &                            harmfreq(j)**2)*                      &
     &                      (V3(i,j,k)*T3(j,k,i)-T3(i,j,k)*T3(k,i,j))
      End If
      End Do
      End If
      End Do
      End Do
!!
!!
! Avoid unused argument warnings
      If (.False.) Call Unused_integer(max_term)
      End
