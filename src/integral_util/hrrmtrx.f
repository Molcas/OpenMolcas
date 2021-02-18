************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1999, Roland Lindh                                     *
************************************************************************
      Subroutine HrrMtrx(HMtrx,np,la,lb,A,B,
     &                   Sph_a,CS_a,nSph_a,Sph_b,Cs_b,nSph_b)
************************************************************************
*                                                                      *
*     Object: to compute the matrix which corresponds to the transfer  *
*             equation.                                                *
*                                                                      *
*     Author: Roland Lindh                                             *
*             Dept of Chem. Phys.                                      *
*             Univ. of Lund, Sweden                                    *
*             February 1999                                            *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "binom.fh"
#include "ican.fh"
#include "real.fh"
      Real*8 A(3), B(3), HMtrx(np,nSph_a,nSph_b), AB(3,0:iTabMx),
     &       CS_a((la+1)*(la+2)/2,nSph_a), CS_b((lb+1)*(lb+2)/2,nSph_b)
      Logical EQ, Sph_a, Sph_b
*
      iOff(ixyz) = ixyz*(ixyz+1)*(ixyz+2)/6
      jCan(ix,iy,iz) = iOff(ix+iy+iz) + (iy+iz)*(iy+iz+1)/2 + iz + 1
*
      iPrint=5
      If (iPrint.ge.99) Then
         Call RecPrt('A',' ',A,1,3)
         Call RecPrt('B',' ',B,1,3)
         Call RecPrt('CS_a',' ',CS_a,(la+1)*(la+2)/2,nSph_a)
         Call RecPrt('CS_b',' ',CS_b,(lb+1)*(lb+2)/2,nSph_b)
         Write (6,*) 'np=',np
      End If
      Call FZero(HMtrx,np*nSph_a*nSph_b)
*
      AB(1,0)=One
      AB(2,0)=One
      AB(3,0)=One
      If (la.ge.lb) Then
         AB(1,1)=A(1)-B(1)
         AB(2,1)=A(2)-B(2)
         AB(3,1)=A(3)-B(3)
      Else
         AB(1,1)=B(1)-A(1)
         AB(2,1)=B(2)-A(2)
         AB(3,1)=B(3)-A(3)
      End If
      Do i = 2, Min(la,lb)
         AB(1,i)=AB(1,i-1)*AB(1,1)
         AB(2,i)=AB(2,i-1)*AB(2,1)
         AB(3,i)=AB(3,i-1)*AB(3,1)
      End Do
*
      If (la.ge.lb) Then
*
         If (Sph_a.and.Sph_b) Then
*
         Do iSph_a = 1, nSph_a
         Do ipa = 1, (la+1)*(la+2)/2
            C_a=CS_a(ipa,iSph_a)
            If (C_a.eq.Zero) Go To 100
            ix=iCan(1,iOff(la)+ipa)
            iy=iCan(2,iOff(la)+ipa)
            iz=iCan(3,iOff(la)+ipa)
            Do iSph_b = 1, nSph_b
            Do ipb = 1, (lb+1)*(lb+2)/2
               C_b=CS_b(ipb,iSph_b)
               If (C_b.eq.Zero) Go To 200
               jx=iCan(1,iOff(lb)+ipb)
               jy=iCan(2,iOff(lb)+ipb)
               jz=iCan(3,iOff(lb)+ipb)
*
               If (EQ(A,B)) Then
                  ixLow=ix+jx
                  iyLow=iy+jy
                  izLow=iz+jz
                  jOff=iOff(la+lb)
               Else
                  ixLow=ix
                  iyLow=iy
                  izLow=iz
                  jOff=iOff(la)
               End If
               Do kx = ixLow, ix+jx
               Do ky = iyLow, iy+jy
               Do kz = izLow, iz+jz
                  ipe=jCan(kx,ky,kz)-jOff
*
                     ABx = AB(1,ix+jx-kx) * Binom(jx,kx-ix)
                     ABy = AB(2,iy+jy-ky) * Binom(jy,ky-iy)
                     ABz = AB(3,iz+jz-kz) * Binom(jz,kz-iz)
                     HMtrx(ipe,iSph_a,iSph_b) =
     &                  HMtrx(ipe,iSph_a,iSph_b) + ABx * ABy * ABz
     &                    * C_a * C_b
               End Do
               End Do
               End Do
*
 200           Continue
            End Do
            End Do
 100        Continue
         End Do
         End Do
*
         Else If (Sph_a) Then
*
         Do iSph_a = 1, nSph_a
         Do ipa = 1, (la+1)*(la+2)/2
            C_a=CS_a(ipa,iSph_a)
            If (C_a.eq.Zero) Go To 101
            ix=iCan(1,iOff(la)+ipa)
            iy=iCan(2,iOff(la)+ipa)
            iz=iCan(3,iOff(la)+ipa)
            Do ipb = 1, (lb+1)*(lb+2)/2
               jx=iCan(1,iOff(lb)+ipb)
               jy=iCan(2,iOff(lb)+ipb)
               jz=iCan(3,iOff(lb)+ipb)
*
               If (EQ(A,B)) Then
                  ixLow=ix+jx
                  iyLow=iy+jy
                  izLow=iz+jz
                  jOff=iOff(la+lb)
               Else
                  ixLow=ix
                  iyLow=iy
                  izLow=iz
                  jOff=iOff(la)
               End If
               Do kx = ixLow, ix+jx
               Do ky = iyLow, iy+jy
               Do kz = izLow, iz+jz
                  ipe=jCan(kx,ky,kz)-jOff
*
                     ABx = AB(1,ix+jx-kx) * Binom(jx,kx-ix)
                     ABy = AB(2,iy+jy-ky) * Binom(jy,ky-iy)
                     ABz = AB(3,iz+jz-kz) * Binom(jz,kz-iz)
                     HMtrx(ipe,iSph_a,ipb) =
     &                  HMtrx(ipe,iSph_a,ipb) + ABx * ABy * ABz
     &                    * C_a
               End Do
               End Do
               End Do
*
            End Do
 101        Continue
         End Do
         End Do
*
         Else If (Sph_b) Then
*
         Do ipa = 1, (la+1)*(la+2)/2
            ix=iCan(1,iOff(la)+ipa)
            iy=iCan(2,iOff(la)+ipa)
            iz=iCan(3,iOff(la)+ipa)
            Do iSph_b = 1, nSph_b
            Do ipb = 1, (lb+1)*(lb+2)/2
               C_b=CS_b(ipb,iSph_b)
               If (C_b.eq.Zero) Go To 201
               jx=iCan(1,iOff(lb)+ipb)
               jy=iCan(2,iOff(lb)+ipb)
               jz=iCan(3,iOff(lb)+ipb)
*
               If (EQ(A,B)) Then
                  ixLow=ix+jx
                  iyLow=iy+jy
                  izLow=iz+jz
                  jOff=iOff(la+lb)
               Else
                  ixLow=ix
                  iyLow=iy
                  izLow=iz
                  jOff=iOff(la)
               End If
               Do kx = ixLow, ix+jx
               Do ky = iyLow, iy+jy
               Do kz = izLow, iz+jz
                  ipe=jCan(kx,ky,kz)-jOff
*
                     ABx = AB(1,ix+jx-kx) * Binom(jx,kx-ix)
                     ABy = AB(2,iy+jy-ky) * Binom(jy,ky-iy)
                     ABz = AB(3,iz+jz-kz) * Binom(jz,kz-iz)
                     HMtrx(ipe,ipa,iSph_b) =
     &                  HMtrx(ipe,ipa,iSph_b) + ABx * ABy * ABz
     &                    * C_b
               End Do
               End Do
               End Do
*
 201           Continue
            End Do
            End Do
         End Do
*
         Else
*
         Do ipa = 1, (la+1)*(la+2)/2
            ix=iCan(1,iOff(la)+ipa)
            iy=iCan(2,iOff(la)+ipa)
            iz=iCan(3,iOff(la)+ipa)
            Do ipb = 1, (lb+1)*(lb+2)/2
               jx=iCan(1,iOff(lb)+ipb)
               jy=iCan(2,iOff(lb)+ipb)
               jz=iCan(3,iOff(lb)+ipb)
*
               If (EQ(A,B)) Then
                  ixLow=ix+jx
                  iyLow=iy+jy
                  izLow=iz+jz
                  jOff=iOff(la+lb)
               Else
                  ixLow=ix
                  iyLow=iy
                  izLow=iz
                  jOff=iOff(la)
               End If
               Do kx = ixLow, ix+jx
               Do ky = iyLow, iy+jy
               Do kz = izLow, iz+jz
                  ipe=jCan(kx,ky,kz)-jOff
*
                     ABx = AB(1,ix+jx-kx) * Binom(jx,kx-ix)
                     ABy = AB(2,iy+jy-ky) * Binom(jy,ky-iy)
                     ABz = AB(3,iz+jz-kz) * Binom(jz,kz-iz)
                     HMtrx(ipe,ipa,ipb) =
     &                  HMtrx(ipe,ipa,ipb) + ABx * ABy * ABz
               End Do
               End Do
               End Do
*
            End Do
         End Do
*
         End If
      Else
*
         If (Sph_a.and.Sph_b) Then
*
         Do iSph_a = 1, nSph_a
         Do ipa = 1, (la+1)*(la+2)/2
            C_a=CS_a(ipa,iSph_a)
            If (C_a.eq.Zero) Go To 300
            ix=iCan(1,iOff(la)+ipa)
            iy=iCan(2,iOff(la)+ipa)
            iz=iCan(3,iOff(la)+ipa)
            Do iSph_b = 1, nSph_b
            Do ipb = 1, (lb+1)*(lb+2)/2
               C_b=CS_b(ipb,iSph_b)
               If (C_b.eq.Zero) Go To 400
               jx=iCan(1,iOff(lb)+ipb)
               jy=iCan(2,iOff(lb)+ipb)
               jz=iCan(3,iOff(lb)+ipb)
*
               If (EQ(A,B)) Then
                  jxLow=ix+jx
                  jyLow=iy+jy
                  jzLow=iz+jz
                  jOff=iOff(la+lb)
               Else
                  jxLow=jx
                  jyLow=jy
                  jzLow=jz
                  jOff=iOff(lb)
               End If
               Do kx = jxLow, ix+jx
               Do ky = jyLow, iy+jy
               Do kz = jzLow, iz+jz
                  ipe=jCan(kx,ky,kz)-jOff
*
                  ABx = AB(1,ix+jx-kx) * Binom(ix,kx-jx)
                  ABy = AB(2,iy+jy-ky) * Binom(iy,ky-jy)
                  ABz = AB(3,iz+jz-kz) * Binom(iz,kz-jz)
                  HMtrx(ipe,iSph_a,iSph_b) =
     &               HMtrx(ipe,iSph_a,iSph_b) + ABx * ABy * ABz
     &                 * C_a * C_b
               End Do
               End Do
               End Do
*
 400           Continue
            End Do
            End Do
 300        Continue
         End Do
         End Do
*
         Else If (Sph_a) Then
*
         Do iSph_a = 1, nSph_a
         Do ipa = 1, (la+1)*(la+2)/2
            C_a=CS_a(ipa,iSph_a)
            If (C_a.eq.Zero) Go To 301
            ix=iCan(1,iOff(la)+ipa)
            iy=iCan(2,iOff(la)+ipa)
            iz=iCan(3,iOff(la)+ipa)
            Do ipb = 1, (lb+1)*(lb+2)/2
               jx=iCan(1,iOff(lb)+ipb)
               jy=iCan(2,iOff(lb)+ipb)
               jz=iCan(3,iOff(lb)+ipb)
*
               If (EQ(A,B)) Then
                  jxLow=ix+jx
                  jyLow=iy+jy
                  jzLow=iz+jz
                  jOff=iOff(la+lb)
               Else
                  jxLow=jx
                  jyLow=jy
                  jzLow=jz
                  jOff=iOff(lb)
               End If
               Do kx = jxLow, ix+jx
               Do ky = jyLow, iy+jy
               Do kz = jzLow, iz+jz
                  ipe=jCan(kx,ky,kz)-jOff
*
                  ABx = AB(1,ix+jx-kx) * Binom(ix,kx-jx)
                  ABy = AB(2,iy+jy-ky) * Binom(iy,ky-jy)
                  ABz = AB(3,iz+jz-kz) * Binom(iz,kz-jz)
                  HMtrx(ipe,iSph_a,ipb) =
     &               HMtrx(ipe,iSph_a,ipb) + ABx * ABy * ABz
     &                 * C_a
               End Do
               End Do
               End Do
*
            End Do
 301        Continue
         End Do
         End Do
*
         Else If (Sph_b) Then
*
         Do ipa = 1, (la+1)*(la+2)/2
            ix=iCan(1,iOff(la)+ipa)
            iy=iCan(2,iOff(la)+ipa)
            iz=iCan(3,iOff(la)+ipa)
            Do iSph_b = 1, nSph_b
            Do ipb = 1, (lb+1)*(lb+2)/2
               C_b=CS_b(ipb,iSph_b)
               If (C_b.eq.Zero) Go To 401
               jx=iCan(1,iOff(lb)+ipb)
               jy=iCan(2,iOff(lb)+ipb)
               jz=iCan(3,iOff(lb)+ipb)
*
               If (EQ(A,B)) Then
                  jxLow=ix+jx
                  jyLow=iy+jy
                  jzLow=iz+jz
                  jOff=iOff(la+lb)
               Else
                  jxLow=jx
                  jyLow=jy
                  jzLow=jz
                  jOff=iOff(lb)
               End If
               Do kx = jxLow, ix+jx
               Do ky = jyLow, iy+jy
               Do kz = jzLow, iz+jz
                  ipe=jCan(kx,ky,kz)-jOff
*
                  ABx = AB(1,ix+jx-kx) * Binom(ix,kx-jx)
                  ABy = AB(2,iy+jy-ky) * Binom(iy,ky-jy)
                  ABz = AB(3,iz+jz-kz) * Binom(iz,kz-jz)
                  HMtrx(ipe,ipa,iSph_b) =
     &               HMtrx(ipe,ipa,iSph_b) + ABx * ABy * ABz
     &                 * C_b
               End Do
               End Do
               End Do
*
 401           Continue
            End Do
            End Do
         End Do
*
         Else
*
         Do ipa = 1, (la+1)*(la+2)/2
            ix=iCan(1,iOff(la)+ipa)
            iy=iCan(2,iOff(la)+ipa)
            iz=iCan(3,iOff(la)+ipa)
            Do ipb = 1, (lb+1)*(lb+2)/2
               jx=iCan(1,iOff(lb)+ipb)
               jy=iCan(2,iOff(lb)+ipb)
               jz=iCan(3,iOff(lb)+ipb)
*
               If (EQ(A,B)) Then
                  jxLow=ix+jx
                  jyLow=iy+jy
                  jzLow=iz+jz
                  jOff=iOff(la+lb)
               Else
                  jxLow=jx
                  jyLow=jy
                  jzLow=jz
                  jOff=iOff(lb)
               End If
               Do kx = jxLow, ix+jx
               Do ky = jyLow, iy+jy
               Do kz = jzLow, iz+jz
                  ipe=jCan(kx,ky,kz)-jOff
*
                  ABx = AB(1,ix+jx-kx) * Binom(ix,kx-jx)
                  ABy = AB(2,iy+jy-ky) * Binom(iy,ky-jy)
                  ABz = AB(3,iz+jz-kz) * Binom(iz,kz-jz)
                  HMtrx(ipe,ipa,ipb) =
     &               HMtrx(ipe,ipa,ipb) + ABx * ABy * ABz
               End Do
               End Do
               End Do
*
            End Do
         End Do
*
         End If
      End If
      If (iPrint.ge.99) Then
         Call RecPrt('HMat ( np x (nSph_a*nSph_b) )','(30F4.1)',HMtrx,
     &             np,nSph_a*nSph_b)
         Write (6,*) DDot_(np*nSph_a*nSph_b,HMtrx,1,HMtrx,1)
      End If
*
      Return
      End
