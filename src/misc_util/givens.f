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
* Copyright (C) 2005, Per-Olof Widmark                                 *
************************************************************************
      Subroutine Givens(H,U,n,nv)
************************************************************************
*                                                                      *
* This routine transforms a symmetric matrix to tridiagonal form using *
* Givens rotations. The matrix is stored in lower triangular form.     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund university, Sweden                                     *
* Written September 2005                                               *
*                                                                      *
************************************************************************
      Implicit None
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
* n   - Dimension of matrix                                            *
* nv  - Length of eigenvectors nv>=n                                   *
* H   - Matrix to be diagonalized                                      *
* U   - Eigenvectors                                                   *
*----------------------------------------------------------------------*
      Integer n
      Integer nv
      Real*8  H(*)
      Real*8  U(nv,n)
*----------------------------------------------------------------------*
* Parameters                                                           *
* zThr - Threshold for when an element is regarded as zero             *
*----------------------------------------------------------------------*
      Real*8  zThr
      Parameter ( zThr = 1.0d-16 )
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Real*8  p,q
      Real*8  Hii,Hjj,Hij
      Real*8  tmp
      Integer iSkip
      Integer ii,ij,jj,ik,jk,im,jm
      Integer i,j,k,m
*----------------------------------------------------------------------*
*                                                                      *
* i,j - rotate around these indices                                    *
* i,k - eliminate this element (k=j-1)                                 *
*                                                                      *
*----------------------------------------------------------------------*
      Do j=2,n-1
         Do i=j+1,n
            k=j-1
            ii=i*(i+1)/2
            jj=j*(j+1)/2
            ij=j+i*(i-1)/2
            ik=k+i*(i-1)/2
            jk=k+j*(j-1)/2
            Hii=H(ii)
            Hjj=H(jj)
            Hij=H(ij)

            iSkip=0
            If(Abs(H(ik)).lt.zThr) Then
               p=1.0d0
               q=0.0d0
               iSkip=1
            Else If(Abs(H(jk)).lt.zThr) Then
               p=0.0d0
               q=1.0d0
            Else If(Abs(H(jk)).lt.Abs(H(ik))) Then
               tmp=H(jk)/H(ik)
               p=tmp/Sqrt(1.0d0+tmp*tmp)
               q=Sqrt(1.0d0-p*p)
               If(p.lt.0.0d0) Then
                  p=-p
                  q=-q
               End If
            Else
               tmp=H(ik)/H(jk)
               q=tmp/Sqrt(1.0d0+tmp*tmp)
               p=Sqrt(1.0d0-q*q)
            End If
            If(iSkip.eq.1) Goto 101

            Do m=1,n
               If(m.lt.j) Then
                  im=m+i*(i-1)/2
                  jm=m+j*(j-1)/2
               Else If(m.lt.i) Then
                  im=m+i*(i-1)/2
                  jm=j+m*(m-1)/2
               Else
                  im=i+m*(m-1)/2
                  jm=j+m*(m-1)/2
               End If
               tmp  =p*H(im)-q*H(jm)
               H(jm)=q*H(im)+p*H(jm)
               H(im)=tmp
            End Do

            H(ii)=p*p*Hii+q*q*Hjj-2.0d0*p*q*Hij
            H(jj)=q*q*Hii+p*p*Hjj+2.0d0*p*q*Hij
            H(ij)=(p*p-q*q)*Hij+p*q*(Hii-Hjj)
            H(ik)=0.0d0
            Do m=1,nv
               tmp   =p*U(m,i)-q*U(m,j)
               U(m,j)=q*U(m,i)+p*U(m,j)
               U(m,i)=tmp
            End Do
101         Continue
         End Do
      End Do
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
