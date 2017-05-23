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
      Subroutine QLdiag(H,U,n,nv,irc)
************************************************************************
*                                                                      *
* This routine diagonalize a tridiagonal symmetric matrix using the    *
* QL algorithm. The matrix is stored in lower triangular form.         *
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
* irc - Return code, 0 = OK                                            *
*                    1 = Not converged                                 *
*                    2 = Too large system.                             *
*----------------------------------------------------------------------*
      Integer n
      Integer nv
      Integer irc
      Real*8  H(*)
      Real*8  U(nv,n)
*----------------------------------------------------------------------*
* Parameters                                                           *
* MxDim - Largest case that can be handled                             *
* zThr  - Threshold for when an element is regarded as zero            *
*----------------------------------------------------------------------*
      Integer MxDim
      Parameter ( MxDim = 5000 )
      Real*8  zThr
      Parameter ( zThr = 1.0d-16 )
      Real*8  qThr
      Parameter ( qThr = 1.0d-20 )
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Real*8  d(MxDim),e(MxDim)
      Real*8  g,r,c,s,p,f,b
      Integer i,j,k,l,m
      Integer iter,maxiter
*----------------------------------------------------------------------*
* Setup                                                                *
*----------------------------------------------------------------------*
      irc=0
      If(n.ge.MxDim) Then
         irc=1
c         Write(6,*) 'QLdiag: system too large!'
         Return
      End If
      If(n.lt.1) Then
         Write(6,*) 'QLdiag: zero size system!'
         Call Abend()
      End If
*----------------------------------------------------------------------*
* Make local copies of diagonal and off-diagonal                       *
*----------------------------------------------------------------------*
      j=1
      Do i=1,n
         d(i)=H(j)
         j=j+i+1
      End Do
      j=2
      Do i=1,n-1
         e(i)=H(j)
         j=j+i+2
      End Do
      e(n)=0.0d0
*----------------------------------------------------------------------*
* Solve it                                                             *
*----------------------------------------------------------------------*
      maxiter=0
      Do l=1,n
         iter=0
1        Continue
         Do j=l,n-1
            If(Abs(e(j)).lt.zThr) Then
               m=j
               Goto 2
            End If
         End Do
         m=n
2        Continue
         If(m.ne.l) Then
            If(iter.eq.25) Then
               irc=1
c               Write(6,*) 'QLdiag: ran out of iterations'
               Goto 900
            End If
            iter=iter+1
            maxiter=Max(maxiter,iter)
            g=(d(l+1)-d(l))/(2.0d0*e(l))
            r=Sqrt(1.0d0+g*g)
            g=d(m)-d(l)+e(l)/(g+Sign(r,g))
            s=1.0d0
            c=1.0d0
            p=0.0d0
            Do i=m-1,l,-1
               f=s*e(i)
               b=c*e(i)
               r=Sqrt(f*f+g*g)
               e(i+1)=r
               If(Abs(r).le.qThr) Then
                  d(i+1)=d(i+1)-p
                  e(m)=0.0d0
                  Goto 1
               End If
               s=f/r
               c=g/r
               g=d(i+1)-p
               r=(d(i)-g)*s+2.0d0*c*b
               p=s*r
               d(i+1)=g+p
               g=c*r-b
               Do k=1,nv
                  f=U(k,i+1)
                  U(k,i+1)=s*U(k,i)+c*f
                  U(k,i)=c*U(k,i)-s*f
               End Do
            End Do
            d(l)=d(l)-p
            e(l)=g
            e(m)=0.0d0
            Goto 1
         End If
      End Do
*----------------------------------------------------------------------*
* Copy back local copy                                                 *
*----------------------------------------------------------------------*
900   Continue
*     Write(6,*) 'QLdiag: maxiter ',maxiter
      j=1
      Do i=1,n
         H(j)=d(i)
         j=j+i+1
      End Do
      j=2
      Do i=1,n-1
         H(j)=e(i)
         j=j+i+2
      End Do
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
