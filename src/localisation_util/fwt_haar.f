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
* Copyright (C) 2009, Francesco Aquilante                              *
************************************************************************
      SUBROUTINE FWT_Haar(n,m,B,X)
************************************************************************
*
*     Author :  F. Aquilante  (Geneva, March 2009)
*
*
*     Purpose : performs a Fast Wavelet Transform (FWT)
*               of 2^m vectors X_1(n),X_2(n),...,X_2^m(n)
*               using m-th order Haar wavelet.
*
*
*     n : number of components of each vector (input)
*     m : order of Haar wavelet = Log_2 of the number of vectors (input)
*     B(n,2^m-1) : scratch array (input)
*     X(n,2^m) : array of initial/transformed vectors (input/output)
*
************************************************************************
*
*     Mallat's pyramid algorithm (IEEE Trans PAMI, 11, 674-693, 1989) :
*
*        (X_1,X_2,X_3,X_4,...,X_2^m) ---> (A[0],B[0],B[1],...,B[m-1])
*
*        j=m,m-1,...,1
*                       B[j-1] = H[j] * A[j]
*                       A[j-1] = L[j] * A[j]
*
*
*        The matrices H[j] and L[j] have:
*
*        2^(j-1) orthonormal rows, and
*
*        2^j cols pair-orthonormal (merging col k and k+1, the result is
*                 2^(j-1) orthonormal cols each containing 2^j elements)
*
*        and very few nonzero entries.
*
*        High Pass (differencing) filters:
*
*                   [  1 -1                          ]
*                   [        1 -1                    ]
*        H[j] = 1/2 [              ....              ]
*                   [                    ....        ]
*                   [                          1 -1  ]
*
*        Low Pass (averaging) filters:
*
*                   [  1  1                          ]
*                   [        1  1                    ]
*        L[j] = 1/2 [              ....              ]
*                   [                    ....        ]
*                   [                          1  1  ]
*
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Integer n, m
      Real*8  B(n,*), X(n,*)

      Parameter ( hlf=0.5d0 )
*
*
*
      If (m .le. 0) Then

         Write(6,*) ' FWT_Haar: Illegal value of m = ',m
         Call Abend

      ElseIf (n .le. 0) Then

         Write(6,*) ' FWT_Haar: Illegal value of n = ',n
         Call Abend

      ElseIf (n .gt. 50) Then  ! use BLAS

         Call FWT_Haar_(n,m,B,X)

      Else  ! do not use BLAS

         fac=sqrt(hlf)
         nv=2**m
         Do j=m,1,-1  ! A[m]:=X
            nv=nv/2
            kB=nv
            Do k=1,nv
               kA=2*k
               lA=kA-1
               Do i=1,n
                  B(i,kB) = fac*(X(i,lA)-X(i,kA))
                  X(i,k) = fac*(X(i,lA)+X(i,kA))
               End Do
               kB=kB+1
            End Do
         End Do
         Do j=1,2**m-1
            k=j+1  ! A[0] is already in place
            Do i=1,n
               X(i,k)=B(i,j)
            End Do
         End Do

      EndIf
*
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      SUBROUTINE FWT_Haar_(n,m,B,X)
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Integer n, m
      Real*8  B(n,*), X(n,*)

      Parameter ( hlf=0.5d0, one=1.0d0, xone=-1.0d0 )

      fac=sqrt(hlf)
      nv=2**m
      mv=n*nv
      Do j=m,1,-1  ! A[m]:=X
         nv=nv/2
         kB=nv
         call dcopy_(n,X(1,1),1,B(1,kB),1)
         Call daxpy_(n,xone,X(1,2),1,B(1,kB),1)
         call dscal_(n,fac,B(1,kB),1)
         Call daxpy_(n,one,X(1,2),1,X(1,1),1)
         call dscal_(n,fac,X(1,1),1)
         Do k=2,nv
            kB=kB+1
            jB=2*k
            lB=jB-1
            call dcopy_(n,X(1,lB),1,B(1,kB),1)
            Call daxpy_(n,xone,X(1,jB),1,B(1,kB),1)
            call dscal_(n,fac,B(1,kB),1)
            Call daxpy_(n,one,X(1,jB),1,X(1,lB),1)
            call dcopy_(n,X(1,lB),1,X(1,k),1)
            call dscal_(n,fac,X(1,k),1)
         End Do
      End Do
      call dcopy_(mv-n,B(1,1),1,X(1,2),1) ! A[0] is already in place
*
      Return
      End
