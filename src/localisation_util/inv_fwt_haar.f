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
      SUBROUTINE Inv_FWT_Haar(n,m,B,X)
************************************************************************
*
*     Author :  F. Aquilante  (Geneva, April 2009)
*
*
*     Purpose : performs an inverse Fast Wavelet Transform (FWT)
*               of 2^m vectors Y_1(n),Y_2(n),...,Y_2^m(n)
*               using m-th order Haar wavelet.
*
*
*     n : number of components of each vector (input)
*     m : order of Haar wavelet = Log_2 of the number of vectors (input)
*     B(n,2^m) : scratch array (input)
*     X(n,2^m) : array of initial/transformed vectors (input/output)
*
************************************************************************
*
*     Mallat's pyramid algorithm (IEEE Trans PAMI, 11, 674-693, 1989) :
*
*      Y=(A[0],B[0],B[1],...,B[m-1]) ---> (X_1,X_2,X_3,X_4,...,X_2^m)
*
*        j=1,2,...,m
*                       A[j] = L[j]^t * A[j-1] + H[j]^t * B[j-1]
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
*
*
*
      If (m .le. 0) Then

         Write(6,*) ' Inv_FWT_Haar: Illegal value of m = ',m
         Call Abend

      ElseIf (n .le. 0) Then

         Write(6,*) ' Inv_FWT_Haar: Illegal value of n = ',n
         Call Abend

      Else  ! do not use BLAS

         fac=sqrt(0.5d0)
         nv=1
         Do j=1,m
            kB=-1
            Do k=1,nv
               kB=kB+2
               lB=kB+1
               kv=k+nv
               Do i=1,n
                  B(i,kB) = fac*(X(i,k)+X(i,kv))
                  B(i,lB) = fac*(X(i,k)-X(i,kv))
               End Do
            End Do
            nv=nv*2
            call dcopy_(n*nv,B,1,X,1)
         End Do

      EndIf
*
      Return
      End
