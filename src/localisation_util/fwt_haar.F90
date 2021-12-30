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
! Copyright (C) 2009, Francesco Aquilante                              *
!***********************************************************************

subroutine FWT_Haar(n,m,B,X)
!***********************************************************************
!
! Author :  F. Aquilante  (Geneva, March 2009)
!
!
! Purpose : performs a Fast Wavelet Transform (FWT)
!           of 2^m vectors X_1(n),X_2(n),...,X_2^m(n)
!           using m-th order Haar wavelet.
!
!
! n : number of components of each vector (input)
! m : order of Haar wavelet = Log_2 of the number of vectors (input)
! B(n,2^m-1) : scratch array (input)
! X(n,2^m) : array of initial/transformed vectors (input/output)
!
!***********************************************************************
!
! Mallat's pyramid algorithm (IEEE Trans PAMI, 11, 674-693, 1989) :
!
!    (X_1,X_2,X_3,X_4,...,X_2^m) ---> (A[0],B[0],B[1],...,B[m-1])
!
!    j=m,m-1,...,1
!                   B[j-1] = H[j] * A[j]
!                   A[j-1] = L[j] * A[j]
!
!
!    The matrices H[j] and L[j] have:
!
!    2^(j-1) orthonormal rows, and
!
!    2^j cols pair-orthonormal (merging col k and k+1, the result is
!             2^(j-1) orthonormal cols each containing 2^j elements)
!
!    and very few nonzero entries.
!
!    High Pass (differencing) filters:
!
!               [  1 -1                          ]
!               [        1 -1                    ]
!    H[j] = 1/2 [              ....              ]
!               [                    ....        ]
!               [                          1 -1  ]
!
!    Low Pass (averaging) filters:
!
!               [  1  1                          ]
!               [        1  1                    ]
!    L[j] = 1/2 [              ....              ]
!               [                    ....        ]
!               [                          1  1  ]
!
!***********************************************************************

implicit real*8(a-h,o-z)
integer n, m
real*8 B(n,*), X(n,*)
parameter(hlf=0.5d0)

if (m <= 0) then

  write(6,*) ' FWT_Haar: Illegal value of m = ',m
  call Abend()

elseif (n <= 0) then

  write(6,*) ' FWT_Haar: Illegal value of n = ',n
  call Abend()

elseif (n > 50) then  ! use BLAS

  call FWT_Haar_(n,m,B,X)

else  ! do not use BLAS

  fac = sqrt(hlf)
  nv = 2**m
  do j=m,1,-1  ! A[m]:=X
    nv = nv/2
    kB = nv
    do k=1,nv
      kA = 2*k
      lA = kA-1
      do i=1,n
        B(i,kB) = fac*(X(i,lA)-X(i,kA))
        X(i,k) = fac*(X(i,lA)+X(i,kA))
      end do
      kB = kB+1
    end do
  end do
  do j=1,2**m-1
    k = j+1  ! A[0] is already in place
    do i=1,n
      X(i,k) = B(i,j)
    end do
  end do

end if

return

end subroutine FWT_Haar
