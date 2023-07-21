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
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!               2006, Francesco Aquilante                              *
!***********************************************************************

subroutine CD_InCore_1p(X,n,Vec,MxVec,NumCho,Thr,ThrNeg,ThrFail,iD,irc)
!
! Author:   F. Aquilante, March 2006
!
! Snicked from CD_InCore_1 :
!
! Thomas Bondo Pedersen, October 2004.
!
! Purpose: Cholesky decompose the n-by-n matrix X.
!          Vectors are returned in Vec array.
!
! Return code (irc):
!
!    0 -- decomposition success.
!  101 -- negative diagonal encountered (i.e. < ThrFail)
!  102 -- number of vectors needed exceeds max. allowed (MxVec)
!
! Note: the algorithm is designed for incomplete Cholesky
! decomposition, i.e. for semi-definitive matrices, and thus makes
! use of level-1 BLAS only.
!
! Feature: In/Out argument iD(MxVec).
!          iD(k) on exit contains the index of the diagonal
!                from which the k-th Cholesky vector was generated

implicit none
integer n, MxVec, NumCho, irc
real*8 X(n,n), Vec(n,MxVec)
real*8 Thr, ThrNeg, ThrFail
integer iD(MxVec)
integer i, imax, j
integer iPass
real*8 Xmax, Factor, xFac, Acc

irc = 0

NumCho = 0
Acc = min(1.0D-12,thr*1.0D-2)
xFac = 0.0d0  ! dummy set
do iPass=1,n

  if (X(1,1) < ThrNeg) then
    if (X(1,1) < ThrFail) then
      irc = 101
      return
    else
      do j=1,n
        X(j,1) = 0.0d0
        X(1,j) = 0.0d0
      end do
    end if
  end if
  Xmax = 0.0d0
  imax = 0
  do i=1,n
    if (X(i,i) < ThrNeg) then
      if (X(i,i) < ThrFail) then
        irc = 101
        return
      else
        do j=1,n
          X(j,i) = 0.0d0
          X(i,j) = 0.0d0
        end do
      end if
    end if
    if (X(i,i) > Xmax+Acc) then
      Xmax = X(i,i)
      imax = i
    end if
  end do

  if (Xmax <= Thr) return ! converged
  xFac = 1.0d0/sqrt(abs(Xmax))

  if (NumCho == MxVec) then ! too many vectors
    irc = 102
    return
  end if

  do j=1,NumCho
    Factor = -Vec(imax,j)
    call dAXPY_(n,Factor,Vec(1,j),1,X(1,imax),1)
  end do
  X(imax,imax) = Xmax

  NumCho = NumCho+1
  do i=1,n
    if (X(i,i) == 0.0d0) then
      Vec(i,NumCho) = 0.0d0
    else
      Vec(i,NumCho) = xFac*X(i,imax)
    end if
  end do

  do i=1,n
    X(i,i) = X(i,i)-Vec(i,NumCho)**2
  end do
  X(imax,imax) = 0.0d0
  iD(NumCho) = imax

end do

end subroutine CD_InCore_1p
