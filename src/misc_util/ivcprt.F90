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
! Copyright (C) 1992, Markus P. Fuelscher                              *
!***********************************************************************

subroutine IVcPrt(Title,FmtIn,X,N)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Print a vector of integer numbers                                *
!                                                                      *
!     calling arguments                                                *
!     Title  : character string containing a title                     *
!              If the string is empty now title will be printed        *
!     X      : vector of integers                                      *
!     N      : dimension of vector X                                   *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M. P. Fuelscher                                                  *
!     University of Lund, Sweden, 1992                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

implicit integer(A-Z)
character*(*) Title
character*(*) FmtIn
dimension X(N)
integer StrnLn
parameter(lPaper=120)
character*(lPaper) Line
real*8 Temp, Pmax, Pmin
character*20 FMT

!----------------------------------------------------------------------*
! print the title                                                      *
!----------------------------------------------------------------------*
lTitle = StrnLn(Title)
if (lTitle > 0) then
  do i=1,lPaper
    Line(i:i) = ' '
  end do
  lLeft = 1
  do i=lTitle,1,-1
    if (Title(i:i) /= ' ') lLeft = i
  end do
  lLeft = lLeft-1
  do i=1,lPaper
    if (i+lLeft <= lTitle) Line(i:i) = Title(i+lLeft:i+lLeft)
  end do
  write(6,*)
  write(6,'(2X,A)') Line
  do i=1,StrnLn(Line)
    Line(i:i) = '-'
  end do
  write(6,'(2X,A)') Line
  write(6,'(2X,A,I6)') 'vec. size = ',N
end if
!----------------------------------------------------------------------*
! determine the printing format                                        *
!----------------------------------------------------------------------*
lFmt = StrnLn(FmtIn)
if (lFmt /= 0) then
  FMT = FmtIn
else
  Xmax = X(1)
  Xmin = X(1)
  do i=1,N
    Xmax = max(Xmax,X(i))
    Xmin = min(Xmin,X(i))
  end do
  Temp = dble(Xmax)
  Pmax = 0d0
  if (abs(Temp) > 1.0D-72) Pmax = log10(abs(Temp))
  iPmax = int(1d0+Pmax)
  iPmax = max(1,iPmax)
  if (Xmax < 0) iPmax = iPmax+1
  Temp = dble(Xmin)
  Pmin = 0d0
  if (abs(Temp) > 1.0D-72) Pmin = log10(abs(Temp))
  iPmin = int(1d0+Pmin)
  iPmin = max(1,iPmin)
  if (Xmin < 0) iPmin = iPmin+1
  lNumbr = max(iPmax,iPmin)+1
  if (50*lNumbr <= lPaper) then
    nCols = 50
  else if (20*lNumbr <= lPaper) then
    nCols = 20
  else if (10*lNumbr <= lPaper) then
    nCols = 10
  else
    nCols = 5
  end if
  lItem = lPaper/nCols
  write(FMT,'(A, I2.2,  A, I2.2,  A)') '(2X,',nCols,'I',lItem,')'
end if
!----------------------------------------------------------------------*
! print the data                                                       *
!----------------------------------------------------------------------*
write(6,*)
write(6,FMT) (X(i),i=1,N)

!----------------------------------------------------------------------*
! End procedure                                                        *
!----------------------------------------------------------------------*
return

end subroutine IVcPrt
