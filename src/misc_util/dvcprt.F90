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

subroutine DVcPrt(Title,FmtIn,X,N)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Print a vector of double precision real numbers                  *
!                                                                      *
!     calling arguments                                                *
!     Title  : character string containing a title                     *
!              If the string is empty now title will be printed        *
!     X      : vector of double precision reals                        *
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

use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: Title, FmtIn
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: X(N)
integer(kind=iwp), parameter :: lPaper = 120
integer(kind=iwp) :: i, iPmax, iPmin, lFmt, lItem, lLeft, lLine, lNumbr, lTitle, nCols, nDecim, nDigit
real(kind=wp) :: Pmax, Pmin, Xmax, Xmin
character(len=lPaper) :: Line
character(len=20) :: FRMT

!----------------------------------------------------------------------*
! print the title                                                      *
!----------------------------------------------------------------------*
lTitle = len_trim(Title)
if (lTitle > 0) then
  Line = ''
  lLeft = 0
  do i=1,lTitle
    if (Title(i:i) /= ' ') then
      lLeft = i-1
      exit
    end if
  end do
  do i=1,lPaper
    if (i+lLeft <= lTitle) Line(i:i) = Title(i+lLeft:i+lLeft)
  end do
  write(u6,*)
  write(u6,'(2X,A)') Line
  do i=1,len_trim(Line)
    Line(i:i) = '-'
  end do
  write(u6,'(2X,A)') Line
  write(u6,'(2X,A,I6)') 'vec. size = ',N
end if
!----------------------------------------------------------------------*
! determine the printing format                                        *
!----------------------------------------------------------------------*
lFmt = len_trim(FmtIn)
if (lFmt /= 0) then
  FRMT = FmtIn
else
  Xmax = -huge(Xmax)
  Xmin = huge(Xmin)
  do i=1,N
    Xmax = max(Xmax,X(i))
    Xmin = min(Xmin,X(i))
  end do
  Pmax = 0
  if (abs(Xmax) > 1.0e-72_wp) Pmax = log10(abs(Xmax))
  iPmax = 1+int(Pmax)
  iPmax = max(1,iPmax)
  Pmin = 0
  if (abs(Xmin) > 1.0e-72_wp) Pmin = log10(abs(Xmin))
  iPmin = 1+int(Pmin)
  iPmin = max(1,iPmin)
  nDigit = 14
  nDecim = min(8,nDigit-max(iPmin,iPmax))
  if (Xmax < 0) iPmax = iPmax+1
  if (Xmin < 0) iPmin = iPmin+1
  lNumbr = max(iPmin,iPmax)+nDecim+1
  nCols = 10
  lLine = nCols*lNumbr
  if (lLine > lPaper) then
    if ((lLine <= lPaper+nCols) .and. (nDecim > 1)) then
      nDecim = nDecim-1
      lNumbr = max(iPmin,iPmax)+nDecim
      lItem = max(lNumbr,lPaper/nCols)
    else
      nCols = 5
      lItem = max(lNumbr,lPaper/nCols)
    end if
  else
    lItem = lNumbr
  end if
  write(FRMT,'(A,I2.2,A,I2.2,A,I2.2,A)') '(2X,',nCols,'F',lItem,'.',nDecim,')'
end if
!----------------------------------------------------------------------*
! print the data                                                       *
!----------------------------------------------------------------------*
write(u6,*)
write(u6,FRMT) X(:)

!----------------------------------------------------------------------*
! End procedure                                                        *
!----------------------------------------------------------------------*
return

end subroutine DVcPrt
