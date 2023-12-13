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

!#define _DEBUGPRINT_
subroutine TrcPrt(Title,FmtIn,A,nRow,nCol)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Print the row and column norms of a rectangular matrix           *
!                                                                      *
!     calling arguments                                                *
!     Title  : character string containing a title                     *
!              If the string is empty now title will be printed        *
!     A      : retangular matrix of double precision reals             *
!     nRow   : row dimension of matrix A                               *
!     nCol   : column dimension of matrix A                            *
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

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
character(len=*), intent(in) :: Title, FmtIn
integer(kind=iwp), intent(in) :: nRow, nCol
real(kind=wp), intent(in) :: A(nRow,nCol)
#include "standard_iounits.fh"
integer(kind=iwp), parameter :: lPaper = 70
integer(kind=iwp) :: i, iPmax, iPmin, j, lFmt, lItem, lLeft, lLine, lNumbr, lTitle, nCols, nDecim, nDigit
real(kind=wp) :: Amax, Amin, Pmax, Pmin, Scal
character(len=lPaper) :: Line
character(len=20) :: FRMT
real(kind=wp), external :: DDot_

!----------------------------------------------------------------------*
! print the title                                                      *
!----------------------------------------------------------------------*
lTitle = len_trim(Title)
if (lTitle > 0) then
  Line = ''
  lLeft = 1
  do i=1,lTitle
    if (Title(i:i) /= ' ') then
      lLeft = i-1
      exit
    end if
  end do
  do i=1,lPaper
    if (i+lLeft <= lTitle) Line(i:i) = Title(i+lLeft:i+lLeft)
  end do
  write(LuWr,*)
  write(LuWr,'(2X,A)') Line
  do i=1,len_trim(Line)
    Line(i:i) = '-'
  end do
  write(LuWr,'(2X,A)') Line
  write(LuWr,'(2X,A,I4,A,I4)') 'mat. size = ',nRow,'x',nCol
end if
!----------------------------------------------------------------------*
! determine the printing format                                        *
!----------------------------------------------------------------------*
lFmt = len_trim(FmtIn)
if (lFmt /= 0) then
  FRMT = FmtIn
else
  Amax = -huge(Amax)
  Amin = huge(Amin)
  do j=1,nCol
    do i=1,nRow
      Amax = max(Amax,A(i,j))
      Amin = min(Amin,A(i,j))
    end do
  end do
  Scal = real(max(nRow,nCol),kind=wp)
  Amax = Amax*Amax*Scal
  Amin = Amin*Amin*Scal
  Pmax = Zero
  if (abs(Amax) > 1.0e-72_wp) Pmax = log10(abs(Amax))
  iPmax = int(One+Pmax)
  iPmax = max(1,iPmax)
  Pmin = Zero
  if (abs(Amin) > 1.0e-72_wp) Pmin = log10(abs(Amin))
  iPmin = int(One+Pmin)
  iPmin = max(1,iPmin)
  nDigit = 14
  nDecim = min(8,nDigit-max(iPmin,iPmax))
  if (Amax < Zero) iPmax = iPmax+1
  if (Amin < Zero) iPmin = iPmin+1
  lNumbr = max(iPmin,iPmax)+nDecim+2
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
  write(FRMT,'(A,I4.4,A,I4.4,A,I4.4,A)') '(2X,',nCols,'F',lItem,'.',nDecim,')'
end if
!----------------------------------------------------------------------*
! print the data                                                       *
!----------------------------------------------------------------------*
#ifdef _DEBUGPRINT_
write(LuWr,*)
write(LuWr,'(ES24.17)') DDot_(nCol*nRow,A,1,A,1),DDot_(nCol*nRow,A,1,[One],0)
#else
write(LuWr,*)
write(LuWr,'(2X,A)') 'row norms'
write(LuWr,FRMT) (DDot_(nCol,A(i,1),nRow,A(i,1),nRow),i=1,nRow)
write(LuWr,'(2X,A)') 'column norms'
write(LuWr,FRMT) (DDot_(nRow,A(1,i),1,A(1,i),1),i=1,nCol)
#endif

!----------------------------------------------------------------------*
! End procedure                                                        *
!----------------------------------------------------------------------*
return

end subroutine TrcPrt
