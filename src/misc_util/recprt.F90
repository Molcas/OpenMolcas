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
! RecPrt
!
!> @brief
!>   Write out a matrix on standard output
!> @author M. P. F&uuml;lscher, Lund, 1992
!>
!> @details
!> The matrix \p A of dimension \p nRow &times; \p nCol is printed
!> in output preceded by the character line \p Title. Format of the numerical
!> output is given by \p FmtIn. If \p FmtIn = ``''`` the utility will decide on format
!> for optimal output.
!>
!> @param[in] Title   Title card
!> @param[in] FmtIn   Format statement
!> @param[in] A       A matrix
!> @param[in] nRow    number of rows of \p A
!> @param[in] nCol    number of columns of \p A
!***********************************************************************

subroutine RecPrt(Title,FmtIn,A,nRow,nCol)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
character(len=*), intent(in) :: Title, FmtIn
integer(kind=iwp), intent(in) :: nRow, nCol
real(kind=wp), intent(in) :: A(nRow,nCol)
#include "standard_iounits.fh"
integer(kind=iwp), parameter :: lMaxTitle = 60, lPaper = 120
integer(kind=iwp) :: i, iPmax, iPmin, j, lFmt, lItem, lLeft, lLine, lNumbr, lTitle, nCols, nDecim, nDigit
real(kind=wp) :: Amax, Amin, Pmax, Pmin
character(len=lMaxTitle) :: Line
character(len=20) :: FRMT

!----------------------------------------------------------------------*
if (nRow*nCol == 0) return
#ifdef _DEBUGPRINT_
call TrcPrt(Title,FmtIn,A,nRow,nCol)
return
#endif
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
  do i=1,lMaxTitle
    if (i+lLeft <= lTitle) Line(i:i) = Title(i+lLeft:i+lLeft)
  end do
  write(LuWr,*)
  write(LuWr,'(2X,A)') Line
  !do i=1,len_trim(Line)
  !  Line(i:i) = '-'
  !end do
  !write(LuWr,'(2X,A)') Line
  write(LuWr,'(2X,A,I5,A,I5)') 'mat. size = ',nRow,'x',nCol
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
  Pmax = Zero
  if (abs(Amax) > 1.0e-72_wp) Pmax = log10(abs(Amax))
  iPmax = 1+int(Pmax)
  iPmax = max(1,iPmax)
  Pmin = Zero
  if (abs(Amin) > 1.0e-72_wp) Pmin = log10(abs(Amin))
  iPmin = 1+int(Pmin)
  iPmin = max(1,iPmin)
  nDigit = 24
  nDecim = min(16,nDigit-max(iPmin,iPmax))
  nDecim = max(nDecim,1)
  if (Amax < Zero) iPmax = iPmax+1
  if (Amin < Zero) iPmin = iPmin+1
  lNumbr = max(iPmin,iPmax)+nDecim+2
  nCols = 9
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
!write(LuWr,*)
do i=1,nRow
  write(LuWr,FRMT) A(i,1:nCol)
end do
!----------------------------------------------------------------------*
! End procedure                                                        *
!----------------------------------------------------------------------*
Call xFlush(LuWr)
return

end subroutine RecPrt
