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
! TriRecPrt
!
!> @brief
!>   Print a triangular matrix
!> @author M. P. F&uuml;lscher, Lund, 1992
!>
!> @details
!> Print a square matrix stored in packed, lower triangular storage mode.
!>
!> @param[in] Title   String containing a title
!> @param[in] FmtIn   String containing a format
!> @param[in] A       Triangular matrix to be printed
!> @param[in] N       Dimension of matrix \p A
!***********************************************************************

subroutine TriPrt(Title,FmtIn,A,N)

use Index_Functions, only: nTri_Elem
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
character(len=*), intent(in) :: Title, FmtIn
integer(kind=iwp), intent(in) :: N
real(kind=wp), intent(in) :: A(nTri_Elem(N))
#include "standard_iounits.fh"
integer(kind=iwp), parameter :: lPaper = 120
integer(kind=iwp) :: i, iPmax, iPmin, jEnd, jStart, lFmt, lItem, lLeft, lLine, lNumbr, lTitle, nCols, nDecim, nDigit
real(kind=wp) :: Amax, Amin, Pmax, Pmin
character(len=lPaper) :: Line
character(len=20) :: FRMT

!----------------------------------------------------------------------*
if (N <= 0) return
!----------------------------------------------------------------------*
#ifdef _DEBUGPRINT_
call TrcPrt(Title,FmtIn,A,1,nTri_Elem(N))
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
  do i=1,lPaper
    if (i+lLeft <= lTitle) Line(i:i) = Title(i+lLeft:i+lLeft)
  end do
  write(LuWr,*)
  write(LuWr,'(2X,A)') Line
  !do i=1,len_trim(Line)
  !  Line(i:i) = '-'
  !end do
  !write(LuWr,'(2X,A)') Line
  write(LuWr,'(2X,A,I5,A,I5)') 'mat. size = ',N,'x',N
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
  do i=1,nTri_Elem(N)
    Amax = max(Amax,A(i))
    Amin = min(Amin,A(i))
  end do
  if (Amax /= Zero) then
    Pmax = log10(abs(Amax))
    iPmax = int(One+Pmax)
    iPmax = max(1,iPmax)
  else
    iPmax = 1
  end if
  if (Amin /= Zero) then
    Pmin = log10(abs(Amin))
    iPmin = int(One+Pmin)
    iPmin = max(1,iPmin)
  else
    iPmin = 1
  end if
  nDigit = 24
  nDecim = min(16,abs(nDigit-max(iPmin,iPmax)))
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
write(LuWr,*)
jEnd = 0
do i=1,N
  jStart = jEnd+1
  jEnd = jEnd+i
  write(LuWr,FRMT) A(jStart:jEnd)
end do

!----------------------------------------------------------------------*
! End procedure                                                        *
!----------------------------------------------------------------------*
return

end subroutine TriPrt
