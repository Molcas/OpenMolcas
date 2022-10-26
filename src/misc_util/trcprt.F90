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

implicit real*8(A-H,O-Z)
#include "standard_iounits.fh"
#include "real.fh"
character*(*) Title
character*(*) FmtIn
dimension A(nRow,nCol)
integer StrnLn
parameter(lPaper=70)
character*(lPaper) Line
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
  write(LuWr,*)
  write(LuWr,'(2X,A)') Line
  do i=1,StrnLn(Line)
    Line(i:i) = '-'
  end do
  write(LuWr,'(2X,A)') Line
  write(LuWr,'(2X,A,I4,A,I4)') 'mat. size = ',nRow,'x',nCol
end if
!----------------------------------------------------------------------*
! determine the printing format                                        *
!----------------------------------------------------------------------*
lFmt = StrnLn(FmtIn)
if (lFmt /= 0) then
  FMT = FmtIn
else
  Amax = A(1,1)
  Amin = A(1,1)
  do j=1,nCol
    do i=1,nRow
      Amax = max(Amax,A(i,j))
      Amin = min(Amin,A(i,j))
    end do
  end do
  Scal = dble(max(nRow,nCol))
  Amax = Amax*Amax*Scal
  Amin = Amin*Amin*Scal
  Pmax = 0d0
  if (abs(Amax) > 1.0D-72) Pmax = log10(abs(Amax))
  iPmax = int(1d0+Pmax)
  iPmax = max(1,iPmax)
  Pmin = 0d0
  if (abs(Amin) > 1.0D-72) Pmin = log10(abs(Amin))
  iPmin = int(1d0+Pmin)
  iPmin = max(1,iPmin)
  nDigit = 14
  nDecim = min(8,nDigit-max(iPmin,iPmax))
  if (Amax < 0d0) iPmax = iPmax+1
  if (Amin < 0d0) iPmin = iPmin+1
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
  write(FMT,'(A,   I4.4,  A, I4.4,  A, I4.4,   A)') '(2X,',nCols,'F',lItem,'.',nDecim,')'
end if
!----------------------------------------------------------------------*
! print the data                                                       *
!----------------------------------------------------------------------*
#ifdef _DEBUGPRINT_
write(LuWr,*)
write(LuWr,'(E24.17)') DDot_(nCol*nRow,A,1,A,1),DDot_(nCol*nRow,A,1,[One],0)
#else
write(LuWr,*)
write(LuWr,'(2X,A)') 'row norms'
write(LuWr,FMT) (DDot_(nCol,A(i,1),nRow,A(i,1),nRow),i=1,nRow)
write(LuWr,'(2X,A)') 'column norms'
write(LuWr,FMT) (DDot_(nRow,A(1,i),1,A(1,i),1),i=1,nCol)
#endif

!----------------------------------------------------------------------*
! End procedure                                                        *
!----------------------------------------------------------------------*
return

end subroutine TrcPrt
