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
!***********************************************************************
! CRecPrt
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
!> @param[in] Type
!***********************************************************************

subroutine CRecPrt(Title,FmtIn,A,nRow,nCol,type)

implicit real*8(A-H,O-Z)
#include "standard_iounits.fh"
character*(*) Title
character*(*) FmtIn
character*1 type
complex*16 A(nRow,nCol)
integer StrnLn
parameter(lPaper=120,lMaxTitle=60)
character*(lMaxTitle) Line
character*20 FMT

!----------------------------------------------------------------------*
if (nRow*nCol == 0) return
#ifdef _DEBUGPRINT_
!call CTrcPrt(Title,FmtIn,A,nRow,nCol,'R')
!call CTrcPrt(Title,FmtIn,A,nRow,nCol,'I')
return
#endif
!----------------------------------------------------------------------*
! print the title                                                      *
!----------------------------------------------------------------------*
lTitle = StrnLn(Title)
if (lTitle > 0) then
  do i=1,lMaxTitle
    Line(i:i) = ' '
  end do
  lLeft = 1
  do i=lTitle,1,-1
    if (Title(i:i) /= ' ') lLeft = i
  end do
  lLeft = lLeft-1
  do i=1,lMaxTitle
    if (i+lLeft <= lTitle) Line(i:i) = Title(i+lLeft:i+lLeft)
  end do
  write(LuWr,*)
  if (type == 'R') then
    write(LuWr,'(2X,A,A)') Line,' Real Component'
  else
    write(LuWr,'(2X,A,A)') Line,' Imaginary Component'
  end if
  !do i=1,StrnLn(Line)
  !  Line(i:i) = '-'
  !end do
  !write(LuWr,'(2X,A)') Line
  write(LuWr,'(2X,A,I5,A,I5)') 'mat. size = ',nRow,'x',nCol
end if
!----------------------------------------------------------------------*
! determine the printing format                                        *
!----------------------------------------------------------------------*
lFmt = Strnln(FmtIn)
if (lFmt /= 0) then
  FMT = FmtIn
else
  if (type == 'R') then
    Amax = dble(A(1,1))
    Amin = dble(A(1,1))
    do j=1,nCol
      do i=1,nRow
        Amax = max(Amax,dble(A(i,j)))
        Amin = min(Amin,dble(A(i,j)))
      end do
    end do
  else
    Amax = DIMAG(A(1,1))
    Amin = DIMAG(A(1,1))
    do j=1,nCol
      do i=1,nRow
        Amax = max(Amax,DIMAG(A(i,j)))
        Amin = min(Amin,DIMAG(A(i,j)))
      end do
    end do
  end if
  Pmax = 0.0d0
  if (abs(Amax) > 1.0D-72) Pmax = log10(abs(Amax))
  iPmax = 1+int(Pmax)
  iPmax = max(1,iPmax)
  Pmin = 0.0d0
  if (abs(Amin) > 1.0D-72) Pmin = log10(abs(Amin))
  iPmin = 1+int(Pmin)
  iPmin = max(1,iPmin)
  nDigit = 15
  nDecim = min(9,nDigit-max(iPmin,iPmax))
  nDecim = max(nDecim,1)
  if (Amax < 0.0d0) iPmax = iPmax+1
  if (Amin < 0.0d0) iPmin = iPmin+1
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
  write(FMT,'(A,   I4.4,  A, I4.4,  A, I4.4,   A)') '(2X,',nCols,'F',lItem,'.',nDecim,')'
end if
!----------------------------------------------------------------------*
! print the data                                                       *
!----------------------------------------------------------------------*
!write(LuWr,*)
if (type == 'R') then
  do i=1,nRow
    write(LuWr,FMT) (dble(A(i,j)),j=1,nCol)
  end do
else
  do i=1,nRow
    write(LuWr,FMT) (DIMAG(A(i,j)),j=1,nCol)
  end do
end if

!----------------------------------------------------------------------*
! End procedure                                                        *
!----------------------------------------------------------------------*
return

end subroutine CRecPrt
