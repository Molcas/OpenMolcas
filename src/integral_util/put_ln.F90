!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine Put_ln(In_Line)
!***********************************************************************
! This function replaces function getln                                *
!                                                                      *
! It reads, broadcasts, and parses an input line                       *
!                                                                      *
! Blank lines or lines containing star (*) in column 1 are skipped     *
! Lines staring with exclaimation (!) in column 1 are skipped          *
!                                                                      *
! After this routine has been called, data can be retrieved using      *
! the subroutines                                                      *
!                                                                      *
!   Get_F(icol,array,n)  (for floating point values)                   *
!   Get_I(icol,iarry,n)  (for integer values)                          *
!   Get_S(icol,strgs,n)  (for character strings)                       *
!                                                                      *
! where icol is the first non-blank work to be taken, and n is the     *
! number of data.                                                      *
!                                                                      *
!***********************************************************************
      use getline_mod, only: Line, nCol, iStrt, iEnd
      implicit None
      Character(LEN=*) In_line

      Integer i, j, l, icom

      Line=In_Line
      l=len(line)
      Do i = 1, l
         if(ichar(line(i:i)).eq.9) line(i:i)=' '
         If (line(i:i).eq.';') line(i:l)=' '
      End Do
      ncol=0
      j=1
      Do
         icom=0
         do i=j,l
           if(line(i:i).eq.',') then
             icom=icom+1
             if(icom.ne.1) goto 20
           else
             if(line(i:i).ne.' ') goto 20
           end if
         end do
         Exit

 20      do j=i,l
           if(line(j:j).eq.' '.or.line(j:j).eq.',') goto 30
         end do
         j=l+1
 30      ncol=ncol+1
         istrt(ncol)=i
         iend(ncol)=j-1
      End Do

      End Subroutine Put_ln
