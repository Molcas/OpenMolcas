************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1992, Markus P. Fuelscher                              *
************************************************************************
      Subroutine TrcPrt(Title,FmtIn,A,nRow,nCol)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Print the row and column norms of a rectangular matrix           *
*                                                                      *
*     calling arguments                                                *
*     Title  : character string containing a title                     *
*              If the string is empty now title will be printed        *
*     A      : retangular matrix of double precision reals             *
*     nRow   : row dimension of matrix A                               *
*     nCol   : column dimension of matrix A                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. P. Fuelscher                                                  *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "standard_iounits.fh"
#include "real.fh"
      Character*(*) Title
      Character*(*) FmtIn
      Dimension A(nRow,nCol)
      Integer StrnLn
      Parameter ( lPaper = 70 )
      Character*(lPaper) Line
      Character*20 FMT
*----------------------------------------------------------------------*
*     print the title                                                  *
*----------------------------------------------------------------------*
      lTitle=StrnLn(Title)
      If ( lTitle.gt.0 ) then
         Do 10 i=1,lPaper
             Line(i:i)=' '
10       Continue
         lLeft=1
         Do 20 i=lTitle,1,-1
            If ( Title(i:i).ne.' ' ) lLeft=i
20       Continue
         lLeft=lLeft-1
         Do 25 i=1,lPaper
            If ( i+lLeft.le.lTitle ) Line(i:i)=Title(i+lLeft:i+lLeft)
25       Continue
         Write(LuWr,*)
         Write(LuWr,'(2X,A)') Line
         Do 30 i=1,StrnLn(Line)
            Line(i:i)='-'
30       Continue
         Write(LuWr,'(2X,A)') Line
         Write(LuWr,'(2X,A,I4,A,I4)') 'mat. size = ',nRow,'x',nCol
      End If
*----------------------------------------------------------------------*
*     determine the printing format                                    *
*----------------------------------------------------------------------*
      lFmt=StrnLn(FmtIn)
      If ( lFmt.ne.0 ) then
         FMT=FmtIn
      Else
         Amax=A(1,1)
         Amin=A(1,1)
         Do 40 j=1,nCol
            Do 50 i=1,nRow
               Amax=Max(Amax,A(i,j))
               Amin=Min(Amin,A(i,j))
50          Continue
40       Continue
         Scal=Dble(Max(nRow,nCol))
         Amax=Amax*Amax*Scal
         Amin=Amin*Amin*Scal
         Pmax=0d0
         if ( Abs(Amax).gt.1.0D-72 ) Pmax=Log10(Abs(Amax))
         iPmax=Int(1d0+Pmax)
         iPmax=Max(1,iPmax)
         Pmin=0d0
         If ( Abs(Amin).gt.1.0D-72 ) Pmin=Log10(Abs(Amin))
         iPmin=Int(1d0+Pmin)
         iPmin=Max(1,iPmin)
         nDigit=14
         nDecim=Min(8,nDigit-Max(iPmin,iPmax))
         If ( Amax.lt.0d0 ) iPmax=iPmax+1
         If ( Amin.lt.0d0 ) iPmin=iPmin+1
         lNumbr=Max(iPmin,iPmax)+nDecim+2
         nCols=10
         lLine=nCols*lNumbr
         If ( lLine.gt.lPaper ) then
            If ( lLine.le.lPaper+nCols .and. nDecim.gt.1 ) then
               nDecim=nDecim-1
               lNumbr=Max(iPmin,iPmax)+nDecim
               lItem=Max(lNumbr,lPaper/nCols)
            Else
               nCols=5
               lItem=Max(lNumbr,lPaper/nCols)
            End If
         Else
            lItem=lNumbr
         End If
         Write(FMT,'(A,   I4.4,  A, I4.4,  A, I4.4,   A)')
     &             '(2X,',nCols,'F',lItem,'.',nDecim,')'
      End if
*----------------------------------------------------------------------*
*     print the data                                                   *
*----------------------------------------------------------------------*
#ifdef _DEBUG_
       Write(LuWr,*)
       Write(LuWr,'(E24.17)') DDot_(nCol*nRow,A,1,A,1),
     &                        DDot_(nCol*nRow,A,1,[One],0)
#else
       Write(LuWr,*)
       Write(LuWr,'(2X,A)') 'row norms'
       Write(LuWr,FMT)(DDot_(nCol,A(i,1),nRow,A(i,1),nRow),i=1,nRow)
       Write(LuWr,'(2X,A)') 'column norms'
       Write(LuWr,FMT)(DDot_(nRow,A(1,i),1,A(1,i),1),i=1,nCol)
#endif
*----------------------------------------------------------------------*
*     End procedure                                                    *
*----------------------------------------------------------------------*
      Return
      End
