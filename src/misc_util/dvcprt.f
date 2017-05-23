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
      Subroutine DVcPrt(Title,FmtIn,X,N)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Print a vector of double precision real numbers                  *
*                                                                      *
*     calling arguments                                                *
*     Title  : character string containing a title                     *
*              If the string is empty now title will be printed        *
*     X      : vector of double precision reals                        *
*     N      : dimension of vector X                                   *
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
      Character*(*) Title
      Character*(*) FmtIn
      Dimension X(N)
      Integer StrnLn
      Parameter (lPaper=120)
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
         Write(6,*)
         Write(6,'(2X,A)') Line
         Do 30 i=1,StrnLn(Line)
            Line(i:i)='-'
30       Continue
         Write(6,'(2X,A)') Line
         Write(6,'(2X,A,I6)') 'vec. size = ',N
      End If
*----------------------------------------------------------------------*
*     determine the printing format                                    *
*----------------------------------------------------------------------*
      lFmt=StrnLn(FmtIn)
      If ( lFmt.ne.0 ) then
         FMT=FmtIn
      Else
         Xmax=X(1)
         Xmin=X(1)
         Do 40 i=1,N
            Xmax=Max(Xmax,X(i))
            Xmin=Min(Xmin,X(i))
40       Continue
         Pmax=0
         If ( Abs(Xmax).gt.1.0D-72 ) Pmax=Log10(Abs(Xmax))
         iPmax=1+Int(Pmax)
         iPmax=Max(1,iPmax)
         Pmin=0
         If ( Abs(Xmin).gt.1.0D-72 ) Pmin=Log10(Abs(Xmin))
         iPmin=1+Int(Pmin)
         iPmin=Max(1,iPmin)
         nDigit=14
         nDecim=Min(8,nDigit-Max(iPmin,iPmax))
         If ( Xmax.lt.0 ) iPmax=iPmax+1
         If ( Xmin.lt.0 ) iPmin=iPmin+1
         lNumbr=Max(iPmin,iPmax)+nDecim+1
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
         Write(FMT,'(A,   I2.2,  A, I2.2,  A, I2.2,   A)')
     &             '(2X,',nCols,'F',lItem,'.',nDecim,')'
      End if
*----------------------------------------------------------------------*
*     print the data                                                   *
*----------------------------------------------------------------------*
      Write(6,*)
      Write(6,FMT)(X(i),i=1,N)
*----------------------------------------------------------------------*
*     End procedure                                                    *
*----------------------------------------------------------------------*
      Return
      End
