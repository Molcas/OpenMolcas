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
* SubRecPrt
*
*> @brief
*>   Print a triangular matrix
*> @author M. P. F&uuml;lscher, Lund, 1992
*>
*> @details
*> Print a square matrix stored in packed, lower triangular storage mode.
*>
*> @param[in] Title   String containing a title
*> @param[in] FmtIn   String containing a format
*> @param[in] A       Triangular matrix to be printed
*> @param[in] N       Dimension of matrix \p A
************************************************************************
      Subroutine TriPrt(Title,FmtIn,A,N)
      Implicit Real*8 (A-H,O-Z)
#include "standard_iounits.fh"
      Character*(*) Title
      Character*(*) FmtIn
      Dimension A(N*(N+1)/2)
      Integer StrnLn
      Parameter (lPaper=120)
      Character*(lPaper) Line
      Character*20 FMT
*----------------------------------------------------------------------*
      If (N.le.0) Return
*----------------------------------------------------------------------*
#ifdef _DEBUGPRINT_
      Call TrcPrt(Title,FmtIn,A,1,N*(N+1)/2)
      Return
#endif
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
c         Do 30 i=1,StrnLn(Line)
c            Line(i:i)='-'
c30       Continue
c         Write(LuWr,'(2X,A)') Line
         Write(LuWr,'(2X,A,I5,A,I5)') 'mat. size = ',N,'x',N
      End If
*----------------------------------------------------------------------*
*     determine the printing format                                    *
*----------------------------------------------------------------------*
      lFmt=StrnLn(FmtIn)
      If ( lFmt.ne.0 ) then
         FMT=FmtIn
      Else
         Amax=A(1)
         Amin=A(1)
         Do 40 i=1,N*(N+1)/2
            Amax=Max(Amax,A(i))
            Amin=Min(Amin,A(i))
40       Continue
         If (Amax.ne.0.0D0) Then
           Pmax=Log10(Abs(Amax))
           iPmax=Int(1d0+Pmax)
           iPmax=Max(1,iPmax)
         Else
           iPmax=1
         End If
         If (Amin.ne.0.0D0) Then
           Pmin=Log10(Abs(Amin))
           iPmin=Int(1d0+Pmin)
           iPmin=Max(1,iPmin)
         Else
           iPmin=1
         End If
         nDigit=24
         nDecim=Min(16,ABS(nDigit-Max(iPmin,iPmax)))
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
      Write(LuWr,*)
      jEnd=0
      Do 60 i=1,N
         jStart=jEnd+1
         jEnd=jEnd+i
         Write(LuWr,FMT)(A(j),j=jStart,jEnd)
60    Continue
*----------------------------------------------------------------------*
*     End procedure                                                    *
*----------------------------------------------------------------------*
      Return
      End
