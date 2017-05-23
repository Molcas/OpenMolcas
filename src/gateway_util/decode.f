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
* Copyright (C) Bjorn O. Roos                                          *
*               1991, Roland Lindh                                     *
************************************************************************
      Subroutine Decode(LBL,string,N,Hit)
************************************************************************
* Object: to find the character string 'string' between                *
*         dots N-1 and N in the label LBL of length lLBL               *
*         blanks in string are removed                                 *
* Called from: Rdbsl                                                   *
* Subroutine calls: none                                               *
*                                                                      *
* Author: Bjoern Roos, University of Lund, Sweden                      *
************************************************************************
      Character*(*) LBL,string
      Character*80 xstring
      Character*1 dot
      Logical Hit
      Data dot/'.'/
*
*     Call qEnter('Decode')
*     write(6,'(1x,a)') LBL
      i1=1
      idot=0
      lstring=0
      lLBL=LEN(LBL)
      Do 10 i=1,lLBL
       if(LBL(i:i).ne.dot) go to 10
       idot=idot+1
       if(idot.eq.N-1) i1=i+1
       if(idot.eq.N) then
        xstring=' '
*       Write (*,'(1x,A,/,1X,A)') ' xstring=',xstring
        if(i.gt.i1) xstring=LBL(i1:i-1)
*       Write (*,'(1x,A,/,1X,A)') ' xstring=',xstring
        lstring=i-i1
        go to 20
       Endif
   10 Continue
      If (Hit) Then
         Call WarningMessage(2,'Decode: error in basis set label')
         Write (6,'(A,A)')'LBL=',LBL
         Call Abend()
      Else
*        Call qExit('Decode')
         Return
      End If
*
*     Pack the string
*
   20 Continue
      Hit=.True.
      i1=0
      string=' '
      Do 30 i=1,lstring
       If(xstring(i:i).eq.' ') go to 30
       i1=i1+1
       string(i1:i1) =xstring(i:i)
   30 continue
*     Call qExit('Decode')
      Return
      End
