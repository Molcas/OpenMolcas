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
* Copyright (C) 2000-2016, Valera Veryazov                             *
************************************************************************
******************************************************************************
*                                                                            *
* Author:   Valera Veryazov 2000-2016                                        *
*           Theoretical Chemistry                                            *
*           Lund University                                                  *
*           Sweden                                                           *
*                                                                            *
******************************************************************************
       Subroutine StdIn_Name(Name)
       Character*(*) Name
       Character*132 Line
*
       nName=LEN(Name)
       If (nName.ne.16) Then
          Write (6,*) 'StdIn_Name: Wrong length of character Name'
          Call Abend()
       End If
C      Write (*,*) 'nName=',nName
*
       Name='Stdin.  '
       Call GetEnvf('EMIL_RC2',Line)
       Read (Line,'(I132.132)') Iter
c      Write (*,*) 'Line=',Line
c      Write (*,*) 'Iter=',Iter
       Iter = Iter+1
       If (Line(1:1).eq.' ') Then
          Name(7:7)='2'
       Else If (Iter.le.9) Then
          Write (Name(7:7),'(I1)') Iter
       Else If (Iter.le.99) Then
          Write (Name(7:8),'(I2)') Iter
       Else
          Write (6,*) 'StdIn_Name: Error in Line!'
          Call Abend()
       End If
       Line=' '
       Call GetEnvf('EMIL_InLoop',Line)
       ib=-1
       ie=-1
       i=1
10     continue
        if(Line(i:i).ne.' '.and. ib.eq.-1) then
           ib=i
           goto 20
        endif
        if(Line(i:i).eq.' '.and.ib.gt.0) then
           ie=i
           goto 30
        endif
20      i=i+1
        goto 10
30      Name(index(Name,' '):)='.'//Line(ib:ie)
c        write(*,*) '>',Line,'<', ib, ie
c       Write (6,*) 'StdIn=',Name
*
       Return
       End
