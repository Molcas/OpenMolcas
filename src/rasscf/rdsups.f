************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine RdSupS(LuInput,n,iBuff)
************************************************************************
*                                                                      *
*     Purpose:                                                         *
*     Read supersymmetry input.                                        *
*                                                                      *
************************************************************************
#include "output_ras.fh"
      Parameter (ROUTINE='RDSUPS  ')
      Integer iBuff(*)
      Integer is(288),ie(288)
      Character*288 Line
*----------------------------------------------------------------------*
*     Start procedure, initialize data counter                         *
*----------------------------------------------------------------------*
      n=0
      k=-1
*---  Read next line as a chatacter string  ---------------------------*
100   Read(LuInput,'(A)',End=900) Line
*---  Left adjust line  -----------------------------------------------*
      Call LeftAd(Line)
      If ( Line(1:1).eq.' ' ) Goto 100
      If ( Line(1:1).eq.'*' ) Goto 100
*---  Remove multiple intervening blanks  -----------------------------*
      Do i=1,287
        nRepeat=0
        Do While ( Line(i:i+1).eq.'  ' .and. nRepeat.lt.288 )
          nRepeat=nRepeat+1
          Line(i:287)=Line(i+1:288)
          Line(288:288)=' '
        End Do
      End Do
*---  Insert commas as the only valid separators  ---------------------*
      Do i=2,287
        If ( Line(i:i).eq.' ' ) then
          If ( Line(i-1:i-1).ne.' ' .and. Line(i-1:i-1).ne.',') then
            If ( Line(i+1:i+1).ne.' ' .and. Line(i+1:i+1).ne.',') then
              Line(i:i)=','
            End If
          End If
        End If
      End Do
*---  Get the last noblank character  ---------------------------------*
      iLast=0
      Do i=1,288
        If ( Line(i:i).ne.' ' ) iLast=i
      End Do
*---  Initialize markers  ---------------------------------------------*
      Do i=1,288
        is(i)=0
        ie(i)=0
      End Do
*---  Divide the line into substrings  --------------------------------*
      m=1
      is(m)=0
      Do i=1,iLast
        If ( Line(i:i).eq.',' ) then
          m=m+1
          is(m)=i
        End If
      End Do
      m=0
      Do i=1,iLast
        If ( Line(i:i).eq.',' ) then
          m=m+1
          ie(m)=i
        End If
      End Do
      If ( Line(iLast:iLast).ne.',' ) m=m+1
      ie(m)=iLast+1
*---  Read by substrings  ---------------------------------------------*
      Do i=1,m
        l=ie(i)-is(i)
        iz=0
        If ( l.gt.2 ) then
          Read(Line(is(i)+1:ie(i)-1),*,Err=910)iz
        Else If ( l.eq.2 .and. Line(is(i)+1:ie(i)-1).ne.' ' ) then
          Read(Line(is(i)+1:ie(i)-1),*,Err=910)iz
        End If
        If ( k.eq.-1 ) then
          k=k+1
          n=iz
        Else if ( k.lt.n ) then
          k=k+1
          iBuff(k)=iz
        End If
      End Do
*---  If necessary continue by next line  -----------------------------*
      If ( k.lt.n ) Goto 100
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
      Return
*----------------------------------------------------------------------*
*     Error exit                                                       *
*----------------------------------------------------------------------*
900   Write(LF,*)
      Write(LF,'(6X,A)') ' RASSCF was reading supersymmetry input from'
      Write(LF,'(6X,A)') 'the input file, when some error occurred,'
      Write(LF,'(6X,A)') 'probably end of file.'
      Write(LF,*)
      Call Quit_OnUserError()
910   Write(LF,*)
      Write(LF,'(6X,A)') ' RASSCF was reading supersymmetry input from'
      Write(LF,'(6X,A)') 'the input file, when some error occurred.'
      Write(LF,'(6X,A)') 'Some of the input data seems to be in error.'
      Write(LF,*)
      Call Quit_OnUserError()
      End
