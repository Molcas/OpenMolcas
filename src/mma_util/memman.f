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
      Subroutine MemMan(ipOut,Ln_,Oper,ir,Label)
************************************************************
*
*   <DOC>
*     <Name>MemMan</Name>
*     <Syntax>Call MemMan(ipOut,Ln\_,Oper,ir,Label)</Syntax>
*     <Arguments>
*       \Argument{ipOut}{Position}{Integer}{out}
*       \Argument{Ln\_}{Length}{Integer}{in}
*       \Argument{Oper}{MARK $|$ ADDL $|$ ADDS $|$ FLUSM $|$ FREE}{Char*(5)}{in}
*       \Argument{ir}{1 $|$ 2}{Integer}{in}
*       \Argument{Label}{Arbitrary}{Char*(6)}{inout}
*     </Arguments>
*     <Purpose>
* Handle calls for memory from LUCIA, by passing the either to
* GETMEM or to the MEMMAN\_LUCIA routine.
*     </Purpose>
*     <Dependencies></Dependencies>
*     <Author></Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*     </Description>
*    </DOC>
*
************************************************************
      Implicit Real*8 (a-h,o-z)
      Character*(*) Oper
      Character*5 OP2
      Character*6 Label
      CHARACTER*4 IRL
      Integer ip(40)
      Logical Mark
      Common /Memilimem/ num,ip,Mark
      INTRINSIC LEN
#include "memman.fh"
*
      L=MIN(LEN(OPER),5)
      OP2=' '
      OP2(1:L)=Oper(1:L)

      Ln = Ln_
      If (ln.eq.0) ln=1
      If (ir.eq.1) IRL='REAL'
      If (ir.eq.2) IRL='REAL'
      If (OP2(1:4).eq.'MARK') Then
         Num_Marks = Num_Marks + 1
         If (Num_Marks .gt. MxpMarks) Then
            Write(6,*) 'MemMan: There is too many active marks.'
            Write(6,*) 'Please either increase MxpMarks or remove'
            Write(6,*) 'some of the existring marks.'
            Write(6,*)
            Write(6,*) 'MxpMarks, Num_Marks = ',MxpMarks, Num_Marks
            Write(6,*)
            write(6,*) 'Num       Label      iPos'
            Write(6,*) '=============================='
            Do i = 1, Num_Marks
               Write(6,'(1X,I3,5X,A8,5X,I9)')
     &               i,label_marks(i),ipos_marks(i)
            End Do
            Call Abend()
         End If
         Call GetMem(Label,'Allo','Real',iPos,1)
         Label_Marks(Num_Marks) = Label(1:6)
         iPos_Marks(Num_Marks)  = iPos
      End If
      If (OP2(1:4).eq.'ADDL' .or. OP2(1:4).eq.'ADDS') Then
         Call GetMem(label,'ALLO',IRL,ipout,Ln)
      End If
      If (OP2(1:5).eq.'FLUSM') Then
         iPos  = 0
         i_Save = 0
         Do i = Num_Marks, 1, -1
            If (iPos .eq. 0 .AND. Label_Marks(i) .eq. Label) Then
               iPos   = iPos_Marks(i)
               i_Save = i
            End If
         End Do
         If (iPos .eq. 0) Then
            Write(6,*) 'MemMan: You have specified a non-existing'
            Write(6,*) 'label. I can''t work with that.'
            Write(6,*)
            Write(6,*) 'Label = ',Label
            Call Abend()
         End If
         Do i = i_Save+1, Num_Marks
            iPos_Marks(i-1)  = iPos_Marks(i)
            Label_Marks(i-1) = Label_Marks(i)
         End Do
         Call GetMem(Label,'Flush','Real',iPos,1)
         Call GetMem(Label,'Free','Real',iPos,1)
         Num_Marks = Num_Marks - 1
      End If
      If (OP2(1:4).eq.'FREE') Then
         Write(6,*) 'Memman: Operation FREE not implemented yet.'
         call abend()
      End If
      Return
      End
