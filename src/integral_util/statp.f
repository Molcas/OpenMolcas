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
      Subroutine  StatP(Index)
      Implicit Real*8 (a-h,o-z)
#include "pstat.fh"
#include "print.fh"
*
      iRout = 10
      iPrint = nPrint(iRout)
      If (Index.eq.0) Then
         Call GetMem('PSOAO0','MAX','Real',iDum,MaxMem)
         MinXtr = MaxMem
      Else
         If (iPrint.ge.6) Then
            Write (6,*)
            Write (6,'(21X,A)') '******* Partitioning Ratios *******'
            Write (6,'(21X,A)') '* Index  i     j     k     l      *'
            Write (6,'(21X,A7,4F6.3,A4)') '* Cont.',
     &                r1/DBLE(iTotal),r2/DBLE(iTotal),
     &                r3/DBLE(iTotal),r4/DBLE(iTotal),'   *'
            Write (6,'(21X,A7,4F6.3,A4)') '* Prim.',
     &                q1/DBLE(iTotal),q2/DBLE(iTotal),
     &                q3/DBLE(iTotal),q4/DBLE(iTotal),'   *'
            Write (6,'(21X,A)') '***********************************'
            Write (6,*)
            Write (6,'(21X,A,I8)') ' Largest Memory Deficiency:',MaxReq
            Write (6,'(21X,A,I8)') ' Least Overflow of Memory :',MinXtr
            Write (6,'(21X,A,I8)') ' Max Available Memory     :',MaxMem
         End If
      End If
*
      Return
      End
