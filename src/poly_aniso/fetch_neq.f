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
      Subroutine fetch_neq(nneq,neq,nexch)
      Implicit None

      Integer, intent(inout) :: nneq
      Integer, intent(inout):: neq(nneq), nexch(nneq)
      !local variables:
      Integer       :: i, Input, LineNr
      Character(72) :: LINE
      Logical       :: ab_initio_all
      Logical       :: DBG
      DBG=.false.
      If(DBG) Write(6,'(A)') 'Enter fetch_neq'
      If(DBG) Write(6,'(A,i3)') 'fetch_neq:  nneq=',nneq

      neq=0
      nexch=0
C=========== End of default settings====================================
      Input = 5
      Rewind(Input)
50    READ(Input,'(A72)',End=998) LINE
      If(DBG) write(6,'(A)') LINE
      Call NORMAL(LINE)
      If(LINE(1:5).ne.'&POLY') Go To 50
      LINENR=0

100   Call xFlush(6)
      READ(Input,'(A72)',End=998) LINE
      LINENR=LINENR+1
      Call NORMAL(LINE)
      If(LINE(1:1).eq.'*')   Go To 100
      If(LINE     .eq.' ')   Go To 100
      If(LINE(1:3).eq.'END') Go To 210 !End

      If (LINE(1:4).eq.'NNEQ') Then
!        number of non-equivalent centers; type of all centers
         READ(Input,*,ERR=997)  NNEQ, ab_initio_all
         If(DBG) Write(6,'(A,i4,A,L2)') 'NNEQ=', NNEQ,
     &                          ' ab_initio_all=',ab_initio_all
!        number of equivalent centers of type "i"
         READ(Input,*,ERR=997) (NEQ(i)  ,i=1,Nneq)
         If(DBG) Write(6,'(A,100I4)') 'NEQ(I)=',(NEQ(i),i=1,nneq)
!        number of RASSI wf for exchange
         READ(Input,*,ERR=997) (Nexch(i),i=1,Nneq)
         If(DBG) Write(6,'(A,100I4)') 'NExch(I)=',(NExch(i),i=1,nneq)
         Go To 210
      End If
      Go To 100

210   Continue



      GoTo 999 ! end

c-----------------------------------------------------------
997   continue
      Write(6,*)' READIN: Error reading "poly_aniso.input" '
      Write(6,*)' near line nr.',LINENR+1
      Go To 999
998   continue
      Write(6,*)' READIN: Unexpected End of input file.'
c-----------------------------------------------------------
999   Continue
      If(DBG) Write(6,'(A)') 'Exit fetch_neq'
      Return
      End Subroutine fetch_neq
