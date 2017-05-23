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
      Subroutine Merge_Constraints(FileIn1,FileIn2,FileOut,
     &                             nLambda,iRow_c)
      Implicit None
      Character(Len=*), Intent(in) :: FileIn1, FileIn2, FileOut
      Integer, Intent(inout) :: nLambda, iRow_c
      Integer :: Lu1, Lu2, Lu3, isFreeUnit, AixRm
      Character(Len=180) :: Line, Get_Ln
      Character(Len=6) :: Tag
      Logical :: Found, AuxFile
      External :: Get_Ln, isFreeUnit, AixRm
      Integer :: i, j, Lu, iErr

* Open the input files if they exist
      Lu1=0
      If (FileIn1.ne.'') Then
        Call f_Inquire(FileIn1,Found)
        If (Found) Then
          Lu1=isFreeUnit(30)
          Call Molcas_Open(Lu1,FileIn1)
        End If
      End If

      Lu2=0
      If ((FileIn2.ne.'').and.(FileIn2.ne.FileIn1)) Then
        Call f_Inquire(FileIn2,Found)
        If (Found) Then
          Lu2=isFreeUnit(30)
          Call Molcas_Open(Lu2,FileIn2)
        End If
      End If

* Set the output file as a third file, or a temporary one
      If (((FileOut.eq.FileIn1).and.(Lu1.ne.0)).or.
     &    ((FileOut.eq.FileIn2).and.(Lu2.ne.0))) Then
        AuxFile=.True.
      Else
        AuxFile=.False.
      End If

      nLambda=0
      iRow_c=0

* If there are no input files, only delete the output file if it exists
      If ((Lu1.eq.0).and.(Lu2.eq.0)) Then
        If (.not.AuxFile) Then
          Call f_Inquire(FileOut,Found)
          If (Found) iErr = AixRm(FileOut)
        End If
        Return
* Otherwise, open it
      Else
        Lu3=isFreeUnit(30)
        If (AuxFile) Then
          Call Molcas_Open(Lu3,'purge')
        Else
          Call Molcas_Open(Lu3,FileOut)
        End If
      End If

* Copy the constraints to the output file, counting the lines
* Copy only from existing files
      Do i=1,2
        If (i.eq.1) Tag='VALUES'
        If (i.eq.2) Tag='END'

        Do j=1,2
          If (j.eq.1) Lu=Lu1
          If (j.eq.2) Lu=Lu2
          If (Lu.eq.0) Cycle
          Line=AdjustL(Get_Ln(Lu))
          Call UpCase(Line)
          Do While (Line(1:4).ne.Tag(1:4))
            If (Index(Line,'&').eq.0) Then
              iRow_c=iRow_c+1
              If (i.eq.2) nLambda=nLambda+1
              If (Index(Line,'=').eq.0) Call FixEqualSign(Line,Lu)
            End If
            Write(Lu3,100) Trim(Line)
            Line=AdjustL(Get_Ln(Lu))
            Call UpCase(Line)
          End Do
        End Do

        Write(Lu3,100) Trim(Tag)
      End Do
      iRow_c=iRow_c+1

* Close files, and copy temporary if used
      If (Lu1.ne.0) Close (Lu1)
      If (Lu2.ne.0) Close (Lu2)
      If (Lu3.ne.0) Close (Lu3)
      If (AuxFile) Call fCopy('purge',FileOut,iErr)

100   Format (A)
      End Subroutine
