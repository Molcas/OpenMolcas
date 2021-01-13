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
      SubRoutine Cho_ParConf(Fake)
C
C     Purpose: set parallel info used in Cholesky decomposition of
C              two-electron integrals. If Fake=.True., run the
C              decomposition in serial and distribute vectors over nodes
C              at the end (no distribution is actually done, of course;
C              the relevant vectors are simply kept on the relevant
C              nodes).
C
      Use Para_Info, Only: Is_Real_Par
      Implicit None
      Logical Fake
#include "cho_para_info.fh"

      Cho_Real_Par = Is_Real_Par() .and. (.not.Fake)

      End
      SubRoutine Cho_GASync()
      Implicit None
#include "cho_para_info.fh"

      If (Cho_Real_Par) Call GASync

      End
      SubRoutine Cho_Init_Tsk(ID,n)
      Implicit None
      Integer ID, n
#include "cho_para_info.fh"

      If (Cho_Real_Par) Then
         Call Init_Tsk(ID,n)
      Else
         Call Cho_Init_Tsk_(ID,n)
      End If

      End
      SubRoutine Cho_Init_Tsk_(ID,n)
      Use Para_Info, Only: Is_Real_Par, Set_Do_Parallel
      Implicit None
      Integer ID, n
#include "cho_para_info.fh"
      Logical Parallel_Mode

      Parallel_Mode = Is_Real_Par()

      Call Set_Do_Parallel(.False.)
      Call Init_Tsk(ID,n)
      Call Set_Do_Parallel(Parallel_Mode)

      End
      SubRoutine Cho_Free_Tsk(ID)
      Implicit None
      Integer ID
#include "cho_para_info.fh"

      If (Cho_Real_Par) Then
         Call Free_Tsk(ID)
      Else
         Call Cho_Free_Tsk_(ID)
      End If

      End
      SubRoutine Cho_Free_Tsk_(ID)
      Use Para_Info, Only: Is_Real_Par, Set_Do_Parallel
      Implicit None
      Integer ID
#include "cho_para_info.fh"
      Logical Parallel_Mode

      Parallel_Mode = Is_Real_Par()

      Call Set_Do_Parallel(.False.)
      Call Free_Tsk(ID)
      Call Set_Do_Parallel(Parallel_Mode)

      End
      Logical Function Cho_Rsv_Tsk(ID,i)
      Implicit None
      Integer  ID, i
#include "cho_para_info.fh"
      Logical  Rsv_Tsk
      External Rsv_Tsk
      Logical  Cho_Rsv_Tsk_
      External Cho_Rsv_Tsk_

      If (Cho_Real_Par) Then
         Cho_Rsv_Tsk = Rsv_Tsk(ID,i)
      Else
         Cho_Rsv_Tsk = Cho_Rsv_Tsk_(ID,i)
      End If

      End
      Logical Function Cho_Rsv_Tsk_(ID,i)
      Use Para_Info, Only: Is_Real_Par, Set_Do_Parallel
      Implicit None
      Integer ID, i
#include "cho_para_info.fh"
      Logical Parallel_Mode
      Logical  Rsv_Tsk
      External Rsv_Tsk

      Parallel_Mode = Is_Real_Par()

      Call Set_Do_Parallel(.False.)
      Cho_Rsv_Tsk_ = Rsv_Tsk(ID,i)
      Call Set_Do_Parallel(Parallel_Mode)

      End
      SubRoutine Cho_GAdGOp(X,n,Op)
      Implicit None
      Integer n, iv, kv
      Real*8  X(n)
      Character*(*) Op
#include "cho_para_info.fh"

      If (Cho_Real_Par) Then
         iv=0
         do while (iv .lt. n)
            kv=Min(n-iv,32000000)
            Call GAdGOp(X(iv+1),kv,Op)
            iv=iv+kv
         end do
      EndIf

      End
      SubRoutine Cho_GAiGOp(X,n,Op)
      Implicit None
      Integer n, iv, kv
      Integer X(n)
      Character*(*) Op
#include "cho_para_info.fh"

      If (Cho_Real_Par) Then
         iv=0
         do while (iv .lt. n)
            kv=Min(n-iv,32000000)
            Call GAiGOp(X(iv+1),kv,Op)
            iv=iv+kv
         end do
      EndIf

      End
