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
* Copyright (C) Jesper Wisborg Krogh                                   *
************************************************************************
*  Set_Do_Parallel
*
*> @brief
*>   Set the variable \c Do_Real_Par (in the \c do_sync common block) to the value of \p Status
*> @author Jesper Wisborg Krogh
*>
*> @details
*> Set the variable \c Do_Real_Par (in the \c do_sync common block) to the value of \p Status.
*> Setting this variable to ``.False.`` will ensure that the code afterwards will be run
*> in serial mode even if the job is started with ``$MOLCAS_NPROCS > 1``, i.e. all tasks are done
*> on each node, and no synchronization. Setting \c Do_Real_Par to ``.True.`` will cause the
*> code to run in parallel again. The variable is initially set in the subroutine
*> ::start, where the value is ``.True.``.
*>
*> @param[in] Status Variable with the new value for \c Do_Real_Par
************************************************************************
      Subroutine Set_Do_Parallel(Status)
      Implicit None
      Logical Status
#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#include "mpp_info.fh"
      mpp_workshare = Status .and. (mpp_nprocs.gt.1)
      If (mpp_workshare) Then
        myRank = mpp_procid
        nProcs = mpp_nprocs
      Else
        myRank = 0
        nProcs = 1
      End If
#else
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(Status)
#endif
      Return
      End
