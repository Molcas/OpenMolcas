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
* Copyright (C) 2006, Thomas Bondo Pedersen                            *
************************************************************************
      Logical Function ChoMP2_ChkPar()
C
C     Thomas Bondo Pedersen, Nov. 2006.
C
C     Purpose: return true if parallel run (nProcs>1), else return
C              false.
C
#if defined (_MOLCAS_MPP_)
      Use Para_Info, Only: nProcs, Is_Real_Par
#endif
      Implicit None
#if defined (_MOLCAS_MPP_)
      ChoMP2_ChkPar=nProcs.gt.1 .and. Is_Real_Par()
#else
      ChoMP2_ChkPar=.False.
#endif
      End
