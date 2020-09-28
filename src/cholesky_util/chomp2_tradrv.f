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
* Copyright (C) 2004, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine ChoMP2_TraDrv(irc,CMO,Diag,DoDiag)
C
C     Thomas Bondo Pedersen, Dec. 2004.
C
C     Purpose: AO-to-MO (ai) transformation of Cholesky vectors
C              performed directly in reduced sets. This assumes
C              that the MP2 program has been appropriately initialized.
C
#include "implicit.fh"
      Real*8  CMO(*), Diag(*)
      Logical DoDiag
#include "cholesky.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"

      Character*6  ThisNm
      Character*13 SecNam
      Parameter (SecNam = 'ChoMP2_TraDrv', ThisNm = 'TraDrv')

      irc = 0

C     Reorder MO coefficients.
C     ------------------------

      l_COcc = nT1AOT(1)
      l_CVir = nAOVir(1)
      Call GetMem('COcc','Allo','Real',ip_COcc,l_COcc)
      Call GetMem('CVir','Allo','Real',ip_CVir,l_CVir)
      Call ChoMP2_MOReOrd(CMO,Work(ip_COcc),Work(ip_CVir))

C     Transform vectors.
C     ------------------

      Call ChoMP2_Tra(Work(ip_COcc),Work(ip_CVir),Diag,DoDiag)

C     Deallocate reordered MO coefficients.
C     -------------------------------------

      Call GetMem('CVir','Free','Real',ip_CVir,l_CVir)
      Call GetMem('COcc','Free','Real',ip_COcc,l_COcc)

      End
