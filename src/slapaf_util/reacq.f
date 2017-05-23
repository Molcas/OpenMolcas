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
* Copyright (C) 2009, Roland Lindh                                     *
************************************************************************
      Subroutine ReacQ(V_X,nX,V_Q,nQ)
************************************************************************
*                                                                      *
*     Objective: Transform the "reaction vector" from cartesian        *
*                coordinates to Internal.                              *
*                                                                      *
*     Roland Lindh, Dept. of Theor. Chem., Lund University, Sweden     *
*     2009                                                             *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "info_slapaf.fh"
      Real*8 V_Q(nQ), V_X(nX)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      Call RecPrt('BMx',' ',Work(ipB),nX,nQ)
      Call RecPrt('V_X',' ',V_X,nX,1)
#endif
*
      M = nX
      N = nQ
      NRHS=1
      Call Eq_Solver('T',M,N,NRHS,Work(ipB),.FALSE.,Degen,V_X,V_Q)
*
#ifdef _DEBUG_
      Call RecPrt('V_Q',' ',V_Q,nQ,1)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
