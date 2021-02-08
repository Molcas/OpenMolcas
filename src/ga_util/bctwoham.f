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
* Copyright (C) 1995, Martin Schuetz                                   *
*               1998, Roland Lindh                                     *
************************************************************************
************************************************************************
* This Module contains subroutines and functions which interface calls *
* to the Global Array Tools (GA)                                       *
*  DISTRIBUTED DATA PARALLEL VERSION for SCF                           *
*                                                                      *
* SubRoutine BCTwoHam(TwoHam,nDens,TCPU,TWall)                         *
* -> partial sum of Fock matrices                                      *
*    TwoHam: TwoHam(nDens): actual Fock matrix                         *
*----------------------------------------------------------------------*
*     written by:                                                      *
*     M. Schuetz                                                       *
*     University of Lund, Sweden, 1995                                 *
*                                                                      *
*     modified by:                                                     *
*     R. Lindh                                                         *
*     University of Lund, Sweden, 1998                                 *
************************************************************************
*
*
*----------------------------------------------------------------------*
      SubRoutine BCTwoHam(TwoHam,nDens,TCPU,TWall)
#ifdef _MOLCAS_MPP_
      Use Para_Info, Only: Is_Real_Par
#endif
      Implicit Real*8 (a-h,o-z)
      Real*8 TwoHam(nDens)
#ifdef _MOLCAS_MPP_
#include "real.fh"
*
#include "global.fh"
*
*
      If (.Not. Is_Real_Par()) Return
      Call CWTime(TCpu1,TWall1)
      Call GADGOP(TwoHam,nDens,'+')
      Call CWTime(TCpu2,TWall2)
      TCPU=TCpu2-TCpu1
      TWall=TWall2-TWall1
*
#else
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(TwoHam)
         Call Unused_real(TCPU)
         Call Unused_real(Twall)
      End If
#endif
      Return
      End
