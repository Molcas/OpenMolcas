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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine IntFcm(lOld,lOld_Implicit)
************************************************************************
*                                                                      *
* Object: to initialize the Hessian matrix for the first iteration.    *
*                                                                      *
*     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
*             University of Lund, SWEDEN                               *
*             May 1991                                                 *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8 rDum(1)
      Logical lOld, lOld_Implicit, Hess_Found, IRC
      Real*8, Allocatable:: Hess(:)
*                                                                      *
************************************************************************
*                                                                      *
      Interface
        Subroutine OldFcm(Hess,nQQ,Lbl)
        Real*8, Allocatable:: Hess(:)
        Integer nQQ
        Character*(*) Lbl
        End Subroutine OldFcm
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
*
*     Read force constant matrix from old interphase
*
      If (lOld) Then
*
*        Explicit request to use an old force constant matrix stored
*        on an old runfile.

         Call OldFcm(Hess,nQQ,'RUNOLD')

      Else
*
*        If this is not an IRC calculation explore if the runfile
*        contains a Hessian. If so, pull it off the runfile.
*
         Call qpg_iScalar('IRC',IRC)

         If (.Not.IRC) Then
            Call qpg_dArray('Hess',Hess_Found,nHess)

            If (Hess_Found.And.(nHess.gt.0)) Then
               lOld_Implicit=.True.
               Call OldFcm(Hess,nQQ,'RUNFILE')
            End If

         End If

      End If

      If (.Not.lOld.and.lOld_Implicit) lOld=.True.

      If (lOld) Then
#ifdef _DEBUGPRINT_
         Call RecPrt('IntFcm: Final Hessian',' ',Hess,nQQ,nQQ)
#endif
         Call Put_dArray('Hss_Q',Hess,nQQ**2)
         Call Put_dArray('Hss_upd',rDum,0)
         Call mma_deallocate(Hess)
      End If

*
      Return
      End
