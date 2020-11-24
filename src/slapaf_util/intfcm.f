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
      SubRoutine IntFcm(ipH,nQQ,lOld,lOld_Implicit)
************************************************************************
*                                                                      *
* Object: to initialize the Hessian matrix for the first iteration.    *
*                                                                      *
*     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
*             University of Lund, SWEDEN                               *
*             May '91                                                  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
      Logical lOld, lOld_Implicit, Hess_Found, IRC
*
      iRout = 51
      iPrint = nPrint(iRout)
*
*     Read force constant matrix from old interphase
*
      If (lOld) Then
*
*        Expl

         Call OldFcm(ipH,nQQ,'RUNOLD')

      Else

         Call qpg_iScalar('IRC',IRC)

         If (.Not.IRC) Then
            Call qpg_dArray('Hess',Hess_Found,nHess)

            If (Hess_Found.And.(nHess.gt.0)) Then
               lOld_Implicit=.True.
               Call OldFcm(ipH,nQQ,'RUNFILE')
            Else
               ipH =  ip_Dummy
            End If

         End If

      End If

      If (.Not.lOld.and.lOld_Implicit) lOld=.True.
      If (iPrint.ge.99.and.lOld) Call RecPrt('IntFcm: Final Hessian',
     &                       ' ',Work(ipH),nQQ,nQQ)
*
      Return
      End
