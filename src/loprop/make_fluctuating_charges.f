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
      Subroutine Make_Fluctuating_Charges(nAtoms,iANr,nij,nPert,rMP,
     &                                    nElem,EC,Alpha)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Real*8 rMP(nij,0:nElem-1,0:nPert-1), EC(3,nij)
      Integer iANr(nAtoms)
*                                                                      *
************************************************************************
*                                                                      *
*     For the A-matrix
*
      Call Allocate_Work(ip_AInv,nAtoms**2)
      Call Allocate_Work(ip_A,nAtoms**2)
*
      Call Build_AMatrix(nAtoms,iANr,Work(ip_A),Work(ip_AInv),EC,nij,
     &     Alpha)
*
      Call Free_Work(ip_A)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the fluctuating charges
*
      Call Allocate_Work(ip_lambda,nAtoms)
      Call Allocate_Work(ip_dq    ,nAtoms)
*
      Call Fluctuating(Work(ip_AInv),nAtoms,Work(ip_lambda),Work(ip_dq),
     &                 nij,nPert,iANr,rMP,nElem,EC,Alpha)
*
      Call Free_Work(ip_lambda)
      Call Free_Work(ip_dq)
      Call Free_Work(ip_AInv)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
