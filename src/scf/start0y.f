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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
************************************************************************
      SubRoutine Start0y(CMO,mBB,nD,EOr,mmB)
************************************************************************
*                                                                      *
* This routine reads old SCf orbitals as start orbitals.               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*                                                                      *
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "mxdm.fh"
#include "infscf.fh"
*
      Real*8 CMO(mBB,nD), EOr(mmB,nD)
      Logical found
*
*----------------------------------------------------------------------*
* Get start orbitals.                                                  *
*----------------------------------------------------------------------*
*
      iRC=0
      Call qpg_darray('SCF orbitals',found,ndata)
      If (Found) Then
         Call get_darray('SCF orbitals',CMO(1,1),ndata)
      Else
         iRC=-1
      End If
      Call qpg_darray('OrbE',found,ndata)
      If (Found) Then
         Call get_darray('OrbE',EOr(1,1),ndata)
      End If
*
      If (nD.eq.2) Then
*
         Call DCopy_(mBB,CMO(1,1),1,CMO(1,2),1)
         Call DCopy_(mmB,EOr(1,1),1,EOr(1,2),1)
*
         Call qpg_darray('SCF orbitals_ab',found,ndata)
         If(Found) Then
            Call get_darray('SCF orbitals_ab',CMO(1,2),ndata)
         Else
            iRC=-1
         End If
         Call qpg_darray('OrbE_ab',found,ndata)
         If(found) Then
            Call get_darray('OrbE_ab',EOr(1,2),ndata)
         End If
      End If
*
      Call qpg_iarray('nDel',Found,ndata)
      nSum=0
      If(Found) Then
         Call Get_iArray('nDel',nDel,ndata)
         Do iSym=1,nSym
            nSum=nSum+nDel(iSym)
         End Do
      End If
      If (nSum.gt.0) Then
         Do iSym=1,nSym
            nOrb(iSym)=nBas(iSym)-nDel(iSym)
         End Do
         Do iD = 1, nD
            Call TrimCMO(CMO(1,iD),CMO(1,iD),nSym,nBas,nOrb)
            Call TrimEor(Eor(1,iD),Eor(1,iD),nSym,nBas,nOrb)
         End Do
      End If
*----------------------------------------------------------------------*
* Done                                                                 *
*----------------------------------------------------------------------*
      Return
      End
