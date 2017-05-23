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
*               2017, Roland Lindh                                     *
************************************************************************
      SubRoutine Start3(CMO,TrM,mBB,nD,OneHam,Ovrlp,mBT,Dens)
************************************************************************
*                                                                      *
*     purpose: Get starting orbitals from density matrix read as input.*
*                                                                      *
*     called from: SOrb                                                *
*                                                                      *
*     calls to: TrGen                                                  *
*                                                                      *
*----------------------------------------------------------------------*
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
#include "real.fh"
#include "file.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "WrkSpc.fh"
      Real*8 CMO(mBB,nD), TrM(mBB,nD), OneHam(mBT), Ovrlp(mBT),
     &       Dens(mBT,nD)
*
*---- Define local variables
      Integer nBasX(8)
      Character*8 Location


#include "SysDef.fh"

*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
       Location='Start3'
*
*---- Compute transformation matrix
      Do iD = 1, nD
         Call TrGen(TrM(1,iD),nBB,Ovrlp,OneHam,nBT)
         Call DCopy_(nBO,TrM(1,iD),1,CMO(1,iD),1)
      End Do
*
*...  read old system definitions
*...  check for compatibility of the calculations
*
*        Call get_iScalar('nSym',nSymX)
         Call Peek_iScalar('nSym',nSymX)
         If (nSymX.ne.nSym) Then
         Call SysWarnMsg(Location,
     &         'Error inconsistent number of Irreps',' ')
         Call SysCondMsg('nSymX=nSym',nSymX,'<>',nSym)
         End If
         Call Get_iArray('nBas',nBasX,nSymX)
         Do iSym=1,nSym
            If (nBasX(iSym).ne.nBas(iSym)) Then
         Call SysWarnMsg(Location,
     &         'Error inconsistent nBas',' ')
         Call SysCondMsg('nBasX(iSym)=nBas (iSym)',
     &         nBasX(iSym),'<>',nBas(iSym))
            End If
         End Do
*
*...  read old density matrix
         Call Get_D1AO(ipD1AO,Length)
         If (Length.ne.nBT) Then
         Call SysWarnMsg(Location,
     &         'Error Reading D1AO',' ')
         Call SysCondMsg('Length.ne.nBT',Length,'<>',nBT)
         End If
         call dcopy_(nBT,Work(ipD1AO),1,Dens(1,1),1)
         Call Free_Work(ipD1AO)
         if (iUHF.eq.1) then
            Call Get_D1sAO(ipD1AO,Length)
            If (Length.ne.nBT) Then
            Call SysWarnMsg(Location,
     &            'Error Reading D1SAO',' ')
            Call SysCondMsg('Length.ne.nBT',Length,'<>',nBT)
         End If
         call dcopy_(nBT,Work(ipD1AO),1,Dens(1,2),1)
         Call Free_Work(ipD1AO)
         Call Quit()
c now we need to fix interface - actually we read a+b,a-b
         Do i=1,nBT
            ra=Half*(Dens(i,1)+Dens(i,2))
            rb=Half*(Dens(i,1)-Dens(i,2))
            Dens(i,1)=ra
            Dens(i,2)=rb
        End Do

         endif
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
