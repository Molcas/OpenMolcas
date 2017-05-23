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
      SubRoutine Start0(CMO,TrM,mBB,nD,OneHam,Ovrlp,mBT,Eor,mmB)
************************************************************************
*                                                                      *
*     purpose: Get starting orbitals from diagonalization of the core. *
*              The result is stored in CMO. Those orbitals are         *
*              then optimized in SubRoutine WfCtl. A set of orthogonal,*
*              symmetry adapted AOs is stored in TrM.                  *
*                                                                      *
*     called from: SOrb                                                *
*                                                                      *
*     calls to: TrGen, DCore                                           *
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
#include "mxdm.fh"
#include "infscf.fh"
      Real*8 CMO(mBB,nD), TrM(mBB,nD), OneHam(mBT), Ovrlp(mBT),
     &       EOr(mmB,nD)
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*
*---- Form transformation matrix
      Call TrGen(TrM(1,1),mBB,Ovrlp,OneHam,nBT)
      if(nD.eq.2) Call dCopy_(mBB,TrM(1,1),1,TrM(1,2),1)
*
*---- Diagonalize core
      Do iD = 1, nD
         Call DCore(OneHam,nBT,CMO(1,iD),TrM(1,iD),nBO,
     &              EOr(1,iD),mmB,nOcc(1,iD),Ovrlp)
      End Do
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
