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
      SubRoutine TrGen(TrMat,nTrMat,Ovlp,OneHam,mBT)
************************************************************************
*                                                                      *
*     purpose: Generate transformation matrix from AO's in which       *
*              integrals are computed to orthogonal, symmetry adapted  *
*              (and spherical, if desired) functions - near linear     *
*              dependencies are also removed.                          *
*                                                                      *
*     input:                                                           *
*       Ovlp    : overlap matrix in AO basis of length nOvlp           *
*                                                                      *
*     output:                                                          *
*       TrMat   : Transformation matrix of length nTrMat               *
*                                                                      *
*     called from: Start1, Start3                                      *
*                                                                      *
*     calls to: Freeze, OvlDel, SetUp, Ortho                           *
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
*
      Implicit Real*8 (a-h,o-z)
*
      Real*8 TrMat(nTrMat),Ovlp(mBT),OneHam(mBT)
*
#include "real.fh"

#include "mxdm.fh"
#include "infscf.fh"
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
#ifdef _DEBUG_
#endif
*
      ind=0
      Do iSym = 1, nSym
         Do i = 1, nBas(iSym)
            Do j = 1, nBas(iSym)
               ind = ind + 1
               TrMat(ind) = 0.0d+00
               If (I.eq.J) TrMat(ind) = 1.0d+00
            End Do
         End Do
      End Do
*
*---- Set up certain parameters (nOrb(i) may be changed)
      Call SetUp
*
*---- Move frozen atomic orbitals to the begining
      If (nnFr.gt.0) Then
         Call Freeze(TrMat,nBO,OneHam,mBT)
         Call SetUp
      End If
*
*---- Remove near linear dependencies from basis set
      If (DelThr.ne.0.0d+00) Then
         Call OvlDel(Ovlp,nBT,TrMat,nBO)
         Call SetUp
      End If
*
*---- Orthogonalize final orbitals
      Call Ortho(TrMat,nBO,Ovlp,nBT)
*
#ifdef _DEBUG_
#endif
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
