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
      SubRoutine UnFold(A,nAA,B,nBB,nSym,nBas)
************************************************************************
*                                                                      *
*     purpose: Expand the density matrix to full storrage and scale    *
*              off diagonal elements                                   *
*                                                                      *
*                                                                      *
*     input:                                                           *
*       A       : triangular matrix of length nAA                      *
*       nSym    : number of blocks to be expanded                      *
*       nBas    : leading dimensions in each block                     *
*                                                                      *
*     output:                                                          *
*       B       : expanded matrix of length nBB                        *
*                                                                      *
*     called from: PMat                                                *
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
#include "real.fh"
*
      Real*8 A(nAA),B(nBB)
      Integer nBas(nSym)
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
#ifdef _DEBUG_
#endif
*
      k1 = 0
      k2 = 0
      Factor=Half
      Do 10 iSym = 1, nSym
        nb = nBas(iSym)
        Do 20 ib = 1, nb
          Do 30 jb = 1, (ib - 1)
             l1    = jb + (ib - 1)*nb + k2
             l2    = ib + (jb - 1)*nb + k2
             l3    = jb + ib*(ib - 1)/2 + k1
             B(l1) = Factor*A(l3)
             B(l2) = Factor*A(l3)
30        Continue
          l2    = ib + (ib - 1)*nb + k2
          l3    = ib + ib*(ib - 1)/2 + k1
          B(l2) = A(l3)
20      Continue
        k1 = k1 + nb*(nb + 1)/2
        k2 = k2 + nb*nb
10    Continue
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
