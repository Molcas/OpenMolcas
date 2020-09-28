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
      SubRoutine Orthox(S,C,nOrb,nBas)
************************************************************************
*                                                                      *
*     purpose: Perform Gram-Schmidt orthogonalization scheme           *
*                                                                      *
*     input:                                                           *
*       S       : overlap in non-orthonormal basis (square storage)    *
*       C       : matrix transforming non-orthonormal basis            *
*      nBas,nOrb: dimensions                                           *
*                                                                      *
*     output:                                                          *
*       C       : matrix transforming orthonormal basis                *
*                                                                      *
*     called from: Ortho                                               *
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
      Real*8 C(nBas,nOrb),S(nOrb,nOrb)
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
#ifdef _DEBUG_
#endif
*
      Do 100 iOrb = 1, nOrb
         F = Zero
         If (S(iOrb,iOrb).gt.Zero) F=One/SqRt(S(iOrb,iOrb))
         Do 110 iBas = 1, nBas
            C(iBas,iOrb) = F*C(iBas,iOrb)
  110    Continue
         Do 120 jOrb = 1, nOrb
            S(iOrb,jOrb) = F*S(iOrb,jOrb)
            S(jOrb,iOrb) = F*S(jOrb,iOrb)
  120    Continue
         S(iOrb,iOrb) = One
         Do 130 jOrb = iOrb + 1, nOrb
            A = S(iOrb,jOrb)
            Do 131 iBas = 1, nBas
               C(iBas,jOrb) = C(iBas,jOrb) - A*C(iBas,iOrb)
  131       Continue
            Do 132 kOrb = 1, nOrb
               S(jOrb,kOrb) = S(jOrb,kOrb) - A*S(iOrb,kOrb)
  132       Continue
            Do 133 kOrb = 1, nOrb
               S(kOrb,jOrb) = S(kOrb,jOrb) - A*S(kOrb,iOrb)
  133       Continue
  130    Continue
  100 Continue
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
