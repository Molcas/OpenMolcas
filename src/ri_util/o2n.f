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
      Subroutine O2N(AA,AB,BB,Temp,nA,nB,Error)
      Implicit Real*8 (a-h,o-z)
      Real*8 AA(nA,nA), AB(nA,nB), BB(nB,nB), Temp(nA,nB)
*                                                                      *
************************************************************************
*                                                                      *
*     (1) Project the new basis on to the old space and compute the
*         matrix elements.
*
      Call DGEMM_('N','N',
     &             nA,nB,nA,
     &             1.0D0,AA,nA,
     &                   AB,nA,
     &             0.0D0,Temp,nA)
      Call DGEMM_('T','N',
     &            nB,nB,nA,
     &            1.0D0,AB,nA,
     &                  Temp,nA,
     &            0.0D0,BB,nB)
*                                                                      *
************************************************************************
*                                                                      *
*     (2) Do diagnostics on how poorly or good the
*         new basis spans the old space
*
      Error = DDot_(nA,A,nA+1,A,nA+1)
     &      - DDot_(nB,B,nB+1,B,nB+1)
*
      Return
      End
