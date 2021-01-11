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
* Copyright (C) 2020, Roland Lindh                                     *
************************************************************************
      Subroutine Trans_K(X,Y,nInter,nIter)
      Use kriging_mod, Only: layer_U
      Implicit None
      Integer nInter, nIter
      Real*8 X(nInter,nIter), Y(nInter,nIter)
*
*     Call RecPrt('layer_U',' ',layer_U,nInter,nInter)
*     Call RecPrt('X',' ',X,nInter,nIter)
      Call DGEMM_('T','N',nInter,nIter,nInter,
     &            1.0D0,layer_U,nInter,
     &                  X,nInter,
     &            0.0D0,Y,nInter)
*     Call RecPrt('Y',' ',Y,nInter,nIter)
*
      Return
      End Subroutine Trans_K
