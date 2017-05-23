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
* Copyright (C) 1990, IBM                                              *
************************************************************************
      Integer Function iPrmt(jOper,iChct)
************************************************************************
*     Returns the phase factor of a basis function under a symmetry    *
*     operation, jOper. iChct contains the information about the       *
*     character of the basis function.                                 *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#ifdef _OLD_
#include "itmax.fh"
#include "info.fh"
#else
      Integer is_nSym, nSym, is_iOper, iOper(0:7)
      Common /iPrmt_Common/ is_nSym, nSym, is_iOper, iOper
*
      if(is_nSym.eq.0) then
        Call get_iScalar('nSym',nSym)
        is_nSym=1
      endif
      if(is_iOper.eq.0) then
        Call Get_iArray('Symmetry operations',iOper,nSym)
        is_iOper=1
      endif
#endif
      iPrmt = 1
      iCom= iAnd(iOper(jOper),iChct)
      Do 10 i = 1, 3
         If (iAnd(iCom,2**(i-1)).ne.0) iPrmt = iPrmt*(-1)
 10   Continue
      Return
      End
      Subroutine iPrmt_Init()
      Implicit Real*8 (a-h,o-z)
      Integer is_nSym, nSym, is_iOper, iOper(0:7)
      Common /iPrmt_Common/ is_nSym, nSym, is_iOper, iOper
*
      is_nSym=0
      is_iOper=0
*
      Return
      End
