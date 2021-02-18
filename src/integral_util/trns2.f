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
* Copyright (C) 1990, Roland Lindh                                     *
************************************************************************
      SubRoutine Trns2(Win,Wout,nvec,nc)
************************************************************************
*                                                                      *
* Object: utility routine to transform a AO batch in case of redun-    *
*         dancy of type aA=bB or cC=dD.                                *
*                                                                      *
* Called from: Trns1                                                   *
*                                                                      *
* Calling    : DCopy  (ESSL)                                           *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             June '90                                                 *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
      Real*8 Win(nvec,nc,nc), Wout(nvec,nc,nc)
*
*     Write (*,*) ' In Trns2, nvec,nc=',nVec,nc
      Do 10 ic = 1, nc
         Do 20 id = 1, nc
            call dcopy_(nvec,Win(1,ic,id),1,Wout(1,id,ic),1)
 20      Continue
 10   Continue
*     Call GetMem(' Exit Trns2','CHECK','REAL',iDum,iDum)
      Return
      End
