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
*               1990, IBM                                              *
************************************************************************
      SubRoutine Trns1(Win,Wout,na,nb,nvec,nc)
************************************************************************
*                                                                      *
* Object: utility routine to transform a AO batch in case of redun-    *
*         dancy of type aA=bB or cC=dD.                                *
*                                                                      *
* Called from: TwoEl                                                   *
*                                                                      *
* Calling    : Trns2                                                   *
*              DGeTMO  (ESSL)                                          *
*              DCopy   (ESSL)                                          *
*                                                                      *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             June '90                                                 *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
      Real*8 Win(na,nb), Wout(nb,na)
*
*     Write (*,*) ' In Trns1: na, nb, nVec, nc=',na,nb,nvec,nc
*     Call RecPrt(' Win',' ',Win,na,nb)
      If (nc.eq.1) Then
         call dcopy_(nvec,Win,1,Wout,1)
         Return
      End If
      If (na.eq.1 .or. nb.eq.1) Then
         Call Trns2(Win,Wout,nvec,nc)
      Else
         Call DGeTMO(Win,na,na,nb,Wout,nb)
*        Call RecPrt(' After first DGeTMO',' ',Wout,nb,na)
         Call Trns2(Wout,Win,nvec,nc)
*        Call RecPrt(' After Trns2',' ',Win,nb,na)
         Call DGeTMO(Win,nb,nb,na,Wout,na)
*        Call RecPrt(' After second DGeTMO',' ',Wout,na,nb)
      End If
*     Call GetMem(' Exit Trns1','CHECK','REAL',iDum,iDum)
      Return
      End
