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
      SubRoutine CrSph_mck(Win,nijx,nab,Coeff1,n1,Tr1,Pr1,
     &                  Wout,mab)
************************************************************************
*                                                                      *
*  Object: to transform the one electron integrals from cartesian      *
*          basis to spherical basis.                                   *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : RecPrt                                                  *
*              DGEMM_   (ESSL)                                         *
*              DGeTMO   (ESSL)                                         *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             February '90                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
c#include "print.fh"
      Real*8 Win(nab*nijx),
     &       Coeff1((n1+1)*(n1+2)/2,(n1+1)*(n1+2)/2),
     &       Wout(mab*nijx)
      Logical Tr1, Pr1
*
c     iRout = 26
c     iPrint = nPrint(iRout)
      l1=(n1+1)*(n1+2)/2
      k1=l1
      If (Pr1) k1 = 2*n1 + 1
*
      If (Tr1) Then
*
*        Starting with a,bIJx transforming to bIJx,A
*
         Call DGEMM_('T','N',
     &               nijx,k1,l1,
     &               1.0d0,Win,l1,
     &               Coeff1,l1,
     &               0.0d0,Wout,nijx)
*
      Else
*
*        Transpose from ab,IJ,x to b,IJ,x,a
*
         Call DGeTmO(Win,l1,l1,nijx,Wout,nijx)
*
*        Start transforming b,IJ,x,a to IJ,x,aB
*
      End If
*
*     Call GetMem('CarSph','CHEC','REAL',iDum,iDum)
      Return
      End
