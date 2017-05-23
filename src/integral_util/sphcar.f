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
* Copyright (C) 1990,1991, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine SphCar(Win,nab,nijx,Scrt,nScrt,Coeff1,n1,Tr1,Pr1,
     &                  Coeff2,n2,Tr2,Pr2,Wout,mab)
************************************************************************
*                                                                      *
*  Object: to project an one-electrom matrix from spherical harmonics  *
*          to cartesians.                                              *
*                                                                      *
*          Matrix on input  AB,ij                                      *
*          Matrix on output ij,ab                                      *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : RecPrt                                                  *
*              DGEMM_   (ESSL)                                         *
*              DGeTMO   (ESSL)                                         *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             February '90                                             *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN                                          *
*             Modified spherical harmonics to cartesians October '91.  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
      Real*8 Win(nab*nijx), Scrt(nScrt),
     &       Coeff1((n1+1)*(n1+2)/2,(n1+1)*(n1+2)/2),
     &       Coeff2((n2+1)*(n2+2)/2,(n2+1)*(n2+2)/2),
     &       Wout(mab*nijx)
      Logical Tr1, Pr1, Tr2, Pr2
*
      iRout = 26
      iPrint = nPrint(iRout)
      l1=(n1+1)*(n1+2)/2
      k1=l1
      If (Pr1) k1 = 2*n1 + 1
      l2=(n2+1)*(n2+2)/2
      k2 = l2
      If (Pr2) k2 = 2*n2 + 1
      if (iprint.ge.99) then
        call recprt(' Win',' ',Win,nab,nijx)
        call recprt('Coeff1',' ',Coeff1,l1,l1)
        call recprt('Coeff2',' ',Coeff2,l2,l2)
      end if
*
      If (Tr1.and.Tr2) Then
*
*        Starting with A,Bij transforming to Bij,a
*
         Call DGEMM_('T','T',
     &               k2*nijx,l1,k1,
     &               1.0d0,Win,k1,
     &               Coeff1,l1,
     &               0.0d0,Scrt,k2*nijx)
*
*        Transform B,ija to ij,ab
*
         Call DGEMM_('T','T',
     &               nijx*l1,l2,k2,
     &               1.0d0,Scrt,k2,
     &               Coeff2,l2,
     &               0.0d0,Wout,nijx*l1)
*
      Else If(Tr2) Then
*
*        Transpose from aB,ij to B,ija
*
         Call DGeTmO(Win,l1,l1,k2*nijx,Scrt,k2*nijx)
*
*        Start transforming B,ija to ij,ab
*
         Call DGEMM_('T','T',
     &               nijx*l1,l2,k2,
     &               1.0d0,Scrt,k2,
     &               Coeff2,l2,
     &               0.0d0,Wout,nijx*l1)
      Else
*
*        Starting with A,bij transforming to a,bij
*
         Call DGEMM_('N','N',
     &               l1,l2*nijx,k1,
     &               1.0d0,Coeff1,l1,
     &               Win,k1,
     &               0.0d0,Scrt,l1)
*
*        Transpose to ij,ab
*
         Call DGeTmO(Scrt,l1*l2,l1*l2,nijx,Wout,nijx)
      End If
*
*     Call GetMem('SphCar','CHEC','REAL',iDum,iDum)
      Return
      End
