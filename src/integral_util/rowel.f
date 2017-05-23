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
* Copyright (C) 1992, Gunnar Karlstrom                                 *
************************************************************************
      Subroutine Rowel(nZeta,r0,Beta,K,alpha,P,a,gri,grin,jsum)
************************************************************************
* 1992                                                                 *
* Gunnar Karlstrom                                                     *
* Department of Theoretical Chemistry                                  *
* University of Lund                                                   *
* Lund                                                                 *
* Sweden                                                               *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "welcom.fh"
      Real*8 gri(nZeta*jsum), alpha(nZeta), a(nZeta),
     &       grin((k+1)*(k/2+1)*(k/4+1)*nZeta), P(nZeta,3)
*
      Call QEnter('Rowel')
*
      Call poti(k,ipot3)
*     Call IecPrt(' ipot3(0:k+1)',ipot3,k+2,1)
      isum=ipot3(k+1)
      Call bino(k+6)
      Call fiin(k+1)
      Call tetin(k+1)
      Call ylmnor(k+1)
*
      Do 200 iZeta = 1, nZeta
         a(iZeta) = Sqrt(P(iZeta,1)**2 +
     &                   P(iZeta,2)**2 +
     &                   P(iZeta,3)**2 )
 200  Continue
*     Call RecPrt(' In Rowel: Distances',' ',a,nZeta,1)
*
      fac(0)=One
      Do 99 i=1,k+2
         fac(i)=fac(i-1)*DBLE(i)
 99   Continue
      iss=k
      itt=k
      ind=1
      Call priwel(k,alpha,beta,r0,a,gri,nZeta,isum,grin)
*     Call RecPrt('Internal well integrals',' ',gri,nZeta,isum)
      Call QExit('Rowel')
      Return
      End
