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
* Copyright (C) 1996, Anders Bernhardsson                              *
************************************************************************
      SubRoutine Precibb_td(ib,is,js,nd,rout,nba,
     &                  Temp1,Scr,Temp2,
     &                  fockii,fockai,
     &                  focki,focka,sign)
************************************************************************
*                                        [2]
*   Calculates the diagonal submatrix of E    that couple
*
*   kappa           with   kappa                for a
*        kinactive,virtual        kinactive,virtual
*
*   single inactive index.
*   Used for preconditioner.
*
*   See Olsen,Yeager, Joergensen:
*    "Optimization and characterization of an MCSCF state"
*
*   Called by prec
*
*   ib,is       :       inactive index for the submatrix
*   js          :       symmetry of virtual,virtual
*   rOut        :       Submatrix
*
************************************************************************
      Implicit Real*8(a-h,o-z)
#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
      Real*8 rout(*)
      Real*8 Temp1(nBa,nBa),Focki(nBa,nBa),Focka(nBa,nBa)
      Real*8 Temp2(*), Scr(*)
*                                                                      *
************************************************************************
*                                                                      *
      iTri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
      iTri1(i,j)=nTri-itri(nd-Min(i,j)+1,nd-Min(i,j)+1)
     &          +Max(i,j)-Min(i,j)+1
*                                                                      *
************************************************************************
*                                                                      *
      nTri=itri(nd,nd)
*
      jVert=nBas(js)-nAsh(js)-nIsh(js)
      if (jvert.eq.0) Return
*
      i1=nD-jVert+1
      ip=itri1(i1,i1)
      ra=4.0d0*sign*(Fockii+Fockai)
      Call COUL(jS,jS,iS,iS,iB,iB,Temp2,Scr)
      Call Dyax(nba**2,-sign*4.0d0,Temp2,1,Temp1,1)
      Call EXCH(js,is,js,is,ib,ib,Temp2,Scr)
      Call DaXpY_(nBa**2,sign*12.0d0,Temp2,1,Temp1,1)
      i=ip-1
      Do kB=nIsh(jS)+nAsh(jS)+1,nBas(jS)
         rOut(i+1)=rout(i+1)-ra
         Do lB=kb,nBAS(JS)
            i=i+1
            rOut(i)=rout(i)+Temp1(kb,lb)+
     &              sign*4.0d0*Focki(kb,lb)+
     &              sign*4.0d0*Focka(kb,lb)
          End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
