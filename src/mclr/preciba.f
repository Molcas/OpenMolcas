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
      SubRoutine Preciba(iB,iS,jS,nd,rOut,nba,
     &                   focki,focka,fock,sign,
     &                   A_J,A_K,Scr,nScr)
************************************************************************
*                                          [2]                         *
*     Calculates the diagonal submatrix of E    that couple            *
*                                                                      *
*     kappa           with   kappa                for a                *
*          kinactive,virtual        kinactive,active                   *
*                                                                      *
*     single inactive index.                                           *
*     Used for preconditioner.                                         *
*                                                                      *
*     See Olsen,Yeager, Joergensen:                                    *
*      "Optimization and characterization of an MCSCF state"           *
*                                                                      *
*     Called by prec                                                   *
*                                                                      *
*     ib,is       :       inactive index for the submatrix             *
*     js          :       symmetry of virtual,active                   *
*     rOut        :       Submatrix                                    *
*                                                                      *
************************************************************************
      Implicit Real*8(a-h,o-z)
#include "Input.fh"
#include "WrkSpc.fh"
#include "Pointers.fh"
      Real*8  Fock(nba,nba),Focki(nba,nba),FockA(nba,nba)
      Real*8 rOut(*), A_J(nScr), A_K(nScr), Scr(nScr)
*                                                                      *
************************************************************************
*                                                                      *
       iTri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
       iTri1(i,j)=nTri-itri(nd-Min(i,j)+1,nd-Min(i,j)+1)
     &           +Max(i,j)-Min(i,j)+1
*                                                                      *
************************************************************************
*                                                                      *
      nTri=itri(nd,nd)
      jVert=nOrb(jS)-nAsh(jS)-nIsh(jS)
      nO=nAsh(jS)+nIsh(jS)
*                                                                      *
************************************************************************
*                                                                      *
*     Get block J^(iB,iB)
      Call Coul(jS,jS,iS,iS,iB,iB,A_J,Scr)
*     Get block K^(iB,iB)
      Call Exch(jS,iS,jS,iS,iB,iB,A_K,Scr)
*
      Do jA=1,nAsh(jS)
         jAA=jA+nIsh(jS)
         ip=itri1(ja,nd-jVert+1)
         Do jB=1,nAsh(jS)
            jBB=jB+nIsh(jS)
*           Get D_(ja,jb)
            rDens=-sign*Work(ipg1-1+(iTri(jA+nA(jS),jB+nA(jS))))
            If (jA.eq.jB) rDens=rdens+sign*2.0d0
*
            ivB=(jBB-1)*nBas(jS) + nO + 1
            Call DaXpY_(jVert,6.0d0*rDens,
     &                 A_K(ivB),1,
     &                 rOut(ip),1) ! ????
            Call DaXpY_(jVert,-2.0d0*rDens,
     &                 A_J(ivB),1,
     &                 rOut(ip),1)
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Do jA=1,nAsh(js)
         ip=iTri1(ja,nAsh(js)+1)
         Call DaXpY_(jVert,sign*4.0d0,Focki(nO+1,ja+nIsh(js)),1,
     &               rout(ip),1)
         Call DaXpY_(jVert,sign*4.0d0,FockA(nO+1,ja+nIsh(js)),1,
     &               rout(ip),1)
         Call DaXpY_(jVert,-sign,     Fock (nO+1,ja+nIsh(js)),1,
     &               rout(ip),1)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
