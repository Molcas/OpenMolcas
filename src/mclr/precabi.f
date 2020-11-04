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
      Subroutine Precabi(ib,is,js,ir,nd,rOut,nba,
     &                   focki,focka,fock,sign,
     &                   A_J,A_K,Scr,nScr)
************************************************************************
*                                          [2]                         *
*     Calculates the diagonal submatrix of E    that couple            *
*                                                                      *
*     kappa           with   kappa                for a                *
*          kactive,virtual        kactive,inactive                     *
*                                                                      *
*     single active index.                                             *
*     Used for preconditioner.                                         *
*                                                                      *
*     See Olsen,Yeager, Joergensen:                                    *
*      "Optimization and characterization of an MCSCF state"           *
*                                                                      *
*     Called by prec                                                   *
*                                                                      *
*     ib,is         :        active index for the submatrix            *
*     js           :        symmetry of inact,virtual                  *
*     rOut          :         Submatrix                                *
*                                                                      *
************************************************************************
      use Arrays, only: G1t
      Implicit Real*8(a-h,o-z)
#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
      Real*8 Fock(nba,nba),Focki(nba,nba),FockA(nba,nba),rOut(*),
     &       A_J(nScr),A_K(nScr),Scr(nScr)
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
      jVert=nOrb(js)-nIsh(js)-nAsh(js)
      If (jVert.eq.0) Return
*
      nO=nAsh(js)+nIsh(js)
      iAA=nA(is)+ib
      itAA=itri(iAA,iAA)
*                                                                      *
************************************************************************
*                                                                      *
      Do kS=1,nSym
         Do kA=1,nAsh(kS)
            kAA=kA+nA(ks)
            kkA=kA+nIsh(ks)
            Do lA=1,nAsh(kS)
               lAA=lA+nA(ks)
               llA=lA+nIsh(ks)
*
               Call Coul(jS,jS,kS,kS,kkA,llA,A_J,Scr)
               Call Exch(jS,kS,jS,kS,kkA,llA,A_K,Scr)
*
               Do jB=1,nIsh(jS)
                  ip=itri1(jB,nd-jVert+1)
*
                  Fact1=-2.0d0*Work(ipG2-1+itri(itAA,itri(kAA,lAA)))
                  Fact2=-4.0d0*
     &                 Work(ipG2-1+(itri(itri(iAA,kAA),itri(iAA,lAA))))
*
                  If (kaa.eq.iaa)
     &               Fact2=Fact2+8.0d0*G1t(itri(iAA,lAA))
                  If (laa.eq.iaa)
     &               Fact1=Fact1-2.0d0*G1t(itri(iAA,kAA))
                  If (laa.eq.iaa)
     &               Fact2=Fact2-2.0d0*G1t(itri(iAA,kAA))
*
                  ivj=(jB-1)*nBas(jS)+no+1
                  Call DaXpY_(jVert,Sign*Fact1,
     &                       A_J(ivj),1,
     &                       rout(ip),1) ! ????
                  Call DaXpY_(jVert,Sign*Fact2,
     &                       A_K(ivj),1,
     &                       rout(ip),1)
*
               End Do
*
            End Do
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Do jB=1,nIsh(jS)
         ip=itri1(jB,nd-jVert+1)
         Fact=(2.0d0-2.0d0*G1t(itAA))
         Call DaxPy_(jVert,Sign*Fact,FockI(nO+1,jB),1,rOut(ip),1)
         Fact=2.0d0
         Call DaxPy_(jVert,Sign*Fact,FockA(nO+1,jB),1,rOut(ip),1)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_integer(ir)
        Call Unused_real_array(fock)
      End If
      End
