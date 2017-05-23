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
      SubRoutine Preciaa(iB,iS,jS,nd,rOut,nbai,nbaj,
     &                   fockii,fockai,fockti,
     &                   focki,focka,fock,sign,
     &                   A_J,A_K,Scr,nScr)
************************************************************************
*     Change Fock(i) ne Fock(j)
*                                          [2]
*     Calculates the diagonal submatrix of E    that couple
*
*     kappa               with   kappa                for a
*          kinactive,active           kinactive,active
*
*     single inactive index.
*     Used for preconditioner.
*
*     See Olsen,Yeager, Joergensen:
*      "Optimization and characterization of an MCSCF state"
*
*     Called by prec
*
*     ib,is       :       inactive index for the submatrix
*     js          :       symmetry of active,active
*     rOut        :       Submatrix
*
************************************************************************
      Implicit Real*8(a-h,o-z)
#include "Input.fh"
#include "Pointers.fh"
#include "standard_iounits.fh"
#include "WrkSpc.fh"
      Real*8 focki(nbaj,nbaj),fock(nbaj,nbaj),focka(nbaj,nbaj),
     &       rout(*), A_J(nScr), A_K(nScr), Scr(nScr)
*                                                                      *
************************************************************************
*                                                                      *
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
      itri1(i,j)=nTri-itri(nd-Min(i,j)+1,nd-Min(i,j)+1)
     &          +Max(i,j)-Min(i,j)+1
*                                                                      *
************************************************************************
*                                                                      *
      nTri=itri(nd,nd)
*                                                                      *
************************************************************************
*                                                                      *
      Do kS=1,nSym
*
         If (nAsh(kS).ne.0)
     &       Call Coul(kS,kS,iS,iS,iB,iB,A_J,Scr)
*
         Do jC=1,nAsh(kS)
            jjC=JC+nA(kS)
            jCC=jC+nIsh(kS)
*
            Call Coul(kS,iS,kS,iS,jCC,iB,A_K,Scr)
*
            Do jD=1,nAsh(kS)
               jjD=JD+nA(kS)
               jDD=jD+nIsh(kS)
*
               ip1=(jDD-1)*nBas(kS)+jCC
               ip2=( iB-1)*nBas(kS)+jDD
               aabb=A_J(ip1)
               abab=A_K(ip2)
*
               Do jA=1,nAsh(jS)
                  jjA=jA+nA(jS)
                  Do jB=1,jA
                     jjB=jB+nA(jS)
                     i=itri1(jA,jB)
*
                     rDens1=sign*work(ipg2-1+
     &                          (itri(itri(jjC,jjD),itri(jjB,jjA))))
                     rDens2=sign*work(ipg2-1+
     &                          (itri(itri(jjB,jjD),itri(jjC,jjA))))
*
                     rout(i)=rout(i)+2.0d0*rDens1*aabb
     &                              +4.0d0*rDens2*abab
*
                  End Do
               End Do
            End Do
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      If (nAsh(jS).ne.0) Then
         Call Exch(jS,iS,jS,iS,iB,iB,A_K,Scr)
         Call Coul(jS,jS,iS,iS,iB,iB,A_J,Scr)
      End If
*
      Do jA=1,nAsh(jS)
         jAA=jA+nA(jS)
         jjA=jA+nIsh(jS)
*
         Do jB=1,jA
            jBB=jB+nA(jS)
            jjB=jB+nIsh(jS)
            i=itri1(jA,jB)
*
            Do jC=1,nAsh(jS)
               jCC=jC+nA(jS)
               jjC=jC+nIsh(jS)
               iBC=(jjC-1)*nBas(jS)+jjB
               BCbb=A_J(iBC)
               BbCb=A_K(iBC)
               iAC=(jjC-1)*nBas(jS)+jjA
               ACbb=A_J(iAC)
               AbCb=A_K(iAC)
*
               rDens1=-sign*Work(ipg1-1+itri(jAA,jCC))
               rDens2=-sign*Work(ipg1-1+itri(jBB,jCC))
               If (jAA.eq.jCC) rDens1=rdens1+sign
               If (jBB.eq.jCC) rDens2=rdens2+sign
*
               rout(i)=rout(i)
     &                +2.0d0*rdens1*(3.0d0*BbCb-BCbb)
     &                +2.0d0*rdens2*(3.0d0*AbCb-ACbb)
*
            End Do
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      rFock=sign*(Fockii+FockAi)
      i=0 ! dummy initialize
      Do jA=1,nAsh(jS)
         jAA=jA+nA(jS)
         jjA=jA+nIsh(js)
         Do jB=1,JA
            jBB=jB+nA(jS)
            jjB=jB+nIsh(js)
            i=itri1(jA,jB)
            rDens=Work(ipG1+itri(jbb,jAA)-1)
            rout(i)=rout(i)+Sign*(2.0d0*rdens*Fockii + ! (ib,ib)+
     &              2.0d0*(2.0d0*Focki(jjA,jjB)+
     &              2.0d0*FockA(jjA,jjB)-Fock(jjB,jjA)))
         End Do
         rout(i) = rout(i) - 4.0d0*rFock
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(nbai)
         Call Unused_real(fockti)
      End If
      End
