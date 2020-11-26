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
      SubRoutine Precaii(iB,is,js,nd,ir,rOut,nbai,nbaj,
     &                   fockii,fockai,fockti,
     &                   focki,focka,fock,sign,
     &                   A_J,A_K,Scr,nScr)
************************************************************************
*                                                                      *
*                                        [2]                           *
*   Calculates the diagonal submatrix of E    that couple              *
*                                                                      *
*   kappa                with   kappa                for a             *
*        kactive,inactive            kactive,inactive                  *
*                                                                      *
*   single active index.                                               *
*   Used for preconditioner.                                           *
*                                                                      *
*   See Olsen,Yeager, Joergensen:                                      *
*    "Optimization and characterization of an MCSCF state"             *
*                                                                      *
*     Eq. C.12a                                                        *
*                                                                      *
*   Called by prec                                                     *
*                                                                      *
*   ib,is       :       active index for the submatrix                 *
*   js          :       symmetry of inactive,inactive                  *
*   rOut        :       Submatrix                                      *
*                                                                      *
************************************************************************
      use Arrays, only: G1t, G2t
      Implicit Real*8(a-h,o-z)
#include "Input.fh"
#include "Pointers.fh"
      Real*8 Fock(nbaj,nbaj),FockA(nBaj,nBaj),Focki(nbaj,nbaj)
      Real*8 rout(nd*(nd+1)/2), A_J(nScr), A_K(nScr), Scr(nScr)
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
      iBB=ib+nA(is)
      iiB=ib+nish(is)
*                                                                      *
************************************************************************
*                                                                      *
      Do jA=1,nIsh(jS)
         Do jB=1,jA
*
            i=itri1(ja,jb)
*
            Do kS=1,nSym
*
               Call Coul(kS,kS,jS,jS,jA,jB,A_J,Scr)
               Call Exch(kS,jS,kS,jS,jA,jB,A_K,Scr)
*
               Do jC=1,nAsh(ks)
                  jCC=jC+nA(ks)
                  jjC=jC+nIsh(ks)
                  Do jD=1,nAsh(kS)
                     jDD=jD+nA(ks)
                     jjD=jD+nIsh(ks)
*
*                    gamma(cdbb)=gamma(bbcd)
*
                     rDens1=sign*G2t(itri(itri(jCC,jDD),
     &                           itri(iBB,iBB)))
*
*                    gamma(bdcb)
*
                     rDens2=sign*G2t(itri(itri(iBB,jDD),
     &                           itri(jCC,iBB)))
*
*                    (cd|ij)
*
                     icd=(jjD-1)*nBas(kS)+jjC
                     cdij=A_J(icd)
                     cidj=A_K(icd)
                     rout(i) = rout(i) + 2.0d0*(rDens2*cidj
     &                                 + 2.0d0* rDens1*cdij)
                  End Do
               End Do
            End Do
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Do iC=1,nAsh(is)
         iCC=ic+nA(is)
         iiC=iC+nIsh(is)
*
*        2*(delta(bc)-D(bc))
*
         rDens=sign*(-G1t(itri(iCC,iBB)))
         If (iCC.eq.iBB) rdens=rdens+sign
         rDens=2.0D0*rDens
*
         Call Coul(jS,jS,iS,iS,iiB,iiC,A_J,Scr)
         Call Exch(jS,iS,jS,iS,iiB,iiC,A_K,Scr)
*
         Do jA=1,nIsh(jS)
            Do jB=1,jA
               i=itri1(jA,jB)
*
*              (ci|bj)
*              (bi|cj)
*
               icb = (jB-1)*nBas(jS) + jA
               ibc = (jA-1)*nBas(jS) + jB
               cibj=A_K(icb)
               bicj=A_K(ibc)
*
*              (bc|ij)
*
               bcij=A_J(ibc)
*
               rout(i)=rout(i)+rdens*(       7.0d0*cibj
     &                                 -sign*      bicj
     &                                 -sign*2.0D0*bcij )
*
            End Do
         End Do
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      rFock = sign*2.0d0*Fockii + sign*2.0d0*Fockai - sign*Fockti
      rdens=sign*2.0d0*G1t(itri(ibb ,ibb))
      i=0 ! dummy initialize
*
      Do jA=1,nIsh(jS)
         Do jB=1,jA
*
            i=itri1(ja,jb)
*
            rout(i) = rout(i) - sign*4.0d0*( Focka(jA,jB)
     &                                      +Focki(jA,jB) )
     &                        + rdens*Focki(ja,jb)
         End Do
         rout(i) = rout(i) + 2.0d0*rfock
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(ir)
         Call Unused_integer(nbai)
         Call Unused_real_array(fock)
      End If
      End
