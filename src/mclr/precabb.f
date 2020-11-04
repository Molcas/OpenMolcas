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
      SubRoutine Precabb(ib,is,js,nd,nba,rout,Temp1,ntemp,
     &                   Scr,Temp2,fockti,
     &                   focki,focka,fock,sign)
************************************************************************
*                                        [2]
*   Calculates the diagonal submatrix of E    that couple
*
*   kappa           with   kappa                for a
*        kactive,virtual        kactive,virtual
*
*   single active index.
*   Used for preconditioner.
*
*   See Olsen,Yeager, Joergensen:
*    "Optimization and characterization of an MCSCF state"
*
*   Called by prec
*
*   ib,is       :       active index for the submatrix
*   js          :       symmetry of virtual,virtual
*   rOut        :       Submatrix
*
************************************************************************
      use Arrays, only: G1t
      Implicit Real*8(a-h,o-z)
#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
      Real*8 rout(*)
      Real*8 Temp2(nBa,nBa),Focki(nBa,nBa),Focka(nBa,nBa),
     &       Fock(nBa,nBa)
      Real*8 Temp1(nTemp),Scr(nTemp)
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
      iib=ib+nA(is)
      jVert=nBas(js)-nAsh(js)-nIsh(js)
      If (jvert.eq.0) Return
*
      i2=nD-jVert+1
      i1=(nAsh(js)+nish(js))*nbas(js)
      ip=iTri1(i2,i2)
      rF=sign*Fockti
      call dcopy_(nBa**2,[0.0d0],0,Temp2,1)
*
      Do kS=1,nSym
         If (nBas(js)*nash(ks).gt.0) Then
            Do kBB=nish(ks)+1,nB(kS)
               Do kCC=nish(ks)+1,kBB
                  Call COUL(jS,jS,kS,kS,kbb,kcc,Temp1,Scr)
                  If (kBB.gt.nish(ks).and.kCC.gt.nish(ks)) Then
                     kkB=kBB+nA(ks)-nish(ks)
                     kkC=kCC+nA(ks)-Nish(ks)
                     rDens1=sign*2.0d0*Work(ipG2-1+
     &                      itri(itri(iib,iib),itri(kkb,kkc)))
*
                     If (kbb.ne.kcc) rdens1=rdens1*2.0d0
                     Call DaxPy_(nBa**2,rdens1,Temp1,1,Temp2,1)
                  End If
               End Do
            End Do
         End If
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Do kS=1,nsym
         If (nBas(jS)*nAsh(kS).ne.0) Then
*
            Do lB=nish(kS)+1,nB(kS)
               kkC=nA(kS)+lB-nIsh(kS)
               Do jB=nIsh(kS)+1,nB(kS)
                  kkb=nA(kS)+jB-nIsh(kS)
                  Call EXCH(js,ks,js,ks,jb,lb,Temp1,Scr)
                  If (lB.gt.nIsh(kS).and.jB.gt.nIsh(kS)) Then
                     rDens2=sign*4.0d0*Work(ipG2-1+
     &                      itri(itri(iib,kkc),itri(kkb,iib)))
                    Call DaXpY_(nBa**2,rDens2,Temp1,1,Temp2,1)
                  End If
               End Do
            End Do
         End If
      End Do
*
      rho=sign*2.0d0*G1t(itri(iib,iib))
      Do iI=nAsh(js)+nIsh(js)+1,nBas(js)
         rOut(ip)=rout(ip)-2.0d0*rF+Rho*FockI(iI,ii)+Temp2(ii,ii)
         ip=ip+1
         Do iJ=iI+1,Nbas(js)
            rOut(ip)=Rho*FockI(iI,iJ)+Temp2(ii,ij)
            ip=ip+1
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(focka)
      If (.False.) Call Unused_real_array(fock)
      End
