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
      SubRoutine Precabb_2(ib,is,js,nd,nba,no,rout,Temp1,ntemp,
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
      use Arrays, only: G1t, G2t
      Implicit Real*8(a-h,o-z)
#include "Input.fh"
#include "Pointers.fh"
      Real*8 rout(*)
      Real*8 Focki(no,no),Focka(no,no),
     &       Fock(no,no)
      Real*8 Temp1(nTemp),Temp2(nO,nO), Scr(nTemp)
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
      jVert=nOrb(js)-nAsh(js)-nIsh(js)
      If (jvert.eq.0) Return
*
      i2=nD-jVert+1
      ip=iTri1(i2,i2)
      rF=sign*Fockti
      call dcopy_(nBa**2,[0.0d0],0,Temp2,1)
*
      Do kS=1,nSym
        iOpt=1
        ijB=1
        ijBas=0
        ijBB=0
        If (nOrb(js)*nash(ks).gt.0) Then

        Do kBB=nish(ks)+1,nB(kS)
         Do kCC=nish(ks)+1,kBB
              Call COUL(jS,jS,kS,kS,kbb,kcc,Temp1,Scr)
              ipT=1

         If (kBB.gt.nish(ks).and.kCC.gt.nish(ks)) Then
           kkB=kBB+nA(ks)-nish(ks)
           kkC=kCC+nA(ks)-Nish(ks)
           rDens1=sign*2.0d0*G2t(
     &            itri(itri(iib,iib),itri(kkb,kkc)))
*
           If (kbb.ne.kcc) rdens1=rdens1*2.0d0

           Call DaxPy_(nO**2,rdens1,Temp1,1,Temp2,1)

         End If
        End Do
       End Do
       End If
      End Do
*
      Do Ks=1,nsym
       iOpt=1
       JLB=1
       JLBas=0
       ijkl=nOrb(js)*nash(ks)
       If (ijkl.ne.0) Then
*

        jlBB=0
        Do LB=nish(ks)+1,nB(KS)
         kkc=nA(ks)+lb-nish(ks)
         Do JB=nish(ks)+1,nB(KS)
          kkb=nA(ks)+jb-nish(ks)
          Call EXCH(js,ks,js,ks,jb,lb,Temp1,Scr)
          ipT=1
          If (LB.gt.nISH(ks).and.jb.gt.nish(ks)) Then
           rDens2=sign*4.0d0*G2t(
     &         itri(itri(iib,kkc),itri(kkb,iib)))
           Call DaXpY_(nO**2,rDens2,Temp1(ipT),1,Temp2,1)
          End If
         End Do
        End Do
       End If
      End Do

      iu=ip

      rho=sign*2.0d0*G1t(itri(iib,iib))
      Do iI=nAsh(js)+nIsh(js)+1,nOrb(js)
       rOut(ip)=rout(ip)-2.0d0*rF+Rho*FockI(iI,ii)+Temp2(ii,ii)
       ip=ip+1
       Do iJ=iI+1,NOrb(js)
        rOut(ip)=Rho*FockI(iI,iJ)+Temp2(ii,ij)
        ip=ip+1
       End Do
      End Do
      return
c Avoid unused argument warnings
      If (.False.) Then
       Call Unused_real_array(focka)
       Call Unused_real_array(fock)
      End If
      end
