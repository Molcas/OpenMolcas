************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine ClasClas(iCNum,iCStart,ncParm,Coord,iFP,iGP,iDT,iFI
     &                ,iDist,iDistIm,Elene,Edisp,Exrep,E2Die,ExDie)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "qmcom.fh"
#include "WrkSpc.fh"

      Dimension iFP(3),iGP(3),iDT(3),iFI(3)
      Dimension Coord(MxAt*3)
      Character Memlabel*20,Memlaabe*20,Memlaaab*20,MemLaaaa*20,ChCo*2
      Parameter (ExLim=10)

*----------------------------------------------------------------------*
* Compute the distance matrices between the classical centers and the  *
* classical image centers.                                             *
*----------------------------------------------------------------------*
      nClas=nPart-iCNum
      Adisp=Disp(1,2)
      nSize=(nClas*(nClas-1)/2)*(nCent**2)   !Get memory
      Call GetMem('DistMat','Allo','Real',iDist,nSize)
      nSizeIm=(nClas*nCent)**2
      Call GetMem('DistMatIm','Allo','Real',iDistIm,nSizeIm)
      Ind=0
      Do 201, ii=iCNum+2,nPart
        Do 202, jj=iCNum+1,ii-1
          Do 203, ij=1,nCent
            i=(ii-1)*nCent+ij
            Do 204, k=1,nCent
              Ind=Ind+1
              j=(jj-1)*nCent+k
              r=0
              Do 205, l=1,3
                r=(Cordst(i,l)-Cordst(j,l))**2+r
205           Continue
              Work(iDist+Ind-1)=1/Sqrt(r)
204         Continue
203       Continue
202     Continue
201   Continue
      Jnd=0
      Do 206, i=iCStart,nCent*nPart
        Do 207, j=iCStart,nCent*nPart
          Jnd=Jnd+1
          r=0
          Do 208, k=1,3
            r=(CordIm(i,k)-Cordst(j,k))**2+r
208       Continue
          Work(iDistIm+Jnd-1)=1.0d0/Sqrt(r)
207     Continue
206   Continue
*----------------------------------------------------------------------*
* Compute the pairwise interaction between the solvent. Classical all  *
* the  way... early NEMO all the way.                                  *
*----------------------------------------------------------------------*
      Elene=0
      Edisp=0
      Exrep=0
      aLim=1.0d0/ExLim
      Sum1=0
      Sum2=0
      Sum3=0
      Sum4=0
*The electrostatic part
      Do 2011, i=1,nSize,nCent**2  !This loop ONLY works for the
                     !early Nemo model of water. If the solvent
                   !model is changed this loop must be rewritten.
        Sum1=Sum1+Work(iDist+i-1+6)*Qsta(1)*Qsta(1)!H-H
        Sum2=Sum2+Work(iDist+i-1+7)*Qsta(1)*Qsta(2)!H-H
        Sum3=Sum3+Work(iDist+i-1+8)*Qsta(1)*Qsta(3)!H-V
        Sum4=Sum4+Work(iDist+i-1+9)*Qsta(1)*Qsta(4)!H-V
        Sum1=Sum1+Work(iDist+i-1+11)*Qsta(2)*Qsta(1)!H-H
        Sum2=Sum2+Work(iDist+i-1+12)*Qsta(2)*Qsta(2)!H-H
        Sum3=Sum3+Work(iDist+i-1+13)*Qsta(2)*Qsta(3)!H-V
        Sum4=Sum4+Work(iDist+i-1+14)*Qsta(2)*Qsta(4)!H-V
        Sum1=Sum1+Work(iDist+i-1+18)*Qsta(3)*Qsta(3)!V-V
        Sum2=Sum2+Work(iDist+i-1+19)*Qsta(3)*Qsta(4)!V-V
        Sum3=Sum3+Work(iDist+i-1+16)*Qsta(3)*Qsta(1)!V-H
        Sum4=Sum4+Work(iDist+i-1+17)*Qsta(3)*Qsta(2)!V-H
        Sum1=Sum1+Work(iDist+i-1+23)*Qsta(4)*Qsta(3)!V-V
        Sum2=Sum2+Work(iDist+i-1+24)*Qsta(4)*Qsta(4)!V-V
        Sum3=Sum3+Work(iDist+i-1+21)*Qsta(4)*Qsta(1)!V-H
        Sum4=Sum4+Work(iDist+i-1+22)*Qsta(4)*Qsta(2)!V-H
2011  Continue
      Elene=Sum1+Sum2+Sum3+Sum4

      Sum1=0
      Sum2=0
      Sum3=0
      Sum4=0
      Sum5=0
*The dispersion, now with damping.
      Do 211, i=1,nSize,nCent**2
        DampFunk=1-Exp(-1.0d0/(Work(iDist+i-1)*2.2677d0))**4
        Sum1=Sum1+Work(iDist+i-1)**6*DampFunk
        DampFunk=1-Exp(-1.0d0/(Work(iDist+i-1+1)*2.2677d0))**4
        Sum2=Sum2+Work(iDist+i-1+1)**6*DampFunk
        DampFunk=1-Exp(-1.0d0/(Work(iDist+i-1+2)*2.2677d0))**4
        Sum3=Sum3+Work(iDist+i-1+2)**6*DampFunk
        DampFunk=1-Exp(-1.0d0/(Work(iDist+i-1+11)*2.2677d0))**4
        Sum4=Sum4+Work(iDist+i-1+11)**6*DampFunk
        DampFunk=1-Exp(-1.0d0/(Work(iDist+i-1+7)*2.2677d0))**4
        Sum5=Sum5+Work(iDist+i-1+7)**6*DampFunk
        DampFunk=1-Exp(-1.0d0/(Work(iDist+i-1+5)*2.2677d0))**4
        Sum2=Sum2+Work(iDist+i-1+5)**6*DampFunk
        DampFunk=1-Exp(-1.0d0/(Work(iDist+i-1+10)*2.2677d0))**4
        Sum3=Sum3+Work(iDist+i-1+10)**6*DampFunk
        DampFunk=1-Exp(-1.0d0/(Work(iDist+i-1+6)*2.2677d0))**4
        Sum4=Sum4+Work(iDist+i-1+6)**6*DampFunk
        DampFunk=1-Exp(-1.0d0/(Work(iDist+i-1+12)*2.2677d0))**4
        Sum5=Sum5+Work(iDist+i-1+12)**6*DampFunk
211   Continue
      Edisp=Sum1*Disp(1,1)+(Sum2+Sum3)*Disp(1,2)+(Sum4+Sum5)*Disp(2,2)
*The exchange repulsion
      Do 221, i=1,nSize,nCent**2
        If(Work(iDist+i-1).gt.aLim) Exrep=Exrep+
     &                          ExNemo(1,1,Work(iDist+i-1))
        If(Work(iDist+i-1+1).gt.aLim) Exrep=Exrep+
     &                          ExNemo(1,2,Work(iDist+i-1+1))
        If(Work(iDist+i-1+2).gt.aLim) Exrep=Exrep+
     &                          ExNemo(1,2,Work(iDist+i-1+2))
        If(Work(iDist+i-1+5).gt.aLim) Exrep=Exrep+
     &                          ExNemo(1,2,Work(iDist+i-1+5))
        If(Work(iDist+i-1+6).gt.aLim) Exrep=Exrep+
     &                          ExNemo(2,2,Work(iDist+i-1+6))
        If(Work(iDist+i-1+7).gt.aLim) Exrep=Exrep+
     &                          ExNemo(2,2,Work(iDist+i-1+7))
        If(Work(iDist+i-1+10).gt.aLim) Exrep=Exrep+
     &                          ExNemo(1,2,Work(iDist+i-1+10))
        If(Work(iDist+i-1+11).gt.aLim) Exrep=Exrep+
     &                          ExNemo(2,2,Work(iDist+i-1+11))
        If(Work(iDist+i-1+12).gt.aLim) Exrep=Exrep+
     &                          ExNemo(2,2,Work(iDist+i-1+12))
221   Continue
*----------------------------------------------------------------------*
* Compute pair-wise interaction with image charges.                    *
*----------------------------------------------------------------------*
      Sum1=0.0d0
      Sum2=0.0d0
      Do 231, i=iCNum+1,nPart
        Do 232, j=nCent-nCha+1,nCent  !Only count over
                                   !charged centers.
          Q1=QIm((i-1)*nCent+j)  !The image charge.
          Inc=ncParm*nCent*(i-(iCNum+1))+(j-1)*ncParm !Counting
                                                  !elements.
          Do 233, k=nCent-nCha+1,nCent
            Inc2=Inc+k
            Q2=QSta(k-nCent+nCha)
            Do 234, l=iCNum+1,nPart !Here is the electrostatic
           !interaction computed. Observe the difference with the real
              !charges, since here interaction between ALL
             !real-image charge pair is computed.
            Sum1=Sum1+Q1*Q2*Work(iDistIm+Inc2+(l-(iCnum+1))*nCent-1)
         Sum1=Sum1-Adisp*Work(iDistIm+Inc2+(l-(iCnum+1))*nCent-1)**6
     &       *(1-Exp(-1/(Work(iDistIm+Inc2+(l-(iCnum+1))*nCent-1)
     &       *2.9677d0)**6))
234         Continue
233       Continue
232     Continue
* Include a repulsion with the boundary to prevent the waters to merge
*into the dielectric continuum. Its construction is such that the
*repulsion only is between the particle and the image of the particle,
*no other repulsion over the boundary.
        Sum2=Sum2+ExNemo(1,2,Work(iDistIm-1+(i-(iCNum+1))
     &       *nClas*nCent**2
     &       +(nClas+i-(iCNum+1))*nCent+2))
        Sum2=Sum2+ExNemo(1,2,Work(iDistIm-1+(i-(iCNum+1))
     &       *nClas*nCent**2
     &       +(2*nClas+i-(iCNum+1))*nCent+3))
        Sum2=Sum2+ExNemo(1,1,Work(iDistIm-1+(i-(iCNum+1))
     &       *nClas*nCent**2
     &       +1+(i-(iCNum+1))*nCent))*Exdt1
231   Continue
      E2Die=sum1*0.5d0  !The half is added since what we actually has
          !computed is the interaction between charge and a part
          !of its reaction field (recall:0.5*q*fi_q).
      EXDie=sum2*0.5d0*ExdTal
*----------------------------------------------------------------------*
* Compute the static electric field on the polarizabilities and obtain *
* initial guess of induced dipoles.                                    *
*----------------------------------------------------------------------*
      IndMa=nPol*nPart
      Do 300,i=1,3   !Allocate memory
        Write(ChCo,'(I2.2)')i
        Write(MemLabel,*)'FP'//ChCo
        Write(MemLaabe,*)'GP'//ChCo
        Write(MemLaaab,*)'DT'//ChCo
        Write(MemLaaaa,*)'FI'//ChCo
        !Explanation: iFP-field iGP plus reaction field,
        !      iGP-field from real charges on polarizable centers,
        !      iFi-induced field.
        Call GetMem(MemLabel,'Allo','Real',iFP(i),IndMa)
        Call GetMem(MemLaabe,'Allo','Real',iGP(i),IndMa)
        Call GetMem(MemLaaab,'Allo','Real',iDT(i),IndMa)
        Call GetMem(MemLaaaa,'Allo','Real',iFi(i),IndMa)
300   Continue
      Do 302, j=1,3
        Do 301, i=0,Indma-1  !Set some zeros
          Work(iFI(j)+i)=0.0d0
          Work(iGP(j)+i)=0.0d0
          Work(iDT(j)+i)=0.0d0
          Work(iFP(j)+i)=0.0d0
301     Continue
302   Continue
*Real centers: The field at the polarizabilities - no reaction field.
      Ind=0
      Do 310, ii=iCNum+2,nPart
        Do 311, jj=iCNum+1,ii-1
          Do 312, ij=1,nCent
            i=(ii-1)*nCent+ij
            Do 313, l=1,nCent
              j=(jj-1)*nCent+l
              Ind=Ind+1
              X=Cordst(i,1)-Cordst(j,1)
              Y=Cordst(i,2)-Cordst(j,2)
              Z=Cordst(i,3)-Cordst(j,3)
              ri=Work(iDist+Ind-1)**3
              If(ij.gt.(nCent-nCha).and.l.le.nPol) then
                          !Given that ij is
                          !counting on centers with charges and
                          !l on polarizable centers, then compute
                         !the field from charge ij on center l.
                Q1=Qsta(ij-nCent+nCha)
                Ind1=(jj-1)*nPol+l
                Work(iGP(1)+Ind1-1)=Work(iGP(1)+Ind1-1)+x*Q1*ri
                Work(iGP(2)+Ind1-1)=Work(iGP(2)+Ind1-1)+y*Q1*ri
                Work(iGP(3)+Ind1-1)=Work(iGP(3)+Ind1-1)+z*Q1*ri
              Endif
              If(l.gt.(nCent-nCha).and.ij.le.nPol) then  !If ij is
                       !on center with polarizability and l is on
                       !center with charge, then compute the field
                       !from charge l on center ij.
                Q2=Qsta(l-nCent+nCha)
                Ind1=(ii-1)*nPol+ij
                Work(iGP(1)+Ind1-1)=Work(iGP(1)+Ind1-1)-x*Q2*ri
                Work(iGP(2)+Ind1-1)=Work(iGP(2)+Ind1-1)-y*Q2*ri
                Work(iGP(3)+Ind1-1)=Work(iGP(3)+Ind1-1)-z*Q2*ri
              Endif
313         Continue
312       Continue
311     Continue
310   Continue
      Epoll=0
      Do 320, i=1+nPol*iCNum,IndMa !Compute polarization energy.
                              !This is only for checking, and will
                              !not enter the energy expression.
        k=i-((i-1)/nPol)*nPol
        Work(iFP(1)+i-1)=Work(iGP(1)+i-1)
        Work(iFP(2)+i-1)=Work(iGP(2)+i-1)
        Work(iFP(3)+i-1)=Work(iGP(3)+i-1)
        F=(Work(iFP(1)+i-1)**2+Work(iFP(2)+i-1)**2+Work(iFP(3)+i-1)**2)
     &   *Pol(k)
        Epoll=Epoll+F
320   Continue
      Epoll=-Epoll*0.5d0
*Image centers: The field at the polarizabilities - reaction field to the
*point charges added.
      Do 330, i=iCStart,nCent*nPart
        Q=Qim(i)
        Do 331, k=1,nPol
          indSep=(i-iCStart)*nCent*(nPart-iCNum)+k
          indR=k+iCNum*nCent
          indF=k+iCNum*nPol
          Do 332, j=nPol*(iCNum+1),IndMa,nPol
            x=CordIm(i,1)-Cordst(indR,1)
            y=CordIm(i,2)-Cordst(indR,2)
            z=CordIm(i,3)-Cordst(indR,3)
            r3=Work(iDistIm-1+indSep)**3
            Work(iFP(1)+indF-1)=Work(iFP(1)+indF-1)+x*Q*r3
            Work(iFP(2)+indF-1)=Work(iFP(2)+indF-1)+y*Q*r3
            Work(iFP(3)+indF-1)=Work(iFP(3)+indF-1)+z*Q*r3
            indR=indR+nCent
            indSep=indSep+nCent
            indF=indF+nPol
332       Continue
331     Continue
330   Continue
*We obtain an initial guess of the induced dipoles on the solvent.
      Do 340, i=1+nPol*iCnum,IndMa
        k=i-((i-1)/nPol)*nPol
        Do 341, l=1,3
          Work(iDt(l)+i-1)=Work(iFP(l)+i-1)*Pol(k)
341     Continue
340   Continue

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Coord)
      End
