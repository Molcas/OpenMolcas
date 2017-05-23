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
      Subroutine PolPrep(iDist,iDistIm,xx,yy,zz,rr3,xxi,yyi,zzi,Gri
     &                  ,iCNum,nSize)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "qmcom.fh"
#include "WrkSpc.fh"

      Dimension xx(nSize,nSize),yy(nSize,nSize),zz(nSize,nSize)
      Dimension rr3(nSize,nSize),xxi(nSize,nSize),yyi(nSize,nSize)
      Dimension zzi(nSize,nSize),Gri(nSize,nSize)
*----------------------------------------------------------------------*
* Simply compute some vectors etc. for the ensuing polarization        *
* calculation.                                                         *
*----------------------------------------------------------------------*
      ncParm=ncent*npart-ncent*icNum
      Do 711, i=nPol*iCNum+1,nPol*nPart  !Loop over solvent polar-
        Do 712, j=nPol*iCNum+1,nPol*nPart  !ization sites.
          rr3(i,j)=0
712     Continue
711   Continue
      IndTr1=0
      !If (iCNum+1).eq.nPart this loop will not be run, but that is
      !okey since then we only have one solvent and thus it can not
      !experience any field from other solvent molecules.
      Do 721, i=1,nPol
        Indp1=i+iCnum*nPol
        Indco1=i+iCNum*nCent
        Do 722, j=iCnum+2,nPart
          IndP1=IndP1+nPol
          IndCo1=IndCo1+nCent
          IndTr=((j-(iCnum+2))*(j-(iCNum+1)))/2*nCent**2+(i-1)*nCent
          Do 723, k=1,nPol
            IndP2=k+(iCnum-1)*nPol
            IndCo2=k+(iCNum-1)*nCent
            Do 724, l=iCnum+1,j-1
              IndTr1=Indtr1+1
              Indp2=Indp2+nPol
              Indco2=Indco2+nCent
              IndTri=IndTr+(l-(iCnum+1))*nCent**2+k
              xx(Indp1,Indp2)=(Cordst(indco1,1)-Cordst(indco2,1))
     &                       *Work(iDist+indtri-1)
              yy(Indp1,Indp2)=(Cordst(indco1,2)-Cordst(indco2,2))
     &                       *Work(iDist+indtri-1)
              zz(Indp1,Indp2)=(Cordst(indco1,3)-Cordst(indco2,3))
     &                       *Work(iDist+indtri-1)
              rr3(IndP1,IndP2)=Work(iDist+Indtri-1)**3
       !Why should xx(indp2,indp1)=xx(indp1,indp2), you wonder, should
       !they not be of different sign? The answer is, it does not
       !matter. Recall that the formula for the field from an ideal
       !dipole change sign two times, thus no time in effect, when
       !the sign of the r-vector is changed.
              xx(Indp2,Indp1)=xx(Indp1,Indp2)
              yy(Indp2,Indp1)=yy(Indp1,Indp2)
              zz(Indp2,Indp1)=zz(Indp1,Indp2)
              rr3(Indp2,Indp1)=rr3(Indp1,Indp2)
724         Continue
723       Continue
722     Continue
721   Continue
      Do 725, ii=1,nSize
        Do 726, jj=1,nSize
          Gri(ii,jj)=0
726     Continue
725   Continue
      Do 731, i=1,nPol
        k=i+(iCnum-1)*nCent
        Do 732, i1=iCnum+1,nPart
          k=k+nCent
          imd=(i1-1)*nPol+i
          Do 733, j=1,nPol
            l=j+(iCnum-1)*nCent
            jnd=((i1-(iCnum+1))*nCent+i-1)*ncparm+j-nCent
            Do 734, j1=iCnum+1,nPart
              l=l+nCent
              ild=(j1-1)*nPol+j
              jnd=jnd+nCent
              Gri(imd,ild)=Work(iDistIm+jnd-1)
          xxi(imd,ild)=(Cordim(k,1)-Cordst(l,1))*Work(iDistIm+jnd-1)
          yyi(imd,ild)=(Cordim(k,2)-Cordst(l,2))*Work(iDistIm+jnd-1)
          zzi(imd,ild)=(Cordim(k,3)-Cordst(l,3))*Work(iDistIm+jnd-1)
734         Continue
733       Continue
732     Continue
731   Continue

      Return
      End
