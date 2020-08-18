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
* Copyright (C) 2000, Gunnar Karlstrom                                 *
*               2000, Roland Lindh                                     *
************************************************************************
      Subroutine lattcr(Grid,nGrid_,nGrid_Eff_,PolEff,DipEff,
     &                  cord,maxato,atorad,nPolComp,
     &                  XF,nXF,nOrd_XF,XEle,iXPolType)

************************************************************************
*                                                                      *
*     Object: to compute effective polarizabilities and dipole moments *
*             on the Langevin grid.                                    *
*                                                                      *
*     Authors: G. Karlstroem                                           *
*              Dept. of Theor. Chem., Univ. of Lund, Sweden.           *
*                                                                      *
*              and                                                     *
*                                                                      *
*              R. Lindh                                                *
*              Dept. of Chem. Phys., Univ. of Lund, Sweden.            *
*                                                                      *
*              March 2000                                              *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "rctfld.fh"
*
      Real*8 Grid(3,nGrid_), PolEff(nPolComp,nGrid_), DipEff(nGrid_)
      Real*8 cord(3,maxato), atorad(maxato),XF(*)
      Integer XEle(nXF)


      Real*8 tr(3,3),co(3)
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

*     Calculate number of entries per XFIELD point
      Inc = 3
      Do iOrdOp = 0, nOrd_XF
         Inc = Inc + nElem(iOrdOp)
      End Do
      If(iXPolType.gt.0) Inc = Inc + 6
*
*     Write (*,*) 'lattcr: polsi,dipsi=',polsi,dipsi

*
*     Rotation matrix for the grid
      tr(1,1)=cos(rotGamma)*cos(rotBeta)*cos(rotAlpha)-sin(rotGamma)
     &     *sin(rotAlpha)
      tr(1,2)=cos(rotGamma)*cos(rotBeta)*sin(rotAlpha)+sin(rotGamma)
     &     *cos(rotAlpha)
      tr(1,3)=-cos(rotGamma)*sin(rotBeta)
      tr(2,1)=-sin(rotGamma)*cos(rotBeta)*cos(rotAlpha)-cos(rotGamma)
     &     *sin(rotAlpha)
      tr(2,2)=-sin(rotGamma)*cos(rotBeta)*sin(rotAlpha)+cos(rotGamma)
     &     *cos(rotAlpha)
      tr(2,3)= sin(rotGamma)*sin(rotBeta)
      tr(3,1)=sin(rotBeta)*cos(rotAlpha)
      tr(3,2)=sin(rotBeta)*sin(rotAlpha)
      tr(3,3)=cos(rotBeta)

      Do ii=-(maxa+1),maxa,nSparse
       Do jj=-(maxb+1),maxb,nSparse
        Do kk=-(maxc+1),maxc,nSparse
         ni=min(nSparse,maxa-ii+1)
         nj=min(nSparse,maxb-jj+1)
         nk=min(nSparse,maxc-kk+1)
         If(LSparse.and.(ni.eq.nSparse).and.(nj.eq.nSparse).and.
     &        (nk.eq.nSparse)) Then
            xs=(DBLE(ii)+(DBLE(nSparse-1)*0.5D0))*scala
            ys=(DBLE(jj)+(DBLE(nSparse-1)*0.5D0))*scala
            zs=(DBLE(kk)+(DBLE(nSparse-1)*0.5D0))*scalc
            nGridOld=nGrid_Eff
            Do l=1,latato
               co(1)=xs+cordsi(1,l)
               co(2)=ys+cordsi(2,l)
               co(3)=zs+cordsi(3,l)
               xp=Zero
               yp=Zero
               zp=Zero
               Do m=1,3
                  xp=xp+co(m)*tr(1,m)
                  yp=yp+co(m)*tr(2,m)
                  zp=zp+co(m)*tr(3,m)
               EndDo
               rp2=xp*xp+yp*yp+zp*zp
               If (rp2.gt.radlat**2) go to 13
               If (lRFCav) Then
                  rp=sqrt(rp2)
                  drp=rds-rp
                  If (drp.le.Zero) go to 13
                  ener1=(diedel/drp)**2
               Else
                  ener1=Zero
               EndIf
               ener=Zero
*
*------- Check if the QC system annhilates the grid point
*
               Do m=1,maxato
                  xa=cord(1,m)
                  ya=cord(2,m)
                  za=cord(3,m)
                  rrr=(atorad(m)*rsca)**2
                  dggx=xa-xp
                  dggy=ya-yp
                  dggz=za-zp
                  rpa2=dggx**2+dggy**2+dggz**2
                  if(rpa2.lt.distSparse**2) Goto 13
                  ener=ener+(rrr/rpa2)**nexpo
               End Do
*
*------- Check if the XFIELD multipoles annhilates the grid point
*
               Do iXF=1,nXF
                  xa=XF((iXF-1)*Inc+1)
                  ya=XF((iXF-1)*Inc+2)
                  za=XF((iXF-1)*Inc+3)
                  If(XEle(iXF).le.0) Then
                     atrad=-DBLE(XEle(iXF))/1000.0D0
                  Else
                     atrad=CovRadT(XEle(iXF))
                  EndIf
                  rrr=(atrad*rsca)**2
                  dggx=xa-xp
                  dggy=ya-yp
                  dggz=za-zp
                  rpa2=dggx**2+dggy**2+dggz**2
                  if(rpa2.lt.distSparse**2) Goto 13
                  ener=ener+(rrr/rpa2)**nexpo
c     If(rpa2.lt.6.0D0)
c     &           write(*,*)'DIST',iGrid,iXF,sqrt(rpa2),atrad
               EndDo

               ener=prefac*ener*500.0D0+ener1
               If (ener.gt.12.D0) Then
                  GoTo 13
               EndIf
               fac=exp(-ener)
               nGrid_Eff_=nGrid_Eff_+1
               Grid(1,nGrid_Eff_)=xp
               Grid(2,nGrid_Eff_)=yp
               Grid(3,nGrid_Eff_)=zp
               PolEff(1,nGrid_Eff_)=polsi*DBLE(nSparse)**3.0D0*fac*fac
               DipEff(nGrid_Eff_)=dipsi*DBLE(nSparse)**1.5D0*fac
               write(6,*)'DGRID',xp,yp,zp,fac
            EndDo
*     Every sparse lattice point in this cell is ok, so skip dense grid

            Goto 14

 13         continue
*     Not every sparse lattice point is ok, so delete the sparse and do dense
            nGrid_Eff_=nGridOld
         EndIf
*     Start the normal, dense grid
         Do i=0,ni-1
            Do j=0,nj-1
               Do k=0,nk-1
                  xs=DBLE(ii+i)*scala
                  ys=DBLE(jj+j)*scala
                  zs=DBLE(kk+k)*scalc
                  Do l=1,latato
                     co(1)=xs+cordsi(1,l)
                     co(2)=ys+cordsi(2,l)
                     co(3)=zs+cordsi(3,l)
                     xp=Zero
                     yp=Zero
                     zp=Zero
                     Do m=1,3
                        xp=xp+co(m)*tr(1,m)
                        yp=yp+co(m)*tr(2,m)
                        zp=zp+co(m)*tr(3,m)
                     EndDo
                     rp2=xp*xp+yp*yp+zp*zp
                     If (rp2.gt.radlat**2) go to 11
                     If (lRFCav) Then
                        rp=sqrt(rp2)
                        drp=rds-rp
                        If (drp.le.Zero) go to 11
                        ener1=(diedel/drp)**2
                     Else
                        ener1=Zero
                     EndIf
                     ener=Zero
*
*------- Check if the QC system annhilates the grid point
*
                     Do m=1,maxato
                        xa=cord(1,m)
                        ya=cord(2,m)
                        za=cord(3,m)
                        rrr=(atorad(m)*rsca)**2
                        dggx=xa-xp
                        dggy=ya-yp
                        dggz=za-zp
                        rpa2=dggx**2+dggy**2+dggz**2
                        ener=ener+(rrr/rpa2)**nexpo
c     If(rpa2.lt.6.0D0)
c     &           write(*,*)'DIST0',iGrid,sqrt(rpa2),atorad(m)
                     End Do

*
*------- Check if the XFIELD multipoles annhilates the grid point
*
                     Do iXF=1,nXF
                        xa=XF((iXF-1)*Inc+1)
                        ya=XF((iXF-1)*Inc+2)
                        za=XF((iXF-1)*Inc+3)
                        If(XEle(iXF).le.0) Then
                           atrad=-DBLE(XEle(iXF))/1000.0D0
                        Else
                           atrad=CovRadT(XEle(iXF))
                        EndIf
                        rrr=(atrad*rsca)**2
                        dggx=xa-xp
                        dggy=ya-yp
                        dggz=za-zp
                        rpa2=dggx**2+dggy**2+dggz**2
                        ener=ener+(rrr/rpa2)**nexpo
c     If(rpa2.lt.6.0D0)
c     &           write(*,*)'DIST',iGrid,iXF,sqrt(rpa2),atrad
                     EndDo

c         ener=prefac*ener*tk5+ener1
                     ener=prefac*ener*500.0D0+ener1
c         ener=prefac*ener*1000.0D0+ener1
                     If (ener.gt.12.D0) Then
                        Write(6,*)'REMOVED',xp,yp,zp
                        Go To 11
                     EndIf
*
                     fac=exp(-ener)
*
                     nGrid_Eff_=nGrid_Eff_+1
                     Grid(1,nGrid_Eff_)=xp
                     Grid(2,nGrid_Eff_)=yp
                     Grid(3,nGrid_Eff_)=zp
                     PolEff(1,nGrid_Eff_)=polsi*fac*fac
                     DipEff(nGrid_Eff_)=dipsi*fac
                     write(6,*)'GRID',xp,yp,zp,fac
 11                  Continue
                  EndDo  ! l
               EndDo  ! k
            EndDo  ! j
         End Do  ! i
 14      Continue
        EndDo  ! kk
       EndDo  ! jj
      EndDo  ! ii
*
*     Write (*,*) 'Grid...',Grid(1,1),Grid(2,1),Grid(3,1)
*     Write (*,*) 'DipEff...',DipEff(1),DipEff(101)
*     Write (*,*) 'PolEff...',PolEff(1),PolEff(101)
*     Write (*,*) 'Lattcr: nGrid_Eff=',nGrid_Eff
*     Call GetMem('Lattcr','Check','Real',idum,idum)
      Return
      End
