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
      subroutine intergeo(FileName,Enrg,Crd,Grd,nAtm,nIter)
*---------------------------------*
* Add geometry optimization info  *
*   to the Molden inputfile       *
*---------------------------------*
      implicit real*8 (a-h,o-z)
#include "info_slapaf.fh"
#include "WrkSpc.fh"
#include "angstr.fh"
#include "periodic_table.fh"
#include "stdalloc.fh"
      Real*8 Charge(Mxdc), Crd(3,nAtm,nIter),
     &       Enrg(nIter), Grd(3,nAtm,nIter)
      Integer iPhase(3,0:7),nStab2(Mxdc),
     &        iChCar(3)
      Integer, Allocatable :: icoset2(:,:,:)
      Character*(*) FileName
*
      Call QEnter('InterGeo')
*                                                                      *
************************************************************************
*                                                                      *
*     Pick information for centers and pseudo centers
*
      Call Get_iScalar('Unique atoms',msAtom)
      Call Get_dArray('Nuclear charge',Charge,msAtom)
*
      Call Get_iScalar('Pseudo atoms',msAtom_p)
      If (msAtom_p.gt.0) Then
         Call Get_dArray('Pseudo charge',Charge(msAtom+1),msAtom_p)
         Call GetMem('Coor_p','Allo','Real',ipCx_p,3*msAtom_p)
         Call Get_dArray('Pseudo Coordinates',Work(ipCx_p),3*msAtom_p)
      Else
         ipCx_p=ip_Dummy
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (msAtom.gt.nAtm) Then
         Call WarningMessage(2,'Error in InterGEO')
         Write (6,*) 'msAtom.gt.nAtm'
         Write (6,*) 'msAtom=',msAtom
         Write (6,*) 'nAtm=',nAtm
         Call Abend()
      End If
*
      Lu_Molden=19
      Call molcas_open(Lu_Molden,FileName)
      Write (Lu_Molden,*) '[MOLDEN FORMAT]'
      Write (Lu_Molden,*) '[N_GEO]'
      Write (Lu_Molden,*) nIter
      Write (Lu_Molden,*) '[GEOCONV]'
      Write (Lu_Molden,*) 'energy'
*
      iEner=0
      Do iIter=1,nIter
         Write (Lu_Molden,'(E24.17)') Enrg(iIter)
         iEner=iEner+1
      End Do
*
      Write (Lu_Molden,*) 'max-force'
      iGx=0
      Do iIter=1,nIter
        grmax=0.0D0
        Do ndc=1,msAtom
           grx=abs(Grd(1,ndc,iIter))
           gry=abs(Grd(2,ndc,iIter))
           grz=abs(Grd(3,ndc,iIter))
           If (grx.gt.grmax) grmax=grx
           If (gry.gt.grmax) grmax=gry
           If (grz.gt.grmax) grmax=grz
           iGx=iGx+3
        End Do
        Write (Lu_Molden,'(F12.7)') grmax
      End Do
*
      Write (Lu_Molden,*) 'rms-force'
      iGx=0
      Do iIter=1,nIter
         grtot=0.0D0
         ngrad=0
         Do ndc=1,msAtom
            grx=Grd(1,ndc,iIter)
            gry=Grd(2,ndc,iIter)
            grz=Grd(3,ndc,iIter)
            Do i=0,nSym/nStab(ndc)-1
               grtot=grtot+grx*grx+gry*gry+grz*grz
               ngrad=ngrad+1
            End Do
            iGx=iGx+3
         End Do
         Write (Lu_Molden,'(F12.7)') sqrt(grtot)/DBLE(ngrad)
      End Do
*
*     This part disabled because gv refuses to open the file
#ifdef write_molden_steps
      Write (Lu_Molden,*) 'max-step'
      Do iIter=1,nIter-1
         stepmax=0.0D0
         Do ndc=1,msAtom
            dx=Crd(1,ndc,iIter+1)-Crd(1,ndc,iIter)
            dy=Crd(2,ndc,iIter+1)-Crd(2,ndc,iIter)
            dz=Crd(3,ndc,iIter+1)-Crd(3,ndc,iIter)
            If (dx.gt.stepmax) stepmax=dx
            If (dy.gt.stepmax) stepmax=dy
            If (dz.gt.stepmax) stepmax=dz
         End Do
         Write (Lu_Molden,'(F12.7)') stepmax
      End Do
*
      Write (Lu_Molden,*) 'rms-step'
      Do iIter=1,nIter-1
         step=0.0D0
         Do ndc=1,msAtom
            dx=Crd(1,ndc,iIter+1)-Crd(1,ndc,iIter)
            dy=Crd(2,ndc,iIter+1)-Crd(2,ndc,iIter)
            dz=Crd(3,ndc,iIter+1)-Crd(3,ndc,iIter)
            Do i=0,nSym/nStab(ndc)-1
               step=step+dx*dx+dy*dy+dz*dz
            End Do
         End Do
         Write (Lu_Molden,'(F12.7)') sqrt(step)/DBLE(ngrad)
      End Do
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Set up the desymmetrization of the coordinates
*
      Do 10 j=0,4
         Do 10 i=1,3
10          iPhase(i,j)=1
      iPhase(1,1)=-1
      iPhase(2,2)=-1
      iPhase(1,3)=-1
      iPhase(2,3)=-1
      iPhase(3,4)=-1
      Do 20 j=5,7
         Do 20 i=1,3
20          iPhase(i,j)=-1
      iPhase(2,5)=1
      iPhase(1,6)=1
*
      iSymX = 0
      iSymY = 0
      iSymZ = 0
      Do i = 0, nSym-1
         If (iAnd(iOper(i),1).ne.0) iSymX = 1
         If (iAnd(iOper(i),2).ne.0) iSymY = 2
         If (iAnd(iOper(i),4).ne.0) iSymZ = 4
      End Do
      iChCar(1) = iSymX
      iChCar(2) = iSymY
      iChCar(3) = iSymZ
*
      ixyz   = ipCx
      ixyz_p = ipCx_p
      MaxDCR=0
      Call mma_allocate(icoset2,[0,7],[0,7],[1,msAtom+msAtom_p],
     &                  label='icoset2')
      Do ndc = 1, msAtom + msAtom_p
         If (ndc.le.msAtom) Then
            iChxyz=iChAtm(Work(ixyz  ),iOper,nSym,iChCar)
            ixyz   = ixyz   + 3
         Else
            iChxyz=iChAtm(Work(ixyz_p),iOper,nSym,iChCar)
            ixyz_p = ixyz_p + 3
         End If
         Call Stblz(iChxyz,iOper,nSym,nStab2(ndc),jStab(0,ndc),
     &              MaxDCR,iCoSet2(0,0,ndc))
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      nAt=0
      Write (Lu_Molden,*) '[GEOMETRIES] (XYZ)'
      Do ndc = 1, msAtom+msAtom_p
         Do i=0,nSym/nStab2(ndc)-1
            nAt=nAt+1
         End do
      End do
*
      Do iIter=1,nIter
         Write (Lu_Molden,'(I4)') nAt
         Write (Lu_Molden,*) Enrg(iIter)
         Do ndc = 1, msAtom+msAtom_p
            If (ndc.le.msAtom) Then
               x=Crd(1,ndc,iIter)
               y=Crd(2,ndc,iIter)
               z=Crd(3,ndc,iIter)
            Else
               ixyz_p = ipCx_p + 3*(ndc-msAtom-1)
               x=Work(ixyz_p  )
               y=Work(ixyz_p+1)
               z=Work(ixyz_p+2)
            End If
            Do i=0,nSym/nStab2(ndc)-1
               iFacx=iPhase(1,icoset2(i,0,ndc))
               iFacy=iPhase(2,icoset2(i,0,ndc))
               iFacz=iPhase(3,icoset2(i,0,ndc))
               x1=angstr*x*DBLE(iFacx)
               y1=angstr*y*DBLE(iFacy)
               z1=angstr*z*DBLE(iFacz)
               LbAtom=int(charge(ndc))
               Write (Lu_Molden,102) pTab(LbAtom),x1,y1,z1
 102           Format(A2,3(3x,F12.7))
            End do
         End do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Write (Lu_Molden,*) '[FORCES]'
*
      Do iIter=1,nIter
         Write (Lu_Molden,'(A,1X,I4)') 'point',iIter
         Write (Lu_Molden,'(I4)') nAt
         Do ndc = 1, msAtom+msAtom_p
            If (ndc.le.msAtom) Then
               x=Grd(1,ndc,iIter)
               y=Grd(2,ndc,iIter)
               z=Grd(3,ndc,iIter)
            Else
               x=0.0D0
               y=0.0D0
               z=0.0D0
            End If
            Do i=0,nSym/nStab2(ndc)-1
               iFacx=iPhase(1,icoset2(i,0,ndc))
               iFacy=iPhase(2,icoset2(i,0,ndc))
               iFacz=iPhase(3,icoset2(i,0,ndc))
               x1=x/angstr*DBLE(iFacx)
               y1=y/angstr*DBLE(iFacy)
               z1=z/angstr*DBLE(iFacz)
               Write (Lu_Molden,103) x1,y1,z1
 103           Format(3(3x,F12.7))
            End do
         End do
      End Do
      Call mma_deallocate(icoset2)
*
      Close(Lu_Molden)
*                                                                      *
************************************************************************
*                                                                      *
      If (msAtom_p.gt.0) Then
         Call GetMem('Coor_p','Free','Real',ipCx_p,3*msAtom_p)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('InterGeo')
      Return
      End
