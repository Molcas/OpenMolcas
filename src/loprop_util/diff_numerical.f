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
      Subroutine Diff_Numerical(nAt,nB,ipMP,ipC,nij,EC,iANr
     &                         ,ip_Ttot,ip_Ttot_Inv
     &                         ,lMax,iTP,dLimmo
     &                         ,Thrs1,Thrs2,nThrs,iPrint
     &                         ,ThrsMul,Pot_Expo,Pot_Point,Pot_Fac
     &                         ,Diffed)
      Implicit real*8 (a-h,o-z)

#include "WrkSpc.fh"

      Dimension EC(3,nij),dMullig((lMax*(lMax**2+6*lMax+11)+6)/6)
      Dimension A(2),dLimmo(2),Pot_Expo(nij*2),Pot_Point(nij)
      Dimension Pot_Fac(nij*4)
      Dimension iANr(nAt)

      Logical AboveThr,Diffed(nij*2),AboveMul(2)

      Character*50 UtChar
      Character*10 OneFile

*
*-- Pick up some auxiliary stuff.
*
      Write(OneFile,'(A)')'ONEINT'
      Call Diff_Aux1(nEPP,ipEPCo,nB,OneFile)
      Call GetMem('BasIndCent','Allo','Inte',ip_Center,nB)
      Call Get_iArray('Center Index',iWork(ip_Center),nB)
      Call GetMem('PickPoints','Allo','Inte',ipPick,nEPP)
      Call GetMem('DistPick','Allo','Real',ipDPick,nEPP)
      nB2=nB*(nB+1)/2

*
*-- Do a 'clever' determination of the threshold for the multipole
*   magnitude.
*
*      Call Diff_ThrsMul(ipMP,ThrsMul,ThrsMul_Clever,nAt,nij,lMax)
      ThrsMul_Clever=ThrsMul

*
*-- Run a numerical fit for each active centre.
*
      kauntA=0
      nAbove=0
      Do iAtom=1,nAt
        ii=iAtom*(iAtom+1)/2
        Do jAtom=1,iAtom
          ij=iAtom*(iAtom-1)/2+jAtom
          jj=jAtom*(jAtom+1)/2
*
*-- Pick up the nuclei+core charge.
*
          If(iAtom.eq.jAtom) then
            chPoint=Work(iTP+iAtom-1)
          Else
            chPoint=0.0d0
          Endif
*
*---- Pick out the multipole, the prefactors. If none is above a
*     certain threshold, then it is not meaningful to make them
*     diffuse. Also check which individual multipoles that should
*     be made diffuse.
*
          kaunt=0
          AboveThr=.false.
          Do l=0,lMax
            kComp=(l+1)*(l+2)/2
            dMag=0.0d0
            Do k=1,kComp
              dM=Work(ipMP+nij*kaunt+kauntA)
              dMullig(kaunt+1)=dM
              dMag=dMag+dM**2
              kaunt=kaunt+1
            Enddo
            dMag=sqrt(dMag)
            If((dMag.gt.ThrsMul_Clever).and.l.le.1) then
              AboveThr=.true.
              AboveMul(l+1)=.true.
            Elseif((dMag.le.ThrsMul_Clever).and.l.le.1) then
              AboveMul(l+1)=.false.
            Endif
          Enddo

          If(AboveThr) then
            nAbove=nAbove+1
*
*---- Select the potential points which should be used for this centre.
*
*            BS=0.5d0*(Bragg_Slater(iANr(iAtom))
*     &               +Bragg_Slater(iANr(jAtom)))
            BS=0.5d0*(vdWRad(iANr(iAtom))+vdWRad(iANr(jAtom)))
            Call PickPoints(nPick,ipPick,ipDPick,nEPP,ipEPCo,EC(1,ij)
     &                     ,dLimmo,BS)
*
*---- Compute the true potential from the density assigned to
*     this centre.
*
            Call GetMem('Potential','Allo','Real',iPotte,nPick)
            Call EPotPoint(iPotte,nPick,ipPick,ipDPick,nEPP
     &                    ,ip_Ttot,ip_Ttot_Inv,iANr(iAtom)
     &                    ,nB,iAtom,jAtom,ip_Center)

*
*---- Print the true potential for given centre if requested.
*
            If(iPrint.ge.5) then
              Write(UtChar,'(A,2I3)')'Partial density potential, centre'
     &                               ,iAtom,jAtom
              Call RecPrt(UtChar,' ',Work(iPotte),nPick,1)
            Endif

*
*---- All hail to the chiefs, i.e. Levenberg and Marquardt.
*
            Call LevMarquart(iPotte,nPick,ipPick,nEPP,ipEPCo,EC(1,ij)
     &                      ,dMullig,lMax,A,iAtom,jAtom,chPoint
     &                      ,Thrs1,Thrs2,nThrs,Chi2B,iPrint,AboveMul)
            Call GetMem('Potential','Free','Real',iPotte,nPick)
          Endif
*
*---- Store things for later.
*
          kaunt=0
          Pot_Point(kauntA+1)=chPoint
          Do iDC=1,2
            nK=iDC*(iDC+1)/2
            Do iK=1,nK
              kaunt=kaunt+1
              Pot_Fac(4*kauntA+kaunt)=dMullig(kaunt)
            Enddo
            If(.not.AboveThr) then
              Diffed(2*kauntA+iDC)=.false.
            Else
              If(A(iDC).lt.3.0d0.and.AboveMul(iDC)) then
                Diffed(2*kauntA+iDC)=.true.
                Pot_Expo(2*kauntA+iDC)=A(iDC)
              Else
                Diffed(2*kauntA+iDC)=.false.
                Pot_Expo(2*kauntA+iDC)=1.0d1 !Dummy
              Endif
            Endif
          Enddo
*
*---- Step the atom-pair counter.
*
          kauntA=kauntA+1
        Enddo
      Enddo

*
*-- Deallocations.
*
      Call GetMem('BasIndCent','Free','Inte',ip_Center,nB)
      Call GetMem('PickPoints','Free','Inte',ipPick,nEPP)
      Call GetMem('DistPick','Free','Real',ipDPick,nEPP)
      Call GetMem('PotPointCoord','Free','Real',ipEPCo,3*nEPP)
      irc=-1
      Call ClsOne(irc,0)

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(ipC)
      End
