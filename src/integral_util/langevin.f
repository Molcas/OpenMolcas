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
      SubRoutine Langevin(h1,TwoHam,D,RepNuc,nh1,First,Dff)
      Use Basis_Info
      use Center_Info
      Use Langevin_arrays
      use External_Centers
      use Phase_Info
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
      Real*8 h1(nh1), TwoHam(nh1), D(nh1)
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Logical First, Dff, Exist
      Save nAnisopol,nPolComp
      Real*8, Allocatable:: D1ao(:)
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

      iRout = 1
      iPrint = nPrint(iRout)
*                                                                      *
************************************************************************
*                                                                      *
*---- Generate list of all atoms
*
*     Cord: list of all atoms
*     Atod: associated effective atomic radius
*

      mdc = 0
      MaxAto=0
      Do iCnttp = 1, nCnttp
         nCnt = dbsc(iCnttp)%nCntr
         Do iCnt = 1, nCnt
            mdc = mdc + 1
            MaxAto = MaxAto + nIrrep/dc(mdc)%nStab
         End Do
      End Do
*
      Call GetMem('Cord','Allo','Real',ipCord,3*MaxAto)
      Call GetMem('Chrg','Allo','Real',ipChrg,MaxAto)
      Call GetMem('Atod','Allo','Real',ipAtod,MaxAto)
*
      ndc = 0
      nc = 1
      Do jCnttp = 1, nCnttp
         Z = dbsc(jCnttp)%Charge
         mCnt = dbsc(jCnttp)%nCntr
         If (dbsc(jCnttp)%AtmNr.ge.1) Then
*            Atod = CovRad (dbsc(jCnttp)%AtmNr)
             Atod = CovRadT(dbsc(jCnttp)%AtmNr)
         Else
             Atod = Zero
         End If
         Do jCnt = 1, mCnt
            ndc = ndc + 1
            Do i = 0, nIrrep/dc(ndc)%nStab - 1
               Call OA(dc(ndc)%iCoSet(i,0),dbsc(jCnttp)%Coor(1:3,jCnt),
     &                 Work(ipCord+(nc-1)*3  :
     &                      ipCord+(nc-1)*3+2))
               Work(ipAtod+(nc-1)) = Atod
               Work(ipChrg+(nc-1)) = Z
*              Write (*,*) 'Z=',Z
               nc = nc + 1
            End Do
         End Do
      End Do
*
      If(LGridAverage) Then
         nAv=nGridAverage
         If(nGridSeed.eq.0) Then
            iSeed=0
         ElseIf(nGridSeed.eq.-1) Then
            iSeed=0
*     Use system_clock only for some systems
c            Call System_clock(iSeed,j,k)
         Else
            iSeed=nGridSeed
         EndIf
      Else
         nAv=1
      EndIf
      sumRepNuc=Zero
      sumRepNuc2=Zero

      Do iAv=1,nAv
         If(LGridAverage) Then
            cordsi(1,1)=Random_Molcas(iSeed)
            cordsi(2,1)=Random_Molcas(iSeed)
            cordsi(3,1)=Random_Molcas(iSeed)
            rotAlpha=Random_Molcas(iSeed)*180.0D0
            rotBeta=Random_Molcas(iSeed)*180.0D0
            rotGamma=Random_Molcas(iSeed)*180.0D0
            write(6,'(a,i4,a,6f10.4)')'GRID NR',iAv,': ',cordsi(1,1),
     &           cordsi(2,1),cordsi(3,1),rotAlpha,rotBeta,rotGamma
            Done_Lattice=.False.
            RepNuc = Zero
         EndIf

      If (.Not.Done_Lattice) Then
         lFirstIter=.True.
         Done_Lattice=.True.
         if(iXPolType.eq.2) Then
            nAnisopol = nXF
            nPolComp = 6
         Else
            nAnisopol = 0
            nPolComp = 1
         EndIf
*
*
*------- Compute Effective Polarizabilities on the Langevin grid,
*        flag also if points on the grid are excluded.
*        This information is computed once!
*
*        Grid : The Langevin grid
*        PolEf: Effective Polarizability on grid
*        DipEf: Effective Dipole moment on grid
*
         nGrid_Eff=0
*        Both these subroutines can increase nGrid_Eff
         If(iXPolType.gt.0) Then
            Call lattXPol(Grid,nGrid,nGrid_Eff,PolEf,DipEf,XF,
     &                    nXF,nOrd_XF,nPolComp)
         EndIf
         If(lLangevin) Then
*            Note: Gen_Grid is now a part of the lattcr subroutine
c            Call Gen_Grid(Work(ipGrid+3*nGrid_Eff),nGrid-nGrid_Eff)
            Call lattcr(Grid,nGrid,nGrid_Eff,PolEf,DipEf,
     &                  Work(ipCord),maxato,Work(ipAtod),nPolComp,
     &                  XF,nXF,nOrd_XF,XEle,iXPolType)
         EndIf
c        Write(6,*) 'nGrid,  nGrid_Eff', nGrid,  nGrid_Eff

*
      End If
*
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute Electric Field from the Quantum Chemical System on the
*     Langevin grid.
*
*     cavxyz: MM expansion of the total charge of the QM system
*     ravxyz: scratch
*     dField: EF on Langevin grid due to QM system
*
*     Get the total 1st order AO density matrix
*
      Call mma_allocate(D1ao,nh1,Label='D1ao')
      Call Get_D1ao(D1ao,nh1)
*
*     Save field from permanent multipoles for use in ener
      Call GetMem('pField','Allo','Real',ippField,nGrid_Eff*4)
      Call GetMem('tmpField','Allo','Real',iptmpField,nGrid_Eff*4)

      Call FZero(Work(ippField ),nGrid_Eff*4)




      Call eperm(D1ao,nh1,Ravxyz,Cavxyz,nCavxyz,
     &           dField,Grid,nGrid_Eff,Work(ipCord),
     &           MaxAto,Work(ipChrg),Work(ippField))

*                                                                      *
************************************************************************
*                                                                      *
*---- Save system info, to be used by visualisation program            *

      Lu=21
      Call OpnFl('LANGINFO',Lu,Exist)
      If(.not.Exist) Then
      Write(Lu,*)nc-1
      do i=0,nc-2
         Write(Lu,11)INT(Work(ipChrg+i)),Work(ipAtod+i),
     &        (Work(ipCord+i*3+j),j=0,2)
 11      format(i3,f10.4,3f16.8)
      enddo
      Write(Lu,*)nXF
      Inc = 3
      Do iOrdOp = 0, nOrd_XF
         Inc = Inc + nElem(iOrdOp)
      End Do
      If(iXPolType.gt.0) Inc = Inc + 6
      Do iXF=1,nXF
         xa=XF(1,iXF)
         ya=XF(2,iXF)
         za=XF(3,iXF)
         If(XEle(iXF).le.0) Then
            atrad=-DBLE(XEle(iXF))/1000.0D0
            iele=0
         Else
            iele=XEle(iXF)
            atrad=CovRadT(iele)
         EndIf
         Write(Lu,11)iele,atrad,xa,ya,za
      EndDo
      Write(Lu,*)nGrid_eff,nAnisopol
      do i=0,nGrid_eff-1
         Write(Lu,12)(Grid(j,i+1),j=1,3),
     &        PolEf(:,i+1),
     &        DipEf(i+1),
     &        (dField(j,i+1),j=1,3),(Work(ippField+i*4+j),j=0,2)
 12      format(11f20.10)
      enddo
      Write(Lu,*)polsi,dipsi,scala,One/tK/3.1668D-6
      Write(Lu,*)(cordsi(k,1),k=1,3)
      Write(Lu,*)rotAlpha, rotBeta, rotGamma
      Write(Lu,*)radlat,nSparse,distSparse
      Write(Lu,*)lDamping, dipCutoff
      Endif
      Close(Lu)

      call dcopy_(nGrid_Eff*4,dField,1,Work(iptmpField),1)

      If(lDiprestart .or. lFirstIter) Then
         Field(:,:)=Zero
         Dip(:,:)=Zero
         Davxyz(:)=Zero
      EndIf

*---- Subtract the static MM from the previous iteration
*     from the static MM of this iteration, and save the
*     untouched static MM of this iteration into Davxyz
*     for use in the next iteration. Ravxyz is
*     just a temporary array

      call dcopy_(nCavxyz,Cavxyz,1,Ravxyz,1)
      Call DaXpY_(nCavxyz,-One,Davxyz,1,Cavxyz,1)
      call dcopy_(nCavxyz,Ravxyz,1,Davxyz,1)

*     Ravxyz(:)=Cavxyz(:)
*     Cavxyz(:)=Cavxyz(:)-Davxyz(:)
*     Davxyz(:)=Ravxyz(:)
*                                                                      *
************************************************************************
*                                                                      *
*---- Equation solver: compute the Langevin dipole moments and the
*                      counter charge on the boundary of the cavity.
*
*     Field : total EF of the Langevin grid
*     Dip   : dipole momement on the Langevin grid
*
      Call edip(Ravxyz,Cavxyz,lmax,
     &          Field,Dip,dField,
     &          PolEf,DipEf,
     &          Grid,nGrid_Eff,nPolComp,nAnisopol,
     &          nXF,iXPolType,nXMolnr,XMolnr)

*                                                                      *
************************************************************************
*                                                                      *
*---- Compute contributions to RepNuc, h1, and TwoHam
*
      Call Ener(h1,TwoHam,D,RepNuc,nh1,First,Dff,D1ao,
     &          Grid,nGrid_Eff,Dip, Field,
     &          DipEf,PolEf,Work(ipCord),MaxAto,
     &          Work(ipChrg),nPolComp,nAnisopol,Work(ippField),
     &     Work(iptmpField))


*     Subtract the static field from the self-consistent field
*     This gives the field from the induced dipoles (saved
*     in Field, to be used in the next iteration if
*     not DRES has been requested

      Call DaXpY_(nGrid*4,-One,Work(iptmpField),1,Field,1)
      lFirstIter=.False.


      Call GetMem('pField','Free','Real',ippField,nGrid_Eff*4)
      Call GetMem('tmpField','Free','Real',iptmpField,nGrid_Eff*4)
      If(LGridAverage) Then
         Write(6,'(a,i4,a,f18.10)')'Solvation energy (Grid nr. ',iAv,
     &        '):',RepNuc
         sumRepNuc = sumRepNuc + RepNuc
         sumRepNuc2 = sumRepNuc2 + RepNuc**2
      EndIf
      EndDo
      If(LGridAverage) Then
         Write(6,'(a,f18.10,f18.10)')
     &        'Average solvation energy and stdev: ',
     &        sumRepNuc/DBLE(nAv),
     &        sqrt(sumRepNuc2/DBLE(nAv)-(sumRepNuc/DBLE(nAv))**2)
      EndIf
      Call GetMem('Atod','Free','Real',ipAtod,MaxAto)
      Call mma_deallocate(D1ao)
      Call GetMem('Chrg','Free','Real',ipChrg,MaxAto)
      Call GetMem('Cord','Free','Real',ipCord,3*MaxAto)

*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
