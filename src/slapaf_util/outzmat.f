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
* Copyright (C) 2007, Giovanni Ghigo                                   *
************************************************************************
      Subroutine OutZMAT(nAtoms,XYZCoords,N_ZMAT)
      Implicit Real*8 (a-h,o-z)
      Real*8 XYZCoords(3,nAtoms)
#include "WrkSpc.fh"
#include "SysDef.fh"
*
      nSymbols=(5+2)*N_ZMAT
      Call GetMem('Symbols','Allo','Char',ip_Symbols,nSymbols)
      Call GetMem('NAT','Allo','Inte',ip_NAT,N_ZMAT)
      Call GetMem('ZMATCoors','Allo','Real',ip_ZMATCoords,3*N_ZMAT)
      Call GetMem('ZMAT','Allo','Real',ip_ZMAT,N_ZMAT*3)
*
      Call OutZMAT_(nAtoms,XYZCoords,N_ZMAT,cWork(ip_Symbols),
     &             iWork(ip_NAT),Work(ip_ZMATCoords),Work(ip_ZMAT))
*
      Call GetMem('ZMAT','Free','Real',ip_ZMAT,N_ZMAT*3)
      Call GetMem('ZMATCoors','Free','Real',ip_ZMATCoords,3*N_ZMAT)
      Call GetMem('NAT','Free','Inte',ip_NAT,N_ZMAT)
      Call GetMem('Symbols','Free','Char',ip_Symbols,nSymbols)
*
      Return
      End
      Subroutine OutZMAT_(nAtoms,XYZCoords,N_ZMAT,Symbols,NAT,
     &                    ZMATCoords,ZMAT)
************************************************************************
* Author: Giovanni Ghigo                                               *
*         Torino (Italy)  February-March 2007                          *
*                                                                      *
* This subroutine generates the nuclear coordinates in Z-matrix        *
* format using the connection indices recovered from RunFile.          *
* Remember:                                                            *
*           X dummy atoms - NAT(i)= 0 - are included in SEWARD         *
*           Z ghost atoms - NAT(i)=-1 - are NOT included in SEWARD     *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "angstr.fh"
      Parameter (MaxAtoms=256) ! See src/gateway_util/g_zmatconv.fh
      Character*5 Symbols(N_ZMAT)   ! Obs: Restricted record
      Integer NAT(N_ZMAT)           ! Obs: Restricted record
      Integer iZmat(MaxAtoms,3)     ! Obs: Full record
      Real*8 XYZCoords(3,nAtoms), ZMATCoords(3,N_ZMAT), ZMAT(N_ZMAT,3)
      Real*8 CTX(3,4), Bt(3,4), dBt(3,4,3,4)
      Logical BigTrasl, IfTest
      Character*8 Label
      Parameter (ThrsTrasl=1.0d0) ! Threshold for warning
*
      IfTest=.False.
#ifdef _DEBUGPRINT_
      Call QEnter('OutZMAT')
      IfTest=.True.
#endif
*
      LuWr = 6
      todeg = 45.0d0 / ATan(1.0d0)
      dMaxTrasl = 0.0d0
      BigTrasl = .False.
      Label =' '
      Call dZero(ZMAT,N_ZMAT*3)
      Call dZero(ZMATCoords,N_ZMAT*3)
      Call Get_cArray('Symbol ZMAT',Symbols,N_ZMAT*5)
      Call Get_iArray('Index ZMAT',iZmat,MaxAtoms*3)
      Call Get_iArray('NAT ZMAT',NAT,N_ZMAT)

      If (IfTest) then
        Write(LuWr,*)
        Write(LuWr,*) '------------------------------------------------'
        Write(LuWr,*) 'OutZMat - Z-Matrix Data :'
        Write(LuWr,*) '          N_ZMAT=',N_ZMAT
        Write(LuWr,*) '   Label  NA      i   j   k'
        Do i=1,N_ZMAT
          Write(LuWr,97) i,Symbols(i),NAT(i),(iZmat(i,j),j=1,3)
        EndDo
        Write(LuWr,*)
        Write(LuWr,*) '------------------------------------------------'
        Write(LuWr,*) 'OutZMat - XYZCoords (Angstrom):'
        Do i=1,nAtoms
          Write(LuWr,98) i,(angstr*XYZCoords(j,i),j=1,3)
        EndDo
        Write(LuWr,*)
      EndIf

      If (N_ZMAT.LT.3) Return
      If (N_ZMAT.GT.(nAtoms+3)) then
        Write(LuWr,'(A)') ' ZMAT cannot be defined :'
        Write(LuWr,'(A)') ' Too many ghost Z atoms. '
        Return
      EndIf

*  Z are ghost atoms
*  A are both dummy (X) and real atoms

*  Z  Z  Z
      If (NAT(1).EQ.-1.and.NAT(2).EQ.-1.and.NAT(3).EQ.-1) then
        iBonded = 0
        Do i = 4, N_ZMAT
          If (iZmat(i,1).EQ.2.and.iBonded.EQ.0) iBonded = i
        EndDo
        If (iBonded.EQ.0) then
          ZMATCoords(3,2) = 1.0d0/angstr
        else
          ZMATCoords(3,2) = XYZCoords(3,iBonded-3)
        EndIf
        iBonded = 0
        Do i = 4, N_ZMAT
          If (iZmat(i,1).EQ.3.and.iBonded.EQ.0) iBonded = i
        EndDo
        If (iBonded.EQ.0) then
          ZMATCoords(1,3) = 1.0d0/angstr
        else
          ZMATCoords(1,3) = XYZCoords(1,iBonded-3)
        EndIf
        If (iZmat(3,1).EQ.2) ZMATCoords(3,3) = ZMATCoords(3,2)
        Do iAtom = 4, N_ZMAT
          ZMATCoords(1,iAtom)=XYZCoords(1,iAtom-3)
          ZMATCoords(2,iAtom)=XYZCoords(2,iAtom-3)
          ZMATCoords(3,iAtom)=XYZCoords(3,iAtom-3)
        EndDo
      EndIf

*  Z  Z  A
      If (NAT(1).EQ.-1.and.NAT(2).EQ.-1.and.NAT(3).GE.0) then
        iBonded = 0
        Do i = 3, N_ZMAT
          If (iZmat(i,1).EQ.2.and.iBonded.EQ.0) iBonded = i
        EndDo
        If (iBonded.EQ.0) then
          ZMATCoords(3,2) = 1.0d0/angstr
        else
          ZMATCoords(3,2) = XYZCoords(3,iBonded-2)
        EndIf
        Do iAtom = 3, N_ZMAT
          ZMATCoords(1,iAtom)=XYZCoords(1,iAtom-2)
          ZMATCoords(2,iAtom)=XYZCoords(2,iAtom-2)
          ZMATCoords(3,iAtom)=XYZCoords(3,iAtom-2)
        EndDo
      EndIf

*  Z  A  Z
      If (NAT(1).EQ.-1.and.NAT(3).EQ.-1.and.NAT(2).GE.0) then
        dMaxTrasl = 0.0d0
        dX = XYZCoords(1,1)
        dY = XYZCoords(2,1)
        Do iAtom = 1, N_ZMAT-2
          XYZCoords(1,iAtom)=XYZCoords(1,iAtom)-dX
          XYZCoords(2,iAtom)=XYZCoords(2,iAtom)-dY
        EndDo
        ZMATCoords(1,2)=XYZCoords(1,1)
        ZMATCoords(2,2)=XYZCoords(2,1)
        ZMATCoords(3,2)=XYZCoords(3,1)
        dMaxTrasl = Max(0.0d0,ABS(dX),ABS(dY))
        If (dMaxTrasl.GE.ThrsTrasl) Write (LuWr,'(A)')
     &  '  Warning: MaxTrasl.GE.ThrsTrasl for atom ',Symbols(2)
        If (iZmat(3,1).EQ.2) ZMATCoords(3,3) = ZMATCoords(3,2)
        iBonded = 0
        Do i = 4, N_ZMAT
          If (iZmat(i,1).EQ.3.and.iBonded.EQ.0) iBonded = i
        EndDo
        If (iBonded.EQ.0) then
          ZMATCoords(1,3) = 1.0d0/angstr
        else
          ZMATCoords(1,3) = XYZCoords(1,iBonded-2)
        EndIf
        Do iAtom = 4, N_ZMAT
          ZMATCoords(1,iAtom)=XYZCoords(1,iAtom-2)
          ZMATCoords(2,iAtom)=XYZCoords(2,iAtom-2)
          ZMATCoords(3,iAtom)=XYZCoords(3,iAtom-2)
        EndDo
      EndIf

*  A  Z  Z
      If (NAT(1).GE.0.and.NAT(2).EQ.-1.and.NAT(3).EQ.-1) then
        dX = XYZCoords(1,1)
        dY = XYZCoords(2,1)
        dZ = XYZCoords(3,1)
        Do iAtom = 1, N_ZMAT-2
          XYZCoords(1,iAtom)=XYZCoords(1,iAtom)-dX
          XYZCoords(2,iAtom)=XYZCoords(2,iAtom)-dY
          XYZCoords(3,iAtom)=XYZCoords(3,iAtom)-dZ
        EndDo
        dMaxTrasl = Max(0.0d0,ABS(dX),ABS(dY),ABS(dZ))
        If (dMaxTrasl.GE.ThrsTrasl) Write (LuWr,'(A)')
     &  '  Warning: MaxTrasl.GE.ThrsTrasl for atom ',Symbols(1)
        iBonded = 0
        Do i = 4, N_ZMAT
          If (iZmat(i,1).EQ.2.and.iBonded.EQ.0) iBonded = i
        EndDo
        If (iBonded.EQ.0) then
          ZMATCoords(3,2) = 1.0d0/angstr
        else
          ZMATCoords(3,2) = XYZCoords(3,iBonded-2)
        EndIf
        iBonded = 0
        Do i = 4, N_ZMAT
          If (iZmat(i,1).EQ.3.and.iBonded.EQ.0) iBonded = i
        EndDo
        If (iBonded.EQ.0) then
          ZMATCoords(1,3) = 1.0d0/angstr
        else
          ZMATCoords(1,3) = XYZCoords(1,iBonded-2)
        EndIf
        Do iAtom = 4, N_ZMAT
          ZMATCoords(1,iAtom)=XYZCoords(1,iAtom-2)
          ZMATCoords(2,iAtom)=XYZCoords(2,iAtom-2)
          ZMATCoords(3,iAtom)=XYZCoords(3,iAtom-2)
        EndDo
      EndIf

*  Z  A  B
      If (NAT(1).EQ.-1.and.NAT(2).GE.0.and.NAT(3).GE.0) then
        dX = XYZCoords(1,1)
        dY = XYZCoords(2,1)
        Do iAtom = 2, N_ZMAT
          ZMATCoords(1,iAtom)=XYZCoords(1,iAtom-1)-dX
          ZMATCoords(2,iAtom)=XYZCoords(2,iAtom-1)-dY
          ZMATCoords(3,iAtom)=XYZCoords(3,iAtom-1)
        EndDo
        dMaxTrasl = Max(0.0d0,ABS(dX),ABS(dY))
        If (dMaxTrasl.GE.ThrsTrasl) Write (LuWr,'(A)')
     &  '  Warning: MaxTrasl.GE.ThrsTrasl for atom ',Symbols(2)
      EndIf

*  A  Z  B
      If (NAT(1).GE.0.and.NAT(2).EQ.-1.and.NAT(3).GE.0) then
        dX = XYZCoords(1,1)
        dY = XYZCoords(2,1)
        dZ = XYZCoords(3,1)
        Do iAtom = 1, N_ZMAT-1
          XYZCoords(1,iAtom)=XYZCoords(1,iAtom)-dX
          XYZCoords(2,iAtom)=XYZCoords(2,iAtom)-dY
          XYZCoords(3,iAtom)=XYZCoords(3,iAtom)-dZ
        EndDo
        dMaxTrasl = Max(0.0d0,ABS(dX),ABS(dY),ABS(dZ))
        If (dMaxTrasl.GE.ThrsTrasl) Write (LuWr,'(A)')
     &  '  Warning: MaxTrasl.GE.ThrsTrasl for atom ',Symbols(1)
        iBonded = 0
        Do i = 3, N_ZMAT
          If (iZmat(i,1).EQ.2.and.iBonded.EQ.0) iBonded = i
        EndDo
        If (iBonded.EQ.0) then
          ZMATCoords(3,2) = 1.0d0/angstr
        else
          ZMATCoords(3,2) = XYZCoords(3,iBonded-1)
        EndIf
        Do iAtom = 3, N_ZMAT
          ZMATCoords(1,iAtom)=XYZCoords(1,iAtom-1)
          ZMATCoords(2,iAtom)=XYZCoords(2,iAtom-1)
          ZMATCoords(3,iAtom)=XYZCoords(3,iAtom-1)
        EndDo
      EndIf

*  A  B  Z
      If (NAT(1).GE.0.and.NAT(2).GE.0.and.NAT(3).EQ.-1) then
        dMaxTrasl = 0.0d0
        dX = XYZCoords(1,1)
        dY = XYZCoords(2,1)
        dZ = XYZCoords(3,1)
        Do iAtom = 1, N_ZMAT-1
          XYZCoords(1,iAtom)=XYZCoords(1,iAtom)-dX
          XYZCoords(2,iAtom)=XYZCoords(2,iAtom)-dY
          XYZCoords(3,iAtom)=XYZCoords(3,iAtom)-dZ
        EndDo
        dMaxTrasl = Max(0.0d0,ABS(dX),ABS(dY),ABS(dZ))
        If (dMaxTrasl.GE.ThrsTrasl) Write (LuWr,'(A)')
     &  '  Warning: MaxTrasl.GE.ThrsTrasl for atom ',Symbols(1)
        ZMATCoords(1,2)=XYZCoords(1,2)
        ZMATCoords(2,2)=XYZCoords(2,2)
        ZMATCoords(3,2)=XYZCoords(3,2)
        If (iZmat(3,1).EQ.1) then
          iBonded = 0
          Do i = 4, N_ZMAT
            If (iZmat(i,1).EQ.3.and.iBonded.EQ.0) iBonded = i
          EndDo
          If (iBonded.EQ.0) then
            ZMATCoords(1,3) = 1.0d0/angstr
          else
            ZMATCoords(1,3) = XYZCoords(1,iBonded-2)
          EndIf
        else
          ZMATCoords(1,3) = XYZCoords(1,2) + 1.0d0/angstr
          ZMATCoords(2,3) = XYZCoords(2,2)
          ZMATCoords(3,3) = XYZCoords(3,2)
        EndIf
        Do iAtom = 4, N_ZMAT
          ZMATCoords(1,iAtom)=XYZCoords(1,iAtom-1)
          ZMATCoords(2,iAtom)=XYZCoords(2,iAtom-1)
          ZMATCoords(3,iAtom)=XYZCoords(3,iAtom-1)
        EndDo
      EndIf

*  A  B  C
      If (NAT(1).GE.0.and.NAT(2).GE.0.and.NAT(3).GE.0) then
        Call dCopy_(N_ZMAT*3,XYZCoords,1,ZMATCoords,1)
      EndIf


      If (IfTest) then
        Write(LuWr,*)
        Write(LuWr,*) '------------------------------------------------'
        Write(LuWr,*) 'OutZMat - ZMATCoords : '
        Do i=1,N_ZMAT
          Write(LuWr,99) i,NAT(i),(angstr*ZMATCoords(j,i),j=1,3)
        EndDo
      EndIf

      iAtom = 2
      iRX = iZMAT(iAtom,1)
      dX = ZMATCoords(1,iAtom) - ZMATCoords(1,iRX)
      dY = ZMATCoords(2,iAtom) - ZMATCoords(2,iRX)
      dZ = ZMATCoords(3,iAtom) - ZMATCoords(3,iRX)
      bond = SQRT( dX*dX + dY*dY + dZ*dZ )
      ZMAT(iAtom,1) = bond*angstr

      If (N_ZMAT.GE.3) then
        iAtom = 3
        iRX = iZMAT(iAtom,1)
        iAX = iZMAT(iAtom,2)
        dX = ZMATCoords(1,iAtom) - ZMATCoords(1,iRX)
        dY = ZMATCoords(2,iAtom) - ZMATCoords(2,iRX)
        dZ = ZMATCoords(3,iAtom) - ZMATCoords(3,iRX)
        dX2 = ZMATCoords(1,iAX) - ZMATCoords(1,iRX)
        dY2 = ZMATCoords(2,iAX) - ZMATCoords(2,iRX)
        dZ2 = ZMATCoords(3,iAX) - ZMATCoords(3,iRX)
        bond = SQRT( dX*dX + dY*dY + dZ*dZ )
        ZMAT(iAtom,1) = bond*angstr
        prod = dX*dX2 + dY*dY2 + dZ*dZ2
        dvec2 = SQRT( dX2*dX2 + dY2*dY2 + dZ2*dZ2 )
        arccos = prod / ( bond * dvec2 )
        If (arccos.GT.1.0d0) arccos = Sign(1.0d0,arccos)
        alpha = ACOS(arccos)
        ZMAT(iAtom,2) = alpha*todeg
      EndIf

      If (N_ZMAT.GE.4) then
        Do iAtom = 4, N_ZMAT
          If (NAT(iAtom).EQ.-1) then
            Write(LuWr,'(A)') ' ZMAT cannot be defined :'
            Write(LuWr,'(A)') ' Found X atom.           '
            Return
          EndIf
          iRX = iZMAT(iAtom,1)
          iAX = iZMAT(iAtom,2)
          iTX = iZMAT(iAtom,3)
          dX = ZMATCoords(1,iAtom) - ZMATCoords(1,iRX)
          dY = ZMATCoords(2,iAtom) - ZMATCoords(2,iRX)
          dZ = ZMATCoords(3,iAtom) - ZMATCoords(3,iRX)
          dX2 = ZMATCoords(1,iAX) - ZMATCoords(1,iRX)
          dY2 = ZMATCoords(2,iAX) - ZMATCoords(2,iRX)
          dZ2 = ZMATCoords(3,iAX) - ZMATCoords(3,iRX)
          bond = SQRT( dX*dX + dY*dY + dZ*dZ )
          ZMAT(iAtom,1) = bond*angstr
          prod = dX*dX2 + dY*dY2 + dZ*dZ2
          dvec2 = SQRT( dX2*dX2 + dY2*dY2 + dZ2*dZ2 )
          arccos = prod / ( bond * dvec2 )
          If (arccos.GT.1.0d0) arccos = Sign(1.0d0,arccos)
          alpha = ACOS(arccos)
          ZMAT(iAtom,2) = alpha*todeg
          Do ii=1,3
            CTX(ii,1) = ZMATCoords(ii,iTX)
            CTX(ii,2) = ZMATCoords(ii,iAX)
            CTX(ii,3) = ZMATCoords(ii,iRX)
            CTX(ii,4) = ZMATCoords(ii,iAtom)
          EndDo
          Call Trsn(CTX,4,theta,Bt,.False.,.False.,Label,dBt,.False.)
          ZMAT(iAtom,3) = -theta*todeg
        EndDo
      EndIf

      Write (LuWr,*)
      Write (LuWr,*) '***********************************************'//
     &               '*****************'
      Write (LuWr,*) '* Nuclear coordinates in Z-Matrix format / Angs'//
     &               'trom and Degree *'
      Write (LuWr,*) '-----------------------------------------------'//
     &               '-----------------'
      Do i=1, N_ZMAT
        If (i.EQ.1)
     &    Write(LuWr,201) Symbols(i)
        If (i.EQ.2)
     &    Write(LuWr,202) Symbols(i),iZmat(i,1),ZMAT(i,1)
        If (i.EQ.3)
     &    Write(LuWr,203) Symbols(i),(iZmat(i,j),ZMAT(i,j),j=1,2)
        If (i.GE.4)
     &    Write(LuWr,204) Symbols(i),(iZmat(i,j),ZMAT(i,j),j=1,3)
      EndDo
      Write (LuWr,*)
201   Format (5X,A5)
202   Format (5X,A5,1X,I4,F13.6)
203   Format (5X,A5,1X,2(I4,F13.6))
204   Format (5X,A5,1X,3(I4,F13.6))

97    Format(I3,1X,A,1X,I3,3X,3I4)
98    Format(I3,1X,3(F12.6))
99    Format(I3,1X,I3,1X,3(F12.6))

#ifdef _DEBUGPRINT_
      Call qExit('OutZMAT')
#endif
      Return
      End
