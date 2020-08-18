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
* Copyright (C) 1991,1995, Roland Lindh                                *
************************************************************************
      SubRoutine DrvN0()
************************************************************************
*                                                                      *
* Object: to compute the nuclear contibutions to the nuclear potential *
*                                                                      *
* Called from: Input                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October '91                                              *
*             Modified for various other contributions May 95', RL     *
************************************************************************
      use external_centers
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
      Real*8 A(3), B(3), RB(3)
      Integer iDCRR(0:7), jCoSet(8,8), iStb(0:7), jStb(0:7)
      Logical EQ, NoLoop
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
C     nElem(ixyz) = 2*ixyz+1
*

      iRout = 33
      iPrint = nPrint(iRout)
      Call qEnter('DrvN0')
      NoLoop=.True.
      iDum=0
      r12_Min=0.0D0
*
*----- Nuclear repulsion, in case of some ECP we include the core electronic
*     contribution (pseudo charges). The interaction of pseudo charges is
*     excluded from the energy term.
*
      PotNuc = Zero
      mdc = 0
      ZB=Zero
      Do iCnttp = 1, nCnttp
         ZA = Charge(iCnttp)
         If (dbsc(iCnttp)%Frag) ZA = FragCharge(iCnttp)
         If (ZA.eq.Zero) Go To 101
         Do iCnt = 1, dbsc(iCnttp)%nCntr
            A(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
*
            ndc = 0
            Do jCnttp = 1, iCnttp
               If (pChrg(iCnttp).and.pChrg(jCnttp)) Go To 201
               If (dbsc(iCnttp)%Frag.and.dbsc(jCnttp)%Frag) Go To 201
               ZB = Charge(jCnttp)
               If (dbsc(jCnttp)%Frag) ZB = FragCharge(jCnttp)
               If (ZB.eq.Zero) Go To 201
               ZAZB = ZA * ZB
               jCntMx = dbsc(jCnttp)%nCntr
               If (iCnttp.eq.jCnttp) jCntMx = iCnt
               Do jCnt = 1, jCntMx
*                 Introduce factor to ensure that contributions from
*                 A>B are the only to be accumulated.
                  Fact = One
                  If (iCnttp.eq.jCnttp.and.iCnt.eq.jCnt) Fact = Half
                  B(1:3)=dbsc(jCnttp)%Coor(1:3,jCnt)
*
*                 Find the DCR for the two centers
*
                  Call DCR(LmbdR,iOper,nIrrep,
     &                     jStab(0,mdc+iCnt),nStab(mdc+iCnt),
     &                     jStab(0,ndc+jCnt),nStab(ndc+jCnt),
     &                     iDCRR,nDCRR)
*
                  temp = Zero
                  Do iR = 0, nDCRR-1
                     RB(1) = DBLE(iPhase(1,iDCRR(iR)))*B(1)
                     RB(2) = DBLE(iPhase(2,iDCRR(iR)))*B(2)
                     RB(3) = DBLE(iPhase(3,iDCRR(iR)))*B(3)
*                    The index A=RB is illegal.
                     If (.Not.EQ(A,RB)) Then
                        r12 = Sqrt((A(1)-RB(1))**2 +
     &                             (A(2)-RB(2))**2 +
     &                             (A(3)-RB(3))**2 )
                        If (r12.lt.r12_min .and.
     &                      .NOT.Expert ) Then
                           Call WarningMessage(2,
     &                        ' The distance between two centers'
     &                      //' are found to be unphysically too'
     &                      //' short.;'
     &                      //' If you want to persist operate in'
     &                      //' expert mode.')
                           Call Abend()
                        End If
                        fab = One
                        If (dbsc(iCnttp)%ECP) Then
*--------------------------Add contribution from M1 operator
                           Do iM1xp=1, dbsc(iCnttp)%nM1
                             Gamma = dbsc(iCnttp)%M1xp(iM1xp)
                             CffM1 = dbsc(iCnttp)%M1cf(iM1xp)
                             fab = fab + CffM1 * Exp(-Gamma*r12**2)
                           End Do
*--------------------------Add contribution from M2 operator
                           Do iM2xp=1, dbsc(iCnttp)%nM2
                             Gamma = dbsc(iCnttp)%M2xp(iM2xp)
                             CffM2 = dbsc(iCnttp)%M2cf(iM2xp)
                             fab = fab + CffM2*r12*Exp(-Gamma*r12**2)
                           End Do
                        End If
                        If (dbsc(jCnttp)%ECP) Then
*--------------------------Add contribution from M1 operator
                           Do iM1xp=1, dbsc(jCnttp)%nM1
                             Gamma = dbsc(jCnttp)%M1xp(iM1xp)
                             CffM1 = dbsc(jCnttp)%M1cf(iM1xp)
                             fab = fab + CffM1 * Exp(-Gamma*r12**2)
                           End Do
*--------------------------Add contribution from M2 operator
                           Do iM2xp=1, dbsc(jCnttp)%nM2
                             Gamma = dbsc(jCnttp)%M2xp(iM2xp)
                             CffM2 = dbsc(jCnttp)%M2cf(iM2xp)
                             fab = fab + CffM2*r12*Exp(-Gamma*r12**2)
                           End Do
                        End If
*
                        temp = temp + fab/r12
                     Else
*---------------------------turn off checking if in expert mode, for
*wavefunction overlaps
                        if (.NOT.Expert) then
                        If (iCnttp.ne.jCnttp .or.
     &                      iCnt.ne.jCnt          ) Then
                           Call WarningMessage(2,
     &                                 'You have two charges on top '//
     &                                 'of each other!;'//
     &                                 'Correct the input!')
                           Call Quit_OnUserError()
                        End If
                        End If
                     End If
                  End Do
                  PotNuc = PotNuc +
     &                 (Fact*ZAZB*temp*DBLE(nIrrep))/DBLE(LmbdR)
*
                  jxyz = jxyz + 3
               End Do
 201           Continue
               ndc = ndc + dbsc(jCnttp)%nCntr
            End Do
         End Do
 101     Continue
         mdc = mdc + dbsc(iCnttp)%nCntr
      End Do
*
      If (Show) Then
         Write (6,*)
         Write (6,'(11X,A,F16.8,A)')
     &     ' Nuclear Potential Energy        ',PotNuc,' au'
         Write (6,*)
      End If
*
      If (lXF.and.(nOrd_XF.ge.0)) Then
*
         If (nIrrep.eq.8) Then
            nOper=3
         Else If (nIrrep.eq.4) Then
            nOper=2
         Else If (nIrrep.eq.2) Then
            nOper=1
         Else
            nOper=0
         End If
*
*--------Add contibution for interaction external field and nuclear
*        charges. Here we will have charge-charge, and charge-dipole
*        inteaction.
*
         ZA = Zero
         DAx= Zero
         DAy= Zero
         DAz= Zero
         Qxx= Zero
         Qxy= Zero
         Qxz= Zero
         Qyy= Zero
         Qyz= Zero
         Qzz= Zero
*
         PNX=Zero
         iDum=0
         Do iFd = 1, nXF
            If (nOrd_XF.eq.0) Then
               ZA = XF(4,iFd)
               NoLoop = ZA.eq.Zero
            Else If (nOrd_XF.eq.1) Then
               ZA = XF(4,iFd)
               DAx= XF(5,iFd)
               DAy= XF(6,iFd)
               DAz= XF(7,iFd)
               NoLoop = ZA.eq.Zero  .and.
     &                  DAx.eq.Zero .and.
     &                  DAy.eq.Zero .and.
     &                  DAz.eq.Zero
            Else If (nOrd_XF.eq.2) Then
               ZA = XF(4,iFd)
               DAx= XF(5,iFd)
               DAy= XF(6,iFd)
               DAz= XF(7,iFd)
               Qxx= XF(8,iFd)
               Qxy= XF(9,iFd)
               Qxz= XF(10,iFd)
               Qyy= XF(11,iFd)
               Qyz= XF(12,iFd)
               Qzz= XF(13,iFd)
               NoLoop = ZA.eq.Zero  .and.
     &                  DAx.eq.Zero .and.
     &                  DAy.eq.Zero .and.
     &                  DAz.eq.Zero .and.
     &                  Qxx.eq.Zero .and.
     &                  Qxy.eq.Zero .and.
     &                  Qxz.eq.Zero .and.
     &                  Qyy.eq.Zero .and.
     &                  Qyz.eq.Zero .and.
     &                  Qzz.eq.Zero
            Else
               Call WarningMessage(2,'Option not implemented yet!')
               Call Quit_OnUserError()
            End If
            If (NoLoop) Go To 102
            A(1:3) = XF(1:3,iFd)
            iChxyz=iChAtm(A,iOper,nOper,iChBas(2))
            Call Stblz(iChxyz,iOper,nIrrep,nStb,iStb,iDum,jCoSet)
*
            ndc = 0
            Do jCnttp = 1, nCnttp
               ZB = Charge(jCnttp)
               If (pChrg(jCnttp)) Go To 202
               If (ZB.eq.Zero) Go To 202
               If (dbsc(jCnttp)%Frag) Go To 202
               ZAZB = ZA * ZB
               Do jCnt = 1, dbsc(jCnttp)%nCntr
                  B(1:3)=dbsc(jCnttp)%Coor(1:3,jCnt)
*
*                 Find the DCR for the two centers
*
                  Call DCR(LmbdR,iOper,nIrrep,
     &                     iStb,nStb,
     &                     jStab(0,ndc+jCnt),nStab(ndc+jCnt),
     &                     iDCRR,nDCRR)
*
                  temp0= Zero
                  temp1= Zero
                  temp2= Zero
                  Do iR = 0, nDCRR-1
                     RB(1) = DBLE(iPhase(1,iDCRR(iR)))*B(1)
                     RB(2) = DBLE(iPhase(2,iDCRR(iR)))*B(2)
                     RB(3) = DBLE(iPhase(3,iDCRR(iR)))*B(3)
*                    The index A=RB is illegal.
                     If (.Not.EQ(A,RB)) Then
                        ABx=A(1)-RB(1)
                        ABy=A(2)-RB(2)
                        ABz=A(3)-RB(3)
                        r12 = Sqrt(ABx**2 + ABy**2 + ABz**2)
*
                        fab=One
                        If (dbsc(jCnttp)%ECP) Then
*--------------------------Add contribution from M1 operator
                           Do iM1xp=1, dbsc(jCnttp)%nM1
                             Gamma = dbsc(jCnttp)%M1xp(iM1xp)
                             CffM1 = dbsc(jCnttp)%M1cf(iM1xp)
                             fab = fab + CffM1 * Exp(-Gamma*r12**2)
                           End Do
*--------------------------Add contribution from M2 operator
                           Do iM2xp=1, dbsc(jCnttp)%nM2
                             Gamma = dbsc(jCnttp)%M2xp(iM2xp)
                             CffM2 = dbsc(jCnttp)%M2cf(iM2xp)
                             fab = fab + CffM2*r12*Exp(-Gamma*r12**2)
                           End Do
                        End If
                        temp0=temp0+fab/r12
                        If (nOrd_XF.ge.1)
     &                     temp1=temp1-fab*(DAx*ABx
     &                                     +DAy*ABy
     &                                     +DAz*ABz)/r12**3
                        If (nOrd_XF.ge.2) then

                            temp2=temp2+fab*0.5D0*(3.0D0*(Qxx*ABx**2
     &                                     +2.0D0*Qxy*ABx*ABy
     &                                     +2.0D0*Qxz*ABx*ABz
     &                                     +Qyy*ABy**2
     &                                     +2.0D0*Qyz*ABy*ABz
     &                                     +Qzz*ABz**2 )/r12**5
     &                                   -One/r12**3*(Qxx+Qyy+Qzz))
                        EndIf

                    End If
                  End Do
                  PNX = PNX + ( ( ZAZB*temp0 + ZB*(temp1+temp2))
     &                * DBLE(nIrrep) ) / DBLE(LmbdR)
*
               End Do
 202           Continue
               ndc = ndc + dbsc(jCnttp)%nCntr
            End Do
 102        Continue
         End Do
*
         If (Show) Then
            Write (6,'(19X,A,F16.8,A)')
     &                 ' Nuclear-External Field Potential Energy ',
     &                    PNX,' au'
         End If
*
         PotNuc=PotNuc+PNX
*
*--------Add contibution for self interaction of the external field.
*        Here we will have charge-charge, charge-dipole, and dipole-
*        dipole inteaction. This term will only be computed but not
*        added to the acutual nuclear potential energy.

* LU:    in the case when the number of point charges is larger than
*        10^4, the algorith below is extremely slow, as it implies a
*        double loop over all field points (nXF). The resulting value
*        is stored in the RunFile (label='PC Self Energy'), but never
*        used again in the entire MOLCAS code.
*        As a first attempt to optimize this part of the code, we skip
*        the computation of the self interaction of the external field.

         if (nXF.gt.9999) go to 1254

         PXX=Zero

         DAx=Zero
         DAy=Zero
         DAz=Zero
         QAxx=Zero
         QAxy=Zero
         QAxz=Zero
         QAyy=Zero
         QAyz=Zero
         QAzz=Zero
         DBx=Zero
         DBy=Zero
         DBz=Zero
         QBxx=Zero
         QBxy=Zero
         QBxz=Zero
         QByy=Zero
         QByz=Zero
         QBzz=Zero

         iDum=0
         Do iFd = 1, nXF
            If (nOrd_XF.eq.0) Then
               ZA = XF(4,iFd)
               NoLoop = ZA.eq.Zero
            ElseIf (nOrd_XF.eq.1) Then
               ZA = XF(4,iFd)
               DAx= XF(5,iFd)
               DAy= XF(6,iFd)
               DAz= XF(7,iFd)
               NoLoop = ZA.eq.Zero .and. DAx.eq.Zero .and. DAy.eq.Zero
     &              .and. DAz.eq.Zero
            ElseIf (nOrd_XF.eq.2) Then
               ZA = XF(4,iFd)
               DAx= XF(5,iFd)
               DAy= XF(6,iFd)
               DAz= XF(7,iFd)
               QAxx=XF(8,iFd)
               QAxy=XF(9,iFd)
               QAxz=XF(10,iFd)
               QAyy=XF(11,iFd)
               QAyz=XF(12,iFd)
               QAzz=XF(13,iFd)
               NoLoop = ZA.eq.Zero .and. DAx.eq.Zero .and. DAy.eq.Zero
     &              .and. DAz.eq.Zero .and.
     &              QAxx.eq.Zero .and.
     &              QAxy.eq.Zero .and.
     &              QAxz.eq.Zero .and.
     &              QAyy.eq.Zero .and.
     &              QAyz.eq.Zero .and.
     &              QAzz.eq.Zero
            Else
               Call WarningMessage(2,'Option not implemented yet!')
               Call Quit_OnUserError()
            End If


            If (NoLoop) Go To 103
            A(1:3) = XF(1:3,iFd)
            iChxyz=iChAtm(A,iOper,nOper,iChBas(2))
            Call Stblz(iChxyz,iOper,nIrrep,nStb,iStb,iDum,jCoSet)
*
            Do jFd = 1, iFd
               If (nOrd_XF.eq.0) Then
                  ZB = XF(4,jFd)
                  NoLoop = ZB.eq.Zero
               ElseIf (nOrd_XF.eq.1) Then
                  ZB = XF(4,jFd)
                  DBx= XF(5,jFd)
                  DBy= XF(6,jFd)
                  DBz= XF(7,jFd)
                  NoLoop=ZB.eq.Zero.and.DBx.eq.Zero .and. DBy.eq.Zero
     &                 .and. DBz.eq.Zero
               ElseIf (nOrd_XF.eq.2) Then
                  ZB = XF(4,jFd)
                  DBx= XF(5,jFd)
                  DBy= XF(6,jFd)
                  DBz= XF(7,jFd)
                  QBxx=XF(8,jFd)
                  QBxy=XF(9,jFd)
                  QBxz=XF(10,jFd)
                  QByy=XF(11,jFd)
                  QByz=XF(12,jFd)
                  QBzz=XF(13,jFd)
                  NoLoop=ZB.eq.Zero.and.DBx.eq.Zero .and. DBy.eq.Zero
     &                 .and. DBz.eq.Zero .and.
     &                 QBxx.eq.Zero .and.
     &                 QBxy.eq.Zero .and.
     &                 QBxz.eq.Zero .and.
     &                 QByy.eq.Zero .and.
     &                 QByz.eq.Zero .and.
     &                 QBzz.eq.Zero
               Else
                  Call WarningMessage(2,'Option not implemented yet!')
                  Call Quit_OnUserError()
               End If

               If (NoLoop) Go To 203
               ZAZB = ZA * ZB
               B(1:3) = XF(1:3,jFd)
               iChxyz=iChAtm(B,iOper,nOper,iChBas(2))
               Call Stblz(iChxyz,iOper,nIrrep,mStb,jStb,iDum,jCoSet)
*              Introduce factor to ensure that contributions from
*              A>B are the only to be accumulated.
               Fact = One
               If (iFd.eq.jFd) Fact = Half
*
*              Find the DCR for the two centers
*
               Call DCR(LmbdR,iOper,nIrrep,iStb,nStb,jStb,mStb,
     &                  iDCRR,nDCRR)
*
               temp = Zero
               Do iR = 0, nDCRR-1
                  RB(1) = DBLE(iPhase(1,iDCRR(iR)))*B(1)
                  RB(2) = DBLE(iPhase(2,iDCRR(iR)))*B(2)
                  RB(3) = DBLE(iPhase(3,iDCRR(iR)))*B(3)
                  DRBx  = DBLE(iPhase(1,iDCRR(iR)))*DBx
                  DRBy  = DBLE(iPhase(2,iDCRR(iR)))*DBy
                  DRBz  = DBLE(iPhase(3,iDCRR(iR)))*DBz
                  QRBxx = QBxx
                  QRByy = QByy
                  QRBzz = QBzz
                  QRBxy = DBLE(iPhase(1,iDCRR(iR))
     &                  *iPhase(2,iDCRR(iR)))*QBxy
                  QRBxz = DBLE(iPhase(1,iDCRR(iR))
     &                  *iPhase(3,iDCRR(iR)))*QBxz
                  QRByz = DBLE(iPhase(2,iDCRR(iR))
     &                  *iPhase(3,iDCRR(iR)))*QByz


*                 The index A=RB is illegal.
                  If (.Not.EQ(A,RB)) Then
                     x=A(1)-RB(1)
                     y=A(2)-RB(2)
                     z=A(3)-RB(3)
                     r12 = Sqrt(x**2 + y**2 + z**2 )

                    eZZ=ZAZB*One/r12
                    eDZ=ZB*(-(DAx *x+DAy *y+DAz *z))/r12**3
                    eZD=ZA*(DRBx*x+DRBy*y+DRBz*z)/r12**3
                    eDD=(DAx*DRBx+DAy*DRBy+DAz*DRBz)/r12**3
     &                   -Three*(DAx* x+DAy* y+DAz *z)
     &                   *(DRBx*x+DRBy*y+DRBz*z)/r12**5

                    If(nOrd_XF.lt.2) Then
                       eTOT = eZZ+eZD+eDZ+eDD
                    Else

                    QAsum=(QAxx*x*x+QAyy*y*y+QAzz*z*z+2.0D0*
     &                   (QAxy*x*y+QAxz*x*z+QAyz*y*z))
                    QBsum=(QRBxx*x*x+QRByy*y*y+QRBzz*z*z+2.0D0*
     &                   (QRBxy*x*y+QRBxz*x*z+QRByz*y*z))

*                   Q-Z
                    eZQ=0.5D0*3.0D0/r12**5*ZA*QBsum
     &                   -0.5D0/r12**3*ZA*(QRBxx+QRByy+QRBzz)
                    eQZ=0.5D0*3.0D0/r12**5*ZB*QAsum
     &                   -0.5D0/r12**3*ZB*(QAxx+QAyy+QAzz)

*                   Q-D
                    eDQ=0.5D0*(-15.0D0/r12**7*
     & (DAx*x+DAy*y+DAz*z)*QBsum
     &                   +3.0D0/r12**5*
     &   (3.0D0*DAx*QRBxx*x
     &   +1.0D0*DAy*QRBxx*y
     &   +1.0D0*DAz*QRBxx*z
     &   +2.0D0*DAx*QRBxy*y
     &   +2.0D0*DAy*QRBxy*x
     &   +2.0D0*DAx*QRBxz*z
     &   +2.0D0*DAz*QRBxz*x
     &   +1.0D0*DAx*QRByy*x
     &   +3.0D0*DAy*QRByy*y
     &   +1.0D0*DAz*QRByy*z
     &   +2.0D0*DAy*QRByz*z
     &   +2.0D0*DAz*QRByz*y
     &   +1.0D0*DAx*QRBzz*x
     &   +1.0D0*DAy*QRBzz*y
     &   +3.0D0*DAz*QRBzz*z))

                    eQD=-0.5D0*(-15.0D0/r12**7*
     & (DRBx*x+DRBy*y+DRBz*z)*QAsum
     &                   +3.0D0/r12**5*
     &   (3.0D0*DRBx*QAxx*x
     &   +1.0D0*DRBy*QAxx*y
     &   +1.0D0*DRBz*QAxx*z
     &   +2.0D0*DRBx*QAxy*y
     &   +2.0D0*DRBy*QAxy*x
     &   +2.0D0*DRBx*QAxz*z
     &   +2.0D0*DRBz*QAxz*x
     &   +1.0D0*DRBx*QAyy*x
     &   +3.0D0*DRBy*QAyy*y
     &   +1.0D0*DRBz*QAyy*z
     &   +2.0D0*DRBy*QAyz*z
     &   +2.0D0*DRBz*QAyz*y
     &   +1.0D0*DRBx*QAzz*x
     &   +1.0D0*DRBy*QAzz*y
     &   +3.0D0*DRBz*QAzz*z))

*                   Q-Q (divided into two lines because of Fortran lim)
                    eQQ=0.25D0*(105.0D0/r12**9 * QAsum*QBsum
     &                   -15.0D0/r12**7*
     &  (6.0D0*QAxx*QRBxx*x*x+6.0D0*QAxx*QRBxy*x*y+6.0D0*QAxx*QRBxz*x*z
     &  +1.0D0*QAxx*QRByy*x*x+1.0D0*QAxx*QRByy*y*y+2.0D0*QAxx*QRByz*y*z
     &  +1.0D0*QAxx*QRBzz*x*x+1.0D0*QAxx*QRBzz*z*z+6.0D0*QAxy*QRBxx*x*y
     &  +4.0D0*QAxy*QRBxy*x*x+4.0D0*QAxy*QRBxy*y*y+4.0D0*QAxy*QRBxz*y*z
     &  +6.0D0*QAxy*QRByy*x*y+4.0D0*QAxy*QRByz*x*z+2.0D0*QAxy*QRBzz*x*y
     &  +6.0D0*QAxz*QRBxx*x*z+4.0D0*QAxz*QRBxy*y*z+4.0D0*QAxz*QRBxz*x*x
     &  +4.0D0*QAxz*QRBxz*z*z+2.0D0*QAxz*QRByy*x*z+4.0D0*QAxz*QRByz*x*y
     &  +6.0D0*QAxz*QRBzz*x*z+1.0D0*QAyy*QRBxx*x*x+1.0D0*QAyy*QRBxx*y*y
     &  +6.0D0*QAyy*QRBxy*x*y+2.0D0*QAyy*QRBxz*x*z+6.0D0*QAyy*QRByy*y*y
     &  +6.0D0*QAyy*QRByz*y*z+1.0D0*QAyy*QRBzz*y*y+1.0D0*QAyy*QRBzz*z*z
     &  +2.0D0*QAyz*QRBxx*y*z+4.0D0*QAyz*QRBxy*x*z+4.0D0*QAyz*QRBxz*x*y
     &  +6.0D0*QAyz*QRByy*y*z+4.0D0*QAyz*QRByz*y*y+4.0D0*QAyz*QRByz*z*z
     &  +6.0D0*QAyz*QRBzz*y*z+1.0D0*QAzz*QRBxx*x*x+1.0D0*QAzz*QRBxx*z*z
     &  +2.0D0*QAzz*QRBxy*x*y+6.0D0*QAzz*QRBxz*x*z+1.0D0*QAzz*QRByy*y*y
     & +1.0D0*QAzz*QRByy*z*z+6.0D0*QAzz*QRByz*y*z+6.0D0*QAzz*QRBzz*z*z))

                    eQQ=eQQ+0.25D0*3.0D0/r12**5*
     &   (3.0D0*QAxx*QRBxx
     &   +1.0D0*QAxx*QRByy
     &   +1.0D0*QAxx*QRBzz
     &   +4.0D0*QAxy*QRBxy
     &   +4.0D0*QAxz*QRBxz
     &   +1.0D0*QAyy*QRBxx
     &   +3.0D0*QAyy*QRByy
     &   +1.0D0*QAyy*QRBzz
     &   +4.0D0*QAyz*QRByz
     &   +1.0D0*QAzz*QRBxx
     &   +1.0D0*QAzz*QRByy
     &   +3.0D0*QAzz*QRBzz)

                    eTOT = eZZ+eZD+eDZ+eDD+eZQ+eQZ+eDQ+eQD+eQQ
c                    write(*,*)'eZZ',eZZ
c                    write(*,*)'eDZ',eDZ
c                    write(*,*)'eZD',eZD
c                    write(*,*)'eDD',eDD
c                    write(*,*)'eZQ',eZQ
c                    write(*,*)'eQZ',eQZ
c                    write(*,*)'eDQ',eDQ
c                    write(*,*)'eQD',eQD
c                    write(*,*)'eQQ',eQQ
c                    write(*,*)'eTOT',eTOT

                    End If     ! if quadrupoles are present

                    temp=temp+eTOT

                  End If
               End Do
c               PXX = PXX + ( Fact * ( ZAZB*temp0 + ZB*temp1
c     &             + ZA*temp2 + temp3)
c     &             * DBLE(nIrrep) ) / DBLE(LmbdR)
               PXX = PXX + ( Fact * temp
     &             * DBLE(nIrrep) ) / DBLE(LmbdR)

*
 203           Continue
            End Do
 103        Continue
         End Do
*
         If (Show) Then
            Write (6,'(19X,A,F16.8,A)')
     &                 ' External Field Potential Energy         ',
     &                    PXX,' au '
         End If
         Call Put_dScalar('PC Self Energy',PXX)
*
*        PotNuc=PotNuc+PXX
*

         If (Show) Then
            Write (6,*)
            Write (6,*)
            Write (6,'(11X,A,F16.8,A)')
     &     ' Total Nuclear Potential Energy        ',PotNuc,' au'
            Write (6,*)
         End If

1254  Continue
      End If
*
      Call Put_dScalar('PotNuc',PotNuc)
      If (isstructure().eq.1) then
         Call Add_Info('PotNuc',[PotNuc],1,6)
      Else
         Call Add_Info('PotNuc',[PotNuc],1,12)
      End If
*
      Call qExit('DrvN0')
      Return
      End
