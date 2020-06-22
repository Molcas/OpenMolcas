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
      Subroutine No_ESPF(natom,Forces,DoTinker)
      use external_centers
      Implicit Real*8 (a-h,o-z)
*
#include "espf.fh"
*
      Real*8 A(3), B(3), RB(3)
      Integer iDCRR(0:7), jCoSet(8,8), iStb(0:7)
      Logical EQ,NoLoop,Forces,DoTinker,DoRys
      Character*180 Line,Get_Ln
      External Get_Ln
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      Call qEnter('no_espf')
      iPL = iPL_espf()
      NoLoop = .True.
*
*     Nothing to do in a gradient calc
*
      If (Forces) Return
*
*     Update the nuclear repulsion energy
*
      Call Get_dScalar('PotNuc',RepNuc)
      RepNuc_old = RepNuc
*
*     Read the MM contribution to the total energy and add it
*     to the Nuclear Repulsion term
*
      If (DoTinker) Then
         Tke = Zero
         ITkQMMM = IsFreeUnit(30)
         Call Molcas_Open (ITkQMMM,'QMMM')
         Line = ' '
         Do While (Index(Line,'TheEnd ') .eq. 0)
            Line=Get_Ln(ITkQMMM)
            If (Index(Line,'MMEnergy ').ne.0) Call Get_F1(2,TkE)
         End Do
         Close (ITkQMMM)
         TkE = TkE * ToHartree
         RepNuc = RepNuc + TkE
         If (iPL.ge.3) Write(6,1000) RepNuc_old,TkE,RepNuc
1000     Format(/,' RepNuc + MM = ',F13.8,' + ',F13.8,' = ',F13.8)
      End If
*
      If (lXF.and.(nOrd_XF.ge.0)) Then
         write(6,*) 'Here we are!!'
*
c         Call Seward_Init()
         DoRys=.True.
         nDiff=0
         Call GetInf(Info,nInfo,DoRys,nDiff,1)
         Primitive_Pass=.True.
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
         write(6,*) 'Work(ip + ...) = ',(Work(ip+ii),ii=1,4)
         write(6,*) 'nCnttp = ',nCnttp
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
            A(1:3)=XF(1:3,iFd)
            iChxyz=iChAtm(A,iOper,nOper,iChBas(2))
            Call Stblz(iChxyz,iOper,nIrrep,nStb,iStb,iDum,jCoSet)
*
            ndc = 0
            Do jCnttp = 1, nCnttp
               ZB = Charge(jCnttp)
               If (pChrg(jCnttp)) Go To 202
               If (ZB.eq.Zero) Go To 202
               If (FragCnttp(jCnttp)) Go To 202
               ZAZB = ZA * ZB
               jxyz = ipCntr(jCnttp)
               Do jCnt = 1, nCntr(jCnttp)
                  B(1) = Work(jxyz  )
                  B(2) = Work(jxyz+1)
                  B(3) = Work(jxyz+2)
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
                        If (ECP(jCnttp)) Then
*--------------------------Add contribution from M1 operator
                           Do iM1xp=0, nM1(jCnttp)-1
                             Gamma = Work(ipM1xp(jCnttp)+iM1xp)
                             CffM1 = Work(ipM1cf(jCnttp)+iM1xp)
                             fab = fab + CffM1 * Exp(-Gamma*r12**2)
                           End Do
*--------------------------Add contribution from M2 operator
                           Do iM2xp=0, nM2(jCnttp)-1
                             Gamma = Work(ipM2xp(jCnttp)+iM2xp)
                             CffM2 = Work(ipM2cf(jCnttp)+iM2xp)
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
                  jxyz = jxyz + 3
               End Do
 202           Continue
               ndc = ndc + nCntr(jCnttp)
            End Do
 102        Continue
         End Do
*
         If (iPL.ge.3) Write(6,1100) RepNuc,PNX,RepNuc+PNX
1100     Format(/,' RepNuc + Point charges = ',F13.8,' + ',F13.8,' = ',
     &          F13.8)
*
         PotNuc=PotNuc+PNX
         Call Put_dScalar('PotNuc',RepNuc)
      End If
*
*     Update the 1-e integrals
*
      Call qExit('no_espf')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(natom)
      End
