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
* Copyright (C) 1995, Roland Lindh                                     *
*               2000, Valera Veryazov                                  *
*               2014, Thomas Dresselhaus                               *
************************************************************************
      Subroutine MOEval(MOValue,nMOs,nCoor,CCoor,CMOs,nCMO,mCoor,DoIt,
     &                  nDrv,mAO,Debug)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from: Drv1EL                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*      Author:Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN. November 1995                           *
*                                                                      *
*      Modified: Thomas Dresselhaus, March 2014                        *
*                Added ability to calculate 2nd derivative as well     *
************************************************************************
      use Real_Spherical
      use Basis_Info
      use Center_Info
      use Phase_Info
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 A(3),Ccoor(3,nCoor),RA(3)
      Integer DoIt(nMOs)
      Integer nDrv ! Between 0 and 2. The highest derivative to be calc.
      Integer mAO  ! Memory slots per point and basis functions. Should
                   ! be >=1 for nDrv=0, >=4 for nDrv=1, >=10 for nDrv=2.
      Real*8 MOValue(mAO,nCoor,nMOs),CMOs(nCMO)
      Logical Debug
*
*     Statement functions
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      iRout = 112
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.99) Then
         Write (6,*) ' In MOEval'
      End If
      call dcopy_(mAO*nCoor*nMOs,[Zero],0,MOValue,1)
*
*     Loop over shells.
*
      iSkal=0
      Thr=0.0D0

      Do iAng = iAngMx , 0, -1

         If (MaxPrm(iAng).eq.0) goto 100
         If (MaxBas(iAng).eq.0) goto 100
*
*        Scratch area for contraction step
*
         nScr1 =  MaxPrm(iAng)* nElem(iAng)
         Call GetMem('Scrtch','ALLO','REAL',iScrt1,nScr1)
*
*        Scratch area for the transformation to spherical gaussians
*
         nScr2=MaxPrm(iAng)*nElem(iAng)
         Call GetMem('ScrSph','Allo','Real',iScrt2,nScr2)
*
*        Loop over basis sets. Skip if basis set do not include
*        angular momentum functions as specified above.
*
         iAOttp=0
         mdc = 0
         Do iCnttp = 1, nCnttp

            nTest = dbsc(iCnttp)%nVal
            If (iAng+1.gt.nTest)  Go To 101
            If (dbsc(iCnttp)%Aux) Go To 101
            If (dbsc(iCnttp)%Frag) Go To 101
            nCnt = dbsc(iCnttp)%nCntr
            iShll = dbsc(iCnttp)%iVal + iAng
            iPrim = Shells(iShll)%nExp
            If (iPrim.eq.0) Go To 101
            iBas  = Shells(iShll)%nBasis
            If (iBas.eq.0) Go To 101
            If (Shells(iShll)%Prjct ) Then
               iCmp = 2*iAng+1
            Else
               iCmp  = nElem(iAng)
            End If

            Call OrdExpD2C(iPrim,Shells(iShll)%Exp,iBas,
     &                           Shells(iShll)%pCff)
            kSh=dbsc(iCnttp)%iVal+iAng
*
*           Compute the total number of function for this
*           basis set, summed over all shells.
*
            Do jAng = 0, nTest-1
               jShll = dbsc(iCnttp)%iVal + jAng
               If (Shells(jShll)%Prjct ) Then
                  jCmp = 2*jAng+1
               Else
                  jCmp  = nElem(jAng)
               End If
            End Do
*
*           Loop over unique centers of basis set "iCnttp"
*
            Do iCnt = 1, nCnt

               iAO = iAOttp + (iCnt-1)*dbsc(iCnttp)%lOffAO
     &             + Shells(kSh)%kOffAO
               A(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)

               Do jAng = 0, iAng-1
                  jShll = dbsc(iCnttp)%iVal + jAng
                  If (Shells(jShll)%Prjct ) Then
                     jCmp = 2*jAng+1
                  Else
                     jCmp  = nElem(jAng)
                  End If
               End Do
*
*--------------Allocate memory for SO and AO values
*
               mRad     = nDrv + 1

               nForm    = 0
               Do iDrv  = 0, nDrv
                 nForm = nForm + nElem(iDrv)
               End Do
               nTerm    = 2**nDrv


               nAO=(iCmp*iBas*nCoor)*(mAO)
               nSO=nAO*nIrrep/dc(mdc+iCnt)%nStab
               nDeg=nIrrep/dc(mdc+iCnt)%nStab
               Call GetMem('AOs','Allo','Real',ipAOs,nAO)
               Call GetMem('SOs','Allo','Real',ipSOs,nSO)
               call dcopy_(nSO,[Zero],0,Work(ipSOs),1)
               nxyz=nCoor*3*(iAng+mRad)
               Call GetMem('xyz','Allo','Real',ipxyz,nxyz)
               ntmp=nCoor
               Call GetMem('tmp','Allo','Real',iptmp,ntmp)
               nRadial=iBas*nCoor*mRad
               Call GetMem('Radial','Allo','Real',ipRadial,nRadial)

               nAngular=5*nForm*nTerm
               Call GetMem('Angular','Allo','Inte',ipAng,nAngular)
*
*------------- Loops over symmetry operations operating on the basis
*              set center.
*
               Do iG = 0, nIrrep/dc(mdc+iCnt)%nStab - 1
                  iSkal=iSkal+1
                  Call OA(dc(mdc+iCnt)%iCoSet(iG,0),A,RA)
                  ipx=iPhase(1,dc(mdc+iCnt)%iCoSet(iG,0))
                  ipy=iPhase(2,dc(mdc+iCnt)%iCoSet(iG,0))
                  ipz=iPhase(3,dc(mdc+iCnt)%iCoSet(iG,0))
                  px=DBLE(iPhase(1,dc(mdc+iCnt)%iCoSet(iG,0)))
                  py=DBLE(iPhase(2,dc(mdc+iCnt)%iCoSet(iG,0)))
                  pz=DBLE(iPhase(3,dc(mdc+iCnt)%iCoSet(iG,0)))
                  nOp = NrOpr(dc(mdc+iCnt)%iCoSet(iG,0))
*
*---------------- Evaluate AOs at RA
*
                  call dcopy_(nAO,[Zero],0,Work(ipAOs),1)
                  mTmp=1
                  Call AOEval(iAng,nCoor,CCoor,Work(ipxyz),RA,
     &                        Shells(iShll)%Transf,
     &                        RSph(ipSph(iAng)),nElem(iAng),iCmp,
     &                        iWork(ipAng),nTerm,nForm,Thr,mRad,
     &                        iPrim,iPrim,Shells(iShll)%Exp,
     &                        Work(ipRadial),iBas,
     &                        Shells(iShll)%pCff,
     &                        Work(ipAOs),mAO,px,py,pz,ipx,ipy,ipz)
*
*---------------- Distribute contributions to the SOs
*
                  Call SOAdpt(Work(ipAOs),mAO,nCoor,iBas,iCmp,nOp,
     &                        Work(ipSOs),nDeg,iAO)
*
               End Do ! iG
*
*------------- Distribute contributions to the MOs
*
               Call SODist(Work(ipSOs),mAO,nCoor,iBas,iCmp,nDeg,
     &                     MOValue,nMOs,iAO,CMOs,nCMO,DoIt)
*
               Call GetMem('Radial','Free','Real',ipRadial,nRadial)
               Call GetMem('Angular','Free','Inte',ipAng,nAngular)
               Call GetMem('tmp','Free','Real',iptmp,ntmp)
               Call GetMem('xyz','Free','Real',ipxyz,nxyz)
               Call GetMem('AOs','Free','Real',ipAOs,nAO)
               Call GetMem('SOs','Free','Real',ipSOs,nSO)
*
            End Do ! iCnt
 101        Continue
            mdc = mdc + dbsc(iCnttp)%nCntr
            iAOttp = iAOttp + dbsc(iCnttp)%lOffAO*dbsc(iCnttp)%nCntr
         End Do
         Call GetMem('ScrSph','Free','Real',iScrt2,nScr2)
         Call GetMem('Scrtch','Free','Real',iScrt1,nScr1)
 100     continue
      End Do
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_integer(mCoor)
        Call Unused_logical(Debug)
      End If
      End


      Subroutine MOEvalDel(MOValueD,nMOs,
     &                     nCoor,CCoor,CMOs,nCMO,mCoor,DoIt,Debug)

      Implicit Real*8 (A-H,O-Z)
      Real*8 Ccoor(3,nCoor),MOValueD(4*nCoor*nMOs),CMOs(nCMO)
      Integer DoIt(nMOs),mAO
      Logical Debug
      integer nDrv

      mAO  = 4
      nDrv = 1

      Call MOEval(MOValueD, nMOs, nCoor, CCoor, CMOs, nCMO, mCoor,
     &                     DoIt, nDrv, mAO, DEBUG)

c        IJ1=1+(I-1)*4
c        IJ2=2+(I-1)*4
c        IJ3=3+(I-1)*4
c        IJ4=4+(I-1)*4
c
c        MOValue(I)=Work(MOTmp-1+IJ1)
c        MOValueDX(I)=Work(MOTmp-1+IJ2)
c        MOValueDY(I)=Work(MOTmp-1+IJ3)
c        MOValueDZ(I)=Work(MOTmp-1+IJ4)
c      END DO

      Return
      End

      Subroutine MOEvalDer(MOValue,iDir,nMOs,
     &                     nCoor,CCoor,CMOs,nCMO,mCoor,DoIt,Debug)

      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
      Real*8 Ccoor(3,nCoor),MOValue(nCoor*nMOs),CMOs(nCMO)
      Integer DoIt(nMOs),mAO
      Logical Debug
      integer nDrv

      mAO  = 4
      nDrv = 1

      Call GetMem('MOTMP','Allo','Real',iMoTmp,4*nCoor*nMOs)

      Call MOEval(work(iMoTmp), nMOs, nCoor, CCoor, CMOs, nCMO, mCoor,
     &                     DoIt, nDrv, mAO, DEBUG)

c iDir = 1 then do dX
c iDir = 2 then do dY
c iDir = 3 then do dZ
      write(6,*) "iDir:",iDir
      if(iDir.gt.0.and.iDir.lt.4) then
        DO I=1,nCoor*nMOs
          IJ=iDir+1+(I-1)*4
          MOValue(I)=Work(iMoTmp-1+IJ)
        END DO
      else ! do gradient
        DO I=1,nCoor*nMOs
          IJX=2+(I-1)*4
          IJY=3+(I-1)*4
          IJZ=4+(I-1)*4
          MOValue(I)=Work(iMoTmp-1+IJX)+
     &               Work(iMoTmp-1+IJY)+
     &               Work(iMoTmp-1+IJZ)
        END DO
      end if
      Call GetMem('MOTMP','Free','Real',iMoTmp,4*nCoor*nMOs)

      Return
      End
