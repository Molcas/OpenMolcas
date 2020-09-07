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
      Subroutine RF_Coord(
     &                 nq,nAtoms,iIter,nIter,Cx,iOper,nSym,jStab,
     &                 nStab,nDim,Smmtrc,Process,Value,
     &                 nB,iANr,qLbl,iRef,fconst,
     &                 rMult,LuIC,Indq,dMass,iCoSet,
     &                 Proc_dB,mB_Tot,mdB_Tot,
     &                 BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,nqB)
      use Phase_Info
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "stdalloc.fh"
#include "real.fh"
#include "sbs.fh"
#include "print.fh"
      Real*8 Cx(3,nAtoms,nIter),
     &       dMass(nAtoms), fconst(nB), Value(nB,nIter), rMult(nB),
     &       Trans(3), RotVec(3), RotMat(3,3),
     &       BM(nB_Tot), dBM(ndB_Tot)
      Integer   nStab(nAtoms), iOper(0:nSym-1), iCoSet(0:7,nAtoms),
     &          jStab(0:7,nAtoms), nqB(nB),
     &          iANr(nAtoms), Indq(3,nB), iBM(nB_Tot), idBM(2,ndB_Tot)
      Logical Smmtrc(3,nAtoms), Process, PSPrint,
     &        TransVar, RotVar, Proc_dB, Invariant
      Character*3 TR_type(6)
      Character*14 Label, qLbl(nB)
#include "ddvdt_RF.fh"
      Real*8, Dimension(:), Allocatable :: xMass
      Real*8, Dimension(:,:), Allocatable :: currXYZ, Ref123, Grad,
     &                                       dRVdxyz, Hess
      Real*8, Dimension(:,:,:), Allocatable :: d2RV
      Integer, Dimension(:), Allocatable :: Ind, iDCR
      Dimension dum(1)
      Data TR_type/'Tx ','Ty ','Tz ','Ryz','Rzx','Rxy'/
*
      iRout=151
      iPrint=nPrint(iRout)
      Call QEnter('RF_Coords')
*
      TransVar=iAnd(iSBS,2**7).eq. 2**7
      RotVar=iAnd(iSBS,2**8).eq. 2**8
      If (.Not.RotVar.and..Not.TransVar) Go To 99
*                                                                      *
************************************************************************
*                                                                      *
      nqRF=0
      PSPrint=.False.
      If (iPrint.ge.99) PSPrint=.True.
      If (PSPrint) Write (6,*) ' Enter RF_Coords.'
*
*---- Find nCent and allocate
*
      nCent=0
      Do iAtom = 1, nAtoms
         nCent=nCent+nSym/nStab(iAtom)
      End Do
      mB = nCent*3
      Call mma_allocate(currXYZ,3,nCent,label='currXYZ')
      Call mma_allocate(Ref123,3,nCent,label='Ref123')
      Call mma_allocate(Grad,3,nCent,label='Grad')
      Call mma_allocate(dRVdxyz,3,3*nCent,label='dRVdxyz')
      Call mma_allocate(xMass,nCent,label='xMass')
      Call mma_allocate(Ind,nCent,label='Ind')
      Call mma_allocate(iDCR,nCent,label='iDCR')
      Call mma_allocate(Hess,mB,mB,label='Hess')
*
*---- Find index of RF center (origin), etc
*
      iCent=0
      Do iAtom = 1, nAtoms
         x=Cx(1,iAtom,iIter)
         y=Cx(2,iAtom,iIter)
         z=Cx(3,iAtom,iIter)
*
         x_ref=Cx(1,iAtom,iRef)
         y_ref=Cx(2,iAtom,iRef)
         z_ref=Cx(3,iAtom,iRef)
         Do i = 0, nSym/nStab(iAtom)-1
            iCent = iCent + 1
            iFacx=iPhase(1,iCoSet(i,iAtom))
            iFacy=iPhase(2,iCoSet(i,iAtom))
            iFacz=iPhase(3,iCoSet(i,iAtom))
*
            currXYZ(1,iCent)=DBLE(iFacx)*x
            currXYZ(2,iCent)=DBLE(iFacy)*y
            currXYZ(3,iCent)=DBLE(iFacz)*z
*
            Ref123(1,iCent)=DBLE(iFacx)*x_ref
            Ref123(2,iCent)=DBLE(iFacy)*y_ref
            Ref123(3,iCent)=DBLE(iFacz)*z_ref
*
            Ind(iCent)=iAtom
            iDCR(iCent)=iCoSet(i,iAtom)
         End Do
      End Do
      nCurrXYZ=iCent
*
      Fact=One
      If (.Not.RotVar) Fact=2.0D-2
*
*     Write (6,*) 'nCent,nDim=',nCent,nDim
*     Write (6,*) (Ind(iCent),iCent=1,nCent)
*
      TMass = Zero
      Do iCent = 1, nCent
         iAtom = Ind(iCent)
         xMass(iCent) = dMass(iAtom)
         TMass = TMass + dMass(iAtom)
      End Do
*---- Loop over cartesian components
*
      Do ixyz = 1, 3
*
         Invariant=.False.
         iTest=2**(ixyz-1)
         Do iSym = 0, nSym-1
            If (iOper(iSym).eq.iTest) Invariant=.True.
         End Do
         If (Invariant) Go To 199
*
*------- Compute total mass and center of mass of the molecule, the center
*        of the RF cavity is at origin with infinite mass. Hence, the latter
*        is ignored!
*
         COM_xyz = Zero
         Do iCent = 1, nCent
            COM_xyz = COM_xyz + currXYZ(ixyz,iCent)*xMass(iCent)
         End Do
         COM_xyz=COM_xyz/TMass
*
         If (.Not.TransVar) Go To 199
         f_Const=One
*        Write (6,*) ' RF Force Constant:',f_Const
*
         iDeg=1
         Deg=Sqrt(DBLE(iDeg))
*
         nq = nq + 1
         If (.Not.Process) mB_Tot = mB_Tot + mB
         If (.Not.Proc_dB) mdB_Tot = mdB_Tot + mB**2
         nqRF = nqRF + 1
         Write (LuIC,'(A,I2.2,2A)')
     &              'TR',nqRF,' = ',TR_type(ixyz)
         Label=' '
         Write (Label,'(A,I2.2)') 'TR',nqRF
*
         Val = COM_xyz
*
*------- Compute the gradient
*
         call dcopy_(mB,[Zero],0,Grad,1)
         Do iCent = 1, nCent
            iAtom=Ind(iCent)
*           Write (6,*) 'iAtom,iCOM=',iAtom,iCOM
            Grad(ixyz,iCent) = dMass(iAtom)/TMass
         End Do
C        Call RecPrt('Grad (Trans)',' ',Grad,3,nCent)
*
*------- Second derivative is trivially zero!
*
         Call FZero(Hess,mB**2)
         If (Process) Then
*
            Indq(1,nq) = -2**(ixyz)/2
            Indq(2,nq) = 0
            Indq(3,nq) = 0
*
C           fconst(nq)=Sqrt(Fact*Trans_Const)
            fconst(nq)=Sqrt(Trans_Const)
            rMult(nq)=Deg
*
            Value(nq,iIter)=Val
            qLbl(nq)=Label
*
*--------   Project the gradient vector
*
            Call ProjSym(nAtoms,nCent,Ind,nStab,jStab,currXYZ,
     &                   iDCR,Grad,Smmtrc,nDim,PSPrint,
     &                   Hess,mB_Tot,mdB_Tot,
     &                   BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,
     &                   Proc_dB,nqB,nB,nq,rMult(nq))
*
         End If
*
 199     Continue
      End Do
*                                                                      *
************************************************************************
*                                                                      *
C     Write (6,*) 'RotVar=',RotVar
      If (.Not.RotVar) Go To 98
*
*     A la Malmqvist
*
      nOrder=2
      nMass=nCent
      Call FZero(Trans,3)
      Call FZero(RotVec,3)
      Call mma_allocate(d2RV,3,3*nCent,3*nCent,label='d2RV')
C     Call RecPrt('xMass',' ',xMass,1,nMass)
      Call RotDer(nMass,xMass,currXYZ,ref123,trans,RotAng,
     &            RotVec,RotMat,nOrder,dRVdXYZ,d2RV,dum,dum)
C     Call RecPrt('RotVec',' ',RotVec,1,3)
C     Call RecPrt('RotMat',' ',RotMat,3,3)
C     Call RecPrt('dRVdXYZ',' ',dRVdXYZ,3,3*nMass)
*
      Do ixyz = 1, 3
*
         Invariant=.False.
         If (ixyz.eq.1) Then
            iTest=6
         Else If (ixyz.eq.2) Then
            iTest=5
         Else
            iTest=3
         End If
         Do iSym = 0, nSym-1
            If (iOper(iSym).eq.iTest) Invariant=.True.
         End Do
         If (Invariant) Go To 299
*
         jxyz = ixyz+1
         If (jxyz.gt.3) jxyz=1
         kxyz = jxyz+1
         If (kxyz.gt.3) kxyz=1
         iDeg=1
         Deg=Sqrt(DBLE(iDeg))
*
         nq = nq + 1
         If (.Not.Process) mB_Tot = mB_Tot + mB
         If (.Not.Proc_dB) mdB_Tot = mdB_Tot + mB**2
         nqRF = nqRF + 1
         Write (LuIC,'(A,I2.2,2A)')
     &              'TR',nqRF,' = ',TR_type(ixyz+3)
         Label=' '
         Write (Label,'(A,I2.2)') 'TR',nqRF
*
         Val = RotVec(ixyz)
*
*------- Compute the gradient
*
         call dcopy_(mB,[Zero],0,Grad,1)
         call dcopy_(mB,dRVdXYZ(ixyz,1),3,Grad,1)
C        Call RecPrt('Grad (Rot)',' ',Grad,3,nCent)
*
*------- Second derivative
*
         Call FZero(Hess,mB**2)
         If (Proc_dB) Call DCopy_(mB**2,d2RV(ixyz,1,1),3,Hess,1)
*
         If (Process) Then
*
            Indq(1,nq) = -( 2**(jxyz)/2  + 2**(kxyz)/2)
            Indq(2,nq) = 0
            Indq(3,nq) = 0
*
            fconst(nq)=Sqrt(Rot_Const)
            rMult(nq)=Deg
*
            Value(nq,iIter)=Val
            qLbl(nq)=Label
*
*--------   Project the gradient vector
*
            Call ProjSym(nAtoms,nCent,Ind,nStab,jStab,currXYZ,
     &                   iDCR,Grad,Smmtrc,nDim,PSPrint,
     &                   Hess,mB_Tot,mdB_Tot,
     &                   BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,
     &                   Proc_dB,nqB,nB,nq,rMult(nq))
*
         End If
*
 299     Continue
      End Do
      Call mma_deallocate(d2RV)
C     Write (6,*) 'nqRF=',nqRF
*                                                                      *
************************************************************************
*                                                                      *
 98   Continue
      Call mma_deallocate(currXYZ)
      Call mma_deallocate(Ref123)
      Call mma_deallocate(Grad)
      Call mma_deallocate(dRVdxyz)
      Call mma_deallocate(xMass)
      Call mma_deallocate(Ind)
      Call mma_deallocate(iDCR)
      Call mma_deallocate(Hess)
 99   Continue
      Call QExit ('RF_Coords')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(iANr)
      End
