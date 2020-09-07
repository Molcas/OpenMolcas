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
      Subroutine MkDisp()
      Use Basis_Info
      use Phase_Info
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "disp.fh"
#include "diff.fh"
#include "WrkSpc.fh"
      Logical TstFnc, Type
C Purpose: Pick out the following data from INFO:
C nUqCnt=Nr of unique centers
C nAlCnt  =Nr of centers (total)
C Coor(3,ic), i=1..nUqCnt, The cartesian coordinates (Unique)
C Coor(3,ic), i=1..nAlCnt,   The cartesian coordinates (All)
C iCntId(Key,i), i=1..nAlCnt, Information about centers:
C iCntId(1,ic)=The center type
C iCntId(2,ic)=An offset,=nr of symmetry-unique centers of earlier type
C iCntId(3,ic)=The symmetry-unique center of this type
C iCntId(4,ic)=(Usually 0) The symmetry op. that relates this center
C             to the corresponding symmetry-unique center
C
C nCntDf=Nr of differentiated centers
C ic=iCntDf(iDf), Which center is being differentiated
C The arrays are dynamic, and created here. Pointers and other data
C is kept in common /DIFF/, see file 'diff.fh'.

      nDiff=0
      Call IniSew(.false.,nDiff)

C Sizes:
      nUqCnt=0
      nAlCnt=0
      ic=0
      Do iCnttp=1,nCnttp
        If(.not.dbsc(iCnttp)%pChrg) Then
          Do iCnt=1,dbsc(iCnttp)%nCntr
            ic=ic+1
            nAlCnt=nAlCnt+nIrrep/nStab(ic)
          End Do
          nUqCnt=nUqCnt+dbsc(iCnttp)%nCntr
        End If
      End Do

C NB. We skip completely the 'pseudo-centers' in the loops below.
C Is this correct?
      Call GetMem('Coor','Allo','Real',ipCoor,3*nAlCnt)
      Call GetMem('iCntId','Allo','Inte',ipCntId,4*nAlCnt)
C Present setup: nCntDf=nUqCnt (see notes below)
      nCntDf=nUqCnt
      Call GetMem('iCntDf','Allo','Inte',ipCntDf,nCntDf)

C The nuclear coordinates: First, the symmetry-unique ones.
      nUqCnt=0
      ic=0
      Do iCnttp=1,nCnttp
        If(.not.dbsc(iCnttp)%pChrg) Then
          Do iCnt=1,dbsc(iCnttp)%nCntr
            ic=ic+1
            iWork(ipCntId+0+4*(ic-1))=iCnttp
            iWork(ipCntId+1+4*(ic-1))=nUqCnt
            iWork(ipCntId+2+4*(ic-1))=iCnt
            iWork(ipCntId+3+4*(ic-1))=0
            Work(0+ipCoor+3*(ic-1))=dbsc(iCnttp)%Coor(1,iCnt)
            Work(1+ipCoor+3*(ic-1))=dbsc(iCnttp)%Coor(2,iCnt)
            Work(2+ipCoor+3*(ic-1))=dbsc(iCnttp)%Coor(3,iCnt)
          End Do
          nUqCnt=nUqCnt+dbsc(iCnttp)%nCntr
        End If
      End Do
C Then, the symmetry related nuclei:
      jc=nUqCnt
      Do ic=1,nUqCnt
        iCnttp=iWork(ipCntId+0+4*(ic-1))
        iOff  =iWork(ipCntId+1+4*(ic-1))
        iCnt  =iWork(ipCntId+2+4*(ic-1))
        Do iOpNr=1,nIrrep/nStab(ic) -1
          jc=jc+1
          iWork(ipCntId+0+4*(jc-1))=iCnttp
          iWork(ipCntId+1+4*(jc-1))=iOff
          iWork(ipCntId+2+4*(jc-1))=iCnt
          iOp=iCoSet(iOpNr,0,ic)
          iWork(ipCntId+3+4*(jc-1))=iOp
          x=DBLE(iPhase(1,iOp))*Work(0+ipCoor+3*(ic-1))
          y=DBLE(iPhase(2,iOp))*Work(1+ipCoor+3*(ic-1))
          z=DBLE(iPhase(3,iOp))*Work(2+ipCoor+3*(ic-1))
          Work(0+ipCoor+3*(jc-1))=x
          Work(1+ipCoor+3*(jc-1))=y
          Work(2+ipCoor+3*(jc-1))=z
        End Do
      End Do
      nAlCnt=jc

C At present, we will assume that the differentiated nuclei are simply
C the symmetry-unique centers, in standard seward order:
      nCntDf=nUqCnt
      Do idf=1,nCntDf
        iWork(ipCntDf+iDf-1)=iDf
      End Do

C With this setup, the differentiated center nr I will simply be the
C same as the symmetry-unique center nr I.
C But this setup may change. It is thus important that externally, the
C differentiated nucleus nr iDf will correspond to the symmetry-unique
C center nr ic=iCntDf(I), which is presumably the number it will have
C on the MCKINT file; unless, of curse, ic is larger than nUqCnt,
C in which case the number of the center in the MCKINT file will be
C actually iCntId(2,ic)+iCntId(3,ic), and integrals with have to
C change sign according to rChTbl(iIrrep,iCntId(4,ic)).

      Write(6,*)' MKDISP: Tables of coordinates and displacements.'
      Write(6,*)' Nr of symmetry-unique centers, nUqCnt=',nUqCnt
      Write(6,*)' Nr of centers (total)        , nAlCnt=',nAlCnt
      Write(6,*)' Table over centers: Center, Center Type, Offset'
      Write(6,*)'   Center within its type, Symmetry op.,'
      Write(6,*)'   and coordinates:'
      Write(6,*)'   iC iCnttp iOff iCnt iOp        x'//
     &            '               y               z'
      Do ic=1,nAlCnt
        iCnttp=iWork(ipCntId+0+4*(ic-1))
        iOff  =iWork(ipCntId+1+4*(ic-1))
        iCnt  =iWork(ipCntId+2+4*(ic-1))
        iOp   =iWork(ipCntId+3+4*(ic-1))
        x=Work(0+ipCoor+3*(ic-1))
        y=Work(1+ipCoor+3*(ic-1))
        z=Work(2+ipCoor+3*(ic-1))
        Write(6,'(1x,5I5,3F16.8)')ic,iCnttp,iOff,iCnt,iOp,x,y,z
      End Do
C-------------------------------------------
CTL2004-start
C-------------------------------------------
      mDisp = 0
      mdc = 0
      Do  iCnttp = 1, nCnttp
         Do  iCnt = 1, dbsc(iCnttp)%nCntr
            mdc = mdc + 1
            mDisp = mDisp + 3*(nIrrep/nStab(mdc))
         End Do
      End Do
      !
      Call ICopy(mxdc*8,[0],0,IndDsp,1)
      Call ICopy(mxdc*3,[0],0,InxDsp,1)
      nDisp = 0
      Do iIrrep = 0, nIrrep-1
         lDisp(iIrrep) = 0
         Type = .True.
*        Loop over basis function definitions
         mdc = 0
         mc  = 1
         Do iCnttp = 1, nCnttp
*           Loop over unique centers associated with this basis set.
            Do iCnt = 1, dbsc(iCnttp)%nCntr
               mdc = mdc + 1
               IndDsp(mdc,iIrrep) = nDisp
*              Loop over the cartesian components
               Do iCar = 0, 2
                  iComp = 2**iCar
                  If ( TstFnc(iOper,nIrrep,iCoSet(0,0,mdc),
     &                 nIrrep/nStab(mdc),iChTbl,iIrrep,
     &                 iComp,nStab(mdc)) ) Then
                     nDisp = nDisp + 1
                     If (nDisp.gt.mDisp) Then
                        Write (6,*) 'nDisp.gt.mDisp'
                        Call AbEnd()
                     End If
                     If (iIrrep.eq.0) InxDsp(mdc,iCar+1) = nDisp
                     lDisp(iIrrep) = lDisp(iIrrep) + 1
                     Type = .False.
                  End If
               End Do
               mc = mc + nIrrep/nStab(mdc)
            End Do
         End Do
      End Do
C-------------------------------------------
CTL2004-end
C-------------------------------------------
      Return
      End
