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
      Subroutine Start_MCLR()
************************************************************************
*                                                                      *
*     Precompute whatever can be before starting the response section  *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit real*8 (a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "WrkSpc.fh"
      Character*5 Fname
*----------------------------------------------------------------------*
*     start                                                            *
*----------------------------------------------------------------------*
      Call QEnter('Start_MCLR')
*
      call setup_MCLR(1)
*
      If ((StepType.ne.'RUN2').and.(iAnd(kPrint,4).eq.4))
     &    Write(6,*) 'Transformation of integrals'
*     For the mp2-gradient calculations we want the transformation
*     routine to produce all integrals of the occupied and virtual
*     orbitals so we tell it that the whole space is inactive and
*     that the active is empty
*
*-----Use the driver from the CASPT2 code (only none-squared).
*
*
      ! re-direct the transformed integrals to the MOTRA file
      ! which is preserved at the end of the calculation.
      ! LuTri1 is deleted.
      If (TwoStep) Then
         LuTri1=LuMOTRA
         FnTri1=FnMOTRA
         Call DaName_MF_wa(LuQDAT,FnQDAT)
      End If

      Call DaName_MF_wa(LuTri1,FnTri1)

      If (newCho) Then
*
**       Compute inverse CMO
*
        lTriDens=0
        lSqrDens=0
        nOrbBas =0
        Do iSym=1,nSym
          lTriDens = lTriDens + nBas(iSym)*(nBas(iSym)+1)/2
          lSqrDens = lSqrDens + nBas(iSym)**2
          nOrbBas  = nOrbBas  + nOrb(iSym)*nBas(iSym)
        End Do
        Call Getmem('OverlapT','Allo','Real', ip_STmat,lTriDens)
        Call Getmem('OverlapS','Allo','Real', ip_Smat,lSqrDens)
*
        iSymlbl=1
        Call RdOne(irc,6,'Mltpl  0',1,Work(ip_STmat),iSymlbl)
*
        index = 0
        iOff = 0
        Do iSym = 1, nSym
           Do i = 1, nBas(iSym)
              Do j = 1, i-1
              Work(ip_Smat + j-1 + nBas(iSym)*(i-1)+iOff) =
     &                        Work(ip_STmat+index)
              Work(ip_Smat + i-1 + nBas(iSym)*(j-1)+iOff) =
     &                        Work(ip_STmat+index)
              index = index + 1
              End Do
              Work(ip_Smat + i-1 + nBas(iSym)*(i-1)+iOff) =
     &                        Work(ip_STmat+index)
              index = index + 1
           End Do
           ioff=ioff+nBas(iSym)**2
        End Do
        Call Getmem('OverlapT','Free','Real', ip_STmat,lTriDens)
*
        Call GetMem('CMO_inv','Allo','Real',ip_CMO_inv,nOrbBas)
        iOff1 = 0
        iOff2 = 0
        Do iSym = 1, nSym
           Call dGemm_('T','N', nOrb(iSym),nBas(iSym),nBas(iSym),
     &                1.0d0,Work(ipCMO+iOff2), nBas(iSym),
     &                      Work(ip_Smat+iOff1), nBas(iSym),
     &                0.0d0,Work(ip_CMO_inv+iOff2), nOrb(iSym))
*
           iOff1 =  iOff1 + nBas(iSym)**2
           iOff2 =  iOff2 + nOrb(iSym)*nBas(iSym)
        End Do
*
        Call Getmem('OverlapS','Free','Real', ip_Smat,lSqrDens)
      EndIf
*
      Call SetUp_CASPT2_Tra(nSym,nBas,nOrb,nIsh,nAsh,
     &                      nFro,nDel,ipCMO,nDens2,
     &                      LuTri1,LuTri2,LuHlf2,LuHlf3)
      iType=3  ! Means that TraCtl is called by MCLR

      If ((.not.newCho).and.(StepType.ne.'RUN2')) Then
        Call TraCtl_Drv(iType,.True.,1)
      End If

      If(TwoStep.and.(StepType.eq.'RUN2')) Then
        ! fetch some data from existing file LuTri1
        ! (from a previous MCLR run)
        ! and make it available to ERI common block intgrl.fh
        ! (LuTRI1=LuMOTRA)
        Call put_temp_data_on_intgrl(LuMOTRA,nSym,nOrb,nIsh,nAsh)
      End If
*
*
*     Init Cholesky informations
      If (NewCho) Then
         BufFrac=0.3D0
         Call Cho_X_Init(irc,BufFrac)
         iSeed=10
         Do i=1,nsym
           LuAChoVec(i) = IsFreeUnit(iSeed)
           iseed=LuAChoVec(i)+1
           Write(Fname,'(A4,I1)') 'CHTA',i
           Call DANAME_MF_WA(LuAChoVec(i),Fname)
         End Do
         Do i=1,nsym
           LuIChoVec(i) = IsFreeUnit(iSeed)
           iSeed=LuIChoVec(i)+1
           Write(Fname,'(A4,I1)') 'CHTI',i
           Call DANAME_MF_WA(LuIChoVec(i),Fname)
         End Do
         LuChoInt(1) = IsFreeUnit(iSeed)
         Write(Fname,'(A4)') 'CHIN'
         Call DANAME_MF_WA(LuChoInt(1),Fname)
         LuChoInt(2) = IsFreeUnit(iSeed)
         Write(Fname,'(A4)') 'CHTU'
         Call DANAME_MF_WA(LuChoInt(2),Fname)
      EndIf
*
      Call DaClos(LuTri2)
      Call DaClos(LuHlf2)
      Call DaClos(LuHlf3)
*
      Call FckMat
      Call StPert
*
**    With Cholesky there is no other choice than computing some
**    integrals used for the preconditioner
*
      If (NewCho) Then
         Call cho_prec_mclr(ipCMO,nIsh,nASh,LuAChoVec,LuChoInt)
      EndIf

*
*----------------------------------------------------------------------*
*     exit                                                             *
*----------------------------------------------------------------------*
      Call QExit('Start_MCLR')
      Return
      End

      Subroutine put_temp_data_on_intgrl(LUINTMZ_, NSYMZ_, NORBZ_,
     &                                   NISHZ_, NASHZ_  )
      Implicit None
      Integer ::       LUINTMZ_, NSYMZ_
      Integer ::       NORBZ_(8), NOSHZ_(8), NISHZ_(8), NASHZ_(8)
      Integer ::       nLength, iAddress, i
#include "intgrl.fh"
      iAddress=0
      IAD2M(1:3,1:36*36)=0
      nLength=3*36*36
      NOSHZ_(1:8)=0
      ! read the address list from the existing file
      Call iDaFile(LUINTMZ_,2,IAD2M,nLength,iAddress)
      NSYMZ=NSYMZ_
      LUINTMZ=LUINTMZ_
      Do i=1,NSYMZ_
        NORBZ(i) = NORBZ_(i)
        NOSHZ(i) = NISHZ_(i) + NASHZ_(i)
      End Do
      Return
      End
