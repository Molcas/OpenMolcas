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
* Copyright (C) 2006, Roland Lindh                                     *
************************************************************************
      Subroutine Gateway(iReturn)
************************************************************************
*                                                                      *
*     In the fall of 2006 Seward, Alaska and McKinley were in need     *
*     of a sibling. The new code would be the gateway into the         *
*     MOLCAS world for our future GUI. After some thoughts we          *
*     desided to call the gateway Gateway after the small village      *
*     Gateway, Alaska.                                                 *
*                                                                      *
*     Roland Lindh, 26th September 2006.                               *
*                                                                      *
************************************************************************
      use Period
      use GeoList
      use MpmC
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
      Integer AixRm
      External Get_Cho_1Center,AixRm
#include "itmax.fh"
#include "info.fh"
#include "status.fh"
#include "gateway.fh"
#include "rctfld.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "vrsn_gateway.fh"
      Character xLblCnt(MxAtom)*(LENIN)
      Character KWord*80, Header(2)*72
      Parameter (nMamn=MaxBfn+MaxBfn_Aux)
      Character*(LENIN8), Allocatable :: Mamn(:)
      Logical lOPTO, Pseudo, Do_OneEl
      Logical Get_Cho_1Center, Cho_1Center
CVV      LOGICAL GA_USES_MA,GA_MEMORY_LIMITED
C-SVC: identify runfile with a fingerprint
      Character cDNA*256
      Logical IsBorn, Found
      Real*8, Allocatable :: DCo(:,:), DCh(:), DCh_Eff(:)
*                                                                      *
************************************************************************
*                                                                      *
C     Call Gateway_banner()
      iReturn = 0
*
*     If Gateway is running the Run_Mode on the runfile should always
*     be G_Mode.
*
      Run_Mode=G_Mode
      Call MkRun(iRC,0)
      Call Put_iScalar('Run_Mode',Run_Mode)
*
*     Determine and save the fingerprint of the runfile in a field with
*     label 'BirthCertificate' if it is empty.  This allows us to
*     uniquely identify the runfile and any later associated files.
*
      Call qpg_cArray('BirthCertificate',IsBorn,nDNA)
      IF (.NOT.IsBorn) THEN
        Call Get_Genome(cDNA,nDNA)
        Call Put_cArray('BirthCertificate',cDNA,nDNA)
      END IF
*                                                                      *
************************************************************************
*                                                                      *
*     Get the memory size available
*
      Call SetMem('Clear=Off')
*                                                                      *
************************************************************************
*                                                                      *
*     Initialize common blocks and start from scratch
*
      Call Seward_Init()
      Call Funi_Init()
      Call NQGrid_Init()
      Call Basis_Info_Init()
*                                                                      *
************************************************************************
*                                                                      *
*     Spool the input
*
      LuSpool=21
      Call SpoolInp(LuSpool)
*                                                                      *
************************************************************************
*                                                                      *
*     Remove possible leftover files
*
      Call f_Inquire('UDC.Gateway',Found)
      If (Found) iRC = AixRm('UDC.Gateway')
      Call f_Inquire('UDC.NG',Found)
      If (Found) iRC = AixRm('UDC.NG')
*                                                                      *
************************************************************************
*                                                                      *
*     Read the input.
*
      lOPTO = .False.
      Call RdCtl_Seward(LuSpool,lOPTO,Do_OneEl)
*
*     Write the Basis_Info data to file. Release the arrays and read
*     them back from the runfile now allocating them to the proper
*     size.
*
      Call Basis_Info_Dmp()
      Call Basis_Info_Free()
      Call Basis_Info_Get()
*                                                                      *
************************************************************************
*                                                                      *
*     Close Spool file
*
      Call Close_LuSpool(LuSpool)
*                                                                      *
************************************************************************
*                                                                      *
*     Print out section
*
      Call Print_Symmetry()
      Call Flip_Flop(.False.)
      Call Print_Basis(lOPTO)
      Call Print_Geometry(0)
      Call Print_Isotopes()
      If (nPrint(2).gt.0) nPrint(117)=6
      Call RigRot(Centr,Mass,kCentr)
      Call Print_Basis2()
      Call Print_OpInfo()
*                                                                      *
************************************************************************
*                                                                      *
*     Genetate the SO/AO basis set here.
*
      Primitive_Pass=.False.
      Call Flip_Flop(Primitive_Pass)
      Call mma_allocate(Mamn,nMamn,label='Mamn')
      Call SOCtl_Seward(Mamn,nMamn)
*                                                                      *
************************************************************************
*                                                                      *
      Call DmpInf()
*                                                                      *
************************************************************************
*                                                                      *
*     Produce minimal set of entries on the runfile to facilitate
*     Grid_It's and ExpBas's needs.
*
      Call Drvn0()
*
      Call Datimx(KWord)
      Header(1)=Title(1)(5:76)
      Write (Header(2),'(4A)')
     &          ' Integrals generated by ',
     &            Vrsn,', ',KWord(1:24)
*
      Call Put_cArray('Seward Title',Header(1),144)
*
      Call Put_iScalar('NSYM',nIrrep)
      Call Put_iArray('Symmetry operations',iOper,nIrrep)
      Call Put_iScalar('Rotational Symmetry Number',iSigma)
      Call Put_cArray('Irreps',lIrrep(0),24)
      Call Put_cArray('Unique Basis Names',Mamn(1),(LENIN8)*nDim)
      Call Put_iArray('NBAS',nBas,nIrrep)
      call basis2run()
      Call mma_deallocate(Mamn)
*
*     Generate list of unique atoms
*
      nNuc = 0
      Do iCnttp = 1, nCnttp
         If (.Not.pChrg(iCnttp).and.
     &       .Not.dbsc(iCnttp)%Frag.and.
     &       .Not.dbsc(iCnttp)%Aux) nNuc = nNuc + dbsc(iCnttp)%nCntr
      End Do
*
      Call mma_allocate(DCo,3,nNuc)
      Call mma_allocate(DCh,nNuc)
      Call mma_allocate(DCh_Eff,nNuc)
      iDCo = 1
      iDCh = 1
      iDChE= 1
      mdc = 0
      iNuc = 0
      Do iCnttp = 1, nCnttp
         If (.Not.pChrg(iCnttp).and.
     &       .Not.dbsc(iCnttp)%Frag.and.
     &       .Not.dbsc(iCnttp)%Aux) Then
            Do iCnt = 1, dbsc(iCnttp)%nCntr
               mdc = mdc + 1
               iNuc = iNuc+ 1
               DCo(1:3,iNuc)=dbsc(iCnttp)%Coor(1:3,iCnt)
               DCh_Eff(iNuc)=Charge(iCnttp)
               DCh(iNuc)=DBLE(dbsc(iCnttp)%AtmNr)
               xLblCnt(iNuc)=LblCnt(mdc)(1:LENIN)
            End Do
         Else
            mdc  = mdc + dbsc(iCnttp)%nCntr
         End If
      End Do
      Call Put_iScalar('Unique atoms',nNuc)
      Call Put_dArray('Unique Coordinates',DCo,3*nNuc)
      Call Put_dArray('Nuclear charge',DCh,nNuc)
      Call Put_dArray('Effective nuclear Charge',DCh_Eff,nNuc)
      Call Put_cArray('Unique Atom Names',xLblCnt(1),LENIN*nNuc)
      Call Put_iArray('nStab',nStab,nNuc)
*
      Call mma_deallocate(DCo)
      Call mma_deallocate(DCh)
      Call mma_deallocate(DCh_Eff)
*
*     Manipulate the option flag
*
      iOption=0
      If (DirInt) iOption=iOr(iOption,1)
      If (Expert) iOption=iOr(iOption,2)
      If (lRF) iOption=iOr(iOption,4)
      If (lLangevin.or.iXPolType.gt.0) iOption=iOr(iOption,8)
      If (PCM) Then
         iOption=iOr(iOption,16)
         nPCM_Info=0
         Call Put_iScalar('PCM info length',nPCM_Info)
      End If
      iOption=iOr(iOption,32)
      If (lRF.and..not.PCM) iOption=iOr(iOption,2**7)
      Pseudo=.False.
      Do iCnttp = 1, nCnttp
         Pseudo = Pseudo .or. (pChrg(iCnttp) .and. Fixed(iCnttp))
      End Do
      If (lXF.or.Pseudo) Then
         iOption=iOr(iOption,2**7)
         iOption=iOr(iOption,2**8)
      End If
      If (VarT) iOption=iOr(iOption,2**7)
      If (VarR) iOption=iOr(iOption,2**8)
* 2el-integrals from the Cholesky vectors
      If (Cholesky.or.Do_RI) iOption=iOr(iOption,2**9)
*     RI-Option
      If (Do_RI) iOption=iOr(iOption,2**10)
*     1C-CD
      Cho_1Center=Get_Cho_1Center()
      If (Cholesky.and.Cho_1Center) iOption=iOr(iOption,2**12)
      Cho_OneCenter=Cho_1Center
      Call Put_iScalar('System BitSwitch',iOption)
      iter_S=0
      Call Put_iScalar('Saddle Iter',iter_S)
      iDNG = 0
      If (Do_Numerical_Gradients) iDNG=1
      Call Put_iScalar('DNG',iDNG)
*                                                                      *
************************************************************************
*                                                                      *
      Call DumpSagit()
      Call ClsSew()
      If (Allocated(AdCell)) Call mma_deallocate(AdCell)
      Call mma_deallocate(Coor_MPM)
      Call mma_deallocate(Chrg)
      Call mma_deallocate(Mass)
      Call mma_deallocate(Centr)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
      Logical Function Get_Cho_1Center()
#include "cholesky.fh"
      Get_Cho_1Center=Cho_1Center
      Return
      End
