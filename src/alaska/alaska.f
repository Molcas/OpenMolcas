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
* Copyright (C) 1989-1992, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      Subroutine Alaska(LuSpool,ireturn)
************************************************************************
*                                                                      *
*  Object: Driver for the one and two electron integral gradient       *
*          program ALASKA.                                             *
*                                                                      *
*          Alaska is a derivative code of Seward 3.1.                  *
*                                                                      *
* Called from: None                                                    *
*                                                                      *
* Calling    : QEnter                                                  *
*              XuFlow (IBM)                                            *
*              SetUp0                                                  *
*              GetMem                                                  *
*              GetInf                                                  *
*              Inputg                                                  *
*              DrvN1                                                   *
*              Drvh1                                                   *
*              PrepP                                                   *
*              Drvg1                                                   *
*              CloseP                                                  *
*                                                                      *
*  Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA     *
*          July '89 - May '90                                          *
*                                                                      *
*          Roland Lindh, Dept. of Theoretical Chemistry, University of *
*          Lund, SWEDEN. Modified to gradient calculations September   *
*          1991 - February 1992.                                       *
************************************************************************
      use Real_Spherical
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
      External RF_On
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"
#include "rctfld.fh"
#include "nsd.fh"
#include "setup.fh"
#include "columbus_gamma.fh"
#include "nac.fh"
#include "alaska_root.fh"
#include "para_info.fh"
      Logical OldTst, DoRys, RF_On, Found
      Logical Do_OFemb,KEonly,OFE_first
      Character(Len=180) Label
      COMMON  / OFembed_L / Do_OFemb,KEonly,OFE_first
************ columbus interface ****************************************
        Integer  Columbus, colgradmode
        Integer lcartgrd
        Real*8 Cgrad(3,MxAtom)
        Character CNames(MxAtom)*(LENIN5),Lab*80
        Integer iatom,icen,j
*                                                                      *
************************************************************************
      npelem=3
*                                                                      *
      Call CWTime(TCpu1,TWall1)
*                                                                      *
************************************************************************
*                                                                      *
*     Prologue
*                                                                      *
************************************************************************
*                                                                      *
      iRout=1
*                                                                      *
************************************************************************
*                                                                      *
*     Print program header
*
      Call qEnter('Alaska')
*                                                                      *
************************************************************************
*                                                                      *
*     Get the input information as Seward dumped on INFO.
*
      nDiff=1
      DoRys=.True.
      Call IniSew(Info,DoRys,nDiff)
      If (RF_On()) Then
         If (NonEq_Ref) Then
            Call WarningMessage(2,'Error in Alaska')
            Write (6,*) 'NonEq=.True., invalid option'
            Call Abend()
         End If
         Call Init_RctFld(.False.,iCharge_Ref)
      End If
*
      OldTst = Test
*                                                                      *
************************************************************************
*                                                                      *
*     Input specific for the gradient calculation.
*
      Call Inputg(LuSpool)
*
*---- Since the input has changed some of the shell information
*     regenerate the tabulated shell information.
*
      iPrint=nPrint(iRout)

      If (iPrint.ge.99) Call GetMem(' LIST ','LIST','REAL',iDum,iDum)
      Call GetMem('Grad','Allo','Real',ipGrad,lDisp(0))
      Call GetMem('Temp','Allo','Real',ipTemp,lDisp(0))
      call dcopy_(lDisp(0),[Zero],0,Work(ipGrad),1)
*
*     remove LuSpool
*
      Call Close_LuSpool(LuSpool)
*
*      identify a Columbus calculation
*      Columbus=1
*      colgradmode=0   standard gradient written to GRAD
*      colgradmode=1   standard gradient written to Grad State1
*      colgradmode=2   standard gradient written to Grad State1
*      colgradmode=3   non-adiabatic coupling vector written to NADC
*
       Call Get_iScalar('Columbus',Columbus)
       Call Get_iScalar('colgradmode',colgradmode)
*
*-----Start computing the gradients
*                                                                      *
************************************************************************
*                                                                      *
*     Compute nuclear contributions.
*
      If (king().or.HF_Force) Then
*
*    per default NADC must not have nuclear contributions added
*
       If (NO_NUC .or. (Columbus.eq.1 .and. colgradmode.eq.3)) then
         write(6,*) 'Skipping Nuclear Charge Contribution'
       Else
         Call DrvN1(Work(ipGrad),Work(ipTemp),lDisp(0))
         If (iPrint.ge.15) Then
            Lab=' Total Nuclear Contribution'
            Call PrGrad(Lab,Work(ipGrad),lDisp(0),lIrrep,ChDisp,iPrint)
         End If
       End If
      End If

!      iPrint=16
*
************************************************************************
*                                                                      *
      If (Do_OFemb) Then
* RepNuc term to the Orbital-Free Embedding gradient
         Call DrvN1_EMB(Work(ipGrad),Work(ipTemp),lDisp(0))
      EndIf
*                                                                      *
************************************************************************
*                                                                      *
      If (Test) Go To 999
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute contribution due to the derivative of the one-
*     electron hamiltonian and the contribution due the re-
*     normalization.
*
      If (.Not.(nProcs.gt.1).or.Onenly) Then
         Call Drvh1(Work(ipGrad),Work(ipTemp),lDisp(0))
      End If ! .Not.(nProcs.gt.1).or.Onenly
      If (Do_OFemb) Then
* NucAtt term to the Orbital-Free Embedding gradient
         Call Drvh1_EMB(Work(ipGrad),Work(ipTemp),lDisp(0))
      EndIf
!         Lab='Nuc + One-electron Contribution'
!         Call PrGrad(Lab,Work(ipGrad),lDisp(0),lIrrep,ChDisp,iPrint)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the DFT contribution to the gradient
*
      Call DrvDFTg(Work(ipGrad),Work(ipTemp),lDisp(0))
*                                                                      *
************************************************************************
*                                                                      *
      If (Onenly) Go To 998
*                                                                      *
************************************************************************
*                                                                      *
*     DFT-type Orbital-Free Embedding term to the gradient
*
      If (Do_OFemb) Call DrvEMBg(Work(ipGrad),Work(ipTemp),lDisp(0))
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute contribution due to 2-electron integrals.
*
      Call GetMem('MemHid','ALLO','REAL',idum,MemHid)
*
*                                                                      *
************************************************************************
*                                                                      *
      If (Cholesky.or.Do_RI) Then
         If (Cholesky) Then
            If (iPrint.ge.6) Write (6,*) 'Cholesky-ERI gradients!'
         Else
            If (iPrint.ge.6) Write (6,*) 'RI-ERI gradients!'
         End If
         Call Drvg1_RI (Work(ipGrad),Work(ipTemp),lDisp(0))
      Else
         If (iPrint.ge.6) Write (6,*) 'Conventional ERI gradients!'
         Call Drvg1    (Work(ipGrad),Work(ipTemp),lDisp(0))
      End If
*
      Call DScal_(lDisp(0),Half,Work(ipTemp),1)
      If (iPrint.ge.15) Then
         Lab=' Two-electron Contribution'
         Call PrGrad(Lab,Work(ipTemp),lDisp(0),lIrrep,ChDisp,iPrint)
      End If
*
*-----Accumulate contribution to the gradient
*
      Call GR_DArray(Work(ipGrad),lDisp(0))
      Call DaXpY_(lDisp(0),One,Work(ipTemp),1,
     &                        Work(ipGrad),1)
*
*                                                                      *
************************************************************************
*                                                                      *
*
      Call GetMem('MemHid','Free','REAL',idum,MemHid)
 998  Continue
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(TCpu2,TWall2)
      Call SavTim(7,TCpu2-TCpu1,TWall2-TWall1)
*                                                                      *
************************************************************************
*                                                                      *
*-----Apply the translational and rotational invariance of
*     the energy.
*
      If (TRSymm) Then
         If (iPrint.ge.99) Then
            Call PrGrad(
     &       ' Molecular gradients (no TR) ',
     &           Work(ipGrad),lDisp(0),lIrrep,ChDisp,iPrint)
            Call RecPrt(' The A matrix',' ',Work(ipAm),
     &                  lDisp(0),lDisp(0))
         End If
         call dcopy_(lDisp(0),Work(ipGrad),1,Work(ipTemp),1)

         Call dGeMV_('N',lDisp(0),lDisp(0),
     &              One,Work(ipAm),lDisp(0),
     &              Work(ipTemp),1,
     &              Zero,Work(ipGrad),1)
         Call GetMem('AMtrx','Free','Real',ipAm,lDisp(0)**2)
      End If ! TRSymm
*                                                                      *
************************************************************************
*                                                                      *
*-----Equivalence option
*
      If (lEq) Then
         Do i = 1, lDisp(0)
            If (IndxEq(i).ne.i) Work(ipGrad-1+i) =
     &         Work(ipGrad-1+IndxEq(i))
         End Do
      End If ! lEq
*                                                                      *
************************************************************************
*                                                                      *
 999  Continue
*                                                                      *
************************************************************************
*                                                                      *
      nCnttp_Valence=0
      Do iCnttp = 1, nCnttp
         If (AuxCnttp(iCnttp)) Go To 1999
         nCnttp_Valence = nCnttp_Valence+1
      End Do
 1999 Continue
*
*     f^AB is the "total derivative coupling"
*     h^AB is the "CI derivative coupling"
*     f^AB = <B|dA/dR> = h^AB/(E_A-E_B) + f_CSF^AB = -f^BA
*     h^AB = <B|dH/dR|A> = h^BA
*     f_CSF^AB = -f_CSF^BA
*
*     Note that we store h^AB + f_CSF^AB*(E_A-E-B), or just h^AB if
*     NOCSF was given, to avoid division by (nearly) zero
*
      If (isNAC) Then
        Call PrGrad('CI derivative coupling ',
     &                 Work(ipGrad),lDisp(0),lIrrep,ChDisp,iPrint)
        EDiff_s = Max(One, Ten**(-Floor(Log10(Abs(EDiff)))-4))
        EDiff_f = EDiff*EDiff_s
        If (DoCSF) Then
          Call Allocate_Work(ipCSFG,lDisp(0))
          Call CSFGrad(Work(ipCSFG),lDisp(0))
          Call PrGrad('CSF derivative coupling ',
     &                   Work(ipCSFG),lDisp(0),lIrrep,ChDisp,iPrint)
          Call daxpy_(lDisp(0),EDiff_f,Work(ipCSFG),1,Work(ipGrad),1)
          Call Free_Work(ipCSFG)
        End If
        write(6,'(15X,A,ES13.6)') 'Energy difference: ',EDiff
        Label = ''
        If (EDiff_s.gt.One)
     &      Write(Label,'(A,ES8.1,A)') ' (divided by',EDiff_s,')'
        Label = 'Total derivative coupling'//Trim(Label)
        Call Allocate_Work(ipTmp,lDisp(0))
        call dcopy_(lDisp(0),Work(ipGrad),1,Work(ipTmp),1)
        call dscal_(lDisp(0),One/EDiff_f,Work(ipTmp),1)
        Call PrGrad(Trim(Label),
     &              Work(ipTmp),lDisp(0),lIrrep,ChDisp,iPrint)
        write(6,'(15X,A,F12.4)') 'norm: ',dnrm2_(lDisp(0),Work(ipTmp),1)
        Call Free_Work(ipTmp)
      ElseIf (iPrint.ge.4) then
         If (HF_Force) Then
            Call PrGrad('Hellmann-Feynman Forces ',
     &                 Work(ipGrad),lDisp(0),lIrrep,ChDisp,iPrint)
         Else
            Call PrGrad(' Molecular gradients',
     &                    Work(ipGrad),lDisp(0),lIrrep,ChDisp,iPrint)
         End If
      End If
      If (isNAC) Then
*        For NAC, the sign is undefined (because the wave functions can change sign),
*        check only absolute values
         Call Allocate_Work(ipTmp,lDisp(0))
         Do i=0,lDisp(0)-1
            Work(ipTmp+i)=Abs(Work(ipGrad+i))
         End Do
         Call Add_Info('Grad',Work(ipTmp),lDisp(0),6)
         Call Free_Work(ipTmp)
      Else
         Call Add_Info('Grad',Work(ipGrad),lDisp(0),6)
      End If
*
*---- Molcas format
*
*-----Write gradient to runfile.
*
*     Save the gradient
*

      Call Get_iScalar('Unique atoms',nsAtom)
      l1 = 3*nsAtom
      Call GetMem('RELAX','Allo','Real',ipRlx,l1)
      mdc = 0
      ndc = 0
      Do iCnttp = 1, nCnttp_Valence
*
*        Skip gradients for pseudo atoms
*
         If (pChrg(iCnttp).or.dbsc(iCnttp)%nFragType.gt.0.or.
     &       FragCnttp(iCnttp)) Then
            mdc=mdc+dbsc(iCnttp)%nCntr
         Else
            Do iCnt = 1, dbsc(iCnttp)%nCntr
               mdc=mdc+1
               ndc=ndc+1
               Do iCar = 1, 3
                  If (InxDsp(mdc,iCar).ne.0) Then
                     temp = Work(ipGrad-1+InxDsp(mdc,iCar))
                     Work(ipRlx-1+(ndc-1)*3+iCar) = temp
                  Else
*
*                    Put in explicit zero if gradient is zero
*                    by symmetry.
*
                     Work(ipRlx-1+(ndc-1)*3+iCar) = Zero
                  End If
               End Do
            End Do
         End If
      End Do
*
      If (HF_Force) Then
         Call Put_dArray('HF-forces',Work(ipRlx),l1)
      Elseif (Columbus.eq.1) then
         Call Put_nadc(colgradmode,Work(ipRlx),l1)
      Else
         Call Put_Grad(Work(ipRlx),l1)
      End If
      Call GetMem('RELAX','Free','Real',ipRlx,l1)


************ columbus interface ****************************************
* print full cartesian gradient in Columbus format
*

      if (Columbus.eq.1) then
*     Real*8 Cgrad(3,mxatom)
*     Character CNames(MxAtom)*9
*     Integer lcartgrd, iatom,icen,j
      Call TrGrd_Alaska_(CGrad,CNames,Work(ipGrad),lDisp(0),iCen)
      lcartgrd=60
      lcartgrd=isFreeUnit(lcartgrd)
      Call Molcas_Open(lcartgrd,'cartgrd')
      DO 300 IATOM = 1,iCen
        write (60,1010) (CGrad(j,iatom), j=1,3)
  300 CONTINUE
      close(lcartgrd)
 1010 format (3d15.6)
      endif

*
*-----At the end of the calculation free all memory to check for
*     corruption of the memory.
*
      Call GetMem('Temp','Free','Real',ipTemp,lDisp(0))
      Call GetMem('Grad','Free','Real',ipGrad,lDisp(0))
*
      If (iPrint.ge.99) Then
         Call GetMem(' LIST ','LIST','REAL',iDum,iDum)
      End If
*
*     Restore iRlxRoot if changed as set by the RASSCF module.
*
      Call qpg_iScalar('Relax CASSCF root',Found)
      If (Found) Then
         Call Get_iScalar('Relax CASSCF root',irlxroot1)
         Call qpg_iScalar('Relax Original root',Found)
         If (Found) Then
            Call Get_iScalar('Relax Original root',irlxroot2)
            If (iRlxRoot1.ne.iRlxRoot2) Then
               Call Put_iScalar('Relax CASSCF root',irlxroot2)
               Call Put_iScalar('NumGradRoot',irlxroot2)
            End If
         End If
      End If
*
*     Epilogue
*
      Call ClsSew
      Call qExit('Alaska')
*
      If (iPrint.ge.6) Then
         Call qStat(' ')
         Call FastIO('STATUS')
      End If
*
      If (Test) Then
         ireturn=20
      Else
         ireturn=0
      End If
      Return
      End
