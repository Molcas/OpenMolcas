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
* Copyright (C) 2016,2017, Roland Lindh                                *
************************************************************************
      SubRoutine PrFin0(Dens,Dens_ab,nDT,EOrb,nEO,CMO,nCMO,KntE)
************************************************************************
*                                                                      *
*     purpose: Final printout                                          *
*                                                                      *
*     input:                                                           *
*       Dens    : the last total density                               *
*       EOrb    : orbital energies of length nEO                       *
*       CMO     : molecular orbital coefficients of length nCMO        *
*                                                                      *
*     called from: Final                                               *
*                                                                      *
*     calls to: R1IntB, RelEny, SetUp, Ortho                           *
*                                                                      *
*----------------------------------------------------------------------*
*  cloned from prfin (the part to call ones)                           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
*
      Real*8 Dens(nDT),Dens_ab(nDT), EOrb(nEO),CMO(nCMO), KntE(nDT)
      Integer iSpn, iMult
*
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "ldfscf.fh"
#include "file.fh"
#include "rctfld.fh"
#include "oneswi.fh"
#ifdef _FDE_
#include "embpotdata.fh"
#endif
#include "scfwfn.fh"
#include "ksdft.fh"

      Logical Do_OFemb,KEonly,OFE_first
      COMMON  / OFembed_L / Do_OFemb,KEonly,OFE_first
      Logical Do_Tw
      COMMON  / Tw_corr_L   / Do_Tw
      Character*16  ADDC_KSDFT
      COMMON  / ADDcorr_C   / ADDC_KSDFT
      Logical Do_Addc
      COMMON  / ADDcorr_L   / Do_Addc
      COMMON  / ADDcorr_R   / DE_KSDFT_c
      Real*8 Erest_xc
      COMMON /dCSCF_xc/ Erest_xc
      Real*8 s2CNO
      COMMON /dCSCF_s2/ s2CNO
      Logical Do_SpinAV
      COMMON  / SPAVE_L  / Do_SpinAV

      Integer  Cho_X_GetTol
      External Cho_X_GetTol
*
*---- Define local variables
      Character*80 Lines(6), Fmt*60
      Logical Reduce_Prt
      External Reduce_Prt
      Dimension Dumm1(1)


#include "SysDef.fh"

*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
      jPrint = iPrint
      iPL=iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=0
      If (iPL.le.1) jPrint=1
*
*----------------------------------------------------------------------*
*
*---- Calculate kinetic energy
      If (iUHF.eq.1) then
         do i=1,nDT
            Dens(i)=Dens(i)+Dens_ab(i)
         end do
      End if
      EKin = DDot_(nBT,KntE,1,Dens,1)
*
*---- Print out header to final results
      Lines(1) = 'SCF/KS-DFT Program, Final results'
      Lines(2) = ' '
      Lines(3)=' '
c      if(.false.) Lines(3) = Molcas_revision
      Lines(4) = ' '
      Lines(5) = 'Final Results'
      If (jPrint.ge.2) Then
         Call Banner(Lines,5,lPaper-7)
         Write(6,*)
      End If
*
*---- Compute energies from aufbau density matrix
      If (Abs(EKin).gt.1.0D-6) Then
         Virial = - EneV/EKin
      Else
         Virial = Zero
      End If
*
*---- Compute estimate of FNO correlation energy from Delta_Tw
      If (Do_Tw) Call Tw_corr_drv(EOrb,nEO,CMO,nCMO,Ecorr)

*
*---- Print out final results
      If(WarnCfg) Then
         call WarningMessage(1,
     &     'Warning:; The program may have converged to a solution;'//
     &     'that does not correspond to the lowest energy!')
      End If
      If(WarnPocc) Then
         call WarningMessage(1,
     &     'Warning:; The program may have converged to a solution;'//
     &     'with partial occupation numbers!')
      End If
      If(WarnSlow) Then
         call WarningMessage(1,
     &     'Warning:; The program had convergence problems;'//
     &     'and terminated with looser convergence')
      End If
      Fmt = '(6X,A,T50,F19.10)'
      suhf=-0.5d0+sqrt(0.25d0+s2uhf)
      Call put_dscalar('UHFSPIN',SUHF)
      iTol = min(Cho_X_GetTol(8),8)
      If (doLDF) Then
         ! the non-robust integral representations are numerically
         ! unstable and may give much larger errors, depending on
         ! machine and compiler options. Be more tolerant in those
         ! cases.
         If (LDF_IntegralMode.ne.1) Then
            If (.not.LDF_UseConventionalIntegrals) Then
               iTol=max(iTol-4,2)
            End If
         End If
      End If
      If (jPrint.ge.2) Then
         If (MxConstr.gt.0) Then
            DE_KSDFT_c=0.0d0
            If (Do_Addc) Then
               Call SetUp_iSD()
               Call Get_DEcorr(nBT,Dumm1,iDumm,'SCF ')
               Call Free_iSD()
            End If
            ECNO=EneV+E_nondyn+DE_KSDFT_c
            If (KSDFT.ne.'SCF') ECNO=ECNO+Erest_xc
            call PrintResult(6,FMT, 'Total energy',0,' ',[ECNO],1)
            call PrintResult(6,FMT, 'Nondynamical correlation energy',
     &                       0,' ',[E_nondyn],1)
            If (KSDFT.ne.'SCF') Then
               call PrintResult(6,FMT, 'Energy-restoring term',
     &                          0,' ',[Erest_xc],1)
            EndIf
            If (Do_Addc) Then
              Call PrintResult(6,FMT,
     &        'Added correlation energy ('//ADDC_KSDFT(1:4)//') ',
     &         0,' ',[DE_KSDFT_c],1)
            EndIf
            Call Add_Info('E_CNO',[ECNO],1,iTol)
         EndIf
         If (Do_Tw) Then
            E_Tw=EneV+Ecorr
            call PrintResult(6,FMT, 'Total energy',0,' ',[E_Tw],1)
            call PrintResult(6,FMT, 'Delta_Tw correlation energy',
     &                       0,' ',[Ecorr],1)
            Call Add_Info('E_Tw',[E_Tw],1,iTol)
         EndIf
         If (KSDFT.eq.'SCF') Then
            call PrintResult(6,FMT, 'Total SCF energy',0,' ',[EneV],1)
c            Write(6,Fmt)'Total SCF energy',EneV
         Else
            call PrintResult(6,FMT, 'Total KS-DFT energy',0,' ',
     &          [EneV],1)
c            Write(6,Fmt)'Total KS-DFT energy',EneV
         End If
         Write(6,Fmt)'One-electron energy',E1V
#ifdef _FDE_
         ! Embedding
         if (embPot) then
          Write(6,Fmt)'E from embedding potential(<Psi|v_emb|Psi>)',Eemb
         end if
#endif
         Write(6,Fmt)'Two-electron energy',E2V
         Write(6,Fmt)'Nuclear repulsion energy',PotNuc
         Write(6,Fmt)'Kinetic energy (interpolated)',EKin
         Write(6,Fmt)'Virial theorem',Virial
         If (.not.Do_SpinAV) Then
            Write(6,Fmt)'Total spin, S(S+1)',s2uhf
            Write(6,Fmt)'Total spin, S',suhf
         EndIf
         If (MxConstr.gt.0) Write(6,Fmt)'Spin deviation',s2uhf-s2CNO
      End If
      iSpn = INT(suhf+0.5d0)
      iMult = 2*iSpn+1
      Call Put_iScalar('Multiplicity',iMult)
      Call Add_Info('E_SCF',[EneV],1,iTol)
#ifdef _HDF5_
      call mh5_put_dset(wfn_energy,EneV)
#endif
*
      If (nIter(nIterP).gt.0.and.jPrint.ge.2) Then
         Write(6,Fmt)'Max non-diagonal density matrix element',DMOMax
         Write(6,Fmt)'Max non-diagonal Fock matrix element',FMOMax
      End If
      if (CoefX.ne.1.0.or.CoefR.ne.1.0) Then
         Write(6,Fmt)'Exchange scaling factor',CoefX
         Write(6,Fmt)'Correlation scaling factor',CoefR
      End If
      If (jPrint.ge.2) Write(6,*)
c      If (jPrint.ge.2 .and. Do_OFemb) Call OFE_print(EneV)
      If (Do_OFemb) Call OFE_print(EneV)
*
* xml tagging
*
      If(KSDFT.eq.'SCF') Then
         Call xml_dDump('energy','Total SCF energy','a.u.',1,[EneV],1,1)
      Else
         Call xml_dDump('energy','Total KS-DFT energy','a.u.',1,
     &          [EneV],1,1)
      End If
      Call xml_dDump('kinetic','Kinetic energy','a.u.',2,[Ekin],1,1)
      Call xml_dDump('virial','Virial coefficient','a.u.',2,
     &         [Virial],1,1)
      Call xml_dDump('spin','UHF spin','',1,[suhf],1,1)
      Call xml_dDump('potnuc','Nuclear repulsion energy','a.u.',
     &         1,[potnuc],1,1)
      Call xml_dDump('energy1el','One electron energy','a.u.',1,
     &         [E1V],1,1)
      Call xml_dDump('energy2el','Two electron energy','a.u.',1,
     &         [E2V],1,1)
      Call xml_iDump('nsym','Number of irreps','',1,[nSym],1,1)
      Call xml_iDump('nbas','Number of basis functions','',
     &         1,nBas,nSym,1)
      Call xml_iDump('norb','Number of orbitals','',1,nOrb,nSym,1)
      If(iUHF.eq.0) Then
         Call xml_iDump('nocc','Number of occupied orbitals','',
     &                  1,nOcc(1,1),nSym,1)
      Else
         Call xml_iDump('nocc_a','Number of occupied alpha orbitals','',
     &                  1,nOcc(1,1),nSym,1)
         Call xml_iDump('nocc_b','Number of occupied beta orbitals','',
     &                  1,nOcc(1,2),nSym,1)
      End If
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
