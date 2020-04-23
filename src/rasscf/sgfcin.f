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
* Copyright (C) 1990, Markus P. Fuelscher                              *
************************************************************************

*>  @brief
*>    Generate the Fock-Matrix for frozen and inactive orbitals in
*>      the basis of active orbitals.
*>
*>  @author
*>    Markus P. Fuelscher
*>
*>  @details
*>  Generate the Fock-matrix for the frozen and inactive orbitals.
*>  Compute also the core energy and write to global variable EMY.
*>  Finally, transform the generated Fock-matrix
*>  into the basis of the active orbitals.
*>  Look into chapters 10.8.3 and 10.8.4 of \cite purple_book.
*>  The one body density matrices are required for e.g. reaction field
*>  or DFT calculations. In this case they are used to create a modified
*>  Fock Matrix.
*>
*>  @param[in] CMO The MO-coefficients
*>  @param[out] F The inactive Fock matrix in the basis of the active MO
*>  @param[inout] FI The inactive Fock matrix in AO-space
*>    \f[\sum_{\sigma\rho} D^I_{\sigma\rho}(g_{\mu\nu\sigma\rho} - \frac{1}{2} g_{\mu\sigma\rho\nu})\f]
*>    In output FI contains also the core energy added to
*>    the diagonal elements.
*>    \f[\sum_{\sigma\rho} D^I_{\sigma\rho}(g_{\mu\nu\sigma\rho} - \frac{1}{2} g_{\mu\sigma\rho\nu}) + \frac{E^{(0)}}{n_{el}} \delta_{\mu\nu} \f]
*>  @param[in] D1I The inactive one-body density matrix in AO-space
*>    \f[D^{\text{AO}, I} = 2 C (C^I)^\dagger \f]
*>    See ::get_D1I_rasscf.
*>  @param[in] D1A The active one-body density matrix in AO-space
*>    \f[ D^{\text{AO}, A} = C^A D^A (C^A)^\dagger \f]
*>    See ::get_D1A_rasscf.
*>  @param[in] D1S The active spin density matrix in AO-space
*>    \f[ D^{\text{AO}, A}_S = C^A (D^A_\alpha - D^A_\beta) (C^A)^\dagger \f]
      Subroutine SGFCIN(CMO, F, FI, D1I, D1A, D1S)
#ifdef _DMRG_
      use qcmaquis_interface_cfg
#endif
      use fcidump, only: DumpOnly
      use fciqmc, only: DoNECI
      use CC_CI_mod, only: Do_CC_CI

      use rasscf_data, only : EMY, KSDFT, dftfock, exfac, nac, nacpar,
     &    noneq, potnuc, rfpert,
     &    tot_charge, tot_el_charge, tot_nuc_charge,
     &    doBlockDMRG
      use general_data, only : iSpin, nActEl, nSym, nTot1,
     &    nBas, nIsh, nAsh, nFro


      implicit none
#include "rasdim.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='SGFCIN  ')
#include "WrkSpc.fh"
#include "rctfld.fh"
#include "pamint.fh"
#include "timers.fh"
#include "SysDef.fh"
*
      real*8, intent(in) :: CMO(*), D1I(*), D1A(*)
      real*8, intent(inout) :: FI(*), D1S(*), F(*)
      Character*8 Label
      Character*8 PAMlbl
      Logical First, Dff, Do_DFT, Found
      Logical Do_ESPF

#ifndef _DMRG_
      logical :: doDMRG = .false.
#endif
*
      Logical Do_OFemb, KEonly, OFE_first
      COMMON  / OFembed_L / Do_OFemb,KEonly,OFE_first
      Character*16 NamRfil
      COMMON  / OFembed_R / Rep_EN,Func_AB,Func_A,Func_B,Energy_NAD,
     &                      V_Nuc_AB,V_Nuc_BA,V_emb
      COMMON  / OFembed_I / ipFMaux, ip_NDSD, l_NDSD
*
      real*8, parameter ::  Zero=0.0d0, One=1.0d0
      real*8 :: CASDFT_Funct, dumm(1), Emyn, energy_nad, Eone,
     &  Erf1, Erf2, Erfhi, Erfx, Etwo, func_a, func_ab, func_b,
     &  potnuc_ref, rep_en, v_emb, v_nuc_ab, v_nuc_ba, dDot_
      integer :: i, iadd, ibas, icharge, iComp,
     &  ioff, iopt, ip_ndsd,
     &  ipam, ipfmaux, iprlev, iptmpfcki, ntmpfck,
     &  irc, iSyLbl,
     &  iSym, iTmp0, iTmp1, iTmp2, iTmp3, iTmp4, iTmp5, iTmp6, iTmp7,
     &  iTmp8, iTmpx, iTmpz, iTu, j, l_ndsd, lx0, lx1, lx2, lx3,
     &  mxna, mxnb, nAt, nst, nt, ntu, nu, nvxc

      Call qEnter(ROUTINE)
C Local print level (if any)
      IPRLEV=IPRLOC(3)
      IPRLEV=0000
      IF(IPRLEV.ge.DEBUG) THEN
         WRITE(LF,*)' Entering ',ROUTINE
      END IF

*
*     Generate molecular charges
      Call GetMem('Ovrlp','Allo','Real',iTmp0,nTot1+4)
      iRc=-1
      iOpt=2
      iComp=1
      iSyLbl=1
      Label='Mltpl  0'
      Call RdOne(iRc,iOpt,Label,iComp,Work(iTmp0),iSyLbl)
      Tot_Nuc_Charge=Work(iTmp0+nTot1+3)
      If ( iRc.ne.0 ) then
         Write(LF,*) 'SGFCIN: iRc from Call RdOne not 0'
         Write(LF,*) 'Label = ',Label
         Write(LF,*) 'iRc = ',iRc
         Call QTrace
         Call Abend
      Endif
      Call GetMem('Ovrlp','Free','Real',iTmp0,nTot1+4)
      Tot_El_Charge=Zero
      Do iSym=1,nSym
         Tot_El_Charge=Tot_El_Charge
     &                -2.0D0*DBLE(nFro(iSym)+nIsh(iSym))
      End Do
      Tot_El_Charge=Tot_El_Charge-DBLE(nActEl)
      Tot_Charge=Tot_Nuc_Charge+Tot_El_Charge
*
*     Load bare nuclei Hamiltonian
      Call GetMem('Fcore','Allo','Real',iTmp1,nTot1)
      iComp  =  1
      iSyLbl =  1
      iRc    = -1
      iOpt   =  6
      Label  = 'OneHam  '
      Call RdOne(iRc,iOpt,Label,iComp,Work(iTmp1),iSyLbl)
      If ( iRc.ne.0 ) then
         Write(LF,*) 'SGFCIN: iRc from Call RdOne not 0'
         Write(LF,*) 'Label = ',Label
         Write(LF,*) 'iRc = ',iRc
         Call QTrace
         Call Abend
      Endif
      If ( IPRLEV.ge.DEBUG ) then
       Write(LF,*)
       Write(LF,*) ' CMO in SGFCIN'
       Write(LF,*) ' ---------------------'
       Write(LF,*)
       ioff=1
       Do iSym = 1,nSym
         iBas = nBas(iSym)
         if(iBas.ne.0) then
           write(6,*) 'Sym =', iSym
           do i= 1,iBas
             write(6,*) (CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
           end do
           iOff = iOff + (iBas*iBas)
         end if
       End Do

       Write(LF,*)
       Write(LF,*) ' D1I in AO basis in SGFCIN'
       Write(LF,*) ' ---------------------'
       Write(LF,*)
       iOff=1
       Do iSym = 1,nSym
         iBas = nBas(iSym)
         call wrtmat(D1I(ioff),iBas,iBas, iBas, iBas)
         iOff = iOff + iBas*iBas
       End Do

       Write(LF,*)
       Write(LF,*) ' D1A in AO basis in SGFCIN'
       Write(LF,*) ' ---------------------'
       Write(LF,*)
       iOff=1
       Do iSym = 1,nSym
         iBas = nBas(iSym)
         call wrtmat(D1A(ioff),iBas,iBas, iBas, iBas)
         iOff = iOff + iBas*iBas
       End Do
      End If
      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' OneHam in AO basis in SGFCIN'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=0
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',Work(iTmp1+iOff),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
        End Do
      End If
*
*     Load the nuclear repulsion energy
*
*     Call Get_PotNuc(PotNuc)
      Call Get_dScalar('PotNuc',potNuc)
*
*     modify the one electron Hamiltonian for reaction
*     field calculations
      ERFX = Zero
      ERFhi = Zero
      iCharge=Int(Tot_Charge)
      Call DecideOnESPF(Do_ESPF)
      If ( Do_ESPF .or. lRF .or. KSDFT.ne.'SCF'
     &                      .or. Do_OFemb) Then
        Call GetMem('DtmpI','Allo','Real',iTmp3,nTot1)
        Call GetMem('DtmpA','Allo','Real',iTmp4,nTot1)
        Call GetMem('DtmpS','Allo','Real',iTmp7,nTot1)
*
*------ Generate total density
*
        Call Fold(nSym,nBas,D1I,Work(iTmp3))
        Call Fold(nSym,nBas,D1A,Work(iTmp4))
        Call Daxpy_(nTot1,1.0D0,Work(iTmp4),1,Work(iTmp3),1)
        Call Put_D1ao(Work(iTmp3),nTot1)
*        Write(LF,*)
*        Write(LF,*) ' D1ao in AO basis in SGFCIN'
*        Write(LF,*) ' ---------------------'
*        Write(LF,*)
*        iOff=0
*        Do iSym = 1,nSym
*          iBas = nBas(iSym)
*          Call TriPrt(' ','(5G17.11)',Work(iTmp3+iOff),iBas)
*          iOff = iOff + (iBas*iBas+iBas)/2
*        End Do
*
*------ Generate spin-density
*

*        do i=1,nBas(1)
*          write(*,*)i,"0000 D1S,Work(iTmp7)",D1S(i),Work(iTmp7+i-1)
*        end do

        Call Fold(nSym,nBas,D1S,Work(iTmp7))
        Call Put_D1Sao(Work(iTmp7),nTot1)

!        do i=1,nTot1
!          Work(iTmp7+i-1)=Zero
!          write(*,*)i,"1111 D1S,Work(iTmp7)",D1S(i),Work(iTmp7+i-1)
!        end do

*
*------ Scratch for one- and two-electron type contributions
*
        Call GetMem('htmp','Allo','Real',iTmp5,nTot1)
        Call GetMem('gtmp','Allo','Real',iTmp6,nTot1)
        Call dCopy_(nTot1,[Zero],0,Work(iTmp5),1)
        Call dCopy_(nTot1,[Zero],0,Work(iTmp6),1)
*
        First=.True.
        Dff=.False.
        Do_DFT=.True.

        Call Timing(Rado_1,Swatch,Swatch,Swatch)
        If(KSDFT(1:3).ne.'SCF' .or. Do_OFemb) Then
           Call Put_iArray('nFro',nFro,nSym)
           Call Put_iArray('nAsh',nAsh,nSym)
           Call Put_iArray('nIsh',nIsh,nSym)
        End If

        Call DrvXV(Work(iTmp5),Work(iTmp6),Work(iTmp3),
     &             PotNuc,nTot1,First,Dff,NonEq,lRF,
     &             KSDFT,ExFac,iCharge,iSpin,D1I,D1A,
     &             nTot1,DFTFOCK,Do_DFT)
      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' Work(iTmp5), h1 (DFT), in AO basis in SGFCIN'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=0
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',Work(iTmp5+iOff),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
        End Do
      End If
        Call Timing(Rado_2,Swatch,Swatch,Swatch)
        Rado_2 = Rado_2 - Rado_1
        Rado_3 = Rado_3 + Rado_2


        ERF1=Zero
        ERF2=dDot_(nTot1,Work(iTmp6),1,Work(iTmp4),1)
        ERFX= ERF1-0.5d0*ERF2
        Call Daxpy_(nTot1,1.0d0,Work(iTmp5),1,Work(iTmp1),1)

!        do i=1,nTot1
!          write(*,*)"   0000  FI",FI(i)
!        end do
        Call Daxpy_(nTot1,1.0d0,Work(iTmp6),1,FI,1)
!        do i=1,nTot
!          write(*,*)"   1111  FI",FI(i)
!        end do
*
*       Insert PAM integrals to one electron Hamiltonian
*
        If(KSDFT(1:3).eq.'PAM') Then
          Call GetMem('gtmp1','Allo','Real',iTmp8,nTot1)
          Do iPAM=1,nPAM
             Write(PAMlbl,'(A,I3.3)') 'PAM  ',ipPam(iPAM)
             Call dCopy_(nTot1,[Zero],0,Work(iTmp8),1)
          iComp=1
             Call RdOne(iRc,iOpt,PAMlbl,iComp,Work(iTmp8),iSyLbl)
             Call Daxpy_(nTot1,CPAM(iPAM),Work(iTmp8),1,Work(iTmp1),1)
          End Do
          Call GetMem('gtmp1','Free','Real',iTmp8,nTot1)
        End If
        Call GetMem('gtmp','Free','Real',iTmp6,nTot1)
        Call GetMem('htmp','Free','Real',iTmp5,nTot1)
        Call GetMem('DtmpS','Free','Real',iTmp7,nTot1)
        Call GetMem('DtmpA','Free','Real',iTmp4,nTot1)
        If(.not.Do_OFemb) Call GetMem('DtmpI','Free','Real',iTmp3,nTot1)
      End If
      If ( RFpert ) then
*
*       Read the reaction field from RunFile or RunOld
*
        Call f_Inquire('RUNOLD',Found)
        If (Found) Call NameRun('RUNOLD')
        Call GetMem('RCTFLD','Allo','Real',iTmpZ,nTot1)
        Call Get_dScalar('RF Self Energy',ERFX)
        Call Get_dArray('Reaction field',Work(iTmpZ),nTot1)
        Call Daxpy_(nTot1,1.0D0,Work(iTmpZ),1,Work(iTmp1),1)
        Call GetMem('RCTFLD','Free','Real',iTmpZ,nTot1)
        If (Found) Call NameRun('RUNFILE')
      End If
      Call GetMem('DoneI','Allo','Real',iTmp2,nTot1)
        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,*)
          Write(LF,*) ' D1I in AO basis in SGFCIN'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            Call TriPrt(' ','(5G17.11)',D1I(ioff),iBas)
            iOff = iOff + (iBas*iBas+iBas)/2
          End Do
        End If
      Call Fold(nSym,nBas,D1I,Work(iTmp2))

*
      If (Do_OFemb) Then
         If (OFE_first) Then
          Call GetMem('FMaux','Allo','Real',ipFMaux,nTot1)
          Call Coul_DMB(.true.,1,Rep_EN,Work(ipFMaux),Work(iTmp3),Dumm,
     &                         nTot1)
          OFE_first=.false.
         Else
          Call Coul_DMB(.false.,1,Rep_EN,Work(ipFMaux),Work(iTmp3),Dumm,
     &                          nTot1)
         EndIf
         Call DaXpY_(nTot1,One,Work(ipFMaux),1,Work(iTmp1),1)
*
         Call Get_NameRun(NamRfil) ! save the old RUNFILE name
         Call NameRun('AUXRFIL')   ! switch the RUNFILE name
         Call Get_dExcdRa(iTmpx,nVxc)
         Call DaXpY_(nTot1,One,Work(iTmpx),1,Work(iTmp1),1)
         If (nVxc.eq.2*nTot1) Then ! Nuc Attr added twice
            Call DaXpY_(nTot1,One,Work(iTmpx+nTot1),1,
     &                           Work(iTmp1),1)
            Call Get_dArray('Nuc Potential',Work(iTmpx),nTot1)
            Call DaXpY_(nTot1,-One,Work(iTmpx),1,Work(iTmp1),1)
         EndIf
         Call Free_Work(iTmpx)
         Call GetMem('DtmpI','Free','Real',iTmp3,nTot1)
         Call NameRun(NamRfil)   ! switch back to old RUNFILE
      End If
*
*     Compute energy contributions
      Eone = dDot_(nTot1,Work(iTmp2),1,Work(iTmp1),1)
*     Call Get_PotNuc(PotNuc_Ref)
      Call Get_dScalar('PotNuc',PotNuc_Ref)
      Eone = Eone + (PotNuc-PotNuc_Ref)
      Etwo = dDot_(nTot1,Work(iTmp2),1,FI,1)
      Call GetMem('DoneI','Free','Real',iTmp2,nTot1)
      EMY  = PotNuc_Ref+Eone+0.5d0*Etwo+ERFX
      CASDFT_Funct = Zero
      If(KSDFT(1:3).ne.'SCF'.and.KSDFT(1:3).ne.'PAM')
     &      Call Get_dScalar('CASDFT energy',CASDFT_Funct)
      If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*) ' Nuclear repulsion 1      :',PotNuc
         Write(LF,*) ' Nuclear repulsion 2 Ref  :',PotNuc_Ref
         Write(LF,*) ' One-electron core energy :',Eone
         Write(LF,*) ' Two-electron core energy :',Etwo
         Write(LF,*) ' Total core energy        :',EMY
         If(KSDFT(1:3).ne.'SCF'.and.KSDFT.ne.'PAM')
     &   Write(LF,*) ' CASDFT Energy            :',CASDFT_Funct
      End If

      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' FI matrix in AO in SGFCIN only 2-electron terms'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=1
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',FI(iOff),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
        End Do
      End If
      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' Work(iTmp1) matrix in SGFCIN, one-electron term'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=iTmp1
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',Work(iOff),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
        End Do
      End If
*
*     Assemble the one-electron (iTmp1) and two-electron contribution to AO Fock matrix
      Call DaXpY_(nTot1,One,Work(iTmp1),1,FI,1)
      Call GetMem('Fcore','Free','Real',iTmp1,nTot1)
*
      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' Inactive Fock matrix in AO basis in SGFCIN'
        Write(LF,*) '(it already contains OneHam and TwoEl contrib.)'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=1
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',FI(iOff),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
        End Do
      End If

*     Transform FI to active orbital basis and move it over to F.
*     Remove also the symmetry blocking.
      MXNB=0
      MXNA=0
      DO ISYM=1,NSYM
        MXNB=MAX(MXNB,NBAS(ISYM))
        MXNA=MAX(MXNA,NASH(ISYM))
      END DO
      CALL GETMEM('XXX0','ALLO','REAL',LX0,NTOT1)
      CALL GETMEM('XXX1','ALLO','REAL',LX1,NTOT1)
      CALL GETMEM('XXX2','ALLO','REAL',LX2,MXNB*MXNB)
      CALL GETMEM('XXX3','ALLO','REAL',LX3,MXNB*MXNA)
      CALL DCOPY_(NTOT1,FI,1,WORK(LX1),1)
      If(KSDFT(1:3).ne.'SCF'.and.KSDFT(1:3).ne.'PAM') Then
         Call Get_dExcdRa(ipTmpFckI,nTmpFck)
         CALL DaXpY_(NTOT1,1.0D0,Work(ipTmpFckI),1,WORK(LX1),1)
        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,*)
          Write(LF,*) ' Exchange correlation in AO basis in SGFCIN'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            Call TriPrt(' ','(5G17.11)',Work(ipTmpFckI+ioff-1),iBas)
            iOff = iOff + (iBas*iBas+iBas)/2
          End Do
        End If
         Call Free_Work(ipTmpFckI)
      End If
      If ( IPRLEV.ge.DEBUG ) then
       Write(LF,*)
       Write(LF,*) ' Modified FI in AO basis in SGFCIN'
       Write(LF,*) ' ---------------------'
       Write(LF,*)
       iOff=1
       Do iSym = 1,nSym
         iBas = nBas(iSym)
         Call TriPrt(' ','(5G17.11)',Work(LX1+ioff-1),iBas)
         iOff = iOff + (iBas*iBas+iBas)/2
       End Do
      End If
      CALL MOTRAC(CMO,WORK(LX1),WORK(LX2),WORK(LX3))
      CALL GETMEM('XXX3','FREE','REAL',LX3,MXNB*MXNA)
      CALL GETMEM('XXX2','FREE','REAL',LX2,MXNB*MXNB)
      CALL DCOPY_(NACPAR,[ZERO],0,F,1)
      NTU=0
      ITU=0
      IADD=0
Cbjp
      IF (NACTEL.NE.0) THEN
         EMYN=EMY/DBLE(NACTEL)
      ELSE
         EMYN=Zero
      ENDIF
      DO NST=1,NSYM
        NAT=NASH(NST)
        IF(NAT.NE.0) THEN
          DO NT=1,NAT
            NTU=NTU+IADD
            DO NU=1,NT
              NTU=NTU+1
              ITU=ITU+1
              F(NTU)=WORK(LX1-1+ITU)
              IF(NT.EQ.NU) F(NTU)=F(NTU)+EMYN
              WORK(LX0-1+ITU) = F(NTU)
            END DO
          END DO
          IADD=IADD+NAT
        ENDIF
      END DO

!Quan: Fix bug, skip Lucia stuff with DMRG
! and other external CI solvers.
      if (.not. any([DoNECI, Do_CC_CI, DumpOnly,
     &              doDMRG, doBlockDMRG])) then
        CALL CP_ONE_INT(WORK(LX0),ITU)
      endif

      CALL GETMEM('XXX1','FREE','REAL',LX1,NTOT1)
      CALL GETMEM('XXX0','FREE','REAL',LX0,NTOT1)

*     print h0
      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' Inactive Fock mat in act MO basis, h0, in SGFCIN'
        Write(LF,*) ' ------------'
        Write(LF,*)
        Call TriPrt(' ',' ',F,NAC)
      End If

      Call qExit('SGFCIN')


      Return
      End
