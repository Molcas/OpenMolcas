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
*>  @param[in,out] FI The inactive Fock matrix in AO-space
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
      use RunFile_procedures, only: Get_dExcdRa
#ifdef _DMRG_
      use qcmaquis_interface_cfg
#endif
      use fcidump, only: DumpOnly
      use fciqmc, only: DoNECI
      use CC_CI_mod, only: Do_CC_CI

      use timers, only: TimeDens
      use lucia_data, only: INT1, INT1O
      use rasscf_global, only : EMY, KSDFT, dftfock, exfac, nac, nacpar,
     &    noneq, potnuc, rfpert,
     &    tot_charge, tot_el_charge, tot_nuc_charge,
     &    doBlockDMRG, doDMRG
      use OneDat, only: sNoNuc, sNoOri
      use general_data, only : iSpin, nActEl, nSym, nTot1,
     &    nBas, nIsh, nAsh, nFro
      use OFEmbed, only: Do_OFemb, OFE_first, FMaux
      use OFEmbed, only: Rep_EN
      use rctfld_module, only:  lRF
      use stdalloc, only: mma_allocate, mma_deallocate
      use PrintLevel, only: DEBUG
      use output_ras, only: LF,IPRLOC

      implicit none
      Character(LEN=16), Parameter :: ROUTINE='SGFCIN  '
*
      real*8, intent(in) :: CMO(*), D1I(*), D1A(*)
      real*8, intent(inout) :: FI(*), D1S(*), F(*)
      Character(LEN=8) Label
      Logical First, Dff, Do_DFT, Found
      Logical Do_ESPF
*
      real*8, parameter ::  Zero=0.0d0, One=1.0d0
      real*8 :: CASDFT_Funct, dum1, dum2, dum3, dumm(1), Emyn, Eone,
     &  Erf1, Erf2, Erfx, Etwo,  potnuc_ref, Time(2)
      integer :: i, iadd, ibas, icharge, iComp,
     &  ioff, iopt, iprlev, ntmpfck,
     &  irc, iSyLbl, iSym, iTu, j,
     &  mxna, mxnb, nAt, nst, nt, ntu, nu, nvxc
      real*8, allocatable :: TmpFckI(:), Tmpx(:)
      real*8, allocatable:: Tmp0(:), Tmp1(:), Tmp2(:), Tmp3(:),
     &                      Tmp4(:), Tmp5(:), Tmp6(:), Tmp7(:),
     &                      Tmpz(:), X0(:), X1(:), X2(:),
     &                      X3(:)
      real*8, external :: dDot_

C Local print level (if any)
      IPRLEV=IPRLOC(3)
      IPRLEV=0000
      IF(IPRLEV.ge.DEBUG) THEN
         WRITE(LF,*)' Entering ',ROUTINE
      END IF

*
*     Generate molecular charges
      Call mma_allocate(Tmp0,nTot1+4,Label='Tmp0')
      iRc=-1
      iOpt=ibset(0,sNoOri)
      iComp=1
      iSyLbl=1
      Label='Mltpl  0'
      Call RdOne(iRc,iOpt,Label,iComp,Tmp0,iSyLbl)
      Tot_Nuc_Charge=Tmp0(nTot1+4)
      If ( iRc.ne.0 ) then
         Write(LF,*) 'SGFCIN: iRc from Call RdOne not 0'
         Write(LF,*) 'Label = ',Label
         Write(LF,*) 'iRc = ',iRc
         Call Abend
      Endif
      Call mma_deallocate(Tmp0)
      Tot_El_Charge=Zero
      Do iSym=1,nSym
         Tot_El_Charge=Tot_El_Charge
     &                -2.0D0*DBLE(nFro(iSym)+nIsh(iSym))
      End Do
      Tot_El_Charge=Tot_El_Charge-DBLE(nActEl)
      Tot_Charge=Tot_Nuc_Charge+Tot_El_Charge
*
*     Load bare nuclei Hamiltonian
      Call mma_allocate(Tmp1,nTot1,Label='Tmp1')
      iComp  =  1
      iSyLbl =  1
      iRc    = -1
      iOpt   =  ibset(ibset(0,sNoOri),sNoNuc)
      Label  = 'OneHam  '
      Call RdOne(iRc,iOpt,Label,iComp,Tmp1,iSyLbl)
      If ( iRc.ne.0 ) then
         Write(LF,*) 'SGFCIN: iRc from Call RdOne not 0'
         Write(LF,*) 'Label = ',Label
         Write(LF,*) 'iRc = ',iRc
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
        iOff=1
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',Tmp1(iOff),iBas)
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
      iCharge=Int(Tot_Charge)

      Call mma_allocate(Tmp3,nTot1,Label='Tmp3')
      Call mma_allocate(Tmp4,nTot1,Label='Tmp4')
      Call mma_allocate(Tmp7,nTot1,Label='Tmp7')

      Call DecideOnESPF(Do_ESPF)
*
*---- Generate total density
*
      Call Fold(nSym,nBas,D1I,Tmp3)
      Call Fold(nSym,nBas,D1A,Tmp4)
      Call Daxpy_(nTot1,1.0D0,Tmp4,1,Tmp3,1)
      Call Put_dArray('D1ao',Tmp3,nTot1)
*     Write(LF,*)
*     Write(LF,*) ' D1ao in AO basis in SGFCIN'
*     Write(LF,*) ' ---------------------'
*     Write(LF,*)
*     iOff=1
*     Do iSym = 1,nSym
*       iBas = nBas(iSym)
*       Call TriPrt(' ','(5G17.11)',Tmp3(iOff),iBas)
*       iOff = iOff + (iBas*iBas+iBas)/2
*     End Do
*
*---- Generate spin-density
*
      Call Fold(nSym,nBas,D1S,Tmp7)
      Call Put_dArray('D1sao',Tmp7,nTot1)

      If(KSDFT(1:3).ne.'SCF' .or. Do_OFemb) Then
        Call Put_iArray('nFro',nFro,nSym)
        Call Put_iArray('nAsh',nAsh,nSym)
        Call Put_iArray('nIsh',nIsh,nSym)
      End If

      If ( Do_ESPF .or. lRF .or. KSDFT.ne.'SCF'
     &                      .or. Do_OFemb ) Then
*
*------ Scratch for one- and two-electron type contributions
*
        Call mma_allocate(Tmp5,nTot1,Label='Tmp5')
        Tmp5(:)=0.0D0
        Call mma_allocate(Tmp6,nTot1,Label='Tmp6')
        Tmp6(:)=0.0D0
*
        First=.True.
        Dff=.False.
        Do_DFT=.True.

        Call Timing(Time(1),dum1,dum2,dum3)

        Call DrvXV(Tmp5,Tmp6,Tmp3,
     &             PotNuc,nTot1,First,Dff,NonEq,lRF,
     &             KSDFT,ExFac,iCharge,iSpin,
     &             DFTFOCK,Do_DFT)
        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,*)
          Write(LF,*) ' Tmp5, h1 (DFT), in AO basis in SGFCIN'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            Call TriPrt(' ','(5G17.11)',Tmp5(iOff),iBas)
            iOff = iOff + (iBas*iBas+iBas)/2
          End Do
        End If
        Call Timing(Time(2),dum1,dum2,dum3)
        TimeDens = TimeDens + Time(2) - Time(1)


        ERF1=Zero
        ERF2=dDot_(nTot1,Tmp6,1,Tmp4,1)
        ERFX= ERF1-0.5d0*ERF2
        Call Daxpy_(nTot1,1.0d0,Tmp5,1,Tmp1,1)

        Call Daxpy_(nTot1,1.0d0,Tmp6,1,FI,1)
*
        Call mma_deallocate(Tmp6)
        Call mma_deallocate(Tmp5)
      End If

      Call mma_deallocate(Tmp7)
      Call mma_deallocate(Tmp4)
      If(.not.Do_OFemb) Call mma_deallocate(Tmp3)
*
      If ( RFpert ) then
*
*       Read the reaction field from RunFile or RunOld
*
        Call f_Inquire('RUNOLD',Found)
        If (Found) Call NameRun('RUNOLD')
        Call mma_allocate(TmpZ,nTot1,Label='TmpZ')
        Call Get_dScalar('RF Self Energy',ERFX)
        Call Get_dArray('Reaction field',TmpZ,nTot1)
        Call Daxpy_(nTot1,1.0D0,TmpZ,1,Tmp1,1)
        Call mma_deallocate(TmpZ)
        If (Found) Call NameRun('#Pop')
      End If
      Call mma_allocate(Tmp2,nTot1,Label='Tmp2')
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
      Call Fold(nSym,nBas,D1I,Tmp2)

*
      If (Do_OFemb) Then
         If (OFE_first) Then
          Call mma_allocate(FMaux,nTot1,Label='FMAux')
          Call Coul_DMB(.true.,1,Rep_EN,FMaux,Tmp3,Dumm,nTot1)
          OFE_first=.false.
         Else
          Call Coul_DMB(.false.,1,Rep_EN,FMaux,Tmp3,Dumm,nTot1)
         EndIf
         Call DaXpY_(nTot1,One,FMaux,1,Tmp1,1)
*
         Call NameRun('AUXRFIL') ! switch the RUNFILE name
         Call Get_dExcdRa(Tmpx,nVxc)
         Call DaXpY_(nTot1,One,Tmpx,1,Tmp1,1)
         If (nVxc.eq.2*nTot1) Then ! Nuc Attr added twice
            Call DaXpY_(nTot1,One,Tmpx(1+nTot1),1,Tmp1,1)
            Call Get_dArray('Nuc Potential',Tmpx,nTot1)
            Call DaXpY_(nTot1,-One,Tmpx,1,Tmp1,1)
         EndIf
         Call mma_deallocate(Tmpx)
         Call mma_deallocate(Tmp3)
         Call NameRun('#Pop')    ! switch back to old RUNFILE
      End If
*
*     Compute energy contributions
      Eone = dDot_(nTot1,Tmp2,1,Tmp1,1)
*     Call Get_PotNuc(PotNuc_Ref)
      Call Get_dScalar('PotNuc',PotNuc_Ref)
      Eone = Eone + (PotNuc-PotNuc_Ref)
      Etwo = dDot_(nTot1,Tmp2,1,FI,1)
      Call mma_deallocate(Tmp2)
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
        Write(LF,*) ' Tmp1 matrix in SGFCIN, one-electron term'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=1
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',Tmp1(iOff),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
        End Do
      End If
*
*     Assemble the one-electron Tmp1 and two-electron contribution to AO Fock matrix
      Call DaXpY_(nTot1,One,Tmp1,1,FI,1)
      Call mma_deallocate(Tmp1)
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
      CALL mma_allocate(X0,NTOT1,Label='X0')
      CALL mma_allocate(X1,NTOT1,Label='X1')
      CALL mma_allocate(X2,MXNB*MXNB,Label='X2')
      CALL mma_allocate(X3,MXNB*MXNA,Label='X3')
      CALL DCOPY_(NTOT1,FI,1,X1,1)
      If(KSDFT(1:3).ne.'SCF'.and.KSDFT(1:3).ne.'PAM') Then
         Call Get_dExcdRa(TmpFckI,nTmpFck)
         CALL DaXpY_(NTOT1,1.0D0,TmpFckI,1,X1,1)
        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,*)
          Write(LF,*) ' Exchange correlation in AO basis in SGFCIN'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            Call TriPrt(' ','(5G17.11)',TmpFckI(ioff),iBas)
            iOff = iOff + (iBas*iBas+iBas)/2
          End Do
        End If
        Call mma_deallocate(TmpFckI)
      End If
      If ( IPRLEV.ge.DEBUG ) then
       Write(LF,*)
       Write(LF,*) ' Modified FI in AO basis in SGFCIN'
       Write(LF,*) ' ---------------------'
       Write(LF,*)
       iOff=1
       Do iSym = 1,nSym
         iBas = nBas(iSym)
         Call TriPrt(' ','(5G17.11)',X1(ioff),iBas)
         iOff = iOff + (iBas*iBas+iBas)/2
       End Do
      End If
      CALL MOTRAC(CMO,X1,X2,X3)
      Call mma_deallocate(X3)
      Call mma_deallocate(X2)
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
              F(NTU)=X1(ITU)
              IF(NT.EQ.NU) F(NTU)=F(NTU)+EMYN
              X0(ITU) = F(NTU)
            END DO
          END DO
          IADD=IADD+NAT
        ENDIF
      END DO

!Quan: Fix bug, skip Lucia stuff with DMRG
! and other external CI solvers.
      if (.not. any([DoNECI, Do_CC_CI, DumpOnly,
     &              doDMRG, doBlockDMRG])) then
        INT1(1:ITU) = X0(1:ITU)
        INT1(ITU+1:) = Zero
        INT1O(1:ITU) = X0(1:ITU)
        INT1O(ITU+1:) = Zero
      endif

      Call mma_deallocate(X1)
      Call mma_deallocate(X0)

*     print h0
      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' Inactive Fock mat in act MO basis, h0, in SGFCIN'
        Write(LF,*) ' ------------'
        Write(LF,*)
        Call TriPrt(' ',' ',F,NAC)
      End If



      Return
      End
