!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1990, Markus P. Fuelscher                              *
!               2013, Giovanni Li Manni                                *
!***********************************************************************
      Subroutine CASDFT_terms(CMO,F,FI,D1I,D1A,D1S)
!
!     This routine is a modification of SGFCIN, adapted to a CASDFT
!     implementation in which the CI step of a CASDFT calculation is
!         not corrected by DFT. DFT will play a role only in the Orbital
!         optimization step.
!     Purpose:
!     Generate the Fock-matrix for the frozen and inactive orbitals.
!     Compute also the core energy. Finally, transform the Fock-matrix
!     into the basis of the active orbitals.
!
!     M.P. Fuelscher, Lund, July 1990
!     GLM, Minneapolis,   May 2013
!
      use stdalloc, only: mma_allocate, mma_deallocate
      use OneDat, only: sNoNuc, sNoOri
      use rctfld_module, only: lRF
#ifdef _DMRG_
!     module dependencies
      use qcmaquis_interface_cfg
#endif
      use Constants, only: Zero, One
      use rasscf_global, only: DFTFOCK, Emy, ExFac, KSDFT_temp,         &
     &                         NAC, NACPAR, NONEQ, PotNuc, Tot_Charge,  &
     &                         Tot_El_Charge, Tot_Nuc_Charge
#ifdef _DMRG_
      use lucia_data, only: INT1, INT1O
      use rasscf_global, only: DoDMRG
#endif
      use PrintLevel, only: DEBUG
      use output_ras, only: LF,IPRLOC
      use general_data, only: NSYM,NTOT1,NACTEL,ISPIN,NASH,NBAS,NFRO,   &
     &                        NISH

      Implicit None
!
      Real*8 CMO(*) ,F(*) , FI(*) , D1I(*) , D1A(*), D1S(*)
      Character(LEN=8) Label
      Logical First, Dff, Do_DFT
      Real*8, Allocatable:: X0(:), X1(:), X2(:), X3(:)
      Real*8, Allocatable:: Tmp0(:), Tmp1(:), Tmp2(:), Tmp3(:),         &
     &                      Tmp4(:), Tmp5(:), Tmp6(:), Tmp7(:)
      Real*8 CASDFT_Funct, Emyn, Eone, ETwo, PotNuc_Ref
      Real*8, External:: DDot_
      Integer i, IADD, iBas, iCharge, iCOmp, iOff, iOpt, iPrLev, iRC,   &
     &        iSyLbl, iSym, ITU, j, MXNA, MXNB, NAT, NST, NT, NTU, NU


!**********************************************************
! Local print level (if any)
!**********************************************************
      IPRLEV=IPRLOC(3)
!      IPRLEV=100
      If ( IPRLEV.ge.DEBUG ) then
       Write(LF,*) 'Printing matrices in CASDFT_Terms'
       Write(LF,*)
       Write(LF,*) ' CMO in CASDFT_terms'
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
       Write(LF,*) ' D1I in AO basis in CASDFT_Terms'
       Write(LF,*) ' ---------------------'
       Write(LF,*)
       iOff=1
       Do iSym = 1,nSym
         iBas = nBas(iSym)
         call wrtmat(D1I(ioff),iBas,iBas, iBas, iBas)
         iOff = iOff + iBas*iBas
       End Do

       Write(LF,*)
       Write(LF,*) ' D1S in AO basis in CASDFT_Terms'
       Write(LF,*) ' ---------------------'
       Write(LF,*)
       iOff=1
       Do iSym = 1,nSym
         iBas = nBas(iSym)
         call wrtmat(D1S(ioff),iBas,iBas, iBas, iBas)
         iOff = iOff + iBas*iBas
       End Do


       Write(LF,*)
       Write(LF,*) ' D1A in AO basis in CASDFT_Terms'
       Write(LF,*) ' ---------------------'
       Write(LF,*)
       iOff=1
       Do iSym = 1,nSym
         iBas = nBas(iSym)
         call wrtmat(D1A(ioff),iBas,iBas, iBas, iBas)
         iOff = iOff + iBas*iBas
       End Do

      End If


!
!**********************************************************
! Generate molecular charges
!**********************************************************
      Call mma_allocate(Tmp0,nTot1+4,Label='Tmp0')
      iRc=-1
      iOpt=ibset(0,sNoOri)
      iComp=1
      iSyLbl=1
      Label='Mltpl  0'
      Call RdOne(iRc,iOpt,Label,iComp,Tmp0,iSyLbl)
      Tot_Nuc_Charge=Tmp0(nTot1+4)
      If ( iRc.ne.0 ) then
         Write(LF,*) 'CASDFT_Terms: iRc from Call RdOne not 0'
         Write(LF,*) 'Label = ',Label
         Write(LF,*) 'iRc = ',iRc
         Call Abend
      Endif
      Call mma_deallocate(Tmp0)

      Tot_El_Charge=Zero
      Do iSym=1,nSym
         Tot_El_Charge=Tot_El_Charge                                    &
     &                -2.0D0*DBLE(nFro(iSym)+nIsh(iSym))
      End Do
      Tot_El_Charge=Tot_El_Charge-DBLE(nActEl)
      Tot_Charge=Tot_Nuc_Charge+Tot_El_Charge
!      If ( IPRLEV.ge.DEBUG ) then
!          write(6,*)
!       write(6,*) 'Total Charge :', Tot_Charge
!      end if
!
!**********************************************************
! Load bare nuclei Hamiltonian
!**********************************************************
      Call mma_Allocate(Tmp1,nTot1,Label='Tmp1')
      iComp  =  1
      iSyLbl =  1
      iRc    = -1
      iOpt   =  ibset(ibset(0,sNoOri),sNoNuc)
      Label  = 'OneHam  '
      Call RdOne(iRc,iOpt,Label,iComp,Tmp1,iSyLbl)
      If ( iRc.ne.0 ) then
         Write(LF,*) 'CASDFT_Terms: iRc from Call RdOne not 0'
         Write(LF,*) 'Label = ',Label
         Write(LF,*) 'iRc = ',iRc
         Call Abend
      Endif
      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' OneHam in AO basis in CASDFT_Terms'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=1
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',Tmp1(iOff),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
        End Do
      End IF
!
!**********************************************************
! Load the nuclear repulsion energy
!**********************************************************
      Call Get_dScalar('PotNuc',potNuc)
!      write(6,*)
!      write(6,*) 'PotNuc in casdft_terms.f:', PotNuc

!
!**********************************************************
! Generate total density
!**********************************************************
!
!      If ( IPRLEV.ge.DEBUG ) then
!       Call mma_allocate(Tmp31,nBas*nBas,Label='Tmp31')
!          Tmp31(:)=0.0D0
!          Call Daxpy_(nBas*nBas,1.0D0,D1I,1,Tmp31,1)
!          Call Daxpy_(nBas*nBas,1.0D0,D1A,1,Tmp31,1)
!      Write(LF,*)
!       Write(LF,*) ' DMAT not folded in AO basis in CASDFT_Terms'
!       Write(LF,*) ' ---------------------'
!       Write(LF,*)
!          call wrtmat(Tmp31,nBas,nBas, nBas, nBas)
!       Call mma_deallocate(Tmp3)
!      End IF

      Call mma_allocate(Tmp3,nTot1,Label='Tmp3')
      Call mma_allocate(Tmp4,nTot1,Label='Tmp4')
      Call Fold(nSym,nBas,D1I,Tmp3)
      Call Fold(nSym,nBas,D1A,Tmp4)
      Call Daxpy_(nTot1,1.0D0,Tmp4,1,Tmp3,1)
      Call Put_dArray('D1ao',Tmp3,nTot1)
      call xflush(6)
!**********************************************************
! Generate spin-density
!**********************************************************
        Call mma_allocate(Tmp7,nTot1,Label='Tmp7')
        Call Fold(nSym,nBas,D1S,Tmp7)
        Call Put_dArray('D1sao',Tmp7,nTot1)
        Call mma_deallocate(Tmp7)
!
!**********************************************************
! One- and two-electron type contributions
!**********************************************************
!
      Call mma_allocate(Tmp5,nTot1,Label='Tmp5')
      Tmp5(:)=0.0D0
      Call mma_allocate(Tmp6,nTot1,Label='Tmp6')
      Tmp6(:)=0.0D0
!
      First=.True.
      Dff=.False.
      Do_DFT=.True.

      Call Put_iArray('nFro',nFro,nSym)
      Call Put_iArray('nAsh',nAsh,nSym)
      Call Put_iArray('nIsh',nIsh,nSym)

      iCharge=Int(Tot_Charge)
! Tmp5 and Tmp6 are not updated in DrvXV...
      Call DrvXV(Tmp5,Tmp6,Tmp3,                                        &
     &             PotNuc,nTot1,First,Dff,NonEq,lRF,                    &
     &             KSDFT_TEMP,ExFac,iCharge,iSpin,                      &
     &             DFTFOCK,Do_DFT)

      Call Daxpy_(nTot1,1.0d0,Tmp5,1,Tmp1,1)
      Call Daxpy_(nTot1,1.0d0,Tmp6,1,FI,1)

      Call mma_deallocate(Tmp6)
      Call mma_deallocate(Tmp5)
      Call mma_deallocate(Tmp4)
      Call mma_deallocate(Tmp3)

!**********************************************************
!     Compute energy contributions
!**********************************************************
      Call mma_allocate(Tmp2,nTot1,Label='Tmp2')

      Call Fold(nSym,nBas,D1I,Tmp2)
!
      Eone = dDot_(nTot1,Tmp2,1,Tmp1,1)
      Call Get_dScalar('PotNuc',PotNuc_Ref)
      Eone = Eone + (PotNuc-PotNuc_Ref)
      Etwo = dDot_(nTot1,Tmp2,1,FI,1)
      Call mma_deallocate(Tmp2)
      EMY  = PotNuc_Ref+Eone+0.5d0*Etwo

      CASDFT_Funct = 0.0D0
      Call Get_dScalar('CASDFT energy',CASDFT_Funct)
      If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*) ' Nuclear repulsion energy :',PotNuc
!         Write(LF,*) ' Nuclear repulsion energy Ref :',PotNuc_Ref
         Write(LF,*) ' One-electron core energy :',Eone
         Write(LF,*) ' Two-electron core energy :',Etwo
         Write(LF,*) ' Total core energy        :',EMY
         Write(LF,*) ' CASDFT Energy            :',CASDFT_Funct
      End If

!**********************************************************
! Printing matrices
!**********************************************************
      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' FI matrix in CASDFT_Terms only 2-electron terms'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=1
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',FI(iOff),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
        End Do
      End If

      Call DaXpY_(nTot1,One,Tmp1,1,FI,1)

      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' Inactive Fock matrix in AO basis in CASDFT_terms'
        write(LF,*) '(it already contains OneHam and TwoEl contrib.)'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=1
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',FI(iOff),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
        End Do
      End If

      Call mma_deallocate(Tmp1)

!**********************************************************
!     Transform FI to active orbital basis and move it over to F.
!     Remove also the symmetry blocking.
!**********************************************************
!**********************************************************
! Shall I add here the DFT contribution? Maybe not yet!
! I am commenting off... if needed we can always re-activate.
!**********************************************************
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
      CALL dcopy_(NTOT1,FI,1,X1,1)
!      Call Get_dExcdRa(TmpFckI,nTmpFck)
!      CALL DaXpY_(NTOT1,1.0D0,TmpFckI,1,X1,1)
!     If ( IPRLEV.ge.DEBUG ) then
!         Write(LF,*)
!         Write(LF,*) ' Exchange Corr. in AO basis in CASDFT_Terms'
!         Write(LF,*) ' ---------------------'
!         Write(LF,*)
!         iOff=1
!         Do iSym = 1,nSym
!           iBas = nBas(iSym)
!           Call TriPrt(' ','(5G17.11)',TmpFckI(ioff),iBas)
!           iOff = iOff + (iBas*iBas+iBas)/2
!         End Do
!     End If
!      Call mma_deallocate(TmpFckI)
!      If ( IPRLEV.ge.DEBUG ) then
!        Write(LF,*)
!        Write(LF,*) ' Modified FI in AO basis in CASDFT_Terms'
!        Write(LF,*) ' ---------------------'
!        Write(LF,*)
!        iOff=1
!        Do iSym = 1,nSym
!          iBas = nBas(iSym)
!          Call TriPrt(' ','(5G17.11)',X1(ioff),iBas)
!          iOff = iOff + (iBas*iBas+iBas)/2
!        End Do
!      End If
!
      CALL MOTRAC(CMO,X1,X2,X3)
      CALL mma_deallocate(X3)
      CALL mma_deallocate(X2)
      CALL dcopy_(NACPAR,[ZERO],0,F,1)
      NTU=0
      ITU=0
      IADD=0

      IF (NACTEL.NE.0) THEN
         EMYN=EMY/DBLE(NACTEL)
      ELSE
         EMYN=0.0d0
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
#ifdef _DMRG_
      if(.not.doDMRG)then
        INT1(1:ITU) = X0(1:ITU)
        INT1(ITU+1:) = Zero
        INT1O(1:ITU) = X0(1:ITU)
        INT1O(ITU+1:) = Zero
      end if
#endif
      CALL mma_deallocate(X1)
      CALL mma_deallocate(X0)

!     print h0
      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*)
        Write(LF,*)
        Write(LF,*)
        Write(LF,*) ' Fock matrix in MO basis, h0, in CASDFT_TERMS'
        Write(LF,*) ' ------------'
        Call TriPrt(' ',' ',F,NAC)
      End If

      End Subroutine CASDFT_terms
