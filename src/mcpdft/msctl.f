!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1990, Markus P. Fuelscher                              *
!               2013, Giovanni Li Manni                                *
!               2016, Andrew M. Sand                                   *
!***********************************************************************
      Subroutine MSCtl(CMO,Ref_Ener)
      use definitions,only:iwp,wp,u6
      use constants,only:zero
      use OneDat, only: sNoNuc, sNoOri
      use mcpdft_input, only: mcpdft_options
      Use KSDFT_Info, only: do_pdftpot
      Use hybridpdft, only: E_NoHyb
      use mspdftgrad,only:P2MOT,D1aoMS,DIDA,D1SaoMS
      use mspdft, only: iIntS
      use printlevel, only: debug
      use mcpdft_output, only:iPrLoc
      use rctfld_module, only: lRF
      use stdalloc, only: mma_allocate, mma_deallocate
      use nq_info, only: Tau_a1, Tau_b1, Tau_a2, Tau_b2,
     &                   Lapl_a1, Lapl_b1, Lapl_a2, Lapl_b2
      use libxc_parameters, only: FuncExtParams
      use wadr, only: FockOcc
      use rasscf_global, only: DFTFOCK, nRoots, ExFac,
     &                         IADR15, IPR, lRoots, lSquare,
     &                         NAC, NACPAR, NACPR2, nFint, NonEq,
     &                         nTot4, PotNuc,
     &                         Tot_Nuc_Charge, ISTORP
      use general_data,only:nash,norb,nsym,ntot2,ntot1,jobiph,ispin,
     &              jobold,nactel,nbas,nish,nfro
      implicit none

      real(kind=wp),intent(inout) :: Ref_Ener(*)
      real(kind=wp),intent(in) :: CMO(*)

      character(len=8) Label
      Logical First, Dff, Do_DFT,Found

      real*8, allocatable:: FI_V(:), FA_V(:), FockI(:),
     &                      Tmp3(:), folded_dm1_cas(:),
     &                      dummy1(:), dummy2(:), folded_dm1s_cas(:),
     &                      dm1_core(:), casdm1(:),
     &                      focka(:), dm1_cas(:), dm1s_cas(:),
     &                      casdm1s(:), P2D(:), PUVX(:), P2t(:),
     &                      OnTopT(:), OnTopO(:),
     &                      TUVX_tmp(:),
     &                      P(:), FOCK(:), Q(:), Coul(:)
      real(kind=wp),allocatable :: int1e_ovlp(:), hcore(:)
      integer(kind=iwp) :: IAD19,iJOB,dmDisk, IADR19(1:30)
      integer(kind=iwp) :: jroot,NQ, isym,i, iCharge, iComp
      integer(kind=iwp) :: iOpt,  iPrLev,irc, iSA, iSyLbl
      integer(kind=iwp) :: niaia, tot_el_charge
      real(kind=wp), external :: energy_mcwfn
      real(kind=wp) :: casdft_e, casdft_funct, e_mcscf
      real(kind=wp) :: Energies(nroots)

      IPRLEV=IPRLOC(3)


!**********************************************************
! Generate molecular charges
!**********************************************************
      call mma_allocate(int1e_ovlp,nTot1+4,label="int1e_ovlp")
      iRc=-1
      iOpt=ibset(0,sNoOri)
      iComp=1
      iSyLbl=1
      Label='Mltpl  0'
      Call RdOne(iRc,iOpt,Label,iComp,int1e_ovlp,iSyLbl)
      If ( iRc.ne.0 ) then
        Write(u6,*) 'msctl: iRc from Call RdOne not 0'
        Write(u6,*) 'Label = ',Label
        Write(u6,*) 'iRc = ',iRc
        Call Abend
      Endif
      ! nuclear charge stored in last element
      Tot_Nuc_Charge=int1e_ovlp(size(int1e_ovlp))
      call mma_deallocate(int1e_ovlp)
      Tot_El_Charge=-2*sum(nFro+nIsh)-nActEl
      icharge = nint(Tot_Nuc_Charge)+Tot_El_Charge


!**********************************************************
! Load bare nuclei Hamiltonian
! This is h_pq but in the AO basis (so h_{mu, nu})
!**********************************************************
      call mma_allocate(hcore,nTot1,label='hcore')
      iComp  =  1
      iSyLbl =  1
      iRc    = -1
      iOpt   =  ibset(ibset(0,sNoOri),sNoNuc)
      Label  = 'OneHam  '
      Call RdOne(iRc,iOpt,Label,iComp,hcore,iSyLbl)
      If ( iRc.ne.0 ) then
        Write(u6,*) 'msctl: iRc from Call RdOne not 0'
        Write(u6,*) 'Label = ',Label
        Write(u6,*) 'iRc = ',iRc
        Call Abend()
      Endif

!Here we calculate the D1 Inactive matrix (AO).
      call mma_allocate(dm1_core,ntot2,label="dm1_core")
      Call Get_D1I_RASSCF(CMO,dm1_core)

      iJOB=0
      IAD19=0
      Call f_Inquire('JOBOLD',Found)
      If (.not.found) then
        Call f_Inquire('JOBIPH',Found)
        if (Found) JOBOLD=JOBIPH
      end if
      If (Found) iJOB=1
      If (iJOB.eq.1) Then
         if(JOBOLD.le.0) Then
           JOBOLD=20
           Call DaName(JOBOLD,'JOBOLD')
         end if
      end if
      IADR19(:)=0
      Call IDaFile(JOBOLD,2,IADR19,15,IAD19)
      IADR15 = IADR19
      dmDisk = IADR19(3)

      Call mma_allocate(casdm1,NACPAR,Label='casdm1')
      Call mma_allocate(dm1_cas,NTOT2,Label='dm1_cas')
      Call mma_allocate(casdm1s,NACPAR,Label='casdm1s')
      Call mma_allocate(dm1s_cas,NTOT2,Label='dm1s_cas')

      Call mma_allocate(Tmp3,nTot1,Label='Tmp3')
      Call mma_allocate(folded_dm1_cas,nTot1,Label='folded_dm1_cas')
      Call mma_allocate(P2d,NACPR2,Label='P2D')


      Call mma_allocate(FockI,ntot1,Label='FockI')
      Call mma_allocate(focka,ntot1,Label='focka')
      Call mma_allocate(coul,ntot1,Label='coul')

!This iSA is used to control gradient calculations.  Analytic gradients
!(in ALASKA) will only run if iSA=1, and iSA will only be set to one if
!the on-top potentials are computed as part of this calculation.
      iSA = 99
      Call Put_iScalar('SA ready',iSA)


!We loop over the number of states.  For each state, we read the density
!matrices from the JOBIPH file.
      do jroot=1,lroots
        iIntS=jRoot
        Tau_a1  = Zero
        Tau_b1  = Zero
        Tau_a2  = Zero
        Tau_b2  = Zero
        Lapl_a1 = Zero
        Lapl_b1 = Zero
        Lapl_a2 = Zero
        Lapl_b2 = Zero

!Read in the density matrices for <jroot>.
        casdm1(:)=0.0D0
        dm1_cas(:)=0.0D0
        casdm1s(:)=0.0D0
        dm1s_cas(:)=0.0D0
        Tmp3(:)=0.0D0
        folded_dm1_cas(:)=0.0D0
        P2D(:)=0.0D0

!Get the D1 Active matrix for this state.  These should probably be
!most easily read from the previous JOBIPH file.  Then, convert D1A from
!the MO to the AO basis.

        Call DDaFile(JOBOLD,2,casdm1,NACPAR,dmDisk)

        if(mcpdft_options%grad) then
          if(mcpdft_options%mspdft) then
            Call P2_contraction(casdm1,P2MOt(:,jroot))
          else if (jroot .eq. mcpdft_options%rlxroot) then
            Call mma_allocate(P2t,NACPR2,Label='P2t')
            P2t(:)=zero
            Call P2_contraction(casdm1,P2t)
            Call Put_dArray('P2MOt',P2t,NACPR2)
            Call mma_deallocate(P2t)
          end if
        end if

        Call Put_dArray('D1mo',casdm1,NACPAR)
        Call DDaFile(JOBOLD,2,casdm1s,NACPAR,dmDisk)
        Call DDaFile(JOBOLD,2,P2d,NACPR2,dmDisk)
        Call Put_dArray('P2mo',P2d,NACPR2)

        ! This dummy read is to cycle the dmDisk so next
        ! time we encounter this block, we load the next
        ! states densities
        Call DDaFile(JOBOLD,0,P2d,NACPR2,dmDisk)

!**********************************************************
! Generate total density
!**********************************************************

         If(NASH(1).ne.NAC) Call DBLOCK(casdm1)
         Call Get_D1A_RASSCF(CMO,casdm1,dm1_cas)

         Call Fold(nSym,nBas,dm1_core,Tmp3)
         Call Fold(nSym,nBas,dm1_cas,folded_dm1_cas)

      if(mcpdft_options%grad .and. mcpdft_options%mspdft)then
         Call Dcopy_(nTot1,folded_dm1_cas,1,DIDA(:,iIntS),1)
         if (iIntS.eq.lRoots)
     &   Call Dcopy_(ntot1,Tmp3,1,DIDA(:,lRoots+1),1)
      end if
         Call Daxpy_(nTot1,1.0D0,folded_dm1_cas,1,Tmp3,1)
!Maybe I can write all of these matrices to file, then modify stuff in
!the nq code to read in the needed density.  In other words, I need to
!replace the next call with something that supports multistates.
         Call Put_dArray('D1ao',Tmp3,nTot1)
         IF(mcpdft_options%grad.and.mcpdft_options%mspdft)THEN
          Call DCopy_(nTot1,Tmp3,1,D1AOMS(:,jRoot),1)
         END IF

!Get the spin density matrix for open shell cases
!**********************************************************
! Generate spin-density
!**********************************************************
         if(iSpin.eq.1) then
            dm1s_cas(:) = zero
         end if
         IF ( NASH(1).NE.NAC ) CALL DBLOCK(casdm1s)
         Call Get_D1A_RASSCF(CMO,casdm1s,dm1s_cas)
      Call mma_allocate(folded_dm1s_cas,nTot1,Label='folded_dm1s_cas')
         Call Fold(nSym,nBas,dm1s_cas,folded_dm1s_cas)
         Call Put_dArray('D1sao',folded_dm1s_cas,nTot1)
         IF(iSpin.ne.1.and. mcpdft_options%grad
     &      .and.mcpdft_options%mspdft) THEN
         Call DCopy_(nTot1,folded_dm1s_cas,1,D1SAOMS(:,jRoot),1)
         END IF
         Call mma_deallocate(folded_dm1s_cas)


!**********************************************************
! Calculation of the CASDFT_energy
!**********************************************************
!Perhaps ideally, we should rework how DrvXV (and its children) handles
!the AO to MO transformation on the grid.  It seems like perhaps we are
!doing redundant transformations by retransforming AOs (which may have
!been included in a previous batch) into MOs.
! dummy1 and dummy2 are not updated in DrvXV...
        Call mma_allocate(dummy1,nTot1,Label='dummy1')
        Call mma_allocate(dummy2,nTot1,Label='dummy2')
        dummy1(:)=0.0D0
        dummy2(:)=0.0D0
        First=.True.
        Dff=.False.
        Do_DFT=.True.

        Call Put_iArray('nFro',nFro,nSym)
        Call Put_iArray('nAsh',nAsh,nSym)
        Call Put_iArray('nIsh',nIsh,nSym)


      do_pdftPot=.false.
      if (mcpdft_options%grad .and.
     &    (mcpdft_options%mspdft .or.
     &     (jroot .eq. mcpdft_options%rlxroot))) then

        do_pdftPot=.true.

      end if

        Call DrvXV(dummy1,dummy2,Tmp3,
     &             PotNuc,nTot1,First,Dff,NonEq,lRF,
     &             mcpdft_options%otfnal%otxc,ExFac,iCharge,iSpin,
     &             DFTFOCK,Do_DFT)

        Call mma_deallocate(dummy1)
        Call mma_deallocate(dummy2)

        CASDFT_Funct = 0.0D0
        Call Get_dScalar('CASDFT energy',CASDFT_Funct)


      call mma_allocate(tuvx_tmp,nacpr2,Label='tuvx_tmp')
      call mma_allocate(puvx,nfint,Label='puvx')
      tuvx_tmp(:) = zero
      puvx(:) = zero
      focka(:) = zero
      focki(:) = zero
      coul(:) = zero

    ! This constructs focki and focka for us. Technically,
    ! focki is a constant but, if we have to recalculate it,
    ! why store it in memory.
         CALL TRACTL2(cmo,PUVX,TUVX_tmp,
     &                dm1_core,focki,
     &                dm1_cas,focka,
     &                IPR,lSquare,ExFac)

      call mma_deallocate(tuvx_tmp)

    ! Note that because ExFac should be 0, we get back
    ! just Coulomb integrals!
      coul(:) = focki(:) + focka(:)

      e_mcscf = energy_mcwfn(tmp3,hcore,coul,PotNuc,ntot1)

         CASDFT_E = e_mcscf+CASDFT_Funct

        IF(mcpdft_options%otfnal%is_hybrid()) THEN
            E_NoHyb=CASDFT_E
            CASDFT_E = mcpdft_options%otfnal%lambda*Ref_Ener(jRoot) +
     &            (1.0-mcpdft_options%otfnal%lambda) * E_NoHyb
        END IF

        Call Print_MCPDFT_2(CASDFT_E,PotNuc,e_mcscf,CASDFT_Funct,
     &         jroot,Ref_Ener)


!JB         replacing ref_ener with MC-PDFT energy for MS-PDFT use
          Energies(jroot)=CASDFT_E
          Ref_Ener(jroot)=CASDFT_E

!At this point, the energy calculation is done.  Now I need to build the
!fock matrix if this root corresponds to the relaxation root.

!***********************************************************************
!
!            BUILDING OF THE NEW FOCK MATRIX                           *
!
!***********************************************************************
      if(mcpdft_options%grad) then
        ! Determine size of Q matrix
         NQ=0
         NIAIA=0
         do ISYM=1,NSYM
           NQ = MAX(NQ,NASH(ISYM)*NORB(ISYM))
           NIAIA = NIAIA+(NASH(ISYM)+NISH(ISYM))**2
         end do
         if(NQ.lt.NIAIA) NQ=NIAIA


       if((.not. mcpdft_options%mspdft)
     &   .and. jroot .eq. mcpdft_options%rlxroot) then

         Write(u6,*) 'Calculating potentials for analytic gradients...'

        IF(ISTORP(NSYM+1).GT.0) THEN
           call mma_allocate(P,ISTORP(NSYM+1),Label='P')
        else
          call mma_allocate(P,1,Label='P')
         END IF
!MCLR requires two sets of things:
!1. The effective one-body Fock matrix and the effective two-body fock
!matrix.  These are used to generate the CI gradient inside of MCLR
!2. The effective generalized fock matrix.  This is used to calculate
!the orbital gradient inside of MCLR and is also used in the
!renormalization/effective lagrangian part of the final gradient
!evalutation.

!I think the plan should be to add on the missing pieces (to Fock_occ)
!which come from the one- and two-electron potentials.  These pieces are
!given by
! F_{xy} = \sum_{p} V_{py} D_{px} + \sum_{pqr} 2v_{pqry}d_{pqrx}.


!I will read in the one- and two-electron potentials here

      Call mma_allocate(ONTOPT,nfint,Label='OnTopT')
      Call mma_allocate(ONTOPO,ntot1,Label='OnTopO')
      Call Get_dArray('ONTOPT',OnTopT,NFINT)
      Call Get_dArray('ONTOPO',OnTopO,NTOT1)

!Grab the active-active part of the FI+FA matrix (currently held in the
!FA matrix) and place it in an array of size NACPAR.  Add the oeotp to
!it.  Write to file.

!I think I need to generate FI, which will contain both the one-electron
!potential contribution and the V_kkpu contribution.



      CALL mma_allocate(FI_V,Ntot1,Label='FI_V')
      Call Get_dArray('FI_V',FI_V,NTOT1)

      call ao2mo_1e(cmo,hcore(:)+coul(:),
     &       focki,nsym,nbas,norb,nfro)
      fi_v(:) = fi_v(:) + ontopo(:) + focki(:)

      Call mma_allocate(TUVX_tmp,NACPR2,Label='TUVX_tmp')
      Call Get_TUVX(OnTopT,TUVX_tmp)

      !Add the V_kktu contribution to Fone_tu?
!STILL MUST DO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This should be addressed in the upd_FI routine.

      call put_darray('F1_PDFT         ',fi_v(:),ntot1)
      call put_darray('F2_PDFT         ',tuvx_tmp,nacpr2)

      Call mma_deallocate(TUVX_tmp)

!____________________________________________________________
!This next part is to generate the MC-PDFT generalized fock operator.

!The corrections (from the potentials) to FI and FA are built in the NQ
!part of the code, for efficiency's sake.  It still needs to be
!debugged.
      CALL mma_allocate(FA_V,Ntot1,Label='FA_V')
      Call Get_dArray('FA_V',FA_V,NTOT1)

      If ( IPRLEV.ge.DEBUG ) then
      write(u6,*) "extra terms to update FI"
      do i=1,ntot1
        write(u6,*) FI_V(i)
      end do
      write(u6,*) "extra terms to update FA"
      do i=1,ntot1
        write(u6,*) FA_V(i)
      end do
        end if

!Reordering of the two-body density matrix.

       IF(ISTORP(NSYM+1).GT.0) THEN
         P(:)=zero
         CALL PMAT_RASSCF(P2d,P)
      END IF

!Must add to existing FOCK operator (occ/act). FOCK is not empty.
         CALL mma_allocate(Q,NQ,Label='Q') ! q-matrix(1symmblock)
         call mma_allocate(fock,ntot4,label='Fock')
        fock(:) = zero
         CALL fock_update(fock,fi_v,fa_v,casdm1,P,
     &                    Q,OnTopT,CMO)

         Call Put_dArray('FockOcc',FockOcc,ntot1)
         call put_darray('Fock_PDFT',fock,ntot4)

         Call mma_deallocate(Q)
      Call mma_deallocate(ONTOPO)
      Call mma_deallocate(ONTOPT)
      CALL mma_deallocate(FI_V)
      CALL mma_deallocate(FA_V)
      call mma_deallocate(P)
      call mma_deallocate(fock)


!Put some information on the runfile for possible gradient calculations.
      Call Put_iScalar('Number of roots',nroots)
      Call Put_dArray('Last energies',Energies,nroots)
      Call Put_dScalar('Last energy',Energies(mcpdft_options%RlxRoot))


!Put information needed for geometry optimizations.
!need to do MCLR for gradient runs. (1 to run, 2 to skip)
      iSA = 1
      Call Put_iScalar('SA ready',iSA)
      Call Put_cArray('MCLR Root','****************',16)
      Call Put_iScalar('Relax CASSCF root',mcpdft_options%rlxroot)

      end if

      if (mcpdft_options%mspdft) then
!      doing exactly the same thing as done in the previous chunck
!      starting from 'BUILDING OF THE NEW FOCK MATRIX'
!      Hopefully this code will be neater.
       call savefock_pdft(CMO,hcore,coul,casdm1,NQ,p2d,jroot)
      end if


      endif

    ! this allocation/deallocation could be done outside the root loop
      call mma_deallocate(puvx)

      end do !loop over roots

      if(mcpdft_options%grad) then
        dmDisk = IADR19(3)
        do jroot=1,mcpdft_options%rlxroot-1
          Call DDaFile(JOBOLD,0,casdm1,NACPAR,dmDisk)
          Call DDaFile(JOBOLD,0,casdm1s,NACPAR,dmDisk)
          Call DDaFile(JOBOLD,0,P2d,NACPR2,dmDisk)
          Call DDaFile(JOBOLD,0,P2d,NACPR2,dmDisk)
        end do
        Call DDaFile(JOBOLD,2,casdm1,NACPAR,dmDisk)
!        Andrew added this line to fix heh2plus
        Call DDaFile(JOBOLD,2,casdm1s,NACPAR,dmDisk)
        Call Put_dArray('D1mo',casdm1,NACPAR)
        Call DDaFile(JOBOLD,2,P2d,NACPR2,dmDisk)
        Call Put_dArray('P2mo',P2d,NACPR2)

         If(NASH(1).ne.NAC) Call DBLOCK(casdm1)
         Call Get_D1A_RASSCF(CMO,casdm1,dm1_cas)

         Call Fold(nSym,nBas,dm1_core,Tmp3)
         Call Fold(nSym,nBas,dm1_cas,folded_dm1_cas)
         Call Daxpy_(nTot1,1.0D0,folded_dm1_cas,1,Tmp3,1)
         Call Put_dArray('D1ao',Tmp3,nTot1)

!Get the spin density matrix for open shell cases
!**********************************************************
! Generate spin-density
!**********************************************************
           dm1s_cas(:) = zero
         IF ( NASH(1).NE.NAC ) CALL DBLOCK(casdm1s)
         Call Get_D1A_RASSCF(CMO,casdm1s,dm1s_cas)
      Call mma_allocate(folded_dm1s_cas,nTot1,Label='folded_dm1s_cas')
         Call Fold(nSym,nBas,dm1s_cas,folded_dm1s_cas)
         Call Put_dArray('D1Sao',folded_dm1s_cas,nTot1)
         Call mma_deallocate(folded_dm1s_cas)

        Call DDaFile(JOBOLD,0,P2d,NACPR2,dmDisk)
      end if

!Free up all the memory we can here, eh?
      Call mma_deallocate(folded_dm1_cas)
      Call mma_deallocate(Tmp3)

      Call mma_deallocate(casdm1)
      Call mma_deallocate(dm1_cas)
      Call mma_deallocate(casdm1s)
      Call mma_deallocate(dm1s_cas)
      call mma_deallocate(hcore)
      Call mma_deallocate(FockI)
      Call mma_deallocate(coul)
      Call mma_deallocate(focka)
      Call mma_deallocate(P2D)
      Call mma_deallocate(dm1_core)

      If (Allocated(FuncExtParams)) Call mma_deallocate(FuncExtParams)

      END Subroutine MSCtl

      Subroutine P2_contraction(D1MO,P2MO)
      use definitions, only: iwp,wp
      use constants,only:one,half
      use rasscf_global, only: NAC

      implicit none

      real(kind=wp), dimension(*), intent(in) :: d1mo
      real(kind=wp), dimension(*), intent(out) :: p2mo

      integer(kind=iwp) :: i, j, k, l, ij, kl, ijkl, lmax
      real(kind=wp) :: fact

      ijkl=0
      do i=1,nac
        do j=1,i
          ij = iTrii(i,j)
          do k=1,i
            if(i == k) then
              lmax = j
            else
              lmax = k
            end if
            do l=1,lmax
              kl = iTrii(k,l)
              ijkl = ijkl + 1
              fact=one
              if(k == l) fact=half
              p2MO(ijkl) = fact*D1MO(ij)*D1MO(kl)
            end do
          end do
        end do
      end do
      contains
        integer function iTrii(i,j)
          integer, intent(in) :: i, j
          itrii = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
        end function
      end Subroutine P2_contraction
