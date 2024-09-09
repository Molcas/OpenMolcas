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
*               2013, Giovanni Li Manni                                *
*               2016, Andrew M. Sand                                   *
************************************************************************
      Subroutine MSCtl(CMO,FI,FA,Ref_Ener)
!CMO,F,FI
*
*     This routine is a modification of SGFCIN, adapted to a CASDFT
*     implementation in which the CI step of a CASDFT calculation is
*         not corrected by DFT. DFT will play a role only in the Orbital
*         optimization step.
*     Purpose:
*     Generate the Fock-matrix for the frozen and inactive orbitals.
*     Compute also the core energy. Finally, transform the Fock-matrix
*     into the basis of the active orbitals.
*
*     M.P. Fuelscher, Lund, July 1990
*     GLM, Minneapolis,   May 2013
*     AMS, Minneapolis,   Feb 2016
*
      use OneDat, only: sNoNuc, sNoOri
      use mcpdft_input, only: mcpdft_options
      Use KSDFT_Info, only: do_pdftpot
      Use hybridpdft, only: E_NoHyb
      use mspdft, only: do_rotate, iIntS, DIDA, P2MOt, D1AOMS, D1SAOMS
      use printlevel, only: debug
      use mcpdft_output, only: lf, iPrLoc
      use rctfld_module
      use stdalloc, only: mma_allocate, mma_deallocate
      use wadr, only: FockOcc, TUVX

      Implicit Real*8 (A-H,O-Z)

      Real*8 CMO(*), FI(*), FA(*), Ref_Ener(*)
*
#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
#include "pamint.fh"
#include "timers.fh"
#include "SysDef.fh"
!      Logical TraOnly

*
      Character(LEN=8) Label
      Logical First, Dff, Do_DFT,Found
      Parameter ( Zero=0.0d0 , One=1.0d0 )
      integer IAD19
      integer iJOB,dmDisk
      integer IADR19(1:30)
      integer NQ
      integer jroot
      real*8,dimension(1:nroots) :: Energies
      integer  i_off1,i_off2
      integer isym,iash
      integer, External:: IsFreeUnit

      real*8, allocatable:: FI_V(:), FA_V(:), FockI(:), Tmp0(:),
     &                      Tmp1(:), Tmp2(:), Tmp3(:), Tmp4(:),
     &                      Tmp5(:), Tmp6(:), Tmp7(:), Tmpn(:),
     &                      Tmpk(:), Tmpa(:), D1I(:), D1Act(:),
     &                      FockA(:), D1ActAO(:), D1SpinAO(:),
     &                      D1Spin(:), P2D(:), PUVX(:), P2t(:),
     &                      DtmpA_g(:), OnTopT(:), OnTopO(:),
     &                      D1Act_FA(:), D1ActAO_FA(:), CMO_X(:),
     &                      FockI_Save(:), TUVX_tmp(:), PUVX_tmp(:),
     &                      P(:), P1(:), FOCK(:), Q(:), BM(:),
     &                      FOne(:), FA_t(:)

***********************************************************
C Local print level (if any)
***********************************************************
      IPRLEV=IPRLOC(3)


*TRS
      Call Get_iScalar('Relax CASSCF root',iRlxRoot)

***********************************************************
* Generate molecular charges
***********************************************************
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
        Tot_El_Charge=Tot_El_Charge
     &                -2.0D0*DBLE(nFro(iSym)+nIsh(iSym))
      End Do
      Tot_El_Charge=Tot_El_Charge-DBLE(nActEl)
      Tot_Charge=Tot_Nuc_Charge+Tot_El_Charge


***********************************************************
* Load bare nuclei Hamiltonian
! This is h_pq but in the AO basis (so h_{mu, nu})
***********************************************************
      Call mma_allocate(Tmp1,nTot1,Label='Tmp1')
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
      End If

      Call mma_allocate(Tmpk,nTot1,Label='Tmpk')
c--reads kinetic energy integrals  Tmpk--(Label=Kinetic)----
      iComp  =  1
      iSyLbl =  1
      iRc    = -1
      iOpt   =  ibset(ibset(0,sNoOri),sNoNuc)
      Label  = 'Kinetic '
      Call RdOne(iRc,iOpt,Label,iComp,Tmpk,iSyLbl)
      If ( iRc.ne.0 ) then
        Write(LF,*) 'CASDFT_Terms: iRc from Call RdOne not 0'
        Write(LF,*) 'Label = ',Label
        Write(LF,*) 'iRc = ',iRc
        Call Abend
      Endif
      Call mma_allocate(Tmpn,nTot1,Label='Tmpn')
      iComp  =  1
      iSyLbl =  1
      iRc    = -1
      iOpt   =  ibset(ibset(0,sNoOri),sNoNuc)
      Label  = 'Attract '
      Call RdOne(iRc,iOpt,Label,iComp,Tmpn,iSyLbl)
      If ( iRc.ne.0 ) then
        Write(LF,*) 'CASDFT_Terms: iRc from Call RdOne not 0'
        Write(LF,*) 'Label = ',Label
        Write(LF,*) 'iRc = ',iRc
        Call Abend
      Endif


!Here we calculate the D1 Inactive matrix (AO).
      Call mma_allocate(D1I,NTOT2,Label='D1I')
      Call Get_D1I_RASSCF(CMO,D1I)

      IF(IPRLEV.ge.DEBUG) THEN
        write(6,*) 'iD1inact'
        do i=1,ntot2
        write(6,*) D1i(i)
        end do
      END IF

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

      Call mma_allocate(D1Act,NACPAR,Label='D1Act')
      Call mma_allocate(D1ActAO,NTOT2,Label='D1ActAO')
      Call mma_allocate(D1Spin,NACPAR,Label='D1Spin')
      Call mma_allocate(D1SpinAO,NTOT2,Label='D1SpinAO')

      Call mma_allocate(Tmp3,nTot1,Label='Tmp3')
      Call mma_allocate(Tmp4,nTot1,Label='Tmp4')
      Call mma_allocate(P2d,NACPR2,Label='P2D')


      Call mma_allocate(FockI,ntot1,Label='FockI')
      Call mma_allocate(FockA,ntot1,Label='FockA')

************************************************************************
* load back two-electron integrals (pu!vx)
************************************************************************

      If ( nFint.gt.0) then
        Call mma_allocate(PUVX,nFint,Label='PUVX')
        iDisk = 0
        Call DDaFile(LUINTM,2,PUVX,nFint,iDisk)
      Else
        Call mma_allocate(PUVX,1,Label='PUVX')
      End If

      IF(IPRLEV.ge.DEBUG) THEN
        write(6,*) 'PUVX integrals in msctl'
        call wrtmat(PUVX,1,nFInt,1,nFInt)
      END IF


!This iSA is used to control gradient calculations.  Analytic gradients
!(in ALASKA) will only run if iSA=1, and iSA will only be set to one if
!the on-top potentials are computed as part of this calculation.
      iSA = 99
      Call Put_iScalar('SA ready',iSA)


!We loop over the number of states.  For each state, we read the density
!matrices from the JOBIPH file.
      do jroot=1,lroots
        iIntS=jRoot

       !Load a fresh FockI and FockA
        Call dcopy_(ntot1,FI,1,FockI,1)
        Call dcopy_(ntot1,FA,1,FockA,1)
*
!Read in the density matrices for <jroot>.
        D1Act(:)=0.0D0
        D1ActAO(:)=0.0D0
        D1Spin(:)=0.0D0
        D1SpinAO(:)=0.0D0
        Tmp3(:)=0.0D0
        Tmp4(:)=0.0D0
        P2D(:)=0.0D0

!Get the D1 Active matrix for this state.  These should probably be
!most easily read from the previous JOBIPH file.  Then, convert D1A from
!the MO to the AO basis.

        Call DDaFile(JOBOLD,2,D1Act,NACPAR,dmDisk)

        if(mcpdft_options%grad) then
          if(mcpdft_options%mspdft) then
            Call P2_contraction(D1Act,P2MOt(:,jroot))
          else if (jroot .eq. irlxroot) then
            Call mma_allocate(P2t,NACPR2,Label='P2t')
            P2t(:)=0.0D0
            Call P2_contraction(D1Act,P2t)
            Call Put_dArray('P2MOt',P2t,NACPR2)
            Call mma_deallocate(P2t)
          end if
        end if

        IF(IPRLEV.ge.DEBUG) THEN
          write(6,*) 'D1Act'
          do i=1,NACPAR
            write(6,*) D1Act(i)
          end do
        end if

        Call Put_dArray('D1mo',D1Act,NACPAR)
        Call DDaFile(JOBOLD,2,D1Spin,NACPAR,dmDisk)
        Call DDaFile(JOBOLD,2,P2d,NACPR2,dmDisk)
        Call Put_dArray('P2mo',P2d,NACPR2)

        Call DDaFile(JOBOLD,0,P2d,NACPR2,dmDisk)

        IF(IPRLEV.ge.DEBUG) THEN
          write(6,*) 'D2'
          do i=1,NACPR2
            write(6,*) p2d(i)
          end do
        END IF

***********************************************************
* Generate total density
***********************************************************

         If(NASH(1).ne.NAC) Call DBLOCK(D1Act)
      IF(IPRLEV.ge.DEBUG) THEN
        write(6,*) 'D1Act'
        do i=1,NACPAR
          write(6,*) D1Act(i)
        end do
      end if
         Call Get_D1A_RASSCF(CMO,D1Act,D1ActAO)

!ANDREW _ RIGHT HERE
      if (mcpdft_options%grad .and.
     &     (mcpdft_options%mspdft .or. (jroot .eq. irlxroot))) then
        Call mma_allocate(DtmpA_g,nTot1,Label='DtmpA_g')
        Call Fold_pdft(nSym,nBas,D1ActAO,DtmpA_g)
        Call put_darray('D1Activeao',DtmpA_g,ntot1)
        Call mma_deallocate(DtmpA_g)
      end if
!END _RIGHT HERE
*
         Call Fold(nSym,nBas,D1I,Tmp3)
         Call Fold(nSym,nBas,D1ActAO,Tmp4)
*
      if(mcpdft_options%grad .and. mcpdft_options%mspdft)then
         Call Dcopy_(nTot1,Tmp4,1,DIDA(:,iIntS),1)
         if (iIntS.eq.lRoots)
     &   Call Dcopy_(ntot1,Tmp3,1,DIDA(:,lRoots+1),1)
      end if
         Call Daxpy_(nTot1,1.0D0,Tmp4,1,Tmp3,1)
!Maybe I can write all of these matrices to file, then modify stuff in
!the nq code to read in the needed density.  In other words, I need to
!replace the next call with something that supports multistates.
         Call Put_dArray('D1ao',Tmp3,nTot1)
         IF(mcpdft_options%grad.and.mcpdft_options%mspdft)THEN
          Call DCopy_(nTot1,Tmp3,1,D1AOMS(:,jRoot),1)
         END IF
      IF(IPRLEV.ge.DEBUG) THEN
         write(6,*) 'd1ao'
         do i=1,ntot1
           write(6,*) tmp3(i)
         end do
cPS         call xflush(6)
      end if

!Get the spin density matrix for open shell cases
***********************************************************
* Generate spin-density
***********************************************************
         if(iSpin.eq.1) then
           Call dcopy_(NACPAR,[0.0d0],0,D1SpinAO,1)
         end if
         IF ( NASH(1).NE.NAC ) CALL DBLOCK(D1Spin)
         Call Get_D1A_RASSCF(CMO,D1Spin,D1SpinAO)
         Call mma_allocate(Tmp7,nTot1,Label='Tmp7')
         Call Fold(nSym,nBas,D1SpinAO,Tmp7)
         Call Put_dArray('D1sao',Tmp7,nTot1)
         IF(iSpin.ne.1.and. mcpdft_options%grad
     &      .and.mcpdft_options%mspdft) THEN
         Call DCopy_(nTot1,Tmp7,1,D1SAOMS(:,jRoot),1)
         END IF
      IF(IPRLEV.ge.DEBUG) THEN
         write(6,*) 'd1so'
         do i=1,ntot1
           write(6,*) tmp7(i)
         end do
      end if
         Call mma_deallocate(Tmp7)


***********************************************************
* Calculation of the CASDFT_energy
***********************************************************
!Perhaps ideally, we should rework how DrvXV (and its children) handles
!the AO to MO transformation on the grid.  It seems like perhaps we are
!doing redundant transformations by retransforming AOs (which may have
!been included in a previous batch) into MOs.
*
        Call mma_allocate(Tmp5,nTot1,Label='Tmp5')
        Call mma_allocate(Tmp6,nTot1,Label='Tmp6')
        Tmp5(:)=0.0D0
        Tmp6(:)=0.0D0
*
        First=.True.
        Dff=.False.
        Do_DFT=.True.

        Call Put_iArray('nFro',nFro,nSym)
        Call Put_iArray('nAsh',nAsh,nSym)
        Call Put_iArray('nIsh',nIsh,nSym)

        iCharge=Int(Tot_Charge)

c Tmp5 and Tmp6 are not updated in DrvXV...
                   NTU=0
                   ITU=0
                   IADD=0

                    Do iSym = 1,nSym
                    iAsh = nAsh(iSym)
                    IF(iAsh.NE.0) THEN
                    DO NT=1,iAsh
                    NTU=NTU+IADD
                    DO NU=1,NT
                    NTU=NTU+1
                     ITU=ITU+1
                    ENDDO
                    ENDDO
                    IADD=IADD+iAsh
                    End If
                    End Do

      do_pdftPot=.false.
      if (mcpdft_options%grad .and.
     &    (mcpdft_options%mspdft .or. (jroot .eq. irlxroot))) then

        do_pdftPot=.true.

        CALL mma_allocate(OnTopT,NFINT,Label='OnTopT')
        OnTopT(:)=0.0D0
        CALL mma_allocate(OnTopO,NTOT1,Label='OnTopO')
        OnTopO(:)=0.0D0

        Call Put_dArray('ONTOPO',OnTopO,NTOT1)
        Call Put_dArray('ONTOPT',OnTopT,NFINT)
        Call Put_dArray('FI_V',OnTopO,NTOT1)
        Call Put_dArray('FA_V',OnTopO,NTOT1)

        CALL mma_deallocate(OnTopO)
        CALL mma_deallocate(OnTopT)
      end if

        Call DrvXV(Tmp5,Tmp6,Tmp3,
     &             PotNuc,nTot1,First,Dff,NonEq,lRF,
     &             mcpdft_options%otfnal%otxc,ExFac,iCharge,iSpin,
     &             D1I,D1ActAO,nTot1,DFTFOCK,Do_DFT)


        Call Daxpy_(nTot1,1.0d0,Tmp5,1,Tmp1,1)
        Call Daxpy_(nTot1,1.0d0,Tmp6,1,FockI,1)

        Call mma_deallocate(Tmp6)
        Call mma_deallocate(Tmp5)

!
***********************************************************
*     Compute energy contributions
***********************************************************
      Call mma_allocate(Tmp2,nTot1,Label='Tmp2')

      Call Fold(nSym,nBas,D1I,Tmp2)
c         call xflush(6)

      Call mma_allocate(Tmpa,nTot1,Label='Tmpa')
c         call xflush(6)
      Call Fold(nSym,nBas,D1ActAO,Tmpa)
c         call xflush(6)
*
      Eone = dDot_(nTot1,Tmp2,1,Tmp1,1)
      Call Get_dScalar('PotNuc',PotNuc_Ref)
      Eone = Eone + (PotNuc-PotNuc_Ref)
      Etwo = dDot_(nTot1,Tmp2,1,FockI,1)

!**************Kinetic energy of inactive electrons********
      Ekin = dDot_(nTot1,Tmp2,1,Tmpk,1)

!*****Nuclear electron attraction for inactive electrons******
      Enuc = dDot_(nTot1,Tmp2,1,Tmpn,1)

c**************Kinetic energy of active electrons*********
      EactK = dDot_(nTot1,Tmpk,1,Tmpa,1)

      EactN = dDot_(nTot1,Tmpn,1,Tmpa,1)
c         call xflush(6)
      EMY  = PotNuc_Ref+Eone+0.5d0*Etwo

      CASDFT_Funct = 0.0D0
      Call Get_dScalar('CASDFT energy',CASDFT_Funct)
*TRS
      If ( IPRLEV.ge.DEBUG ) then
       Write(LF,'(4X,A35,F18.8)')
     &  'Nuclear repulsion energy :',PotNuc
       Write(LF,'(4X,A35,F18.8)')
     &  'One-electron kinetic core energy:',Ekin
       Write(LF,'(4X,A35,F18.8)')
     &   'Nuc-elec attraction core energy:',Enuc
       Write(LF,'(4X,A35,F18.8)') 'One-electron core energy:',Eone
       Write(LF,'(4X,A35,F18.8)') 'Two-electron core energy:',Etwo
       Write(LF,'(4X,A35,F18.8)') 'Total core energy:',EMY
       Write(LF,'(4X,A35,F18.8)') 'Active Kinetic energy:',EactK
       Write(LF,'(4X,A35,F18.8)')
     &  'Active nuc-elec attraction energy:',EactN
c       Write(LF,*) ' CASDFT Energy            :',CASDFT_Funct
      End If
c         call xflush(6)
***********************************************************
* Printing matrices
***********************************************************
      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' FI matrix in CASDFT_Terms only 2-electron terms'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=1
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',FockI(iOff),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
        End Do
      End If

      Call DaXpY_(nTot1,One,Tmp1,1,FockI,1)

      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' Inactive Fock matrix in AO basis in CASDFT_terms'
        write(LF,*) '(it already contains OneHam and TwoEl contrib.)'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=1
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',FockI(iOff),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
        End Do
      End If


************************************************************************
* update and transform the Fock matrices FI and FA ----> FMAT routine
************************************************************************

        if (iprlev.ge.debug) then
              write(6,*) 'D1Act before reading in'
              do i=1,nacpar
                write(6,*) D1Act(i)
              end do
        end if
*
        Call mma_allocate(D1Act_FA,nacpar,Label='D1Act_FA')
        Call mma_allocate(D1ActAO_FA,ntot2,Label='D1ActAO_FA')
*
        itsDisk = IADR19(3)
        do i=1,irlxroot-1
          Call DDaFile(JOBOLD,0,D1Act_FA,NACPAR,itsDisk)
          Call DDaFile(JOBOLD,0,D1Spin,NACPAR,itsDisk)
          Call DDaFile(JOBOLD,0,P2d,NACPR2,itsDisk)
          Call DDaFile(JOBOLD,0,P2d,NACPR2,itsDisk)
        end do
        Call DDaFile(JOBOLD,2,D1Act_FA,NACPAR,itsDisk)
        IF ( NASH(1).NE.NAC ) CALL DBLOCK(D1Act_FA)
        Call Get_D1A_RASSCF(CMO,D1Act_FA,D1ActAO_FA)
*****
        Call mma_allocate(CMO_X,ntot2,Label='CMO_X')
        CALL DCOPY_(NTOT2,CMO,1,CMO_X,1)
        if(iprlev.ge.debug) then
            write(6,*) 'cmo before tractl'
            do i=1,ntot2
              write(6,*) CMO_X(i)
            end do
            write(6,*) 'puvx before tractl'
            do i=1,nfint
              write(6,*) puvx(i)
            end do
            write(6,*) 'tuvx after tractl'
            do i=1,nacpr2
              write(6,*) tuvx(i)
            end do
            write(6,*) 'D1Act_FA before tractl'
            do i=1,nacpar
              write(6,*) D1Act_FA(i)
            end do
            write(6,*) 'D1ActAO_FA before tractl'
            do i=1,ntot2
              write(6,*) D1ActAO_FA(i)
            end do
            write(6,*) 'D1Act before tractl'
            do i=1,nacpar
              write(6,*) D1Act(i)
            end do
            write(6,*) 'D1Actao before tractl'
            do i=1,ntot2
              write(6,*) D1Actao(i)
            end do
            write(6,*) 'd1i before tractl'
            do i=1,ntot2
              write(6,*) d1i(i)
            end do
        end if
*
        Call mma_allocate(FockI_Save,ntot1,Label='FockI_Save')
        if(iprlev.ge.debug) then
            write(6,*) 'ifocki before tractl'
            do i=1,ntot1
              write(6,*) FockI(i)
            end do
        end if
*
        call  dcopy_(ntot1,FockI,1,FockI_Save,1)
*
        if (iprlev.ge.debug) then
             write(6,*) 'focki_save before tractl'
             do i=1,ntot1
               write(6,*) focki_save(i)
             end do

             write(6,*) 'FockA before tractl'
             do i=1,ntot1
               write(6,*) FockA(i)
             end do
         end if
*
      Call mma_allocate(tuvx_tmp,nacpr2,Label='TUVX_tmp')
      Call mma_allocate(puvx_tmp,nfint,Label='PUVX_tmp')
      TUVX_tmp(:)=0.0D0
      PUVX_tmp(:)=0.0D0
*
      if (iprlev.ge.debug) then
            write(6,*) 'tuvx_tmp before !!! tractl'
            do i=1, nacpr2
              write(6,*) tuvx_tmp(i)
            end do
      end if
*
      if (iprlev.ge.debug) then
            write(6,*) 'puvx before tractl'
            do i=1,nfint
              write(6,*) puvx_tmp(i)
            end do
      end if
*
*
         CALL TRACTL2(CMO_X,PUVX_tmp,TUVX_tmp,D1I,FockI,D1ActAO,FockA,
     &                IPR,lSquare,ExFac)
*        Call dcopy_(ntot1,FA,1,FockA,1)
        if (iprlev.ge.debug) then
             write(6,*) 'FA tractl msctl'
             call wrtmat(FockA,1,ntot1,1,ntot1)
             write(6,*) 'FI tractl msctl'
             call wrtmat(FockI,1,ntot1,1,ntot1)
        end if
*
      Call mma_deallocate(tuvx_tmp)
      Call mma_deallocate(puvx_tmp)
*
        if(iprlev.ge.debug) then
            write(6,*) 'cmo after tractl'
            do i=1,ntot2
              write(6,*) CMO_X(i)
            end do
            write(6,*) 'puvx after tractl'
            do i=1,nfint
              write(6,*) puvx(i)
            end do
            write(6,*) 'tuvx after tractl'
            do i=1,nacpr2
              write(6,*) tuvx(i)
            end do
           write(6,*) 'D1Act_FA after tractl'
            do i=1,nacpar
              write(6,*) D1Act_FA(i)
            end do
            write(6,*) 'D1ActAO_FA after tractl'
            do i=1,ntot2
              write(6,*) D1ActAO_FA(i)
            end do
            write(6,*) 'D1Act after tractl'
            do i=1,nacpar
              write(6,*) D1Act(i)
            end do
            write(6,*) 'D1Actao after tractl'
            do i=1,ntot2
              write(6,*) D1Actao(i)
            end do
            write(6,*) 'd1i before tractl'
            do i=1,ntot2
              write(6,*) d1i(i)
            end do
        end if

        if(iprlev.ge.debug) then
            write(6,*) 'ifocki after tractl'
            do i=1,ntot1
              write(6,*) FockI(i)
            end do
        end if

        if (iprlev.ge.debug) then
             write(6,*) 'focki_save after tractl'
             do i=1,ntot1
               write(6,*) focki_save(i)
             end do
*
             write(6,*) 'focka after tractl'
             do i=1,ntot1
               write(6,*) FockA(i)
             end do
        end  if
*
*       if (jroot.ne.irlxroot) then
*          Call dcopy_(ntot1,FI,1,FockI,1)
*          Call dcopy_(ntot1,FA,1,FockA,1)
*       end if
*
        Call Fmat_m(CMO,PUVX,D1Act,D1ActAO,FockI_save,FockA)
        call  dcopy_(ntot1,focki_save,1,FockI,1)
*       call  dcopy_(nacpar,D1Act_FA,1,D1Act,1)
        Call mma_deallocate(FockI_Save)
        Call mma_deallocate(CMO_X)
        Call mma_deallocate(D1Act_FA)
        Call mma_deallocate(D1ActAO_FA)
!
!
        if (iprlev.ge.debug) then
             write(6,*) 'FA after fmat 1'
             call wrtmat(FockA,1,ntot1,1,ntot1)
             write(6,*) 'FI after fmat 1'
             call wrtmat(FockI,1,ntot1,1,ntot1)
        end if
*      end if
*TRS
******
*
      if (iprlev.ge.debug) then
           write(6,*) 'D1Act after copy in tractl'
            do i=1,nacpar
              write(6,*) D1Act(i)
            end do
      end if
*
         IF(ISTORP(NSYM+1).GT.0) THEN
           CALL mma_allocate(P,ISTORP(NSYM+1),Label='P')
           CALL DmatDmat(D1Act,P)
           CALL mma_allocate(P1,ISTORP(NSYM+1),Label='P1')
           CALL PMAT_RASSCF(P2d,P1)
         Else
           CALL mma_allocate(P,1,Label='P')
           CALL mma_allocate(P1,1,Label='P1')
         END IF

       if(iprlev.ge.debug) then
         write(6,*) 'dmatdmat'
         do i=1,istorp(nsym+1)
           write(6,*) P(i)
         end do
       end if
         NQ=0
         NSXS=0
         NIAIA=0
         do ISYM=1,NSYM
           NQ = MAX(NQ,NASH(ISYM)*NORB(ISYM))
           NSXS= NSXS+(NISH(ISYM)+NASH(ISYM))*(NASH(ISYM)+NSSH(ISYM))
           NIAIA = NIAIA+(NASH(ISYM)+NISH(ISYM))**2
         end do
         if(NQ.lt.NIAIA) NQ=NIAIA

         CALL mma_allocate(FOCK,NTOT4,Label='FOCK')
         CALL mma_allocate(BM,NSXS,Label='BM')
         CALL mma_allocate(Q,NQ,Label='Q') ! q-matrix(1symmblock)
         IFINAL = 1
         CALL FOCK_m(FOCK,BM,FockI,FockA,D1Act,P,Q,PUVX,IFINAL,CMO)
!TMP TEST
!         Call Put_Darray('fock_tempo',FockOcc,ntot1)
!END TMP TEST


         CASDFT_Funct = 0
         Call Get_dScalar('CASDFT energy',CASDFT_Funct)

         CASDFT_E = ECAS+CASDFT_Funct

        IF(mcpdft_options%otfnal%is_hybrid()) THEN
            E_NoHyb=CASDFT_E
            CASDFT_E = mcpdft_options%otfnal%lambda*Ref_Ener(jRoot) +
     &            (1.0-mcpdft_options%otfnal%lambda) * E_NoHyb
        END IF

        Call Print_MCPDFT_2(CASDFT_E,PotNuc,EMY,ECAS,CASDFT_Funct,
     &         jroot,Ref_Ener)


         IF(Do_Rotate) Then
            Energies(jroot)=CASDFT_E
*JB         replacing ref_ener with MC-PDFT energy for MS-PDFT use
            Ref_Ener(jroot)=CASDFT_E
         ELSE
            Energies(jroot)=CASDFT_E
            ener(jroot,1)=CASDFT_E
         END IF


         Call mma_deallocate(BM)
         Call mma_deallocate(Q)
!At this point, the energy calculation is done.  Now I need to build the
!fock matrix if this root corresponds to the relaxation root.

!***********************************************************************
*
*            BUILDING OF THE NEW FOCK MATRIX                           *
*
************************************************************************
      if(mcpdft_options%grad .and. (.not. mcpdft_options%mspdft)
     &   .and. jroot .eq. irlxroot) then

         Write(LF,*) 'Calculating potentials for analytic gradients...'
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

!
!      write(6,*) 'NACPAR (input fock)',nacpar
!      write(6,*) 'ntot1 (# of V, fock_occ)',ntot1
!      write(6,*) 'nfint (# of v)',nfint
cPS         call xflush(6)

!I will read in the one- and two-electron potentials here

      Call mma_allocate(ONTOPT,nfint,Label='OnTopT')
      OnTopT(:)=0.0D0
      Call mma_allocate(ONTOPO,ntot1,Label='OnTopO')
      OnTopO(:)=0.0D0


      Call Get_dArray('ONTOPT',OnTopT,NFINT)
      Call Get_dArray('ONTOPO',OnTopO,NTOT1)
!
        If ( IPRLEV.ge.DEBUG ) then
        write(6,*) 'One-electron potentials'
        do i=1,ntot1
          write(6,*) OnTopO(i)
        end do
        write(6,*) 'Two-electron potentials'
        do i=1,nfint
          if (abs(puvx(i)).ge.1d-10)then
            write(6,*) OnTopT(i),puvx(i)
          end if
        end do
        end if


!______________________________________________________
!Grab the active-active part of the FI+FA matrix (currently held in the
!FA matrix) and place it in an array of size NACPAR.  Add the oeotp to
!it.  Write to file.
        If ( IPRLEV.ge.DEBUG ) then
      write(6,*) "FA+FI to send to MCLR"
      do i=1,Ntot1
        write(6,*) FockA(i)
      end do
        end if

      Call mma_allocate(Fone,NTOT1,Label='FOne')
      FOne(:)=0.0D0

!I think I need to generate FI, which will contain both the one-electron
!potential contribution and the V_kkpu contribution.



      CALL mma_allocate(FI_V,Ntot1,Label='FI_V')
      Call Get_dArray('FI_V',FI_V,NTOT1)
!         Call Dscal_(nTOT1,4.0d0,FI_V,1)
*         write(6,*) 'fiv after tractl'
*         do i=1,ntot1
*           write(6,*) FI_V(i)
*         end do

      !Call daxpy_(ntot1,0.5d0,FI_V,1,FockA,1)
      Call daxpy_(ntot1,1.0d0,FI_V,1,FockA,1)
      Call daxpy_(ntot1,1.0d0,OnTopO,1,FockA,1)


      i_off1=1
      i_off2=1


      Do iSym = 1,nSym
        iBas = nBas(iSym)
!        iAct = nAsh(iSym)
!        iIsh = nIsh(iSym)

      !FI + FA + V_oe

        do i=1,iBas
          do j=1,i
!            if (i.gt.iIsh.and.j.gt.iIsh) then
!              if(i.le.iAct+iIsh.and.j.le.iAct+iIsh) then
                Fone(i_off1) = fone(i_off1) + FockA(i_off2)
                i_off1 = i_off1 + 1
!              end if
!            end if
            i_off2 = i_off2 + 1
          end do
        end do
      end do
        If ( IPRLEV.ge.DEBUG ) then
      write(6,*) 'F1 to send'
      do i=1,NTOT1
        write(6,*) Fone(i)
      end do
        end if

      !Add the V_kktu contribution to Fone_tu?
!STILL MUST DO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This should be addressed in the upd_FI routine.

      !Write to file.
      LUTMP=87
      LUTMP=IsFreeUnit(LUTMP)
      Call Molcas_Open(LUTMP,'TmpFock')
      do i=1,ntot1
        write(LUTMP,*) Fone(i)
      end do

!Write the TUVX teotp to file.q

!      Call mma_allocate(TUVX_tmp,NACPR2,Label='TUVX_tmp')
!      Call Get_TUVX(PUVX,TUVX_tmp)
!      write(6,*) 'TUVX'
!      do i=1,nacpr2
!      write(6,*) TUVX_tmp(i)
!      end do
!      Call mma_deallocate(TUVX_tmp,NACPR2)


      Call mma_allocate(TUVX_tmp,NACPR2,Label='TUVX_tmp')
      TUVX_tmp(:)=0.0D0
      Call Get_TUVX(OnTopT,TUVX_tmp)
      !Call Get_TUVX(puvx,TUVX_tmp)

      !Unpack TUVX to size
      do i=1,nacpr2
        write(LUTMP,*) TUVX_tmp(i)
      end do
      Close(LUTMP)
      Call mma_deallocate(FOne)
      Call mma_deallocate(TUVX_tmp)

!____________________________________________________________
!This next part is to generate the MC-PDFT generalized fock operator.


!      CALL DCOPY_(nFint,[0.0D0],0,OnTopT,1)
!      CALL DCOPY_(ntot1,[0.0D0],0,OnTopO,1)
!        write(6,*) 'ONTOPT'
!        call wrtmat(OnTopT,1,nFInt,1,nFInt)
!        write(6,*) 'ONTOPO'
!        call wrtmat(OnTopO,1,ntot1,1,ntot1)

!Zero out the matrices.  We will be adding the potential-containing
!terms as a correction to the Focc component already on the runfile.
      CALL DCOPY_(ntot1,[0.0D0],0,FockA,1)
      CALL DCOPY_(ntot1,[0.0D0],0,FockI,1)

!The corrections (from the potentials) to FI and FA are built in the NQ
!part of the code, for efficiency's sake.  It still needs to be
!debugged.
      CALL mma_allocate(FA_V,Ntot1,Label='FA_V')
      Call Get_dArray('FA_V',FA_V,NTOT1)
!         Call Dscal_(nTOT1,4.0d0,FA_V),1)
      !FA_V(2) = 0d0


!         Call Dscal_(ntot1,0.5d0,FI_V,1)
!         Call Dscal_(ntot1,0.5d0,FA_V,1)
      !Call Dscal_(nfint,-0.5d0,OnTopT,1)
        If ( IPRLEV.ge.DEBUG ) then
      write(6,*) "extra terms to update FI"
      do i=1,ntot1
      write(6,*) FI_V(i)
      end do
      !Call TriPrt(' ','(5G18.10)',FI_V,norb(1))
      write(6,*) "extra terms to update FA"
      do i=1,ntot1
      write(6,*) FA_V(i)
      end do
      !Call TriPrt(' ','(5G18.10)',FA_V,norb(1))
        end if

        If ( IPRLEV.ge.DEBUG ) then
      CALL mma_allocate(FA_t,Ntot1,Label='FA_t')
      FA_t(:)=0.0D0
      Call DaXpY_(NTOT1,1.0D0,OnTopO,1,FA_t,1)
      Call daxpy_(NTOT1,1.0D0,FI_V,1,FA_t,1)
      Call daxpy_(NTOT1,1.0D0,FA_V,1,FA_t,1)
      write(6,*) "Total F additions:"
      Call TriPrt(' ','(5G18.10)',FA_t,norb(1))
      CALL mma_deallocate(FA_t)
        end if



      !Add one e potential, too.
      Call DaXpY_(NTOT1,1.0D0,OnTopO,1,FockI,1)
      !Add two e potentials
      Call daxpy_(NTOT1,1.0D0,FI_V,1,FockI,1)
      Call daxpy_(NTOT1,1.0D0,FA_V,1,FockA,1)
        If ( IPRLEV.ge.DEBUG ) then
      write(6,*) "new FI"
      Call TriPrt(' ','(5G18.10)',FockI,norb(1))
      write(6,*) "new FA"
      Call TriPrt(' ','(5G18.10)',FockA,norb(1))
        end if

      CALL mma_deallocate(FI_V)
      CALL mma_deallocate(FA_V)

       IF(ISTORP(NSYM+1).GT.0) THEN
         P(:)=0.0D0
!p = P2d
         CALL PMAT_RASSCF(P2d,P)
      END IF

!Must add to existing FOCK operator (occ/act). FOCK is not empty.
         CALL mma_allocate(BM,NSXS,Label='BM')
         CALL mma_allocate(Q,NQ,Label='Q') ! q-matrix(1symmblock)
         CALL FOCK_update(FOCK,BM,FockI,FockA,D1Act,P,
     &                    Q,OnTopT,IFINAL,CMO)

         Call Put_dArray('FockOcc',FockOcc,ntot1)
        If ( IPRLEV.ge.DEBUG ) then
        write(6,*) 'FOCC_OCC'
        call wrtmat(FockOcc,1,ntot1,1,ntot1)


      write(6,*) 'DONE WITH NEW FOCK OPERATOR'
        end if

         Call mma_deallocate(BM)
         Call mma_deallocate(Q)
      Call mma_deallocate(ONTOPO)
      Call mma_deallocate(ONTOPT)


!Put some information on the runfile for possible gradient calculations.
      Call Put_iScalar('Number of roots',nroots)
      Call Put_dArray('Last energies',Energies,nroots)
      Call Put_cArray('Relax Method','MCPDFT  ',8)
      Call Put_dScalar('Last energy',Energies(iRlxRoot))
      iSA = 1

!      Call mma_allocate(P2dt1,NACPR2,Label='P2dt1')
!      P2dt1(:)=0.0D0
!      !I need the non-symmetry blocked D1Act, hence the read.
!      Call mma_allocate(D1Act1,NACPAR,Label='D1Act1')
!      Call Get_dArray_chk('D1mo',D1Act1,NACPAR)
!        write(6,*) 'D1Act'
!        do i=1,NACPAR
!          write(6,*) D1Act1(i)
!        end do
!      Call P2_contraction(D1Act1,P2dt1)
!      Call Put_dArray('P2MOt',P2dt1,NACPR2)
!      Call mma_deallocate(P2dt1)
!      Call mma_deallocate(D1Act1)

!Put information needed for geometry optimizations.
      !if (jroot.eq.iRlxRoot) then
          iSA = 1 !need to do MCLR for gradient runs. (1 to run, 2 to
*skip)
       !MUST MODIFY THIS.  I need to check that the calculation is not
       !SA, and if it is, set iSA to -1.
      Call Put_iScalar('SA ready',iSA)
      Call Put_cArray('MCLR Root','****************',16)
      Call Put_iScalar('Relax CASSCF root',irlxroot)
      !end if



      end if

      if (mcpdft_options%grad .and. mcpdft_options%mspdft) then
*      doing exactly the same thing as done in the previous chunck
*      starting from 'BUILDING OF THE NEW FOCK MATRIX'
*      Hopefully this code will be neater.
       call savefock_pdft(CMO,FockI,FockA,D1Act,Fock,
     &                    P,NQ,PUVX,p2d,jroot)
      end if



      !if (jroot.eq.iRlxRoot) then
      Call Put_iScalar('Number of roots',nroots)
      Call Put_dArray('Last energies',Energies,nroots)
      Call Put_cArray('Relax Method','MCPDFT  ',8)
      Call Put_dScalar('Last energy',Energies(iRlxRoot))
      !end if

      Call mma_deallocate(Tmp2)
      Call mma_deallocate(Tmpa)
      Call mma_deallocate(FOCK)
      CALL mma_deallocate(P)
      CALL mma_deallocate(P1)
      end do !loop over roots

      if(mcpdft_options%grad) then
        dmDisk = IADR19(3)
        do jroot=1,irlxroot-1
          Call DDaFile(JOBOLD,0,D1Act,NACPAR,dmDisk)
          Call DDaFile(JOBOLD,0,D1Spin,NACPAR,dmDisk)
          Call DDaFile(JOBOLD,0,P2d,NACPR2,dmDisk)
          Call DDaFile(JOBOLD,0,P2d,NACPR2,dmDisk)
        end do
        Call DDaFile(JOBOLD,2,D1Act,NACPAR,dmDisk)
*        Andrew added this line to fix heh2plus
        Call DDaFile(JOBOLD,2,D1Spin,NACPAR,dmDisk)
        Call Put_dArray('D1mo',D1Act,NACPAR)
*        write(6,*) 'd1Spin'
*        do i=1,NACPAR
*          write(6,*) D1spin(i)
*        end do
*TRS commenting out because we already read over this
*        Call DDaFile(JOBOLD,0,D1Spin,NACPAR,dmDisk)
        Call DDaFile(JOBOLD,2,P2d,NACPR2,dmDisk)
        Call Put_dArray('P2mo',P2d,NACPR2)
*        write(6,*) 'D2'
*        do i=1,NACPR2
*          write(6,*) p2d(i)
*        end do

         If(NASH(1).ne.NAC) Call DBLOCK(D1Act)
         Call Get_D1A_RASSCF(CMO,D1Act,D1ActAO)

         Call Fold(nSym,nBas,D1I,Tmp3)
         Call Fold(nSym,nBas,D1ActAO,Tmp4)
         Call Daxpy_(nTot1,1.0D0,Tmp4,1,Tmp3,1)
         Call Put_dArray('D1ao',Tmp3,nTot1)
!         write(6,*) 'd1ao'
!         do i=1,ntot1
!           write(6,*) tmp3(i)
!         end do

!Get the spin density matrix for open shell cases
***********************************************************
* Generate spin-density
***********************************************************
*TRS ams also commented out this if and endif part of this
* statement
         !if(iSpin.eq.1) then
           Call dcopy_(NACPAR,[Zero],0,D1SpinAO,1)
         !end if
         IF ( NASH(1).NE.NAC ) CALL DBLOCK(D1Spin)
         Call Get_D1A_RASSCF(CMO,D1Spin,D1SpinAO)
         Call mma_allocate(Tmp7,nTot1,Label='Tmp7')
         Call Fold(nSym,nBas,D1SpinAO,Tmp7)
         Call Put_dArray('D1Sao',Tmp7,nTot1)
!         write(6,*) 'd1so'
!         do i=1,ntot1
!           write(6,*) tmp7(i)
!         end do
         Call mma_deallocate(Tmp7)




        Call DDaFile(JOBOLD,0,P2d,NACPR2,dmDisk)
      end if

      Call mma_deallocate(PUVX)

!Free up all the memory we can here, eh?
      Call mma_deallocate(Tmp4)
      Call mma_deallocate(Tmp3)

      Call mma_deallocate(D1Act)
      Call mma_deallocate(D1ActAO)
      Call mma_deallocate(D1Spin)
      Call mma_deallocate(D1SpinAO)
      Call mma_deallocate(Tmp1)
      Call mma_deallocate(FockI)
      Call mma_deallocate(FockA)
*      Call DDaFile(JOBOLD,0,P2d,NACPR2,dmDisk)
      Call mma_deallocate(P2D)
      Call mma_deallocate(D1I)
      Call mma_deallocate(Tmpk)
      Call mma_deallocate(Tmpn)
c      call xflush(6)

      END Subroutine MSCtl

      Subroutine P2_contraction(D1MO,P2MO)
      use definitions, only: wp
      implicit none

      real(kind=wp), dimension(*), intent(in) :: d1mo
      real(kind=wp), dimension(*), intent(out) :: p2mo

#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
      integer :: i, j, k, l, ij, kl, ijkl, lmax
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
              fact=1.0d0
              if(k == l) fact=0.5d0
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


      Subroutine Fold_pdft(nSym,nBas,A,B)

      Implicit Real*8 (A-H,O-Z)

      Dimension nBas(*) , A(*) , B(*)

      iOff1 = 0
      iOff2 = 0
      Do iSym = 1, nSym
        mBas = nBas(iSym)
        Do iBas= 1, mBas
          Do jBas = 1 , iBas-1
            B(iOff2+jBas) =   A(iOff1+jBas)
          End Do
          B(iOff2+iBas) =  A(iOff1+iBas)
          iOff1 = iOff1 + mBas
          iOff2 = iOff2 + iBas
        End Do
      End Do

      Return
      end Subroutine Fold_pdft
