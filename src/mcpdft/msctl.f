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
      Subroutine MSCtl(CMO,F,FI,FA,Ref_Ener)
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
      Implicit Real*8 (A-H,O-Z)
      Dimension CMO(*) ,F(*) , FI(*), FA(*), Ref_Ener(*)
*
#include "rasdim.fh"
#include "general.fh"
#include "input_ras.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='MSCTL   ')
#include "rasscf.fh"
#include "WrkSpc.fh"
#include "rctfld.fh"
#include "pamint.fh"
#include "timers.fh"
#include "SysDef.fh"
#include "csfbas.fh"
#include "gugx.fh"
#include "casvb.fh"
#include "wadr.fh"
#include "rasscf_lucia.fh"
#include "raswfn.fh"
#include "ksdft.fh"
      Logical DoActive,DoQmat,DoCholesky
!      Logical TraOnly
      Integer ALGO

      COMMON /CHOTODO /DoActive,DoQmat,ipQmat
      COMMON /CHLCAS /DoCholesky,ALGO

*
      Character*8 Label
      Logical First, Dff, Do_DFT,Found
      Parameter ( Zero=0.0d0 , One=1.0d0 )
      integer iD1I,iD1Act,iD1ActAO,iD1Spin,iD1SpinAO,IAD19
      integer iJOB,dmDisk,iP2d
      integer itmp0,itmp1,itmp2,itmp3,itmp4
      integer itmp5,itmp6,itmp7
      integer ifocki,ifocka
      integer IADR19(1:30)
      integer LP,NQ,LQ,LPUVX
      integer  LOEOTP,NACP,NACP2
      integer vdisk,jroot
      real*8,dimension(1:nroots) :: Energies
      integer count_tmp1,count_tmp2
      integer  i_off1,i_off2,ifone
      integer isym,iorb,iash,iish,jsym
      iTrii(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)


      Call qEnter('MSCTL')
      Call unused_real_array(F)
***********************************************************
C Local print level (if any)
***********************************************************
      IPRLEV=IPRLOC(3)


***********************************************************
* Load the nuclear repulsion energy
***********************************************************

      Call Get_dScalar('PotNuc',potNuc)

***********************************************************
* Generate molecular charges
***********************************************************
      Call GetMem('Ovrlp','Allo','Real',iTmp0,nTot1+4)
      iRc=-1
      iOpt=2
      iComp=1
      iSyLbl=1
      Label='Mltpl  0'
      Call RdOne(iRc,iOpt,Label,iComp,Work(iTmp0),iSyLbl)
      Tot_Nuc_Charge=Work(iTmp0+nTot1+3)
      If ( iRc.ne.0 ) then
         Write(LF,*) 'CASDFT_Terms: iRc from Call RdOne not 0'
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


***********************************************************
* Load bare nuclei Hamiltonian
***********************************************************
      Call GetMem('Fcore','Allo','Real',iTmp1,nTot1)
      iComp  =  1
      iSyLbl =  1
      iRc    = -1
      iOpt   =  6
      Label  = 'OneHam  '
      Call RdOne(iRc,iOpt,Label,iComp,Work(iTmp1),iSyLbl)
      If ( iRc.ne.0 ) then
         Write(LF,*) 'CASDFT_Terms: iRc from Call RdOne not 0'
         Write(LF,*) 'Label = ',Label
         Write(LF,*) 'iRc = ',iRc
         Call QTrace
         Call Abend
      Endif
      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' OneHam in AO basis in CASDFT_Terms'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=0
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',Work(iTmp1+iOff),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
        End Do
      End If

!Do I need this Kincore and NucElcore?
      if(iPrLev.ge.DEBUG) then
      Call GetMem('Kincore','Allo','Real',iTmpk,nTot1)
      iComp  =  1
      iSyLbl =  1
      iRc    = -1
      iOpt   =  6
      Label  = 'Kinetic '
      Call RdOne(iRc,iOpt,Label,iComp,Work(iTmpk),iSyLbl)
      If ( iRc.ne.0 ) then
         Write(LF,*) 'CASDFT_Terms: iRc from Call RdOne not 0'
         Write(LF,*) 'Label = ',Label
         Write(LF,*) 'iRc = ',iRc
         Call QTrace
         Call Abend
      Endif
      Call GetMem('NucElcore','Allo','Real',iTmpn,nTot1)
      iComp  =  1
      iSyLbl =  1
      iRc    = -1
      iOpt   =  6
      Label  = 'Attract '
      Call RdOne(iRc,iOpt,Label,iComp,Work(iTmpn),iSyLbl)
      If ( iRc.ne.0 ) then
         Write(LF,*) 'CASDFT_Terms: iRc from Call RdOne not 0'
         Write(LF,*) 'Label = ',Label
         Write(LF,*) 'iRc = ',iRc
         Call QTrace
         Call Abend
      Endif
      Call GetMem('Kincore','free','Real',iTmpk,nTot1)
      Call GetMem('NucElcore','free','Real',iTmpn,nTot1)
      end if


!Here we calculate the D1 Inactive matrix (AO).
      Call GetMem('D1Inact','Allo','Real',iD1I,NTOT2)
      Call Get_D1I_RASSCF_m(CMO,Work(iD1I))

      IF(IPRLEV.ge.DEBUG) THEN
         write(6,*) 'iD1inact'
         do i=1,ntot2
         write(6,*) Work(iD1i-1+i)
         end do
      END IF

        iJOB=0
        IAD19=0
        Call f_Inquire('JOBOLD',Found)
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
       vDisk =  IADR19(4)
       dmDisk = IADR19(3)
!       Call GetMem('jVector','Allo','Real',ijVec,nConf)
       Call GetMem('D1Active','Allo','Real',iD1Act,NACPAR)
       Call GetMem('D1ActiveAO','Allo','Real',iD1ActAO,NTOT2)
       Call GetMem('D1Spin','Allo','Real',iD1Spin,NACPAR)
       Call GetMem('D1SpinAO','Allo','Real',iD1SpinAO,NTOT2)

       Call GetMem('DtmpI','Allo','Real',iTmp3,nTot1)
       Call GetMem('DtmpA','Allo','Real',iTmp4,nTot1)
       Call GetMem('P2','Allo','Real',iP2d,NACPR2)


       Call GetMem('FockI','ALLO','Real',ifocki,ntot1)
       Call GetMem('FockA','ALLO','Real',ifocka,ntot1)

************************************************************************
* load back two-electron integrals (pu!vx)
************************************************************************
      lPUVX = 1
      If (.not.DoCholesky .or. ALGO.eq.1) Then
        If ( nFint.gt.0) then
          iDisk = 0
          Call GetMem('PUVX','Allo','Real',lPUVX,nFint)
          Call DDaFile(LUINTM,2,Work(lPUVX),nFint,iDisk)
        End If
      End If
      IF(IPRLEV.ge.DEBUG) THEN
        write(6,*) 'PUVX integrals in SXCTL'
        call wrtmat(Work(lPUVX),1,nFInt,1,nFInt)
      END IF


!This iSA is used to control gradient calculations.  Analytic gradients
!(in ALASKA) will only run if iSA=1, and iSA will only be set to one if
!the on-top potentials are computed as part of this calculation.
          iSA = 99
      Call Put_iScalar('SA ready',iSA)


!We loop over the number of states.  For each state, we read the density
!matrices from the JOBIPH file.
       do jroot=1,lroots
!       do d_off=1,1!100!,43!33,33!1,43

       !Load a fresh FockI and FockA
         Call dcopy_(ntot1,FI,1,Work(ifocki),1)
         Call dcopy_(ntot1,FA,1,Work(ifocka),1)


!Read in the density matrices for <jroot>.
!         Call Fzero(Work(ijVec),nConf)
         Call Fzero(Work(iD1Act),NACPAR)
         Call Fzero(Work(iD1ActAO),NTOT2)
         Call Fzero(Work(iD1Spin),NACPAR)
         Call Fzero(Work(iD1SpinAO),NTOT2)
         Call Fzero(Work(iTmp3),nTot1)
         Call Fzero(Work(iTmp4),nTot1)
         Call Fzero(Work(iP2d),NACPR2)
!         Call DDaFile(JOBIPH,2,Work(ijVec),nConf,jDisk)
!Get the D1 Active matrix for this state.  These should probably be
!most easily read from the previous JOBIPH file.  Then, convert D1A from
!the MO to the AO basis.

         Call DDaFile(JOBOLD,2,Work(iD1Act),NACPAR,dmDisk)

      if(.false.) then
         open(unit=90,File='DMs.out',action='read')
         do i=1,NACPAR
           read(90,*) Work(id1Act-1+i)
         end do
      end if


      IF(IPRLEV.ge.DEBUG) THEN
        write(6,*) 'd1act'
        do i=1,NACPAR
          write(6,*) work(iD1Act-1+i)
        end do
      end if

         Call Put_D1MO(Work(iD1Act),NACPAR)
         Call DDaFile(JOBOLD,2,Work(iD1Spin),NACPAR,dmDisk)
         Call DDaFile(JOBOLD,2,Work(iP2d),NACPR2,dmDisk)
         Call Put_P2MO(Work(iP2d),NACPR2)

         Call DDaFile(JOBOLD,0,Work(iP2d),NACPR2,dmDisk)
      IF(IPRLEV.ge.DEBUG) THEN
        write(6,*) 'D2'
        do i=1,NACPR2
          write(6,*) Work(ip2d-1+i)
        end do
      END IF



***********************************************************
* Generate total density
***********************************************************

         If(NASH(1).ne.NAC) Call DBLOCK_m(Work(iD1Act))
      IF(IPRLEV.ge.DEBUG) THEN
        write(6,*) 'd1act'
        do i=1,NACPAR
          write(6,*) work(iD1Act-1+i)
        end do
      end if
         Call Get_D1A_RASSCF_m(CMO,Work(iD1Act),Work(iD1ActAO))
!         write(6,*) "is this it?",ntot2
!         do i=1,ntot2
!           write,*) Work(iD1ActAO-1+i)
!         end do

         Call Fold(nSym,nBas,Work(iD1I),Work(iTmp3))
         Call Fold(nSym,nBas,Work(iD1ActAO),Work(iTmp4))
         Call Daxpy_(nTot1,1.0D0,Work(iTmp4),1,Work(iTmp3),1)
!Maybe I can write all of these matrices to file, then modify stuff in
!the nq code to read in the needed density.  In other words, I need to
!replace the next call with something that supports multistates.
         Call Put_D1ao(Work(iTmp3),nTot1)
      IF(IPRLEV.ge.DEBUG) THEN
         write(6,*) 'd1ao'
         do i=1,ntot1
           write(6,*) work(itmp3-1+i)
         end do
         call xflush(6)
      end if

!Get the spin density matrix for open shell cases
***********************************************************
* Generate spin-density
***********************************************************
         if(iSpin.eq.1) then
           Call dcopy_(NACPAR,0.0d0,0,Work(iD1SpinAO),1)
         end if
         IF ( NASH(1).NE.NAC ) CALL DBLOCK_m(Work(iD1Spin))
         Call Get_D1A_RASSCF_m(CMO,Work(iD1Spin),
     &                      Work(iD1SpinAO))
         Call GetMem('DtmpS','Allo','Real',iTmp7,nTot1)
         Call Fold(nSym,nBas,Work(iD1SpinAO),Work(iTmp7))
         Call Put_D1Sao(Work(iTmp7),nTot1)
         Call GetMem('DtmpS','Free','Real',iTmp7,nTot1)


***********************************************************
* Calculation of the CASDFT_energy
***********************************************************
!Perhaps ideally, we should rework how DrvXV (and its children) handles
!the AO to MO transformation on the grid.  It seems like perhaps we are
!doing redundant transformations by retransforming AOs (which may have
!been included in a previous batch) into MOs.
*
        Call GetMem('htmp','Allo','Real',iTmp5,nTot1)
        Call GetMem('gtmp','Allo','Real',iTmp6,nTot1)
        Call dCopy_(nTot1,0.0d0,0,Work(iTmp5),1)
        Call dCopy_(nTot1,0.0d0,0,Work(iTmp6),1)
*
        First=.True.
        Dff=.False.
        Do_DFT=.True.

        Call Put_iArray('nFro',nFro,nSym)
        Call Put_iArray('nAsh',nAsh,nSym)
        Call Put_iArray('nIsh',nIsh,nSym)

         call xflush(6)

        iCharge=Int(Tot_Charge)

c iTmp5 and iTmp6 are not updated in DrvXV...
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
          NACP=(NAC+NAC**2)/2
          NACP2=(NACP+NACP**2)/2
          NCEH2=1

      do_pdftPot=.false.
      if(DoGradPDFT) then

        do_pdftPot=.true.

        CALL GETMEM('TE_POTG','ALLO','REAL',LTEOTPG,NFINT)
        CALL GETMEM('OE_POT','ALLO','REAL',LOEOTP,NTOT1)
        Call DCopy_(NTOT1,0.0d0,0,Work(LOEOTP),1)
        Call DCopy_(NFINT,0.0d0,0,work(LTEOTPG),1)

      !preallocate the runfile stuff for inact-containing potentials.
        count_tmp1 = 0
        count_tmp2 = 0
        do isym=1,nsym
          do jsym=1,nsym
         count_tmp2 = count_tmp2 + nIsh(isym)*(nIsh(jsym)+nAsh(jsym))**2
          end do
          count_tmp1 = count_tmp1 + nIsh(isym)*(nish(isym)+nAsh(isym))
        end do
        CALL GETMEM('JUNK','ALLO','REAL',ijunk,count_tmp2)
        Call DCopy_(count_tmp2,0.0d0,0,Work(ijunk),1)
!        Call Put_dArray('TEP_I',work(ijunk),count_tmp2)
        CALL GETMEM('JUNK','Free','REAL',ijunk,count_tmp2)

        CALL GETMEM('JUNK','ALLO','REAL',ijunk,count_tmp1)
        Call DCopy_(count_tmp1,0.0d0,0,Work(ijunk),1)
!        Call Put_dArray('OEP_I',work(ijunk),count_tmp1)
        CALL GETMEM('JUNK','Free','REAL',ijunk,count_tmp1)

        Call Put_dArray('ONTOPO',work(LOEOTP),NTOT1)
        Call Put_dArray('ONTOPT',work(LTEOTPG),NFINT)

        Call DCopy_(NTOT1,0.0d0,0,Work(LOEOTP),1)
        Call DCopy_(NFINT,0.0d0,0,work(LTEOTPG),1)

        Call Put_dArray('FI_V',work(LOEOTP),NTOT1)
        Call Put_dArray('FA_V',work(LOEOTP),NTOT1)

        CALL GETMEM('OE_POT','FREE','REAL',LOEOTP,NTOT1)
        CALL GETMEM('TE_POTG','FREE','REAL',LTEOTPG,NFINT)
      end if !DoGradPdft

        Call DrvXV(Work(iTmp5),Work(iTmp6),Work(iTmp3),
     &             PotNuc,nTot1,First,Dff,NonEq,lRF,
     &             KSDFT_TEMP,ExFac,iCharge,iSpin,
     &             Work(iD1I),Work(iD1ActAO),
     &             nTot1,DFTFOCK,Do_DFT)


        Call Daxpy_(nTot1,1.0d0,Work(iTmp5),1,Work(iTmp1),1)
        Call Daxpy_(nTot1,1.0d0,Work(iTmp6),1,Work(iFockI),1)

        Call GetMem('gtmp','Free','Real',iTmp6,nTot1)
        Call GetMem('htmp','Free','Real',iTmp5,nTot1)

!
***********************************************************
*     Compute energy contributions
***********************************************************
      Call GetMem('DoneI','Allo','Real',iTmp2,nTot1)

      Call Fold(nSym,nBas,Work(iD1I),Work(iTmp2))
*
      Eone = dDot_(nTot1,Work(iTmp2),1,Work(iTmp1),1)
      Call Get_dScalar('PotNuc',PotNuc_Ref)
      Eone = Eone + (PotNuc-PotNuc_Ref)
      Etwo = dDot_(nTot1,Work(iTmp2),1,Work(iFockI),1)
      Call GetMem('DoneI','Free','Real',iTmp2,nTot1)
      EMY  = PotNuc_Ref+Eone+0.5d0*Etwo

      CASDFT_Funct = 0.0D0
      Call Get_dScalar('CASDFT energy',CASDFT_Funct)
      If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*) ' Nuclear repulsion energy :',PotNuc
         Write(LF,*) ' One-electron core energy :',Eone
         Write(LF,*) ' Two-electron core energy :',Etwo
         Write(LF,*) ' Total core energy        :',EMY
         Write(LF,*) ' CASDFT Energy            :',CASDFT_Funct
      End If
         call xflush(6)
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
          Call TriPrt(' ','(5G17.11)',Work(IFockI+iOff-1),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
        End Do
      End If

      Call DaXpY_(nTot1,One,Work(iTmp1),1,Work(iFockI),1)

      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' Inactive Fock matrix in AO basis in CASDFT_terms'
        write(LF,*) '(it already contains OneHam and TwoEl contrib.)'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=1
        Do iSym = 1,nSym
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',Work(iFockI+iOff-1),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
        End Do
      End If


************************************************************************
* update and transform the Fock matrices FI and FA ----> FMAT routine
************************************************************************

      If (.not.DoCholesky .or. ALGO.eq.1) Then
         Call Fmat_m(CMO,Work(lPUVX),Work(iD1Act),Work(iD1ActAO),
     &             Work(iFockI),Work(iFockA))
!
!
        If ( IPRLEV.ge.DEBUG ) then
        write(6,*) 'FA_old'
        call wrtmat(Work(ifocka),1,ntot1,1,ntot1)
        write(6,*) 'FI_old'
        call wrtmat(Work(ifocki),1,ntot1,1,ntot1)
        End if
      Else

         Write(LF,*)'SXCTL: Illegal Cholesky parameter ALGO= ',ALGO
         call qtrace()
         call abend()

      EndIf

         IF(ISTORP(NSYM+1).GT.0) THEN
           CALL GETMEM('ISTRP','ALLO','REAL',LP,ISTORP(NSYM+1))
           CALL DmatDmat_m(Work(iD1Act),WORK(LP))
           CALL GETMEM('ISTRP','ALLO','REAL',LP1,ISTORP(NSYM+1))
           CALL PMAT_RASSCF_M(Work(iP2d),WORK(LP1))
         END IF

       if(iprlev.ge.debug) then
         write(6,*) 'dmatdmat'
         do i=1,istorp(nsym+1)
           write(6,*) Work(LP-1+i)
         end do
       end if
!DANGER!
         call xflush(6)
         NQ=0
         NSXS=0
         NIAIA=0
         do ISYM=1,NSYM
           NQ = MAX(NQ,NASH(ISYM)*NORB(ISYM))
           NSXS= NSXS+(NISH(ISYM)+NASH(ISYM))*(NASH(ISYM)+NSSH(ISYM))
           NIAIA = NIAIA+(NASH(ISYM)+NISH(ISYM))**2
         end do
         if(NQ.lt.NIAIA) NQ=NIAIA
         call xflush(6)

         CALL GETMEM('FOCK','ALLO','REAL',LFOCK,NTOT4)
         CALL GETMEM('SXBM','ALLO','REAL',LBM,NSXS)
         CALL GETMEM('SXLQ','ALLO','REAL',LQ,NQ) ! q-matrix(1symmblock)
         IFINAL = 1
         call xflush(6)
         CALL FOCK_m(WORK(LFOCK),WORK(LBM),Work(iFockI),Work(iFockA),
     &         Work(iD1Act),WORK(LP),WORK(LQ),WORK(LPUVX),IFINAL,CMO)
!         CALL GETMEM('FOCK','FREE','REAL',LFOCK,NTOT4)
!TMP TEST
!         Call Put_Darray('fock_tempo',Work(ipFocc),ntot1)
!END TMP TEST


         CASDFT_Funct = 0
         Call Get_dScalar('CASDFT energy',CASDFT_Funct)

         CASDFT_E = ECAS-EVAC+CASDFT_Funct

!         Write(*,*) '**************************************************'
!         write(*,*) 'ENERGY REPORT FOR STATE',jroot
        Call Print_MCPDFT_2(CASDFT_E,PotNuc,EMY,ECAS,CASDFT_Funct,
     &         jroot,Ref_Ener)
         call xflush(6)

         Energies(jroot)=CASDFT_E

         CALL GETMEM('SXBM','FREE','REAL',LBM,NSXS)
         CALL GETMEM('SXLQ','FREE','REAL',LQ,NQ)
!At this point, the energy calculation is done.  Now I need to build the
!fock matrix if this root corresponds to the relaxation root.


!***********************************************************************
*
*            BUILDING OF THE NEW FOCK MATRIX                           *
*
************************************************************************
        if(DoGradPDFT) then

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
         call xflush(6)

!I will read in the one- and two-electron potentials here

      Call GetMem('ONTOPT','ALLO','Real',ipTmpLTEOTP,nfint)
      Call GetMem('ONTOPO','ALLO','Real',ipTmpLOEOTP,ntot1)
      Call FZero(Work(iptmplteotp),Nfint)
      Call FZero(Work(iptmploeotp),ntot1)


      Call Get_dArray('ONTOPT',work(ipTmpLTEOTP),NFINT)
      Call Get_dArray('ONTOPO',work(ipTmpLOEOTP),NTOT1)
!
        If ( IPRLEV.ge.DEBUG ) then
        write(6,*) 'One-electron potentials'
        do i=1,ntot1
          write(6,*) Work(iptmploeotp-1+i)
        end do
        write(6,*) 'Two-electron potentials'
        do i=1,nfint
          if (abs(work(lpuvx-1+i)).ge.1d-10)then
            write(6,*) Work(iptmplteotp-1+i),work(lpuvx-1+i)
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
        write(6,*) Work(ifocka-1+i)
      end do
        end if

      Call GetMem('F_ONE','ALLO','Real',iFone,NTOT1)
      CALL DCOPY_(NTOT1,0.0D0,0,WORK(iFone),1)

!I think I need to generate FI, which will contain both the one-electron
!potential contribution and the V_kkpu contribution.



      CALL GETMEM('FI_V','ALLO','REAL',ifiv,Ntot1)
      Call Get_dArray('FI_V',work(ifiv),NTOT1)
!         Call Dscal_(nTOT1,4.0d0,Work(ifiv),1)

      !Call daxpy_(ntot1,0.5d0,Work(ifiv),1,Work(iFocka),1)
      Call daxpy_(ntot1,1.0d0,Work(ifiv),1,Work(iFocka),1)
      Call daxpy_(ntot1,1.0d0,Work(iptmploeotp),1,Work(iFocka),1)


      i_off1=0
      i_off2=0


      Do iSym = 1,nSym
        iBas = nBas(iSym)
        iOrb = nOrb(iSym)
        iFro = nFro(iSym)
        iAct = nAsh(iSym)
        iIsh = nIsh(iSym)

      !FI + FA + V_oe

        do i=1,iBas
          do j=1,i
!            if (i.gt.iIsh.and.j.gt.iIsh) then
!              if(i.le.iAct+iIsh.and.j.le.iAct+iIsh) then
                Work(iFone+i_off1) = Work(ifone+i_off1) +
     &              Work(ifocka+i_off2)
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
        write(6,*) work(iFone-1+i)
      end do
        end if

      !Add the V_kktu contribution to Fone_tu?
!STILL MUST DO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This should be addressed in the upd_FI routine.

      !Write to file.
      Open(unit=87,file='TmpFock', action='write',iostat=ios)
      if (ios.ne.0) then
        write(6,*) "error opening file!"
      end if
      do i=1,ntot1
        write(87,*) Work(iFone+i-1)
      end do

!Write the TUVX teotp to file.q

!      Call GetMem('ttTUVX2','Allo','Real',ittTUVX2,NACPR2)
!      Call Get_TUVX(Work(lPUVX),Work(ittTUVX2))
!      write(6,*) 'TUVX'
!      do i=1,nacpr2
!      write(6,*) work(ittTUVX2-1+i)
!      end do
!      Call GetMem('ttTUVX2','Free','Real',ittTUVX2,NACPR2)


      Call GetMem('ttTUVX','Allo','Real',ittTUVX,NACPR2)
      CALL DCOPY_(nacpr2,0.0D0,0,WORK(ittTUVX),1)
      Call Get_TUVX(Work(ipTmpLTEOTP),Work(ittTUVX))
      !Call Get_TUVX(Work(lpuvx),Work(ittTUVX))

      !Unpack TUVX to size
      do i=1,nacpr2
        write(87,*) Work(ittTUVX+i-1)
      end do
      Close(87)
      Call GetMem('F_ONE','Free','Real',iFone,NTOT1)
      Call GetMem('ttTUVX','Free','Real',ittTUVX,NACPR2)


!____________________________________________________________
!This next part is to generate the MC-PDFT generalized fock operator.


!      CALL DCOPY_(nFint,0.0D0,0,WORK(ipTmpLTEOTP),1)
!      CALL DCOPY_(ntot1,0.0D0,0,WORK(ipTmpLOEOTP),1)
!        write(6,*) 'ONTOPT'
!        call wrtmat(Work(ipTmpLTEOTP),1,nFInt,1,nFInt)
!        write(6,*) 'ONTOPO'
!        call wrtmat(Work(ipTmpLOEOTP),1,ntot1,1,ntot1)

!Zero out the matrices.  We will be adding the potential-containing
!terms as a correction to the Focc component already on the runfile.
      CALL DCOPY_(ntot1,0.0D0,0,WORK(iFocka),1)
      CALL DCOPY_(ntot1,0.0D0,0,WORK(iFocki),1)

!The corrections (from the potentials) to FI and FA are built in the NQ
!part of the code, for efficiency's sake.  It still needs to be
!debugged.
      CALL GETMEM('FA_V','ALLO','REAL',ifav,Ntot1)
      Call Get_dArray('FA_V',work(ifav),NTOT1)
!         Call Dscal_(nTOT1,4.0d0,Work(ifav),1)
      !work(ifav-1+2) = 0d0


!         Call Dscal_(ntot1,0.5d0,Work(ifiv),1)
!         Call Dscal_(ntot1,0.5d0,Work(ifav),1)
      !Call Dscal_(nfint,-0.5d0,Work(ipTmpLTEOTP),1)
        If ( IPRLEV.ge.DEBUG ) then
      write(6,*) "extra terms to update FI"
      do i=1,ntot1
      write(6,*) Work(ifiv-1+i)
      end do
      !Call TriPrt(' ','(5G18.10)',Work(ifiv),norb(1))
      write(6,*) "extra terms to update FA"
      do i=1,ntot1
      write(6,*) Work(ifav-1+i)
      end do
      !Call TriPrt(' ','(5G18.10)',Work(ifav),norb(1))
        end if

        If ( IPRLEV.ge.DEBUG ) then
      CALL GETMEM('FA_t','ALLO','REAL',ifat,Ntot1)
      Call dcopy_(ntot1,0.0d0,0,work(ifat),1)
      Call DaXpY_(NTOT1,1.0D0,Work(ipTmpLOEOTP),1,Work(ifat),1)
      Call daxpy_(NTOT1,1.0D0,Work(ifiv),1,Work(ifat),1)
      Call daxpy_(NTOT1,1.0D0,Work(ifav),1,Work(ifat),1)
      write(6,*) "Total F additions:"
      Call TriPrt(' ','(5G18.10)',Work(ifat),norb(1))
      CALL GETMEM('FA_t','free','REAL',ifat,Ntot1)
        end if



      !Add one e potential, too.
!test comment
      Call DaXpY_(NTOT1,1.0D0,Work(ipTmpLOEOTP),1,Work(ifocki),1)
      !Add two e potentials
!test comment
      Call daxpy_(NTOT1,1.0D0,Work(ifiv),1,Work(ifocki),1)
!test comment
      Call daxpy_(NTOT1,1.0D0,Work(ifav),1,Work(ifocka),1)
        If ( IPRLEV.ge.DEBUG ) then
      write(6,*) "new FI"
      Call TriPrt(' ','(5G18.10)',Work(ifocki),norb(1))
      write(6,*) "new FA"
      Call TriPrt(' ','(5G18.10)',Work(ifocka),norb(1))
        end if

      CALL GETMEM('FI_V','Free','REAL',ifiv,Ntot1)
      CALL GETMEM('FA_V','Free','REAL',ifav,Ntot1)

!Reordering of the two-body density matrix.
!      Call Getmem('test_p2','allo','real',ip2test,ISTORP(NSYM+1))
!      Call dcopy_(ISTORP(NSYM+1),Work(LP),1,Work(ip2test),1)

       IF(ISTORP(NSYM+1).GT.0) THEN
      CALL DCOPY_(ISTORP(NSYM+1),0.0D0,0,WORK(LP),1)
!p = work(iP2d)
         CALL PMAT_RASSCF_M(Work(iP2d),WORK(LP))
      END IF
!test comment add
!      Call FZero(Work(iptmplteotp),nfint)
!test comment add end

!Must add to existing FOCK operator (occ/act). FOCK is not empty.
         CALL GETMEM('SXBM','ALLO','REAL',LBM,NSXS)
         CALL GETMEM('SXLQ','ALLO','REAL',LQ,NQ) ! q-matrix(1symmblock)
         CALL FOCK_update(WORK(LFOCK),WORK(LBM),Work(iFockI),
     &        Work(iFockA),Work(iD1Act),WORK(LP),
     &        WORK(LQ),WORK(ipTmpLTEOTP),IFINAL,CMO)

         Call Put_Fock_Occ(Work(ipFocc),ntot1)
        If ( IPRLEV.ge.DEBUG ) then
        write(6,*) 'FOCC_OCC'
        call wrtmat(Work(ipFocc),1,ntot1,1,ntot1)


      write(6,*) 'DONE WITH NEW FOCK OPERATOR'
        end if
      call xflush(6)

         CALL GETMEM('FOCK','Free','REAL',LFOCK,NTOT4)
         CALL GETMEM('SXBM','Free','REAL',LBM,NSXS)
         CALL GETMEM('SXLQ','Free','REAL',LQ,NQ) ! q-matrix(1symmblock)
         IF(ISTORP(NSYM+1).GT.0) THEN
           CALL GETMEM('ISTRP','FREE','REAL',LP,ISTORP(NSYM+1))
           CALL GETMEM('ISTRP','FREE','REAL',LP1,ISTORP(NSYM+1))
         END IF
      Call GetMem('ONTOPO','FREE','Real',ipTmpLOEOTP,ntot1)
      Call GetMem('ONTOPT','FREE','Real',ipTmpLTEOTP,nfint)


!Put some information on the runfile for possible gradient calculations.
      Call Put_iScalar('Number of roots',nroots)
      Call Put_dArray('Last energies',Energies,nroots)
      Call Put_cArray('Relax Method','MCPDFT  ',8)
      Call Put_dScalar('Last energy',Energies(iRlxRoot))
          iSA = 1 !need to do MCLR for gradient runs. (1 to run, 2 to skip)
       !MUST MODIFY THIS.  I need to check that the calculation is not
       !SA, and if it is, set iSA to -1.
      Call Put_iScalar('SA ready',iSA)
      Call Put_cArray('MCLR Root','****************',16)


      Call GetMem('P2t','allo','Real',iP2dt1,NACPR2)
      Call FZero(Work(ip2dt1),Nacpr2)
      !I need the non-symmetry blocked d1act, hence the read.
      Call Get_D1MO(iD1Act1,NACPAR)
!        write(*,*) 'd1act'
!        do i=1,NACPAR
!          write(*,*) work(iD1Act1-1+i)
!        end do
      Call P2_contraction(Work(iD1Act1),Work(iP2dt1))
      Call Put_P2MOt(Work(iP2dt1),NACPR2)

      Call GetMem('P2t','free','Real',iP2dt1,NACPR2)
      Call GetMem('Dens','free','Real',iD1Act1,NACPAR)

!Put information needed for geometry optimizations.
      !if (jroot.eq.iRlxRoot) then
          iSA = 1 !need to do MCLR for gradient runs. (1 to run, 2 to skip)
       !MUST MODIFY THIS.  I need to check that the calculation is not
       !SA, and if it is, set iSA to -1.
      Call Put_iScalar('SA ready',iSA)
      Call Put_cArray('MCLR Root','****************',16)
      !end if

      end if !DoGradPDFT

      !if (jroot.eq.iRlxRoot) then
      Call Put_iScalar('Number of roots',nroots)
      Call Put_dArray('Last energies',Energies,nroots)
      Call Put_cArray('Relax Method','MCPDFT  ',8)
      Call Put_dScalar('Last energy',Energies(iRlxRoot))
      !end if

!      end do
      end do !loop over roots






      If (.not.DoCholesky .or. ALGO.eq.1) Then
        if (nFint.gt.0) then
          Call GetMem('PUVX','FREE','Real',lPUVX,nFint)
        end if
      End If

!Free up all the memory we can here, eh?
      Call GetMem('DtmpA','Free','Real',iTmp4,nTot1)
      Call GetMem('DtmpI','Free','Real',iTmp3,nTot1)

      Call GetMem('D1Active','free','Real',iD1Act,NACPAR)
      Call GetMem('D1ActiveAO','free','Real',iD1ActAO,NTOT2)
      Call GetMem('D1Spin','free','Real',iD1Spin,NACPAR)
      Call GetMem('D1SpinAO','free','Real',iD1SpinAO,NTOT2)
      Call GetMem('Fcore','Free','Real',iTmp1,nTot1)
      Call GetMem('FockI','FREE','Real',ifocki,ntot1)
      Call GetMem('FockA','FREE','Real',ifocka,ntot1)
      Call DDaFile(JOBOLD,0,Work(iP2d),NACPR2,dmDisk)
      Call GetMem('P2','Free','Real',iP2d,NACPR2)
      Call GetMem('D1Inact','Free','Real',iD1i,NTOT2)
      call xflush(6)
      Call qExit('MSCTL')
      Return
      END

      Subroutine P2_contraction(D1MO,P2MO)
      Implicit Real*8 (A-H,O-Z)
      Dimension D1MO(*),P2MO(*)
#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
#include "WrkSpc.fh"
      integer :: i,j,k,l,ij,kl,ijkl
      integer :: iD1c
      real*8 :: fact

      iTrii(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)

      Call GetMem('D1copy','Allo','Real',iD1c,NACPAR)
      Call DCOPY_(NACPAR,D1MO,1,Work(iD1c),1)
      Call DSCAL_(NACPAR,0.50d0,Work(iD1c),1)
      do i=1,nac
         Work(iD1c + iTrii(i,i)-1) =
     &   work(iD1c + iTrii(i,i)-1)*2
      end do
      call xflush(6)
       ijkl=0
       do i=1,nac
         do j=1,i
           ij = iTrii(i,j)
           do k=1,i
             if(i.eq.k) then
               lmax = j
             else
               lmax = k
             end if
             do l=1,lmax
               kl = iTrii(k,l)
               ijkl = ijkl + 1
               fact=1.0d0
               if(k.eq.l) fact=0.5d0
               p2MO(ijkl) = fact*Work(iD1c-1+ij)*Work(iD1c-1+kl)
             end do
           end do
         end do
       end do
       Call GetMem('D1copy','FREE','Real',iD1c,NACPAR)
      end subroutine
