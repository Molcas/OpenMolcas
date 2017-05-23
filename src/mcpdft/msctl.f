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
      integer iJOB,dmDisk,iCrap
      integer itmp0,itmp1,itmp2,itmp3,itmp4
      integer itmp5,itmp6,itmp7
      integer ifocki,ifocka
      integer IADR19(1:15)
      integer LP,NQ,LQ,LPUVX
      integer vdisk,jroot
      real*8,dimension(1:nroots) :: Energies
!      integer w,x,y,z,wx,yz,wxyz
!      integer ipd1
      iTrii(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)


      Call qEnter('MSCTL')

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

!Here we calculate the D1 Inactive matrix (AO).
      Call GetMem('D1Inact','Allo','Real',iD1I,NTOT2)
      Call Get_D1I_RASSCF_m(CMO,Work(iD1I))

!         write(*,*) 'iD1inact'
!         do i=1,ntot2
!         write(*,*) Work(iD1i-1+i)
!         end do


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
       vDisk =  IADR19(4)
       dmDisk = IADR19(3)
!       Call GetMem('jVector','Allo','Real',ijVec,nConf)
       Call GetMem('D1Active','Allo','Real',iD1Act,NACPAR)
       Call GetMem('D1ActiveAO','Allo','Real',iD1ActAO,NTOT2)
       Call GetMem('D1Spin','Allo','Real',iD1Spin,NACPAR)
       Call GetMem('D1SpinAO','Allo','Real',iD1SpinAO,NTOT2)

       Call GetMem('DtmpI','Allo','Real',iTmp3,nTot1)
       Call GetMem('DtmpA','Allo','Real',iTmp4,nTot1)
       Call GetMem('Crap','Allo','Real',iCrap,NACPR2)


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


!We loop over the number of states.  For each state, we read the density
!matrices from the JOBIPH file.
       do jroot=1,nroots

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
!         Call DDaFile(JOBIPH,2,Work(ijVec),nConf,jDisk)
!Get the D1 Active matrix for this state.  These should probably be
!most easily read from the previous JOBIPH file.  Then, convert D1A from
!the MO to the AO basis.

         Call DDaFile(JOBOLD,2,Work(iD1Act),NACPAR,dmDisk)
         Call Put_D1MO(Work(iD1Act),NACPAR)
         Call DDaFile(JOBOLD,2,Work(iD1Spin),NACPAR,dmDisk)
         Call DDaFile(JOBOLD,2,Work(iCrap),NACPR2,dmDisk)
         Call Put_P2MO(Work(iCrap),NACPR2)

         Call DDaFile(JOBOLD,0,Work(iCrap),NACPR2,dmDisk)


***********************************************************
* Generate total density
***********************************************************

         If(NASH(1).ne.NAC) Call DBLOCK_m(Work(iD1Act))
         Call Get_D1A_RASSCF_m(CMO,Work(iD1Act),Work(iD1ActAO))
!         write(*,*) "is this it?",ntot2
!         do i=1,ntot2
!           write(*,*) Work(iD1ActAO-1+i)
!         end do

         Call Fold(nSym,nBas,Work(iD1I),Work(iTmp3))
         Call Fold(nSym,nBas,Work(iD1ActAO),Work(iTmp4))
         Call Daxpy_(nTot1,1.0D0,Work(iTmp4),1,Work(iTmp3),1)
!Maybe I can write all of these matrices to file, then modify stuff in
!the nq code to read in the needed density.  In other words, I need to
!replace the next call with something that supports multistates.
         Call Put_D1ao(Work(iTmp3),nTot1)
!         write(*,*) 'end of d1 AO'
         call xflush(6)

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
      Else

         Write(LF,*)'SXCTL: Illegal Cholesky parameter ALGO= ',ALGO
         call qtrace()
         call abend()

      EndIf

         IF(ISTORP(NSYM+1).GT.0) THEN
           CALL GETMEM('ISTRP','ALLO','REAL',LP,ISTORP(NSYM+1))
           CALL DmatDmat_m(Work(iD1Act),WORK(LP))
         END IF
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
     &          D,WORK(LP),WORK(LQ),WORK(LPUVX),IFINAL,CMO)
         CALL GETMEM('FOCK','FREE','REAL',LFOCK,NTOT4)
         CALL GETMEM('SXBM','FREE','REAL',LBM,NSXS)
         CALL GETMEM('SXLQ','FREE','REAL',LQ,NQ)

         CASDFT_Funct = 0
         Call Get_dScalar('CASDFT energy',CASDFT_Funct)

         CASDFT_E = ECAS-EVAC+CASDFT_Funct

!         Write(*,*) '**************************************************'
!         write(*,*) 'ENERGY REPORT FOR STATE',jroot
        Call Print_MCPDFT_2(CASDFT_E,PotNuc,EMY,ECAS,CASDFT_Funct,
     &         jroot,Ref_Ener)
         call xflush(6)

         Energies(jroot)=CASDFT_E

!At this point, the energy calculation is done.  Now I need to build the
!fock matrix if this root corresponds to the relaxation root.
!
!I think the best plan is to reload FI,


!At this point, I think the Fock matrix corresponding to the CASSCF part
!of the calculation is stored in Work(ipFocc).  I think I should add
!this to the DFT part stored on the runfile.
       !Call GetMem('FOCKDFT','Allo','Real',iFDFT,nTot1)





!Put information needed for geometry optimizations.
!      write(*,*) jroot,iRLXROOT
      if (jroot.eq.iRlxRoot) then
!      write(*,*) "THIS IS THE ROOT TO RELAX"
!!!!!!!!!!!!!THIS IS BROKEN...can't put to P2MO (need that, yet!)
      if (.false.) then
      end if
      end if

           CALL GETMEM('ISTRP','FREE','REAL',LP,ISTORP(NSYM+1))

      end do !loop over roots

      If (.not.DoCholesky .or. ALGO.eq.1) Then
        if (nFint.gt.0) then
          Call GetMem('PUVX','FREE','Real',lPUVX,nFint)
        end if
      End If
!Put some information on the runfile for possible gradient calculations.
      Call Put_iScalar('Number of roots',nroots)
      Call Put_dArray('Last energies',Energies,nroots)
      Call Put_cArray('Relax Method','MCPDFT  ',8)
      Call Put_dScalar('Last energy',Energies(iRlxRoot))

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
      Call GetMem('Crap','Free','Real',iCrap,NACPR2)
      Call GetMem('D1Inact','Free','Real',iD1i,NTOT2)



!What is the multistate analogue of the following? Do we only want to
!print this for the IRLXROOT state, maybe?

***********************************************************
*     Transform FI to active orbital basis and move it over to F.
*     Remove also the symmetry blocking.
***********************************************************
***********************************************************
* Shall I add here the DFT contribution? Maybe not yet!
* I am commenting off... if needed we can always re-activate.
***********************************************************
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
!      Call Get_dExcdRa(ipTmpFckI,nTmpFck)
!      CALL DaXpY_(NTOT1,1.0D0,Work(ipTmpFckI),1,WORK(LX1),1)
!     If ( IPRLEV.ge.DEBUG ) then
!         Write(LF,*)
!         Write(LF,*) ' Exchange Corr. in AO basis in CASDFT_Terms'
!         Write(LF,*) ' ---------------------'
!         Write(LF,*)
!         iOff=1
!         Do iSym = 1,nSym
!           iBas = nBas(iSym)
!           Call TriPrt(' ','(5G17.11)',Work(ipTmpFckI+ioff-1),iBas)
!           iOff = iOff + (iBas*iBas+iBas)/2
!         End Do
!     End If
!      Call Free_Work(ipTmpFckI)
!      If ( IPRLEV.ge.DEBUG ) then
!        Write(LF,*)
!        Write(LF,*) ' Modified FI in AO basis in CASDFT_Terms'
!        Write(LF,*) ' ---------------------'
!        Write(LF,*)
!        iOff=1
!        Do iSym = 1,nSym
!          iBas = nBas(iSym)
!          Call TriPrt(' ','(5G17.11)',Work(LX1+ioff-1),iBas)
!          iOff = iOff + (iBas*iBas+iBas)/2
!        End Do
!      End If
*
      CALL MOTRAC(CMO,WORK(LX1),WORK(LX2),WORK(LX3))
      CALL GETMEM('XXX3','FREE','REAL',LX3,MXNB*MXNA)
      CALL GETMEM('XXX2','FREE','REAL',LX2,MXNB*MXNB)
      CALL DCOPY_(NACPAR,0.0D0,0,F,1)
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
              F(NTU)=WORK(LX1-1+ITU)
              IF(NT.EQ.NU) F(NTU)=F(NTU)+EMYN
              WORK(LX0-1+ITU) = F(NTU)
            END DO
          END DO
          IADD=IADD+NAT
        ENDIF
      END DO
      CALL CP_ONE_INT(WORK(LX0),ITU)
      CALL GETMEM('XXX1','FREE','REAL',LX1,NTOT1)
      CALL GETMEM('XXX0','FREE','REAL',LX0,NTOT1)

*     print h0
      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*)
        Write(LF,*)
        Write(LF,*)
        Write(LF,*) ' Fock matrix in MO basis, h0, in CASDFT_TERMS'
        Write(LF,*) ' ------------'
        Call TriPrt(' ',' ',F,NAC)
      End If

      Call qExit('CASDFT_terms')

      Return
      End











