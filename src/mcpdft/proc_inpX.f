************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Proc_InpX(DSCF,iRc)

! module dependencies
      use mspdft_grad, only: dogradmspd
      use mspdft, only: cmsNACstates, doNACMSPD, doMECIMSPD
#ifdef module_DMRG
!     use molcas_dmrg_interface !stknecht: Maquis-DMRG program
#endif
      use Fock_util_global, only: DoCholesky
      use Cholesky, only: ChFracMem
#ifdef _HDF5_
      Use mh5, Only: mh5_is_hdf5, mh5_open_file_r, mh5_exists_attr,
     &               mh5_exists_dset, mh5_fetch_attr, mh5_fetch_dset,
     &               mh5_close_file
      use stdalloc, only: mma_allocate, mma_deallocate
#endif
      use KSDFT_Info, only: CoefR, CoefX
      use hybridpdft, only: Ratio_WF, Do_Hybrid
      use UnixInfo, only: SuperName
      use write_pdft_job, only: iwjob, hasHDF5ref, hasMPSref
      use mcpdft_output, only: terse, debug, insane, lf, iPrLoc
      use lucia_interface, only: lucia_util

      Implicit Real*8 (A-H,O-Z)
#include "SysDef.fh"
#include "rasdim.fh"
#include "warnings.h"
#include "WrkSpc.fh"
#include "gas.fh"
#include "rasscf.fh"
#include "input_ras_mcpdft.fh"
#include "general.fh"
* Lucia-stuff:
#include "ciinfo.fh"
#include "spinfo.fh"
#include "lucia_ini.fh"
*
      Character*180  Line

      logical lExists, RunFile_Exists
* Some strange extra logical variables...
      logical DSCF
      logical RF_On
      logical Langevin_On
      logical PCM_On

      Logical DBG

#ifdef _HDF5_
* Local NBAS_L, NORB_L .. avoid collision with items in common.
      DIMENSION NFRO_L(8),NISH_L(8),NRS1_L(8),NRS2_L(8)
      DIMENSION NRS3_L(8),NSSH_L(8),NDEL_L(8)
      character(len=1), allocatable :: typestring(:)
      DIMENSION NBAS_L(8)
#endif
* TOC on JOBOLD (or JOBIPH)
      DIMENSION IADR19(15)

      Character*180 Get_LN
      External Get_LN
      Real*8   Get_ExFac
      External Get_ExFac
      Character*72 ReadStatus
      Character*72 JobTit(mxTit)
      Character*256 RealName
      Logical, External :: Is_First_Iter
      Character*(LENIN8*mxOrb) lJobH1
      Character*(2*72) lJobH2
      CHARACTER*(80) OriginalKS

      INTEGER :: iDNG,IPRLEV
      Logical :: DNG
      Character*8 emiloop
      Character*8 inGeo
      logical :: keyJOBI
      Intrinsic DBLE
C...Dongxia note for GAS:
C   No changing about read in orbital information from INPORB yet.

      Call StatusLine('MCPDFT:','Processing Input')

      IPRLEV = TERSE

      doGradPDFT = .false.
      doGradMSPD = .false.
      doNOGRad = .false.
      DoGSOR=.false.

! TODO PUT THIS INITIALITION IN MODULE FILE
*     CMS NACs variables
      doNACMSPD = .false.
      doMECIMSPD = .false.
      cmsNACstates(1) = 0
      cmsNACstates(2) = 0

*    GAS flag, means the INPUT was GAS
      iDoGas = .false.

!> default for MC-PDFT: read/write from/to JOBIPH-type files
      keyJOBI = .true.

      iRc=_RC_ALL_IS_WELL_


      If (SuperName(1:6).eq.'mcpdft') Then
* For geometry optimizations use the old CI coefficients.
        If (.Not.Is_First_Iter()) Then
          KeyFILE=.false.
        End If
      Else If (SuperName(1:18).eq.'numerical_gradient') Then
        KeyFILE=.false.
      End If

      !> Local print level in this routine:
      IPRLEV=IPRLOC(1)

      DBG= (IPRLEV >= DEBUG)

* ==== Check if there is any runfile ====
      Call F_Inquire('RUNFILE',RunFile_Exists)
      If (DBG) Write(6,*)' Inquire about RUNFILE.'
      IF (RunFile_Exists) Then
       If (DBG) Write(6,*)' Yes, there is one.'
       NSYM=0
       Call qpg_iScalar('nSym',lExists)
       IF (lExists) Then
        Call Get_iScalar('nSym',nSym)
        Call Get_iArray('nBas',nBas,nSym)
        If (DBG) Then
          Write(6,*)' The following information exists on runfile:'
          Write(6,*)' Nr of symmetries, NSYM:',NSYM
          Write(6,*)' Nr of basis functions/symmetry:'
          Write(6,'(1x,8I5)')(NBAS(I),I=1,NSYM)
          Call XFlush(6)
        End If
       ELSE
        Call WarningMessage(2,'No symmetry info on runfile.')
        Write(6,*)' There seems to be no information about symmetry'
        Write(6,*)' on the runfile! This is an unexpected error.'
        Call Quit(_RC_IO_ERROR_READ_)
       END IF
      ELSE
       Call WarningMessage(2,'Cannot find runfile.')
       Write(6,*)' PROC_INP: Cannot find RUNFILE. This is an'//
     &           ' unexpected error.'
        Call Quit(_RC_IO_ERROR_READ_)
      END IF
* ==== End check if there is any runfile ====

      iOrbData=0
* iOrbData=0: no orbital space data is specified
*         >0: specifications from some orbital file (JOBOLD, JOBIPH, HDF5)
      INVEC=0
* INVEC=0, no source for orbitals (yet)
*       3, take from JOBOLD, or JOBIPH file
*       4, take from an HDF5 file

*---  ==== FILE(ORB) keyword =====
      If (DBG) Write(6,*)' Where to read MOs/CI vectors? '
      StartOrbFile=''

      If (KeyFILE) Then
       If (DBG) Then
         Write(6,*)' Reading file name for start orbitals.'
       End If
       Call SetPos_m(LUInput,'FILE',Line,iRc)
       Line=Get_Ln(LUInput)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       If (DBG) Then
         Write(6,*) ' Calling fileorb with filename='
         Write(6,*) Line
       End If
       call fileorb(Line,StartOrbFile)
#ifdef _HDF5_
       if (mh5_is_hdf5(StartOrbFile)) then
         hasHDF5ref = .true.
       end if
#endif
!> we do not need a JOBIPH file if we have HDF5 - override the default!
       if(hasHDF5ref) keyJOBI = .false.

      End If

*---  ==== JOBI(PH) keyword =====

* The JOBIPH file, following decisions from the Torre Normanna labour camp:
* Default, the file name is 'JOBIPH'.
* However, if keyword IPHNAME was used, then the name was given in input.
* Also, if instead the keyword NEWIPH was given, then a new name will be
* chosen as the first not-already-used name in the sequence
* 'JOBIPH', 'JOBIPH01', 'JOBIPH02', etc.

      if(keyJOBI)then
        IPHNAME='ToBeFoun'
        If (KeyIPHN) Then
          If (DBG) Then
            Write(6,*)' Reading file name for JOBIPH file.'
          End If
          Call SetPos_m(LUInput,'IPHN',Line,iRc)
          If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
          ReadStatus=' Failure reading IPHNAME string.'
          Read(LUInput,*,End=9910,Err=9920) IPHNAME
          ReadStatus=' O.K. after reading IPHNAME string.'
          Call UpCase(IPHNAME)
        End If

        IF(IPHNAME.EQ.'ToBeFoun') IPHNAME='JOBIPH'

        !> check first for JOBOLD
        Call f_Inquire('JOBOLD',lExists)
        If (lExists) Then
          INVEC=3
        else !> check next for JOBIPH
          Call F_Inquire(trim(IPHNAME),lExists)
          If(lExists) Then
            Call PrgmTranslate(IPHNAME,RealName,not_sure)
            If (DBG) Then
                Write(6,*)' A JOBIPH file has been found:'
                write(6,*) RealName
            End If
            INVEC=3
          Else
            Write(LF,*)
            Write(LF,*)'******************************************'
            Write(LF,*)' JOBIPH does not seem to exist,           '
            Write(LF,*)' so the calculation cannot continue.      '
            Write(LF,*)'******************************************'
            Call Abend()
          End If
        end if !> check for JOBOLD

        if(JOBIPH.gt.0) Then
          Call DaClos(JOBIPH)
          JOBIPH=-1
        end if
        JOBIPH=IsFreeUnit(15)
        CALL DANAME(JOBIPH,IPHNAME)
        INVEC=3
      end if !> JOBI(PH) keyword

*---  ==== JOBI(PH) keyword =====

*---  process KSDF command --------------------------------------------*
      If (DBG) Write(6,*) ' Check if KSDFT was requested.'
      If (.not.KeyKSDF) Then
        Call WarningMessage(2,'No KSDFT functional specified')
        Write(LF,*) ' ************* ERROR **************'
        Write(LF,*) ' KSDFT functional type must be     '
        Write(LF,*) ' specified for MCPDFT calculations '
        Write(LF,*) ' **********************************'
        Call Abend()
      End If
      If (DBG) Write(6,*) ' KSDFT command was given.'
      DFTFOCK='ROKS'
      Call SetPos_m(LUInput,'KSDF',Line,iRc)
      If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
      Read(LUInput,*,End=9910,Err=9920) Line
      KSDFT=Line(1:80)
      Call UpCase(KSDFT)
* checking KSDFT input for MC-PDFT
      IF(KSDFT(1:2) == 'T:') THEN
       OriginalKS=KSDFT(3:80)
      ELSE IF(KSDFT(1:3) == 'FT:') THEN
       OriginalKS=KSDFT(4:80)
      ELSE
       Call WarningMessage(2,'Wrong on-top functional for MC-PDFT')
       Write(LF,*) ' ************* ERROR **************'
       Write(LF,*) ' Current on-top functionals are:   '
       Write(LF,*) ' T :  translated functionals       '
       Write(LF,*) ' FT:  fully translated functionals '
       Write(LF,*) ' e.g. T:PBE for tPBE functional    '
       Write(LF,*) ' **********************************'
       Call Abend()
      END IF
      If (DBG) write(LF,*) 'The original KS functional is', OriginalKS
      ExFac=Get_ExFac(OriginalKS)
* Assuming hybrid KS functionals contain more than 10E-6 percent
* Hartree-Fock exchagne.
      IF(Abs(ExFac).gt.1.0d-8) THEN
       Call WarningMessage(2,'Hybrid functionals not supported')
       Write(LF,*) ' ************* ERROR **************'
       Write(LF,*) ' MC-PDFT does not translate hybrid '
       Write(LF,*) ' functionals. If you want to run   '
       Write(LF,*) ' hybrid MC-PDFT, use the LAMBda    '
       Write(LF,*) ' keyword instead.                  '
       Write(LF,*) '                                   '
       Write(LF,*) ' EXAMPLE:                          '
       Write(LF,*) '  tPBE0 = 75% tPBE + 25% MCSCF.    '
       Write(LF,*) ' Usage:                            '
       Write(LF,*) '  KSDFT=T:PBE                      '
       Write(LF,*) '  LAMB =0.25                       '
       Write(LF,*) ' **********************************'
       Call Abend()
      END IF
* End of checking KSDFT input for MC-PDFT

      ExFac=Get_ExFac(KSDFT)

*---  Process DFCF command (S Dong, 2018)--------------------------*
      If (DBG) Write(6,*) ' Check if DFCF was provided.'
      If (KeyDFCF) Then
       If (DBG) Write(6,*) ' DFCF command has been used.'
       Call SetPos_m(LUInput,'DFCF',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       ReadStatus=' Failure after reading DFCF keyword.'
       Read(LUInput,*,End=9910,Err=9920) CoefX,CoefR
       ReadStatus=' O.K. after reading DFCF keyword.'
      End If
*---  Process MSPD command --------------------------------------------*
      If (DBG) Write(6,*) ' Check if Multi-state MC-PDFT case.'
      If (KeyMSPD) Then
       If (DBG) Write(6,*) ' MSPD keyword was used.'
       iMSPDFT=1
       Call SetPos_m(LUInput,'MSPD',Line,iRc)
       if(dogradpdft) then
        dogradmspd=.true.
        dogradpdft=.false.
       end if
      End If
*---  Process WJOB command --------------------------------------------*
      If (DBG) Write(6,*) ' Check if write JOBIPH case.'
      If (KeyWJOB) Then
       If (DBG) Write(6,*) ' WJOB keyword was used.'
       iWJOB=1
       Call SetPos_m(LUInput,'WJOB',Line,iRc)
      End If
*---  Process LAMB command --------------------------------------------*
      If (KeyLAMB) Then
       If (DBG) Write(6,*) 'Check if hybrid PDFT case'
       Call SetPos_m(LUInput,'LAMB',Line,iRc)
       ReadStatus=' Failure reading data following LAMB keyword.'
       Read(LUInput,*,End=9910,Err=9920) Ratio_WF
       ReadStatus=' O.K. reading data following LAMB keyword.'
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       If(Ratio_WF.gt.0.0d0) Then
        Do_Hybrid=.true.
        CALL Put_DScalar('R_WF_HMC',Ratio_WF)
       End If
       If (DBG) Write(6,*) 'Wave Funtion Ratio in hybrid PDFT',Ratio_WF
       If (dogradmspd.or.dogradpdft) Then
        Call WarningMessage(2,'GRAD currently not compatible with HPDF')
        GoTo 9810
       End If
      End If

*---  Process HDF5 file --------------------------------------------*
      If (hasHDF5ref) Then
#ifdef _HDF5_
        mh5id = mh5_open_file_r(StartOrbFile)
*     read basic attributes
        call mh5_fetch_attr(mh5id, 'NSYM', NSYM_L)
        if (nsym.ne.nsym_l) then
          write (LF,*) 'Number of symmetries on HDF5 file does not'
          write (LF,*) 'match the number of symmetries on the'
          write (LF,*) 'RunFile, calculation will stop now.'
          call Quit(_RC_INPUT_ERROR_)
        end if
        call mh5_fetch_attr(mh5id, 'NBAS', NBAS_L)
        ierr=0
        do isym=1,nsym
          if (nbas(isym).ne.nbas_l(isym)) ierr=1
        end do
        if (ierr.eq.1) then
          write (LF,*) 'Number of basis functions on HDF5 file does not'
          write (LF,*) 'match the number of basis functions on the'
          write (LF,*) 'RunFile, calculation will stop now.'
          call Quit(_RC_INPUT_ERROR_)
        end if
*     orbitals available?
        if (mh5_exists_dset(mh5id, 'MO_VECTORS')) then
          inVec=4
        else
          write (LF,*)'The HDF5 ref file does not contain MO vectors.'
          write (LF,*)'Fatal error, the calculation will stop now.'
          call Quit(_RC_INPUT_ERROR_)
        end if
*     typeindex data available?
        if (mh5_exists_dset(mh5id, 'MO_TYPEINDICES')) then
          iOrbData=3
          call mma_allocate(typestring, sum(nbas(1:nsym)))
          call mh5_fetch_dset(mh5id, 'MO_TYPEINDICES', typestring)
          call tpstr2orb(nSym,nbas_l,
     $            typestring,
     $            nFro_L,nISh_L,
     $            NRS1_L,NRS2_L,NRS3_L,
     $            nSSh_L,nDel_L)
          call mma_deallocate(typestring)
        else
          write (LF,*)'The HDF5 ref file does not contain TYPEindices.'
          write (LF,*)'Fatal error, the calculation will stop now.'
          call Quit(_RC_INPUT_ERROR_)
        end if

#ifdef _DMRG_
        !> QCMaquis MPS reference wave function
        if (mh5_exists_dset(mh5id, 'QCMAQUIS_CHECKPOINT')) then
           hasMPSref  = .true.
        end if
#endif
        if(.not.hasMPSref)then
          if (mh5_exists_dset(mh5id, 'CI_VECTORS'))then
            write (LF,*)' CI vectors will be read from HDF5 ref file.'
          else
            write (LF,*)'The HDF5 ref file does not contain CI vectors.'
            write (LF,*)'Fatal error, the calculation will stop now.'
            call Quit(_RC_INPUT_ERROR_)
          end if
        end if
        call mh5_close_file(mh5id)
#else
        write (6,*) 'The format of the start orbital file was'
        write (6,*) 'specified by the user as HDF5, but this'
        write (6,*) 'is not implemented in this installation.'
        call Quit(_RC_INPUT_ERROR_)
#endif
      End If

* =======================================================================
#ifdef _HDF5_
      !> transfer orbital space data read from HDF5 file
      IF(IORBDATA.eq.3) THEN
        DO ISYM=1,NSYM
          NFRO(ISYM)=NFRO_L(ISYM)
          NISH(ISYM)=NISH_L(ISYM)
          NRS1(ISYM)=NRS1_L(ISYM)
          NRS2(ISYM)=NRS2_L(ISYM)
          NRS3(ISYM)=NRS3_L(ISYM)
          NSSH(ISYM)=NSSH_L(ISYM)
          NDEL(ISYM)=NDEL_L(ISYM)
        END DO
      END IF
#endif
* =======================================================================
      iprlev=insane

!> read orbital space data AND CI optimiation parameters from JOBIPH
      IF (IORBDATA.EQ.0) THEN
        IAD19=0
        Call IDaFile(JOBIPH,2,IADR19,10,IAD19)
        iAd19=iAdr19(1)
        CALL WR_RASSCF_Info(JobIPH,2,iAd19,NACTEL,ISPIN,NSYM,STSYM,
     &                      NFRO,NISH,NASH,NDEL,NBAS,
     &                      mxSym,lJobH1,LENIN8*mxOrb,NCONF,
     &                      lJobH2,2*72,JobTit,4*18*mxTit,
     &                      POTNUCDUMMY,LROOTS,NROOTS,IROOT,mxRoot,
     &                      NRS1,NRS2,NRS3,NHOLE1,NELEC3,IPT2,WEIGHT)
      End If  !> IORBDATA

!> read CI optimiation parameters from HDF5 file
      if(hasHDF5ref)then
#ifdef _HDF5_
        mh5id = mh5_open_file_r(StartOrbFile)
        call mh5_fetch_attr (mh5id,'SPINMULT', iSpin)
        call mh5_fetch_attr (mh5id,'NSYM', nSym)
        call mh5_fetch_attr (mh5id,'LSYM', stSym)
        call mh5_fetch_attr (mh5id,'NBAS', nBas)

        call mh5_fetch_attr (mh5id,'NACTEL', nactel)
        call mh5_fetch_attr (mh5id,'NHOLE1', nhole1)
        call mh5_fetch_attr (mh5id,'NELEC3', nelec3)
        call mh5_fetch_attr (mh5id,'NCONF',  nconf)
        call mh5_fetch_attr (mh5id,'NSTATES', lroots)
        If (mh5_exists_attr(mh5id, 'NROOTS')) Then
          call mh5_fetch_attr (mh5id,'NROOTS', nroots)
        Else
          nroots = lroots
        End If
        call mh5_fetch_attr (mh5id,'STATE_WEIGHT', weight)

        call mh5_close_file(mh5id)
#endif
      end if

!AMS - this may be closing either JOBOLD or JOBIPH. Close only JOBOLD.
      IF(JOBOLD>0) then
        IF(JOBOLD.ne.JOBIPH) THEN
            Call DaClos(JOBOLD)
        END IF
      end if


!AMS - make sure we change to a different JOBIPH file - we don't want to
!overwrite any existing JOBIPH file.
!
!Close the old JOBIPH file
      if(JOBIPH.gt.0) Then
        Call DaClos(JOBIPH)
        JOBIPH=-1
      end if
!Rename JOBIPH file, and open it.
      JOBIPH=IsFreeUnit(15)
      CALL DANAME(JOBIPH,IPHNAME)

*---  complete orbital specifications ---------------------------------*
      Do iSym=1,nSym
        nash(isym)=nrs1(isym)+nrs2(isym)+nrs3(isym)
        NORB(ISYM)=NBAS(ISYM)-NFRO(ISYM)-NDEL(ISYM)
        NSSH(ISYM)=NORB(ISYM)-NISH(ISYM)-NASH(ISYM)
      End Do
*---  Related data for sizes, etc.
      NTOT=0
      NTOT1=0
      NTOT2=0
      NO2M=0
      NISHT=0
      NASHT=0
      NDELT=0
      NFROT=0
      NSEC=0
      NORBT=0
      NTOT3=0
      NTOTSP=0
      NTOT4=0
      NRS1T=0 ! for RASSCF
      NRS2T=0
      NRS3T=0
c      Call FZero(NGSSH_tot,ngas)
c      do igas=1,ngas
c        NGSSH_tot(igas) = SUM(NGSSH(IGAS,1:NSYM))
c      end do
      DO ISYM=1,NSYM
         NTOT=NTOT+NBAS(ISYM)
         NTOT1=NTOT1+NBAS(ISYM)*(NBAS(ISYM)+1)/2
         NTOT2=NTOT2+NBAS(ISYM)**2
         NO2M=MAX(NO2M,NBAS(ISYM)**2)
         NRS1T=NRS1T+NRS1(ISYM)  ! for RAS
         NRS2T=NRS2T+NRS2(ISYM)
         NRS3T=NRS3T+NRS3(ISYM)
         NFROT=NFROT+NFRO(ISYM)
         NISHT=NISHT+NISH(ISYM)
         NASHT=NASHT+NASH(ISYM)
         NDELT=NDELT+NDEL(ISYM)
         NSEC=NSEC+NSSH(ISYM)
         NORBT=NORBT+NORB(ISYM)
         NTOT3=NTOT3+(NORB(ISYM)+NORB(ISYM)**2)/2
         NTOTSP=NTOTSP+(NASH(ISYM)*(NASH(ISYM)+1)/2)
         NTOT4=NTOT4+NORB(ISYM)**2
      END DO
      NACPAR=(NASHT+NASHT**2)/2
      NACPR2=(NACPAR+NACPAR**2)/2
* NASHT is called NAC in some places:
      NAC=NASHT
* Same, NISHT, NIN:
      NIN=NISHT
      NFR=NFROT
      If (DBG) Write(6,*)' The iOrbData code is now',iOrbData
* =======================================================================
* Compute effective nuclear charge.
* Identical to nr of protons for conventional basis sets only, not ECP.
      Call Get_iScalar('Unique atoms',nNuc)
      Call GetMem('EffNChrg','Allo','Real',ipENC,nNuc)
      Call Get_dArray('Effective nuclear Charge',Work(ipENC),nNuc)
      TEffNChrg=0.0D0
      Call GetMem('nStab','Allo','Inte',ipStab,nNuc)
      Call Get_iArray('nStab',iWork(ipStab),nNuc)
      do i=1,nNuc
       TEffNChrg=TEffNChrg+Work(ipENC-1+i)*DBLE(nSym/iWork(ipStab-1+i))
      end do
      Call GetMem('nStab','Free','Inte',ipStab,nNuc)
      Call GetMem('EffNChrg','Free','Real',ipENC,nNuc)
      If (DBG) Write(6,*)
     &             ' Effective nuclear charge is TEffNChrg=',TEffNChrg
      TotChrg=0.0D0
      If (DBG) Write(6,*)' Set TotChrg=',TotChrg
*---  Process GRAD command --------------------------------------------*
      If (DBG) Write(6,*) ' Check if GRADient case.'
      If (KeyGRAD) Then
       If (DBG) Write(6,*) ' GRADient keyword was used.'
       DoGradPDFT=.true.
       if(iMSPDFT==1) then
        dogradmspd=.true.
        dogradpdft=.false.
       end if
       Call SetPos_m(LUInput,'GRAD',Line,iRc)
*TRS - Adding else statement to make nograd the default if the grad
*keyword isn't used
       Else
       DoNoGrad=.true.
*TRS
      End If
*---  Process NAC command --------------------------------------------*
      If (DBG) Write(6,*) ' Check if NAC case.'
      If (KeyNAC) Then
       If (DBG) Write(6,*) ' NAC keyword was used.'
       if(iMSPDFT==1 .and. DoGradMSPD .eqv. .true.)  then
           doNACMSPD=.true.
       end if
       if(iMSPDFT==0) then
        Call WarningMessage(2,'NACs implemented only for MS-PDFT')
        Write(LF,*) ' ************* ERROR **************'
        Write(LF,*) ' NACs are only implemented         '
        Write(LF,*) ' for Multistate PDFT               '
        Write(LF,*) ' **********************************'
        Call Abend()
       end if
       if(DoGradMSPD .eqv. .false.) then
        Call WarningMessage(2,'NACs implemented with GRAD Code')
        Write(LF,*) ' ************* ERROR **************'
        Write(LF,*) ' NACs require the GRAD Keyword     '
        Write(LF,*) ' **********************************'
        Call Abend()
       end if
       Call SetPos_m(LUInput,'NAC',Line,iRc)
       If(iRc.ne._RC_ALL_IS_WELL_) GoTo 9810
       Read(LUInput,*,End=9910,Err=9920) cmsNACstates(1),
     & cmsNACstates(2)
      End If
*
*---  Process MECI command --------------------------------------------*
      If (DBG) Write(6,*) ' Check if MECI case.'
      If (KeyMECI) Then
       If (DBG) Write(6,*) ' MECI keyword was used.'
       if(iMSPDFT==1 .and. DoGradMSPD .eqv. .true. .and.
     & doNACMSPD .eqv. .true.)  then
           doMECIMSPD=.true.
       end if
       if(iMSPDFT==0) then
        Call WarningMessage(2,'NACs implemented only for MS-PDFT')
        Write(LF,*) ' ************* ERROR **************'
        Write(LF,*) ' MECI is only implemented          '
        Write(LF,*) ' for Multistate PDFT               '
        Write(LF,*) ' **********************************'
        Call Abend()
       end if
       if(DoGradMSPD .eqv. .false.) then
        Call WarningMessage(2,'NACs implemented with GRAD Code')
        Write(LF,*) ' ************* ERROR **************'
        Write(LF,*) ' MECI requires the GRAD Keyword    '
        Write(LF,*) ' **********************************'
        Call Abend()
       end if
       if(DoNACMSPD .eqv. .false.) then
        Call WarningMessage(2,'NACs implemented with GRAD Code')
        Write(LF,*) ' ************* ERROR **************'
        Write(LF,*) ' MECI requires the NAC Keyword     '
        Write(LF,*) ' **********************************'
        Call Abend()
       end if
       Call SetPos_m(LUInput,'MECI',Line,iRc)
      End If
*
*---  Process GSOR command --------------------------------------------*
      If (DBG) Write(6,*) ' Check if Gram-Schmidt case.'
      If (KeyGSOR) Then
       If (DBG) Write(6,*) ' GSOR keyword was used.'
       DoGSOR=.true.
       Call SetPos_m(LUInput,'GSOR',Line,iRc)
      End If
*
*---  All keywords have been processed ------------------------------*

************************************************************************
* Generate artificial splitting or RAS into GAS for parallel blocking  *
************************************************************************
* SVC: convert CAS/RAS to general GAS description here, then we only
* need to copy it for lucia later, which always uses GAS description.
      NGSSH(1,1:NSYM)=NRS1(1:NSYM)
      NGSSH(2,1:NSYM)=NRS2(1:NSYM)
      NGSSH(3,1:NSYM)=NRS3(1:NSYM)
      IGSOCCX(1,1) = MAX(2*SUM(NRS1(1:NSYM))-NHOLE1,0)
      IGSOCCX(1,2) = 2*SUM(NRS1(1:NSYM))
      IGSOCCX(2,1) = NACTEL - NELEC3
      IGSOCCX(2,2) = NACTEL
      IGSOCCX(3,1) = NACTEL
      IGSOCCX(3,2) = NACTEL
*
!Considerations for gradients/geometry optimizations

*     Numerical gradients requested in GATEWAY
      Call Qpg_iScalar('DNG',DNG)
      If (DNG) Then
         Call Get_iScalar('DNG',iDNG)
         DNG = iDNG.eq.1
      End If
      DNG=DoNoGrad.or.DNG
*
*     Inside LAST_ENERGY we do not need analytical gradients
      If (SuperName(1:11).eq.'last_energy') DNG=.true.
*
*     Inside NUMERICAL_GRADIENT override input!
      If (SuperName(1:18).eq.'numerical_gradient') DNG=.true.
*
*
      If (DNG) Then
         DoGradPDFT=.false.
         if(iMSPDFT==1) then
          dogradmspd=.false.
         end if
      End If
*
*     Check to see if we are in a Do While loop
      Call GetEnvF('EMIL_InLoop',emiloop)
      If (emiloop.eq.' ') emiloop='0'
      Call GetEnvF('MOLCAS_IN_GEO',inGeo)
      If ((emiloop(1:1).ne.'0') .and. inGeo(1:1) .ne. 'Y'
     &    .and. .not.DNG) Then
         DoGradPDFT=.true.
         if(iMSPDFT==1) then
          dogradmspd=.true.
          dogradpdft=.false.
         end if
      End If
*                                                                      *
*---  Compute IZROT. IZROT is a matrix (lower triangular over the -----*
*     active space), which specifies which t,u rotations should be
*     avoided, since the orbitals belong to the same RAS space.
*     This is the only way the RAS concept is explicitly used in the
*     SX section of the program.
      ITU=0
      DO ISYM=1,NSYM
        NAO=NASH(ISYM)
          DO NT=2,NAO
            DO NU=1,NT-1
              ITU=ITU+1
              IZROT(ITU)=0
CSVC: check if NU<NT are included in the same gas space
              NGSSH_LO=0
              DO IGAS=1,NGAS
                NGSSH_HI=NGSSH_LO+NGSSH(IGAS,ISYM)
                IF (NU.GT.NGSSH_LO.AND.NT.LE.NGSSH_HI) THEN
                  IZROT(ITU)=1
                END IF
                NGSSH_LO=NGSSH_HI
              END DO
            END DO
          END DO
      END DO
*
      Call Put_iArray('nIsh',nIsh,nSym)
      Call Put_iArray('nAsh',nAsh,nSym)
      Call Put_iScalar('Multiplicity',ISPIN)
*
*---  Initialize Cholesky information if requested
      if (DoCholesky) then
         Call Cho_X_init(irc,ChFracMem)
         if (irc.ne.0) Go To 9930
      endif

* ===============================================================

*
*     Initialize seward
*
      If (DBG) Write(6,*)' Initialize seward.'
      nDiff = 0
      If (DSCF           .or.
     &    RF_On()        .or.
     &    Langevin_On()  .or.
     &    PCM_On()       .or.
     &    KSDFT.ne.'SCF'     )
     &    Call IniSew(DSCF.or.Langevin_On().or.PCM_On(),nDiff)
* ===============================================================
*
*     Check the input data
*
      If (DBG) Then
        Write(6,*)' Call ChkInp.'
        Call XFlush(6)
      End If
      Call ChkInp_m()
* ===============================================================

      if(.not.DoGSOR) then
        NCONF=1
        GoTo 9000
      end if

      If(hasMPSref) GoTo 9000
* ===============================================================
*
*     Construct the Guga tables
*
      call gugactl()
* ===============================================================
*
*     Construct the determinant tables
*

      If (DBG) Write(6,*)' Construct the determinant tables.'
      MS2 = iSpin - 1
*
* Set variables needed in Lucia_Ini
*
      Call iCopy(mxGAS*mxSym,ngssh,1,ngssh_Molcas,1)
      Call iCopy(mxGAS*2,igsoccx,1,igsoccx_Molcas,1)
      Call iCopy(nSym,norb,1,norb_Molcas,1)
      Call iCopy(nSym,nbas,1,nbas_Molcas,1)
      Call iCopy(nSym,nish,1,nish_Molcas,1)
      potnuc_Molcas    = potnuc
      thre_Molcas      = thre
      nsym_Molcas      = nsym
      nactel_Molcas    = nactel
      ms2_Molcas       = ms2
      ispin_Molcas     = ispin
      lsym_Molcas      = stsym
      NHOLE1_Molcas    = NHOLE1
      NELEC3_Molcas    = NELEC3
      itmax_Molcas     = itmax
      rtoi_Molcas      = rtoi
      nroots_Molcas    = Max(nroots,lRoots)
      ipt2_Molcas      = ipt2
      iprci_molcas     = iprloc(3)
      ngas_molcas      = ngas
      INOCALC_MOLCAS   = INOCALC
      ISAVE_EXP_MOLCAS = ISAVE_EXP
      IEXPAND_MOLCAS   = IEXPAND

*
* And call Lucia_Ini to initialize LUCIA
      CALL Lucia_Util('Ini')
* to get number of CSFs for GAS
      nconf=0
      do i=1,mxsym
        nconf=nconf+ncsasm(i)
      end do
*
      ISCF=0
      IF (ISPIN.EQ.NAC+1.AND.NACTEL.EQ.NAC) ISCF=1
      IF (ISPIN.EQ.1.AND.NACTEL.EQ.2*NAC)   ISCF=1
      IF (ISCF.EQ.1) THEN
         NCONF=1
      END IF

* ===============================================================

      Go to 9000
*
*---  Error exits -----------------------------------------------------*
9810  CONTINUE
      If (IPRLEV.ge.TERSE) Then
       Call WarningMessage(2,'Error in input preprocessing.')
       Write(6,*)' PROC_INP: A keyword was found during prescanning'
       Write(6,*)' the input file, but when later trying to locate'
       Write(6,*)' this input, it could not be found. Something has'
       Write(6,*)' happened to the input file, or else there is some'
       Write(6,*)' strange program error.'
       iRc=_RC_INPUT_ERROR_
      End If
      Go to 9900

9910  CONTINUE
      Call WarningMessage(2,'End of input file during preprocessing.')
      Call WarningMessage(2,ReadStatus)
      If (IPRLEV.ge.TERSE) Write(6,*)' Error exit 9910 from PROC_INP.'
      iRc=_RC_INPUT_ERROR_
      Go to 9900
*
9920  CONTINUE
      Call WarningMessage(2,'Read error during input preprocessing.')
      Call WarningMessage(2,ReadStatus)
      If (IPRLEV.ge.TERSE) Write(6,*)' Error exit 9920 from PROC_INP.'
      iRc=_RC_INPUT_ERROR_
      Go to 9900
*
9930  CONTINUE
      Call WarningMessage(2,'Error during input preprocessing.')
      Call WarningMessage(2,ReadStatus)
      If (IPRLEV.ge.TERSE) Write(6,*)' Error exit 9930 from PROC_INP.'
      iRc=_RC_INPUT_ERROR_
      Go to 9900

*---  Normal exit -----------------------------------------------------*
9000  CONTINUE
      close(989)
      If (DBG) Write(6,*)' Normal exit from PROC_INP.'
      Return
*---  Abnormal exit -----------------------------------------------------*
9900  CONTINUE
      If (DBG) Write(6,*)' Abnormal exit from PROC_INP.'
      Return
      End
