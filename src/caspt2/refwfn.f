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
      module refwfn
      Implicit None
      Logical :: refwfn_active = .False.
      Character(128) :: refwfn_filename
      Integer :: refwfn_id
      Logical :: refwfn_is_h5 = .False.
      Integer :: IADR15(30)
      Save

      Contains

************************************************************************
      Subroutine refwfn_init(Filename)
************************************************************************
      Implicit None
      Character(*) :: Filename
      Integer :: I, IAD15
#ifdef _HDF5_
#  include "mh5.fh"
#endif

      If (refwfn_active) Then
        write(6,*) ' trying to activate refwfn twice, aborting!'
        call abend
      Else
        refwfn_active = .True.
      End If

* if not a standard filename, call fileorb??
      If (FileName.ne.'JOBIPH') Then
        call fileorb(Filename,refwfn_filename)
      Else
        refwfn_filename = 'JOBIPH'
      End If

#ifdef _HDF5_
      If (mh5_is_hdf5(refwfn_filename)) Then
        refwfn_is_h5 = .True.
        write(6,'(1X,A)') 'wavefunction data from HDF5 file:'
        write(6,'(3X,A)') TRIM(refwfn_filename)
        refwfn_id = mh5_open_file_r(refwfn_filename)
      Else
#endif
        refwfn_is_h5 = .False.
* Assume reference wavefunction is stored as JobIph format
        refwfn_id=15
        CALL DANAME(refwfn_id,refwfn_filename)
* Read table of contents into IADR15() array.
* There are two possible different layouts, 15 or 30 integers:
        IAD15=0
        CALL IDAFILE(refwfn_id,2,IADR15,15,IAD15)
        IF(IADR15(15).EQ.-1) THEN
          IAD15=0
          CALL IDAFILE(refwfn_id,2,IADR15,30,IAD15)
        ELSE
          DO I=16,30
            IADR15(I)=0
          END DO
          Call WarningMessage(1,'Old JOBIPH file layout.')
        END IF
#ifdef _HDF5_
      End If
#endif
      End Subroutine

************************************************************************
      Subroutine refwfn_close
************************************************************************
      Implicit None
#ifdef _HDF5_
#  include "mh5.fh"
#endif

#ifdef _HDF5_
      If (refwfn_is_h5) Then
        call mh5_close_file(refwfn_id)
      Else
#endif
        call DaClos(refwfn_id)
#ifdef _HDF5_
      End If
#endif
      refwfn_id = -1
      refwfn_active = .False.
      End SUbroutine

************************************************************************
      Subroutine refwfn_info
************************************************************************
CSVC: initialize the reference wavefunction info
      Implicit None
#include "rasdim.fh"
#include "caspt2.fh"
#include "stdalloc.fh"
#ifdef _HDF5_
#  include "mh5.fh"
      character(1), allocatable :: typestring(:)
#endif
      Integer iSym, ref_nSym, ref_nBas(mxSym)
      Real*8 :: Weight(mxRoot)
      Integer IAD15

      If (.NOT.refwfn_active) Then
        Write(6,*) ' refwfn not yet activated, aborting!'
        call abend
      End If

#ifdef _HDF5_
      If (refwfn_is_h5) Then
*     general wavefunction attributes
*        call mh5_fetch_attr (refwfn_id, 'TITLE', Title)
        call mh5_fetch_attr (refwfn_id, 'SPINMULT', iSpin)
        call mh5_fetch_attr (refwfn_id, 'NSYM', ref_nSym)
        call mh5_fetch_attr (refwfn_id, 'LSYM', lSym)
        call mh5_fetch_attr (refwfn_id, 'NBAS', ref_nBas)

        call mh5_fetch_attr (refwfn_id, 'NACTEL', nActEl)
        call mh5_fetch_attr (refwfn_id, 'NHOLE1', nHole1)
        call mh5_fetch_attr (refwfn_id, 'NELEC3', nEle3)
        call mh5_fetch_attr (refwfn_id, 'NCONF',  nConf)
        call mh5_fetch_attr (refwfn_id, 'NSTATES', nRoots)
        call mh5_fetch_attr (refwfn_id, 'NROOTS', lRoots)
        call mh5_fetch_attr (refwfn_id, 'STATE_ROOTID', iRoot)
        call mh5_fetch_attr (refwfn_id, 'STATE_WEIGHT', Weight)

        call mma_allocate (typestring, sum(ref_nbas(1:nsym)))
        call mh5_fetch_dset (refwfn_id, 'MO_TYPEINDICES', typestring)
        call tpstr2orb (ref_nsym,ref_nbas,typestring,
     $          nfro,nish,nras1,nras2,nras3,nssh,ndel)
        nash = nras1 + nras2 + nras3
        call mma_deallocate (typestring)

        If (.not.mh5_exists_dset(refwfn_id, 'CI_VECTORS')) Then
          Write(6,'(1X,A)') 'The HDF5 file does not contain CI vectors,'
          Write(6,'(1X,A)') 'make sure it was created by rasscf/caspt2.'
          Call AbEnd()
        End If
        If (.not.mh5_exists_dset(refwfn_id, 'MO_VECTORS')) Then
          Write(6,'(1X,A)') 'The HDF5 file does not contain MO vectors,'
          Write(6,'(1X,A)') 'make sure it was created by rasscf/caspt2.'
          Call AbEnd()
        End If
        IFQCAN=0
      Else
#endif
C Sizes in the GSLIST is counted in INTEGERS.
C Note that the title field in the JOBIPH file is not used for anything
C in this program, it is just a dummy read.
C Another title field is read from input a little later, it is called
C TITLE2. That one is printed out in PRINP_CASPT2.
        IAD15=IADR15(1)
        CALL WR_RASSCF_Info(refwfn_id,2,iAd15,
     &                      NACTEL,ISPIN,REF_NSYM,LSYM,
     &                      NFRO,NISH,NASH,NDEL,REF_NBAS,8,
     &                      NAME,LENIN8*MXORB,NCONF,HEADER,144,
     &                      TITLE,4*18*mxTit,POTNUC,
     &                      LROOTS,NROOTS,IROOT,MXROOT,NRAS1,
     &                      NRAS2,NRAS3,NHOLE1,NELE3,IFQCAN,
     &                      Weight)
      nssh = ref_nbas - nfro - nish - nash - ndel
#ifdef _HDF5_
      End If
#endif
      If (nSym.ne.ref_nSym) Then
        write(6,*) ' Number of irreps of the reference wavefunction'
        write(6,*) ' does not match the data on the RunFile, abort!'
        Call AbEnd
      Else
        Do iSym=1,nSym
          If (nBas(iSym).ne.ref_nBas(iSym)) Then
            write(6,*) ' Number of basis functions of the reference'
            write(6,*) ' wavefunction does not match the data on the'
            write(6,*) ' RunFile, abort!'
            Call AbEnd
          End If
        End Do
      End If
      End Subroutine

************************************************************************
      subroutine refwfn_data
************************************************************************
CSVC: initialize the reference wavefunction data
      Implicit None
#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#ifdef _HDF5_
#  include "mh5.fh"
#endif

      Integer :: I, IAD15, II, IDISK, ID
      Integer :: J, JSNUM, LHEFF1

      Real*8 :: Root_Energies(mxRoot)
      Real*8 :: AEMAX, E
      Integer :: IAD, LEJOB, NEJOB, IT, NMAYBE, ISNUM

      If (.NOT.refwfn_active) Then
        Write(6,*) ' refwfn not yet activated, aborting!'
        call abend
      End If

*---  Read the MO coefficients from HDF5/JOBIPH and store on LUONEM
      NCMO=NBSQT
      CALL GETMEM('LCMORAS','ALLO','REAL',LCMORAS,NCMO)
#ifdef _HDF5_
      If (refwfn_is_h5) Then
        call mh5_fetch_dset_array_real(refwfn_id,
     &         'MO_VECTORS', WORK(LCMORAS))
      Else
#endif
        IAD15=IADR15(9)
        IF(IFQCAN.EQ.0) IAD15=IADR15(2)
        CALL DDAFILE(refwfn_id,2,WORK(LCMORAS),NCMO,IAD15)
#ifdef _HDF5_
      End If
#endif
      IEOF1M=0
      IDISK=IEOF1M
      IAD1M(1)=IDISK
      CALL DDAFILE(LUONEM,1,WORK(LCMORAS),NCMO,IDISK)
      CALL GETMEM('LCMORAS','FREE','REAL',LCMORAS,NCMO)
      IEOF1M=IDISK

C IDCIEX: Present EOF on LUCIEX.
      ID=IDCIEX
C Skip when using cumulant reconstruction of (3-,) 4-RDM
      IF((.Not.DoCumulant).AND.(ISCF.EQ.0)) THEN
        CALL GETMEM('LCI','ALLO','REAL',LCI,NCONF)
        DO I=1,NSTATE
          ISNUM=MSTATE(I)
#ifdef _HDF5_
          If (refwfn_is_h5) Then
*---  Read the CI coefficients from the HDF5 file
            call mh5_fetch_dset_array_real(refwfn_id,'CI_VECTORS',
     $             Work(LCI),[nconf,1],[0,ISNUM-1])
          Else
#endif
*---  Read the CI coefficients from the JOBIPH file
            IDISK=IADR15(4)
            DO II=1,ISNUM-1
              CALL DDAFILE(refwfn_id,0,WORK(LCI),NCONF,IDISK)
            END DO
            CALL DDAFILE(refwfn_id,2,WORK(LCI),NCONF,IDISK)
#ifdef _HDF5_
          End If
#endif
C Copy selected vectors to LUCI:
          CALL DDAFILE(LUCIEX,1,WORK(LCI),NCONF,ID)
        END DO
C Disk address = present EOF on LUCIEX.
C IDTCEX = Disk address to transformed CI.
        IF(ORBIN.EQ.'TRANSFOR') THEN
          IDTCEX=ID
C Dummy writes:
          DO II=1,NSTATE
            CALL DDAFILE(LUCIEX,0,WORK(LCI),NCONF,ID)
          END DO
        ELSE
          IDTCEX=IDCIEX
        END IF
        CALL GETMEM('LCI','FREE','REAL',LCI,NCONF)
      ELSE
* If this is Closed-shell or Hi-spin SCF case
* Just in case...
        IF (.Not.DoCumulant .and. (NSTATE.ne.1 .or. NCONF.ne.1)) THEN
          write(6,*)' readin_caspt2: A Closed-shell or Hi-spin SCF'
          write(6,*)' but nr of states is: NSTATE=', NSTATE
          write(6,*)' and nr of CSFs is    NCONF= ', NCONF
          write(6,*)' Program error?? Must stop.'
          CALL ABEND
        END IF
* This should be solved elsewhere in the code...just for the now,
* make a write of a CI vector to LUCIEX, so other routines do not get
* their knickers into a twist:
        NCONF=1
        CALL GETMEM('LCI','ALLO','REAL',LCI,NCONF)
        WORK(LCI)=1.0D0
        CALL DDAFILE(LUCIEX,1,WORK(LCI),NCONF,ID)
        CALL GETMEM('LCI','FREE','REAL',LCI,NCONF)
      END IF
C Now, the selected original CASCI expansions are on LUCIEX
C beginning from disk address 0.

CSVC: read the L2ACT and LEVEL arrays
#ifdef _HDF5_
      If (refwfn_is_h5) Then
        call mh5_fetch_attr (refwfn_id,'L2ACT', L2ACT)
        call mh5_fetch_attr (refwfn_id,'A2LEV', LEVEL)
      Else
#endif
        IAD15=IADR15(18)
        CALL IDAFILE(refwfn_id,2,L2ACT,mxAct,IAD15)
        CALL IDAFILE(refwfn_id,2,LEVEL,mxAct,IAD15)
#ifdef _HDF5_
      End If
#endif

#ifdef _HDF5_
      If (refwfn_is_h5) Then
        call mh5_fetch_dset_array_real(refwfn_id,
     &         'ROOT_ENERGIES', ROOT_ENERGIES)
      Else
#endif
*PAM 2015: We will no longer recompute the RASSCF energies
* but just assume they can be obtained from the JOBIPH file.
* There is table with unknown length with all energies from
* all iterations (!) there.
        NEJOB=MXROOT*MXITER
        CALL GETMEM('EJOB','ALLO','REAL',LEJOB,NEJOB)
        IAD=IADR15(6)
        CALL DDAFILE(refwfn_id,2,WORK(LEJOB),NEJOB,IAD)
C Note that there is no info on nr of iterations
C so we cannot know what energies to pick...
C Let us make a guess: The correct set of energy values in the
C table of energies/iteration is the last one with not all zeroes.
        NMAYBE=0
        DO IT=1,MXITER
          AEMAX=0.0D0
          DO I=1,MXROOT
            E=WORK(LEJOB+MXROOT*(IT-1)+(I-1))
            AEMAX=MAX(AEMAX,ABS(E))
          END DO
          IF(ABS(AEMAX).LE.1.0D-12) GOTO 11
          NMAYBE=IT
        END DO
11      CONTINUE
        IF(NMAYBE.EQ.0) THEN
          WRITE(6,*)' PT2INI tried to read energies from the'
          WRITE(6,*)' JOBIPH file, but could not find any.'
          CALL ABEND()
        END IF
* And then put the energies into the Hamiltonian matrix,
* unless already filled in by the EFFE keyword
        DO I=1,mxRoot
          Root_Energies(I)=WORK(LEJOB+MXROOT*(NMAYBE-1)+(I-1))
        END DO
* No more use for the array EJOB.
        CALL GETMEM('EJOB','FREE','REAL',LEJOB,NEJOB)
#ifdef _HDF5_
      End If
#endif
      DO I=1,NSTATE
        REFENE(I)=ROOT_ENERGIES(MSTATE(I))
      END DO

* Read HEFF1 from JobMix file
      IF (IFSC) THEN
        CALL GETMEM('HEFF1','ALLO','REAL',LHEFF1,NSTATE**2)
        IAD15=IADR15(17)
        CALL DDAFILE(refwfn_id,2,WORK(LHEFF1),NSTATE**2,IAD15)
        DO I=1,NSTATE
          DO J=1,NSTATE
            HEFF1(I,J)=WORK(LHEFF1-1+I+NSTATE*(J-1))
          END DO
        END DO
        CALL GETMEM('HEFF1','FREE','REAL',LHEFF1,NSTATE**2)
      END IF

      End Subroutine

      End Module
