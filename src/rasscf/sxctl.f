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
* Copyright (C) 1989, Bjorn O. Roos                                    *
*               1989, Per Ake Malmqvist                                *
*               1991,1993,1996, Markus P. Fuelscher                    *
************************************************************************
      Subroutine SXCtl(CMO,OCC,D,P,PA,FI,FA,D1A,THMAX,IFINAL)
********************************************************************
*                                                                  *
* super-CI control section                                         *
*                                                                  *
* calling arguments:                                               *
* CMO     : array of real*8                                        *
*           MO-coefficients                                        *
* OCC     : array of real*8                                        *
*           orbital occupation numbers                             *
* D       : array of real*8                                        *
*           averaged one-body density matrix                       *
* P       : array of real*8                                        *
*           averaged two-body density matrix                       *
* PA      : array of real*8                                        *
*           averaged antisymmetric two-body density matrix         *
* FI      : array of real*8                                        *
*           inactive Fock matrix                                   *
* FA      : array of real*8                                        *
*           active Fock matrix                                     *
* D1A     : array of real*8                                        *
*           active one-body density matrix in AO basis             *
* IFINAL  : integer                                                *
*           termination flag                                       *
*                                                                  *
*------------------------------------------------------------------*
*                                                                  *
* written by:                                                      *
* B.O. Roos and P.Aa. Malmqvist                                    *
* University of Lund, Sweden, 1989                                 *
*                                                                  *
*------------------------------------------------------------------*
*                                                                  *
* history:                                                         *
* - updated for MOLCAS version 2                                   *
*   M.P. Fuelscher, University of Lund, Sweden, 1991               *
* - updated for MOLCAS version 3                                   *
*   M.P. Fuelscher, University of Lund, Sweden, 1993               *
* - updated for integral direct and reaction field calculations    *
*   M.P. Fuelscher, University of Lund, Sweden, 1996               *
*                                                                  *
********************************************************************

#ifdef _DMRG_
      use qcmaquis_interface_cfg
#endif
      use fciqmc, only : DoNECI
      Implicit Real*8 (A-H,O-Z)

      Dimension CMO(*),OCC(*),D(*),P(*),PA(*),FI(*),FA(*),D1A(*)

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "input_ras.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='SXCTL   ')
#include "WrkSpc.fh"
#include "wadr.fh"
#include "raswfn.fh"
      Character*4 Word
* PAM 2008 IndType, VecTyp added, see below at call to WrVec
      Integer IndType(56)
      Character*80 VecTyp
      Save nCall
      Logical DoActive,DoQmat,DoCholesky,TraOnly,l_casdft
      Integer ALGO
      Dimension P2act(1),CIDUMMY(1)

      COMMON /CHOTODO /DoActive,DoQmat,ipQmat
      COMMON /CHLCAS /DoCholesky,ALGO
#ifndef _DMRG_
      logical :: doDMRG = .false.
#endif
      ipDMAT=ip_Dummy
      nDMAT = 1

C PAM01 The SXCI part has been slightly modified by P-AA M Jan 15, 2001:
C Changes affect several of the subroutines of this part.
C The program should now behave more gracefully when there are (almost)
C redundant orbital rotation modes.
C For individual 2x2 rotations that are to be excluded, the SX hamiltonian
C diagonal, used in the preconditioning step in DAVCRE, is set to a huge
C value. Corresponding elements of trial vectors and sigma vectors are
C zeroed.
C To take care of 'hidden' redundant linear combinations, a small extra
C term, proportional to the square norm of the orbital rotations, is
C added to the optimization problem: this means in practice that a small
C extra term is added in the SIGVEC routine.
C Also, an extra term is added to the overlap calculation. The
C proportionality factor in SIGVEC should be small, but larger than the
C one used in the COVLP routine.
C Presently I try the following values:
C Huge elements in diagonal values: 1.0D32
C (Tested against 1.0D20 in IF-statements)
C Extra term in overlaps (COVLP, SXHAM): 1.0D-14
C Extra term in SIGVEC:                  1.0D-12


      Call qEnter(ROUTINE)
C Local print level (if any)
      IPRLEV=IPRLOC(4)
c      write(6,*) 'Entering SXCTL!'
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
3333  FORMAT(1X,'WORK SPACE VARIABLES IN SUBR. SXCTL: ',/,
     &       1X,'SUBSECTION: ',A6,/,(1X,12I10,/))

C --- Check for Cholesky ---------------

c      Call DecideOnCholesky(DoCholesky)
C --------------------------------------

* compute constants needed for addressing
      NLX1=0
      NLX2=0
      NQ=0
      NSXS=0
      NIAIA=0
      NAEAE=0
      NAOAE=0
      MNO=0
      MIAAE=0
      NLX=0
      DO 10 ISYM=1,NSYM
       NB=NBAS(ISYM)
       NA=NASH(ISYM)
       NI=NISH(ISYM)
       NE=NSSH(ISYM)
       NO=NORB(ISYM)
       NIA=NI+NA
       NAE=NA+NE
       NSXS=NSXS+NIA*NAE
       NIAIA=NIAIA+NIA**2
       NAEAE=NAEAE+NAE**2
       NAOAE=NAOAE+NA*NAE
       NLX1=MAX(NLX1,NO**2)
       NLX2=MAX(NLX2,NB*NO)
       NQ=MAX(NQ,NA*NO)
       MNO=MAX(MNO,NO)
       MIAAE=MAX(MIAAE,NIA*NAE)
       NLX=MAX(NLX,NIA*NA,NIA*NAE)
10    CONTINUE
      IF(NQ.LT.NIAIA) NQ=NIAIA
      NROOT=1
      NDIMSX=NROOT+NSXS
************************************************************************
* load back two-electron integrals (pu|vx)
************************************************************************
      lPUVX = 1
      If (.not.DoCholesky .or. ALGO.eq.1) Then
        If ( nFint.gt.0) then
          iDisk = 0
          Call GetMem('PUVX','Allo','Real',lPUVX,nFint)
          Call DDaFile(LUINTM,2,Work(lPUVX),nFint,iDisk)
        EndIf
      End If
      IF(IPRLEV.ge.DEBUG) THEN
       write(6,*) 'PUVX integrals in SXCTL'
       call wrtmat(Work(lPUVX),1,nFInt,1,nFInt)
      END IF
*********************************************************************************
* update and transform the Fock matrices FI and FA in MO basis ----> Fmat routine
*********************************************************************************
      If (.not.DoCholesky .or. ALGO.eq.1) Then
         WORD='FMAT'
         IF(IPRLEV.GE.DEBUG) THEN
           Write(LF,3333)WORD,lPUVX
         END IF
         Call Fmat(CMO,Work(lPUVX),D,D1A,FI,FA)

      ElseIf (ALGO.eq.2) Then

*     Inactive-active contribution to ECAS
         Call GetMem('DALT','Allo','Real',iDA,nTot1)
         Call Fold(nSym,nBas,D1A,Work(iDA)) !get the packed DA
         VIA=dDot_(nTot1,FI,1,Work(iDA),1)
         ECAS=EMY+VIA
         If ( IPRLEV.ge.DEBUG ) then
           Write(LF,'(A,E20.10)') ' Total core energy:            ',EMY
           Write(LF,'(A,E20.10)') ' inactive-active interaction:  ',VIA
           Write(LF,'(A,E20.10)') ' CAS energy (core+interaction):',ECAS
         End If
         Call GetMem('DALT','Free','Real',iDA,nTot1)

         TraOnly=.true.
         Call CHO_CAS_DRV(irc,CMO,D,FI,D1A,FA,P,TraOnly)

         if (irc.ne.0) then
         Write(LF,*)'SXCTL: Cho_cas_drv non-zero return code! rc= ',irc
         call qtrace()
         call abend()
         endif

      Else

         Write(LF,*)'SXCTL: Illegal Cholesky parameter ALGO= ',ALGO
         call qtrace()
         call abend()

      EndIf
************************************************************************
* reorder the two-body density matrix P
************************************************************************
      LP=1
      l_casdft = KSDFT(1:5).eq.'TLSDA'   .or.
     &           KSDFT(1:6).eq.'TLSDA5'  .or.
     &           KSDFT(1:5).eq.'TBLYP'   .or.
     &           KSDFT(1:6).eq.'TSSBSW'  .or.
     &           KSDFT(1:5).eq.'TSSBD'   .or.
     &           KSDFT(1:5).eq.'TS12G'   .or.
     &           KSDFT(1:4).eq.'TPBE'    .or.
     &           KSDFT(1:5).eq.'FTPBE'   .or.
     &           KSDFT(1:7).eq.'TREVPBE' .or.
     &           KSDFT(1:8).eq.'FTREVPBE'.or.
     &           KSDFT(1:6).eq.'FTLSDA'  .or.
     &           KSDFT(1:6).eq.'FTBLYP'
      IF(.not.l_casdft) then
* ISTORP(NSYM+1) represents the size of the 2-body density matrix,d(vwxy), with vwxy all active.
* the size is computed as NAP*NAQ*NRS (sum over all symmetries). If Sym_R = Sym_S then triangular
* form over NRS... with R.ge.S, rectanguar otherwise.
       IF(ISTORP(NSYM+1).GT.0) THEN
         WORD='FMAT'
         CALL GETMEM('ISTRP','ALLO','REAL',LP,ISTORP(NSYM+1))
         IF(IPRLEV.GE.DEBUG) THEN
           Write(LF,3333)WORD,LP
         END IF

         CALL PMAT_RASSCF(P,WORK(LP))

         If (ExFac.ne.1.0D0 .and.KSDFT(1:3).ne.'SCF') Then
            CALL GETMEM('P2_reo','ALLO','REAL',
     &                               ipP2reo,ISTORP(NSYM+1))
            Call Get_Temp('nP2Act  ',P2Act,1)
            nP2Act=Int(P2Act(1))
            CALL GETMEM('P2RAW','ALLO','REAL',ipP2_RAW,nP2Act)
            Call Get_Temp('P2_RAW  ',Work(ipP2_RAW),nP2Act)
            CALL PMAT_RASSCF(Work(ipP2_RAW),WORK(ipP2reo))
            CALL GETMEM('P2RAW','FREE','REAL',ipP2_RAW,nP2Act)
            P2reo=DBLE(ISTORP(NSYM+1))
            Call Put_Temp('nP2reo  ',[P2reo],1)
            Call Put_Temp('P2_reo  ',Work(ipP2reo),ISTORP(NSYM+1))
            CALL GETMEM('P2_reo','FREE','REAL',
     &                               ipP2reo,ISTORP(NSYM+1))
         End If
       END IF
      ELSE ! GLM-CASDFT
* ISTORP(NSYM+1) here represents the size of the Dvw*Dxy array (product of one-body
* density matrix,d(vwxy), with vwxy all active. The size is computed as NAP*NAQ*NRS
* (sum over all symmetries). If Sym_R = Sym_S then triangular form over NRS...
* with R.ge.S, rectanguar otherwise. Basically we will use same simmetry as for dvwxy.
       IF(ISTORP(NSYM+1).GT.0) THEN
c         Write(LF,*)
c         Write(LF,*) ' ---------------------'
C         Call Get_D1MO(ipDMAT,nDmat)
c         CALL TRIPRT('Averaged 1-body Dmat D in MO in SXCTL',' ',D,NAC)
         CALL GETMEM('ISTRP','ALLO','REAL',LP,ISTORP(NSYM+1))
         CALL DmatDmat(D,WORK(LP))
C         If(ipDMAT.ne.ip_Dummy) Call Free_Work(ipDMAT)
       END IF
      end if
************************************************************************
* Compute the MCSCF generalized Fock matrix and Brillouin matrix elements
************************************************************************
      WORD='FOCK'
      CALL GETMEM('FOCK','ALLO','REAL',LFOCK,NTOT4)
      CALL GETMEM('SXBM','ALLO','REAL',LBM,NSXS)
      CALL GETMEM('SXLQ','ALLO','REAL',LQ,NQ) ! q-matrix(1symmblock)
      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,3333)WORD,LFOCK,LBM,LP,LQ
      END IF
      CALL FOCK(WORK(LFOCK),WORK(LBM),FI,FA,
     &          D,WORK(LP),WORK(LQ),WORK(LPUVX),IFINAL,CMO)
c Now FA = FI + FA. Original FA has been overwritten in FOCK routine.
      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,*)
        Write(LF,*)'FI+FA in MO-basis in sxctl (overwritten on FA)'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FA(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
      END IF
      CALL GETMEM('FOCK','FREE','REAL',LFOCK,NTOT4)
      CALL GETMEM('SXLQ','FREE','REAL',LQ,NQ)
cGLM      If(KSDFT(1:3).ne.'GLM') then
       IF(ISTORP(NSYM+1).GT.0) THEN
         CALL GETMEM('ISTRP','FREE','REAL',LP,ISTORP(NSYM+1))
       END IF
cGLM      end if

      If (.not.DoCholesky .or. ALGO.eq.1) Then
         If ( nFint.gt.0 ) Then
            Call GetMem('PUVX','Free','Real',lPUVX,nFint)
         EndIf
      EndIf

* PAM 2008: Orbital files should be updated each iteration
* for easy access in case of catastrophic failure.
      IF ( IFINAL.NE.1 ) THEN
       iShift=0
       DO ISYM=1,NSYM
         IndT=0
         IndType(1+iShift)= NFRO(ISYM)
         IndT=IndT+NFRO(ISYM)
         IndType(2+iShift)= NISH(ISYM)
         IndT=IndT+NISH(ISYM)
         IndType(3+iShift)= NRS1(ISYM)
         IndT=IndT+NRS1(ISYM)
         IndType(4+iShift)= NRS2(ISYM)
         IndT=IndT+NRS2(ISYM)
         IndType(5+iShift)= NRS3(ISYM)
         IndT=IndT+NRS3(ISYM)
         IndType(7+iShift)= NDEL(ISYM)
         IndT=IndT+NDEL(ISYM)
         IndType(6+iShift)= NBAS(ISYM)-IndT
         iShift=iShift+7
        EndDo
* Note: This is not the final orbitals, and the orbital energies and
* active occupation numbers may be meaningless.
* There is an array with occupation numbers, so use it, even if
* possibly irrelevant. But put zeroes as orbital energies:
        Call GetMem('EDUM','ALLO','REAL',LEDUM,NTOT)
        call dcopy_(NTOT,[0.0D0],0,WORK(LEDUM),1)

        IF (DoNECI) THEN
          write(6,*)'For NECI orbital energies are approximated'
          write(6,*)'to the diagonal value of the Fock matrix (SXCTL)'
          write(6,*)'These values are going to the temporary RasOrb'
          write(6,*)'together with CMOs before orb rot is performed'
c This is of course not true other than for special cases ... but some might find it beneficial!
          iOff  = 0
          iOff2 = 0
          Do iSym = 1,nSym
           iBas = nBas(iSym)
           if(iBas.gt.0) then
            do iDiag = 1,iBas
             WORK(LEDUM+iDiag-1+iOff) = FA(iOff2+(iDiag*(iDiag+1)/2))
            end do
            write(6,*) 'OrbEn in line for sym = ',iSym
            write(6,*) (Work(LEDUM+i+iOff), i = 0,iBas-1)
             iOff  = iOff + iBas
             iOff2 = iOff2 + (iBas*iBas+iBas)/2
           end if
          End Do
        END IF
        Write(VecTyp,'(A)')
        VecTyp='* RASSCF average (pseudo-natural) orbitals (Not final)'
        LuvvVec=50
        LuvvVec=isfreeunit(LuvvVec)
        Call WrVec('RASORB',LuvvVec,'COE',NSYM,NBAS,NBAS,
     &           CMO, OCC, WORK(LEDUM), INDTYPE,VECTYP)
        Call WrVec('RASORB',LuvvVec,'AI',NSYM,NBAS,NBAS,
     &           CMO, OCC, WORK(LEDUM), INDTYPE,VECTYP)
        Call GetMem('EDUM','FREE','REAL',LEDUM,NTOT)

#ifdef _HDF5_
        call mh5_put_dset(wfn_mocoef,CMO)
        call mh5_put_dset(wfn_occnum,OCC)
#endif
      END IF

      IF ( IFINAL.EQ.1 ) THEN
C If ifinal=1 , this is the last call to SXCTL for calculating
C the MCSCF Fock matrix for occupied orbitals, and new orbitals).
C First generate final orbitals:
C Diagonalize inactive and secondary part of FP = FI + FA.
C         diagonal blocks of the active density matrix
C for RAS1, RAS2, and RAS3. All this is done in NEWORB.
C Finally generate orbitals and Fock matrices for CAS-PT2 runs.
C in FCKPT2

C Memory allocation and call to NEWORB and FCKPT2
C LCMON: New molecular orbitals (NTOT2)
C LFTR:  Temporary area for part of the Fock matrix FP (NTOT1)
C LVEC:  EIGENVECTORS OF FTR (NO2M)
C LSQ and LWO: scratch areas

        WORD='FPT2'
        CALL GETMEM('XCMO','ALLO','REAL',LCMON,NTOT2)
        CALL GETMEM('XFTR','ALLO','REAL',LFTR,NTOT1)
        CALL GETMEM('XVEC','ALLO','REAL',LVEC,NTOT2)
        CALL GETMEM('SXWO','ALLO','REAL',LWO,NTOT2)
        CALL GETMEM('SXSQ','ALLO','REAL',LSQ,NTOT2)
        CALL GETMEM('SXMX','ALLO','REAL',LCMOX,NTOT2)
        IF(IPRLEV.GE.DEBUG) THEN
          Write(LF,3333) WORD,LCMON,LFTR,LVEC,LWO,LSQ,LCMOX
        END IF
        If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*)
         Write(LF,*) ' CMO in SXCTL for IFINAL=1'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         ioff=0
         Do iSym = 1,nSym
          iBas = nBas(iSym)
          if(iBas.ne.0) then
            write(6,*) 'Sym =', iSym
            do i= 1,iBas
              write(6,*) (CMO(ioff+iBas*(i-1)+j),j=1,iBas)
            end do
            iOff = iOff + (iBas*iBas)
          end if
         End Do
        End If

        If ( iCIonly.eq.1 ) then
          IDISK=IADR15(2)
          CALL DDAFILE(JOBIPH,1,CMO,NTOT2,IDISK)
          CALL DDAFILE(JOBIPH,1,OCC,NTOT,IDISK)
#ifdef _HDF5_
          call mh5_put_dset(wfn_mocoef,CMO)
          call mh5_put_dset(wfn_occnum,OCC)
#endif
        Else
!          this part (TRACI) need to be changed to "TRAMPS", not yet ! Yingjin
           CALL NEWORB_RASSCF(CMO,WORK(LCMON),FA,WORK(LFTR),WORK(LVEC),
     &                        WORK(LWO),WORK(LSQ),WORK(LCMOX),D,OCC)
* compute orbital overlap matrix
c           IF (NACTEL.GT.0) THEN
* NN.14 Skip this when DMRG-CASSCF due to CI-vector dependency
           !IF(.NOT.(DoDMRG.or.doBlockDMRG).AND.NACTEL.GT.0) THEN
           IF(NACTEL.GT.0)THEN
             CALL GETMEM('SMAT','ALLO','REAL',LSMAT,NAC*NAC)
             IWAY = 1
             CALL OVLP(IWAY,CMO,WORK(LCMON),WORK(LSMAT))

             if(dodmrg)then
#ifdef _DMRG_
#ifdef BLUBB
               call mpsrot(work(lsmat),nac,nrs2,nsym)
#endif
#endif
             else if(doBlockDMRG .or. DoNECI)then
             else !CI
               IDISK=IADR15(4)
               CALL LUCIA_UTIL('TRACI',IDISK,JOBIPH,WORK(LSMAT))
             end if
             CALL GETMEM('SMAT','FREE','REAL',LSMAT,NAC*NAC)
           ELSE
             CIDUMMY=1.0D0
             IDISK=IADR15(4)
             CALL DDAFILE(JOBIPH,1,CIDUMMY,1,IDISK)
           END IF
        End If

* IPT2 = 1 for OUTO, CANOnical option...
        IF(IPT2.NE.0)
     &  CALL FCKPT2(CMO,WORK(LCMON),FI,FA,
     &              WORK(LFTR),WORK(LVEC),
     &              WORK(LWO),WORK(LSQ),WORK(LCMOX))

        CALL GETMEM('SXMX','FREE','REAL',LCMOX,NTOT2)
        CALL GETMEM('SXSQ','FREE','REAL',LSQ,NTOT2)
        CALL GETMEM('SXWO','FREE','REAL',LWO,NTOT2)
        CALL GETMEM('XVEC','FREE','REAL',LVEC,NTOT2)
        CALL GETMEM('XFTR','FREE','REAL',LFTR,NTOT1)
        CALL GETMEM('XCMO','FREE','REAL',LCMON,NTOT2)

        CALL TIMING(CPTS,CPES,TIOS,TIOES)

        GOTO 9990
      ENDIF

C Memory allocation and calling sequence for SXHAM
C LSXN: Normalization constants for super-CI vector
C LF1 and LF2: parts of the Fock matrix FP
C LDIA: Occupied part of the density matrix (squared)
C LG: The G matrix(used in sigvec)
C LH: The H matrix( "    "   "   )
C LHD: The diagonal of the super-CI Hamiltonian
C LDF: The matrix D*FP
C LDDIA: Diagonal of the density matrix (all elements one symmetry)

      WORD='SXHA'
      LH=1
      CALL GETMEM('XSXN','ALLO','REAL',LSXN,NSXS)
      CALL GETMEM('SXF1','ALLO','REAL',LF1,NIAIA)
      CALL GETMEM('SXF2','ALLO','REAL',LF2,NAEAE)
      CALL GETMEM('XDIA','ALLO','REAL',LDIA,NIAIA)
      CALL GETMEM('SXG1','ALLO','REAL',LG,NIAIA)
      IF(NAOAE.GT.0) CALL GETMEM('SXH1','ALLO','REAL',LH,NAOAE)
      CALL GETMEM('SXHD','ALLO','REAL',LHD,NDIMSX)
      CALL GETMEM('SXDF','ALLO','REAL',LDF,NQ)
      CALL GETMEM('SXDD','ALLO','REAL',LDDIA,MNO)
      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,3333)WORD,LSXN,LF1,LF2,LDIA,LG,
     &                               LH,LHD,LDF,LDDIA
      END IF

c         CALL TRIPRT(' Dmat in MO in SXCTL bf call to SXHAM ',' ',D,NAC)
c         CALL TRIPRT(' Pmat in MO in SXCTL bf call to SXHAM ',
c     &              ' ',P,NACPAR)
c         CALL TRIPRT(' PAmat in MO in SXCTL bf call to SXHAM',
c     &              ' ',PA,NACPAR)
      CALL SXHAM(D,P,PA,FA,WORK(LSXN),
     &               WORK(LF1),WORK(LF2),WORK(LDIA),WORK(LG),
     &               WORK(LH),WORK(LHD),WORK(LDF),WORK(LDDIA))

      CALL GETMEM('SXDD','FREE','REAL',LDDIA,MNO)
      CALL GETMEM('SXDF','FREE','REAL',LDF,NQ)

C PAM01 Removal of certain rotations from the BLB elements.
C Some additional rotations (besides those listed in IZROT) may
C need to be suppressed. These are rotations that are (very nearly)
C redundant -- they hardly affect the wave function at all.
C All suppressed rotations can be identified because the corresponding
C diagonal elements have been set to a huge number in SXHAM.
C Use this criterion to set some BLB elements exactly =0:
      DO I=1,NSXS
       IF(WORK(LHD+NROOT-1+I).GT.1.0D20) THEN
        WORK(LBM-1+I)=0.0D0
       END IF
      END DO

C MEMORY ALLOCATION AND CALLING SEQUENCE FOR SX DIAGONALIZATION

C LCSX: The super-CI vectors
C LSIGMA: The sigma vectors
C LHH:  The Davidson H matrix
C LCC:   "     "     egenvectors
C LENER: "     "     energies
C LSC:   Scratch area
C LQ:    Davidson update vectors
C LQQ:   Norm of update vectors
C LOVL:  Overlap matrix

      WORD='SXDA'
      NCR=NDIMSX*NROOT*ITMAX
      KMAX=ITMAX*NROOT
      NCR1=NDIMSX*NROOT*(ITMAX+1)
      CALL GETMEM('XCSX','ALLO','REAL',LCSX,NCR1)
      CALL GETMEM('XSIG','ALLO','REAL',LSIGMA,NCR)
      NLHH=KMAX**2+KMAX
      NLCC=KMAX**2
      NLQ=NDIMSX*(NROOT+1)
      NLOVL=ITMAX*NROOT**2
      CALL GETMEM('SXHH','ALLO','REAL',LHH,NLHH)
      CALL GETMEM('SXCC','ALLO','REAL',LCC,NLCC)
      CALL GETMEM('ENER','ALLO','REAL',LENER,KMAX)
      CALL GETMEM('SXSC','ALLO','REAL',LSC,NDIMSX)
      CALL GETMEM('SXLQ','ALLO','REAL',LQ,NLQ)
      CALL GETMEM('SXQQ','ALLO','REAL',LQQ,NROOT)
      CALL GETMEM('XOVL','ALLO','REAL',LOVL,NLOVL)
      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,3333)WORD,LCSX,LSIGMA,LHH,LCC,LENER,
     &         LHD,LBM,LSC,LQ,LQQ,LOVL
      END IF

      CALL DAVCRE(WORK(LCSX),WORK(LSIGMA),WORK(LHH),WORK(LCC),
     &            WORK(LENER),WORK(LHD),WORK(LSC),
     &            WORK(LQ),WORK(LQQ),WORK(LOVL),SXSEL,
     &            NROOT,ITMAX,NDIMSX,ITERSX,NSXS)

      ESX=WORK(LENER)
      CALL GETMEM('XSIG','FREE','REAL',LSIGMA,NCR)
      CALL GETMEM('SXHH','FREE','REAL',LHH,NLHH)
      CALL GETMEM('SXCC','FREE','REAL',LCC,NLCC)
      CALL GETMEM('ENER','FREE','REAL',LENER,KMAX)
      CALL GETMEM('SXSC','FREE','REAL',LSC,NDIMSX)
      CALL GETMEM('SXLQ','FREE','REAL',LQ,NLQ)
      CALL GETMEM('SXQQ','FREE','REAL',LQQ,NROOT)
      CALL GETMEM('XOVL','FREE','REAL',LOVL,NLOVL)
      CALL GETMEM('SXF1','FREE','REAL',LF1,NIAIA)
      CALL GETMEM('SXF2','FREE','REAL',LF2,NAEAE)
      CALL GETMEM('SXG1','FREE','REAL',LG,NIAIA)
      IF(NAOAE.GT.0) CALL GETMEM('SXH1','FREE','REAL',LH,NAOAE)
      CALL GETMEM('SXHD','FREE','REAL',LHD,NDIMSX)

C Renormalize the SX-coefficients

      IREF=1
      LCSXI=LCSX+NDIMSX*(IREF-1)
      IC=NROOT+LCSXI-1
      XSXMAX=0.0D0
      DO 54 I=1,NSXS
       WORK(IC+I)=WORK(LSXN+I-1)*WORK(IC+I)/WORK(LCSXI)
       XSXMAX=MAX(XSXMAX,abs(WORK(IC+I)))
54    CONTINUE
      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,*) 'SXCTL after DAVCRE, Renormalized SX coeffs:'
        Write(LF,'(1X,8F14.6)')(WORK(IC+I),I=1,NSXS)
      END IF

* Step size control has been built into qune now.
*C Step length control, just for safety.
*      DO I=1,NSXS
*        VAL=WORK(IC+I)
*        WORK(IC+I)=VAL/(1.0D0+1.7D0*ABS(VAL))
*      END DO

C Intercept XSX and BM, to use (perhaps) Quasi-Newton or Line Search

      WORD='QUNE'
*      IF(ITER.EQ.1) NCALL=0
      IF(ITER.LE.4) NCALL=0
      IF(KeyHEUR.AND.ITER.GT.10.AND.MOD(ITER,10).LT.4) NCALL=0
      if(doDMRG) then
        IF(ITER.LE.2) NCALL=0  ! YM: change 4 -> 2, for saving time
      end if
      IF(XSXMAX.GT.0.5D0) NCALL=0
      IF(NQUNE.NE.0.AND.XSXMAX.LT.0.5D0) THEN
        CALL GETMEM('SXVT','ALLO','REAL',LVT,NSXS)
        CALL GETMEM('SXVL','ALLO','REAL',LVL,NSXS)
        CALL GETMEM('SXQN','ALLO','REAL',LXQN,NSXS)
        CALL GETMEM('SXSC','ALLO','REAL',LSCR,NSXS)
        CALL GETMEM('XV11','ALLO','REAL',LV1,NSXS)
        CALL GETMEM('XV22','ALLO','REAL',LV2,NSXS)
        IF(IPRLEV.GE.DEBUG) THEN
          Write(LF,3333)WORD,LBM,(NROOT+LCSXI),LVL,LVT,
     &                               LXQN,LSCR,LV1,LV2
        END IF
        CASDFT_En=0.0d0
        If(KSDFT.ne.'SCF'.and.KSDFT(1:3).ne.'PAM')
     &      Call Get_dScalar('CASDFT energy',CASDFT_En)
        CASDFT_En=ECAS+CASDFT_En
        CALL QUNE(NCALL,CASDFT_En,WORK(LBM),WORK(NROOT+LCSXI),
     &            WORK(LVL),WORK(LVT),WORK(LXQN),WORK(LSCR),
     &            WORK(LV1),WORK(LV2),NSXS,LUQUNE,
     &            TMIN,QNSTEP,QNUPDT,KSDFT)

        CALL GETMEM('XXXX','FREE','REAL',LVT,NSXS)
        CALL GETMEM('XXXX','FREE','REAL',LVL,NSXS)
        CALL GETMEM('XXXX','FREE','REAL',LXQN,NSXS)
        CALL GETMEM('XXXX','FREE','REAL',LSCR,NSXS)
        CALL GETMEM('XXXX','FREE','REAL',LV1,NSXS)
        CALL GETMEM('XXXX','FREE','REAL',LV2,NSXS)
      ENDIF

C Rotation of orbitals with exp(x) where x is obtained from
C the super-CI coefficients, with a Quasi Newton update (NQUNE=1)


C CMO:  before - old MO's           after - new MO's
C LCMON: intermediate storage for new MO's (moved to CMO in ORTHO)
C LVEC:  eigenvectors of exp(X)
C LX2:  work area, also in ORTHO (AO overlap matrix)
C LWSQ:  "     "     "    "   "
C LY,LA, AND LB WORK AREAS

      WORD='ROTO'
      CALL GETMEM('CMO1','ALLO','REAL',LCMON,NTOT2)
      CALL GETMEM('VEC1','ALLO','REAL',LVEC,NO2M)
      CALL GETMEM('XMAT','ALLO','REAL',LXMAT,NO2M)
      CALL GETMEM('SXX2','ALLO','REAL',LX2,NTOT1)
      CALL GETMEM('SXY2','ALLO','REAL',LY,NO2M)
      CALL GETMEM('SXA1','ALLO','REAL',LA,MNO)
      CALL GETMEM('SXB2','ALLO','REAL',LB,MNO)
      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,3333)WORD,LCMON,LSXN,LCSXI,LXMAT,LX2,
     &                          LY,LVEC,LA,LB
      END IF

      CALL ROTORB(CMO,WORK(LCMON),WORK(LCSXI),WORK(LXMAT),
     &       WORK(LX2),WORK(LY),WORK(LVEC),WORK(LA),WORK(LB),THMAX,FA)

      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,*)
        Write(LF,*)'FI+FA in SXCTL after Unitary transform in ROTORB'
        Write(LF,*) ' --------------'
        Write(LF,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          Call TriPrt(' ',' ',FA(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
      END IF
      CALL GETMEM('CMO1','FREE','REAL',LCMON,NTOT2)
      CALL GETMEM('VEC1','FREE','REAL',LVEC,NO2M)
      CALL GETMEM('XMAT','FREE','REAL',LXMAT,NO2M)
      CALL GETMEM('SXX2','FREE','REAL',LX2,NTOT1)
      CALL GETMEM('SXY2','FREE','REAL',LY,NO2M)
      CALL GETMEM('SXA1','FREE','REAL',LA,MNO)
      CALL GETMEM('SXB2','FREE','REAL',LB,MNO)
      CALL GETMEM('XSXN','FREE','REAL',LSXN,NSXS)
      CALL GETMEM('XDIA','FREE','REAL',LDIA,NIAIA)
      CALL GETMEM('XCSX','FREE','REAL',LCSX,NCR1)

      IDISK=IADR15(2)
      CALL DDAFILE(JOBIPH,1,CMO,NTOT2,IDISK)
      CALL DDAFILE(JOBIPH,1,OCC,NTOT,IDISK)
#ifdef _HDF5_
      call mh5_put_dset(wfn_mocoef,CMO)
      call mh5_put_dset(wfn_occnum,OCC)
#endif
      CALL TIMING(CPTS,CPES,TIOS,TIOES)

9990  CONTINUE
      CALL GETMEM('SXBM','FREE','REAL',LBM,NSXS)
      Call qExit(ROUTINE)
      RETURN
      END
