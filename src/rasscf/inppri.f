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
* Copyright (C) 1993, Markus P. Fuelscher                              *
************************************************************************
      Subroutine InpPri(lOPTO)
************************************************************************
*                                                                      *
*     Echo input                                                       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
#ifdef _DMRG_
!     module dependencies
      use qcmaquis_interface_cfg
      use qcmaquis_interface_utility_routines, only: print_dmrg_info
#endif
      use OneDat, only: sNoOri
      use fcidump, only : DumpOnly
      use fciqmc, only: DoNECI
      use CC_CI_mod, only: Do_CC_CI
      use Fock_util_global, only: DoLocK
      use Functionals, only: Init_Funcs, Print_Info
      use KSDFT_Info, only: CoefR, CoefX
      use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
      use rctfld_module, only: lRF
      use gas_data, only: iDoGAS, NGAS, NGSSH, IGSOCCX
      use Constants, only: Zero
      use rasscf_global, only: KSDFT, DoBlockDMRG, DoDMRG, ICICH,
     &                         ICIRST, iPCMRoot, iRLXRoot, iSupSM,
     &                         ITMAX, l_casdft, lRoots, lSquare, LvShft,
     &                         MAXIT, n_Det, NAC, NFR, NIN, NONEQ,
     &                         NROOTS, NSEC, nTit, RFPert, ThrE, ThrSX,
     &                         ThrTE, Tot_Charge, Tot_El_Charge,
     &                         Tot_Nuc_Charge, Title, Header, iRoot,
     &                         Weight, iCI, cCI, ixSym, NQUNE
#ifdef _ENABLE_DICE_SHCI_
      use rasscf_global, only: dice_eps1, dice_eps2, dice_iter,
     &                         dice_Restart, dice_SampleN,
     &                         Dice_Stoc, nRef_Dice, diceocc
#endif
#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
      use rasscf_global, only: MXDMRG, ChemPS2_blb, ChemPS2_lreStart,
     &                         ChemPS2_Noise, ChemPS2_Restart,
     &                         Davidson_tol, Do3RDM, HFOcc,
     &                         Max_canonical, Max_Sweep
#endif
      use SplitCas_Data, only: DoSPlitCas,MxIterSplit,ThrSplit,
     &                         lRootSplit,NumSplit,iDimBlockA,
     &                         EnerSplit,GapSpli,PerSplit,PerCSpli,
     &                         fOrdSplit
      use printlevel, only: USUAL,SILENT
      use output_ras, only: LF,IPRLOC
      use general_data, only: NACTEL,NHOLE1,NELEC3,ISPIN,STSYM,NSYM,
     &                        NSEL,NTOT1,NASH,NBAS,NDEL,NFRO,NISH,
     &                        NRS1,NRS2,NRS3,NSSH
      use spinfo, only: DoComb,NCNFTP,NCSASM,NDTASM,NDTFTP
      use spinfo, only: I_ELIMINATE_GAS_MOLCAS,NCSF_HEXS

      Implicit None
      Logical lOPTO

#include "rasdim.fh"
      Character(LEN=8)   Fmt1,Fmt2,Label
      Character(LEN=120)  Line,BlLine,StLine
      Character(LEN=3) lIrrep(8)
      Character(LEN=80) KSDFT2
#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
      Character(LEN=3) SNAC
      Integer iHFOcc
#endif
      Logical DoCholesky
#ifdef _DMRG_
      character(len=100) :: dmrg_start_guess
#endif
#ifdef _ENABLE_DICE_SHCI_
      Integer iRef_Dice
#endif
      Real*8, Allocatable:: Tmp0(:)
      Real*8 AvailMB, WillNeedMB
      Integer i, iCharge, iComp, iDoRI, iEnd, iGAS, iOpt, iPrLev, iRC,
     &        iRef, iStart, iSyLbl, iSym, iTemp, left, lLine, lPaper,
     &        MaxRem, n_paired_elec, n_unpaired_elec, nLine

* Print level:
      IPRLEV=IPRLOC(1)
*----------------------------------------------------------------------*
*     Start and define the paper width                                 *
*----------------------------------------------------------------------*
      lPaper=132
*----------------------------------------------------------------------*
*     Initialize blank and header lines                                *
*----------------------------------------------------------------------*
      lLine=Len(Line)
      Do i=1,lLine
        BlLine(i:i)=' '
        StLine(i:i)='*'
      End Do
      lPaper=132
      left=(lPaper-lLine)/2
      Write(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
      Write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
      IF (IPRLEV.EQ.SILENT) GOTO 900
*----------------------------------------------------------------------*
*     Print the project title                                          *
*----------------------------------------------------------------------*
      IF(IPRLEV.GE.USUAL) THEN
       If ( nTit.gt.0 ) then
         Write(LF,*)
         nLine=nTit+5
         Do i=1,nLine
           Line=BlLine
           If ( i.eq.1 .or. i.eq.nLine )
     &     Line=StLine
           If ( i.eq.3 )
     &     Line='Project:'
           If ( i.ge.4 .and. i.le.nLine-2 )
     &     Write(Line,'(A72)')Title(i-3)
           Call Center_Text(Line)
           Write(LF,Fmt1) '*'//Line//'*'
         End Do
       Write(LF,*)
       End If
      END IF
*----------------------------------------------------------------------*
*     Print the ONEINT file identifier                                 *
*----------------------------------------------------------------------*
      IF(IPRLEV.GE.USUAL .AND..NOT.lOPTO) THEN
       Write(LF,*)
       Write(LF,Fmt1) 'Header of the ONEINT file:'
       Write(LF,Fmt1) '--------------------------'
       Write(Line,'(36A2)') (Header(i),i=1,36)
       Write(LF,Fmt1)  trim(adjustl(Line))
       Write(Line,'(36A2)') (Header(i),i=37,72)
       Write(LF,Fmt1)  trim(adjustl(Line))
       Write(LF,*)
*----------------------------------------------------------------------*
*     Print the status of ORDINT                                       *
*----------------------------------------------------------------------*
       Write(LF,*)
       If (lSquare) Then
         Write(LF,Fmt1) 'OrdInt status: squared'
       Else
         Write(LF,Fmt1) 'OrdInt status: non-squared'
       End If
       Write(LF,*)
*----------------------------------------------------------------------*
*     Print cartesian coordinates of the system                        *
*----------------------------------------------------------------------*
       Call PrCoor
      END IF
*----------------------------------------------------------------------*
*     Print orbital and wavefunction specifications                    *
*----------------------------------------------------------------------*
      IF(IPRLEV.GE.USUAL .AND..NOT.lOPTO) THEN
      Write(LF,*)
      Line=' '
      Write(Line(left-2:),'(A)') 'Wave function specifications:'
      Call CollapseOutput(1,Line)
      Write(LF,Fmt1)'-----------------------------'
      Write(LF,*)
      If (NFR.gt.0)
     &Write(LF,Fmt2//'A,T45,I6)')'Number of frozen shell electrons',
     &                           2*NFR
      Write(LF,Fmt2//'A,T45,I6)')'Number of closed shell electrons',
     &                           2*NIN
      Write(LF,Fmt2//'A,T45,I6)')'Number of electrons in active shells',
     &                           NACTEL
C.. for RAS
      if (.not.idogas) then
      Write(LF,Fmt2//'A,T45,I6)')'Max number of holes in RAS1 space',
     &                           NHOLE1
      Write(LF,Fmt2//'A,T45,I6)')'Max nr of electrons in RAS3 space',
     &                           NELEC3
C.. for GAS
      else
        DO IGAS=1,NGAS
          Write(LF,Fmt2//'A,I1,A,T45,2I6)')
     &      'Min/Max nr of electrons up to GAS',IGAS,' sp.',
     &                           igsoccx(igas,1),igsoccx(igas,2)
        END DO
      end if

      If (NFR.gt.0)
     &Write(LF,Fmt2//'A,T45,I6)')'Number of frozen orbitals',
     &                           NFR
      Write(LF,Fmt2//'A,T45,I6)')'Number of inactive orbitals',
     &                           NIN
      Write(LF,Fmt2//'A,T45,I6)')'Number of active orbitals',
     &                           NAC
      Write(LF,Fmt2//'A,T45,I6)')'Number of secondary orbitals',
     &                           NSEC
      Write(LF,Fmt2//'A,T45,F6.1)')'Spin quantum number',
     &                           (DBLE(ISPIN-1))/2.0d0
      Write(LF,Fmt2//'A,T45,I6)')'State symmetry',
     &                           STSYM
      Call CollapseOutput(0,'Wave function specifications:')
*
      Call Get_cArray('Irreps',lIrrep,24)
      Do iSym = 1, nSym
         lIrrep(iSym) = adjustr(lIrrep(iSym))
      End Do
*
      Write(LF,*)
      Line=' '
      Write(Line(left-2:),'(A)') 'Orbital specifications:'
      Call CollapseOutput(1,Line)
      Write(LF,Fmt1)'-----------------------'
      Write(LF,*)
      Write(LF,Fmt2//'A,T47,8I4)') 'Symmetry species',
     &                            (iSym,iSym=1,nSym)
      Write(LF,Fmt2//'A,T47,8(1X,A))') '                ',
     &                            (lIrrep(iSym),iSym=1,nSym)
      Write(LF,Fmt2//'A,T47,8I4)') 'Frozen orbitals',
     &                            (nFro(iSym),iSym=1,nSym)
      Write(LF,Fmt2//'A,T47,8I4)') 'Inactive orbitals',
     &                            (nIsh(iSym),iSym=1,nSym)
      Write(LF,Fmt2//'A,T47,8I4)') 'Active orbitals',
     &                            (nAsh(iSym),iSym=1,nSym)
      IF(.not.iDoGas)then
        Write(LF,Fmt2//'A,T47,8I4)') 'RAS1 orbitals',
     &                            (nRs1(iSym),iSym=1,nSym)
        Write(LF,Fmt2//'A,T47,8I4)') 'RAS2 orbitals',
     &                            (nRs2(iSym),iSym=1,nSym)
        Write(LF,Fmt2//'A,T47,8I4)') 'RAS3 orbitals',
     &                            (nRs3(iSym),iSym=1,nSym)
      Else
        DO IGAS=1,NGAS
          Write(LF,Fmt2//'A,I1,A,T47,8I4)') 'GAS',IGAS,' orbitals',
     &                            (ngssh(igas,iSym),iSym=1,nSym)
        END DO
      End If

      Write(LF,Fmt2//'A,T47,8I4)') 'Secondary orbitals',
     &                            (nSsh(iSym),iSym=1,nSym)
      Write(LF,Fmt2//'A,T47,8I4)') 'Deleted orbitals',
     &                            (nDel(iSym),iSym=1,nSym)
      Write(LF,Fmt2//'A,T47,8I4)') 'Number of basis functions',
     &                            (nBas(iSym),iSym=1,nSym)
      Call CollapseOutput(0,'Orbital specifications:')
      Write(LF,*)

#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
      If(.Not.DoBlockDMRG) GoTo 113

#ifdef _ENABLE_DICE_SHCI_
      Line=' '
      Write(Line(left-2:),'(A)') 'DICE specifications:'
      Call CollapseOutput(1,Line)
      Write(LF,Fmt1)'--------------------------'
      Write(LF,*)
      Write(LF,Fmt2//'A,T70,L6)')'Heat-bath configuration interaction
     &(JCTC, 2017, 13, 1595)', DoBlockDMRG
      Write(LF,Fmt2//'A,T45,L6)')'Semistochastic algorithm',Dice_stoc
      Write(LF,Fmt2//'A,T45,L6)')'Full restart',dice_restart
      Write(LF,Fmt2//'A,T45,I6)')'Max iterations',dice_iter
      Write(LF,Fmt2//'A,T45,ES10.3)')'Epsilon1',
     &                           dice_eps1
      Write(LF,Fmt2//'A,T45,ES10.3)')'Epsilon2',
     &                           dice_eps2
      Write(LF,Fmt2//'A,T45,I6)')'SampleN',
     &                           dice_sampleN
      Write(LF,Fmt2//'A,T45)')'Occupation guess'
      do iref_dice=1,nref_dice
         write(LF,Fmt2//'A)') trim(diceocc(iref_dice))
      enddo
      Call CollapseOutput(0,'DICE specifications:')

*     Skip printing CI specifications in DICE
      GoTo 114
#endif

#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_)
      Line=' '
      Write(Line(left-2:),'(A)') 'DMRG sweep specifications:'
      Call CollapseOutput(1,Line)
      Write(LF,Fmt1)'--------------------------'
      Write(LF,*)
      Write(LF,Fmt2//'A,T45,I6)')'Number of renormalized basis',
     &                           MxDMRG
      Write(LF,Fmt2//'A,T45,I6)')'Number of root(s) required',
     &                           NROOTS

      Write(LF,Fmt2//'A,T45,I6)')'Maximum number of sweeps',
     &                           max_sweep
      Write(LF,Fmt2//'A,T45,I6)')'Maximum number of sweeps in RDM',
     &                           max_canonical
      Write(LF,Fmt2//'A,T45,ES10.3)')'Threshold for restarting',
     &                           chemps2_blb
      Write(LF,Fmt2//'A,T45,ES10.3)')'Minimum Davidson tolerance',
     &                           davidson_tol
      Write(LF,Fmt2//'A,T45,ES10.3)')'DMRG convergence threshold',
     &                           THRE/2.0
      Write(LF,Fmt2//'A,T45,ES10.3)')'Noise prefactor',
     &                           chemps2_noise
      Write(LF,Fmt2//'A,T45,L6)')'Restart from previous calculation',
     &                           chemps2_restart
      Write(LF,Fmt2//'A,T45,L6)')'Calculate 3-RDM and F.4-RDM',
     &                           Do3RDM
      Write(LF,Fmt2//'A,T45,I6)')'Restart scheme in 3-RDM and F.4-RDM',
     &                           chemps2_lrestart
      write(SNAC, '(I3)') NAC
      Write(LF,Fmt2//'A,T45,'//trim(adjustl(SNAC))//'I2)')
     &                           'Occupation guess',
     &                           (HFOCC(ihfocc), ihfocc=1,NAC)

* NN.14 FIXME: haven't yet checked whether geometry opt. works correctly with DMRG
      Write(LF,Fmt2//'A,T45,I6)')'Root chosen for geometry opt.',
     &                           IRLXROOT
      Call CollapseOutput(0,'DMRG sweep specifications:')

*     Skip printing CI specifications in DMRG-CASSCF
      GoTo 114

 113  Continue
#endif
#endif

      Line=' '
      if(doDMRG)then
        Write(Line(left-2:),'(A)') 'DMRG specifications:'
      else
        Write(Line(left-2:),'(A)') 'CI expansion specifications:'
      end if
      Call CollapseOutput(1,Line)
      Write(LF,Fmt1)'----------------------------'
      Write(LF,*)

      if(doDMRG)then  !> Information for QCMaquis-DMRG
#ifdef _DMRG_
       if(dmrg_orbital_space%initial_occ(1,1) > 0)then
         dmrg_start_guess = "Single determinant"
       else
         if(dmrg_warmup%doCIDEAS)then
           dmrg_start_guess = "CI-DEAS"
         else
           dmrg_start_guess = "Random numbers (default)"
         end if
       end if
       call print_dmrg_info(lf,fmt2,1,dmrg_start_guess,nroots,thre)
#endif
      else
        Write(LF,Fmt2//'A,T40,I11)')'Number of CSFs',
     &                           NCSASM(STSYM)
        if (I_ELIMINATE_GAS_MOLCAS > 0) Then
          Write(LF,Fmt2//'A,T40,I11)')'Number of highly excited CSFs',
     &                           nCSF_HEXS
        EndIf
        If (DoComb) Then
          Write(LF,Fmt2//'A,T40,I11)')'Number of spin combinations',
     &                                NDTASM(STSYM)
          Write(LF,Fmt2//'A,T40,I11)')'Number of determinants',
     &                         2*NDTASM(STSYM)-NDTFTP(1)*NCNFTP(1,STSYM)
        Else
          Write(LF,Fmt2//'A,T40,I11)')'Number of determinants',
     &                                NDTASM(STSYM)
        EndIf
      end if
        n_Det=2
        n_unpaired_elec=(iSpin-1)
        n_paired_elec=nActEl-n_unpaired_elec
        If(n_unpaired_elec+n_paired_elec/2.eq.nac.or.
     &     NDTASM(STSYM).eq.1) n_Det = 1
        If(KSDFT.eq.'DIFF')   n_Det = 1
        If(KSDFT.eq.'ROKS')   n_Det = 1

      if(.not.DoSplitCAS) then  ! GLMJ
        Write(LF,Fmt2//'A,T45,I6)')'Number of root(s) required',
     &                             NROOTS
        If (irlxroot.ne.0)
     &  Write(LF,Fmt2//'A,T45,I6)')'Root chosen for geometry opt.',
     &                             IRLXROOT
        If ( ICICH.eq.0 ) then

          If ( nRoots.eq.1 ) then

           if(doDMRG)then
            Write(LF,Fmt2//'A,(T45,10I6))')'DMRG root used',
     &                                  IROOT(1)
           else
            Write(LF,Fmt2//'A,(T45,10I6))')'CI root used',
     &                                  IROOT(1)
           end if
           Write(LF,Fmt1)"   "

          else

           if(doDMRG)then
            Write(LF,Fmt2//'A,(T45,10I6))')'DMRG roots used',
     &                                    (IROOT(i),i=1,nRoots)
           else
            Write(LF,Fmt2//'A,(T45,10I6))')'CI roots used',
     &                                    (IROOT(i),i=1,nRoots)
           end if
              Write(LF,Fmt2//'A,(T45,10F6.3))')'weights',
     &                                    (Weight(i),i=1,nRoots)
          end if

        else
          Do i=1,nRoots
            Write(LF,Fmt2//'A,T45,I6)')'selected root',iRoot(i)
            Write(LF,Fmt2//'A,T45,10I6)')'Reference configurations',
     &                                (iCI(i,iRef),iRef=1,mxRef)
            Write(LF,Fmt2//'A,T45,10F6.3)')'CI-coeff',
     &                                (cCI(i,iRef),iRef=1,mxRef)
          End Do
        End If
        if(.not.doDMRG)then
          Write(LF,Fmt2//'A,T45,I6)')'highest root included in the CI',
     &                           LROOTS
          Write(LF,Fmt2//'A,T45,I6)')'max. size of the explicit '//
     &                          'Hamiltonian',NSEL
        end if
      else
        write(LF,Fmt2//'A,T45,I6)')  'Root required ', lrootSplit
        if (NumSplit)
     &  write(LF,Fmt2//'A,T45,I6)')'A-Block Size in '//
     &                            'SplitCAS',iDimBlockA
        if (EnerSplit)
     &  write(LF,Fmt2//'A,T44,F7.2)')'Energy Gap (eV) in SplitCAS',
     &                              GapSpli
        if (PerSplit)
     &  write(LF,Fmt2//'A,T44,F7.1)')'Percentage sought '//
     &                            'in SplitCAS',PercSpli
        Write(LF,Fmt2//'A,T45,ES10.3)')'Threshold for SplitCAS',
     &                              ThrSplit
*       write(LF,Fmt2//'A,T49,G10.3)')'Thrs over the root to be '//
*     &                        'opt in SplitCAS', ThrSplit
        write(LF,Fmt2//'A,T47,I4)') 'Maximum number of SplitCAS '//
     &                       'iterations', MxIterSplit
        if (FordSplit) then
          write(LF,Fmt2//'A,T47)') 'CI coeff: 1st-order approximation'
        else
          write(LF,Fmt2//'A,T47)')'CI coeff: 0th-order approximation'
        end if
      end if
      Call CollapseOutput(0,'CI expansion specifications:')

#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
 114  Continue
#endif

      END IF
      IF (lOPTO) THEN
        Write(LF,*)
        Line=' '
        Write(Line(left-2:),'(A)') 'RASSCF input specifications:'
        Call CollapseOutput(1,Line)
        Write(LF,Fmt1)'----------------------------'
        if(.not.DoSplitCAS) then
          Write(LF,Fmt2//'A,T45,I6)')'Number of root(s) required',
     &                             NROOTS
          If (irlxroot.ne.0)
     &      Write(LF,Fmt2//'A,T45,I6)')'Root chosen for geometry opt.',
     &                                 IRLXROOT
          If ( ICICH.eq.0 ) then
            If ( nRoots.eq.1 ) then
              Write(LF,Fmt2//'A,(T45,10I6))')'CI root used',
     &                                    IROOT(1)
            Else
              Write(LF,Fmt2//'A,(T45,10I6))')'CI roots used',
     &                                    (IROOT(i),i=1,nRoots)
              Write(LF,Fmt2//'A,(T45,10F6.3))')'weights',
     &                                     (Weight(i),i=1,nRoots)
            End If
          Else
           Do i=1,nRoots
              Write(LF,Fmt2//'A,T45,I6)')'selected root',iRoot(i)
              Write(LF,Fmt2//'A,T45,10I6)')'Reference configurations',
     &                                  (iCI(i,iRef),iRef=1,mxRef)
              Write(LF,Fmt2//'A,T45,10F6.3)')'CI-coeff',
     &                                  (cCI(i,iRef),iRef=1,mxRef)
            End Do
          End If
        end if
        Call CollapseOutput(0,'RASSCF input specifications:')
      ENDIF
* Check that the user doesn't try to calculate more roots than it's possible
* NN.14 FIXME: in DMRG-CASSCF, skip this check for the time
*              since Block DMRG code will check this internally
*     If (NROOTS .GT. NCSASM(STSYM)) Then
      If (.not. any([DoNECI, Do_CC_CI, DumpOnly, doDMRG, doBlockDMRG])
     &    .and. NROOTS > NCSASM(STSYM)) Then
         Write(LF,*) '************ ERROR ***********'
         Write(LF,*) ' You can''t ask for more roots'
         Write(LF,*) ' than there are configurations '
         Write(LF,*) '******************************'
         Write(LF,*)
         Call Quit_OnUserError()
      End If
* If the calculation will be too big:
      Call mma_MaxDBLE(MaxRem)
      WillNeedMB=(8.0D0*1.50D0*6.0D0*NDTASM(STSYM)/1.048D6)
      AvailMB=(8.0D0*MaxRem/1.048D6)
      if (WillNeedMB .gt. AvailMB) then
        write(6,*)
        write(6,*)' *************************************************'
        write(6,*)' Sorry, but your calculation will probably be too'
        write(6,*)' large for the available memory.'
        write(6,*)' The number of determinants is ',NDTASM(STSYM)
        write(6,*)' During CI equation solution, there will be'
        write(6,*)' up to six vectors of this size in memory.'
        write(6,*)' We estimate an additional 50% for other stuff.'
        write(6,*)
        write(6,'(A,F9.1,A)')' This alone will need at least ',
     &                                            WillNeedMB,' MB.'
        write(6,'(A,F9.1,A)')' Available at this point is ',
     &                                               AvailMB,' MB.'
        write(6,*)' Please increase MOLCAS_MEM, and try again.'
        write(6,*)' *************************************************'
        write(6,*)
        Call Quit_OnUserError()
      end if

      IF(IPRLEV.GE.USUAL .AND..NOT.lOPTO) THEN
       Write(LF,*)
       Line=' '
       Write(Line(left-2:),'(A)') 'Optimization specifications:'
       Call CollapseOutput(1,Line)
       Write(LF,Fmt1)'----------------------------'
       Write(LF,*)
       call DecideOnCholesky(DoCholesky)
       If (DoCholesky) Then
        Call Get_iScalar('System BitSwitch',iDoRI)
        if (Iand(iDoRI,1024).Eq.1024) then
           if (DoLocK) then
             Write(LF,Fmt2//'A,T45,I6)')'RASSCF algorithm: LK RI/DF'
           else
             Write(LF,Fmt2//'A,T45,I6)')'RASSCF algorithm: RI/DF'
           endif
        else
           if (DoLocK) then
             Write(LF,Fmt2//'A,T45,I6)')'RASSCF algorithm: LK Cholesky'
           else
             Write(LF,Fmt2//'A,T45,I6)')'RASSCF algorithm: Cholesky'
           endif
        endif
       Else
        Write(LF,Fmt2//'A,T45,I6)')'RASSCF algorithm: Conventional'
       EndIf
************************************************************************
* Some printout for mcpdft method
************************************************************************
       KSDFT2 = KSDFT
       IF(l_casdft) then
          KSDFT2 = KSDFT(index(KSDFT,'T:')+2:)
          Write(LF,Fmt2//'A)') 'This is a MC-PDFT calculation '//
     &                         'with functional: '//KSDFT
          Write(LF,Fmt2//'A,T45,ES10.3)')'Exchange scaling factor',CoefX
          Write(LF,Fmt2//'A,T45,ES10.3)')'Correlation scaling factor',
     &                                 CoefR
       end if
************************************************************************

       Write(LF,Fmt2//'A,T45,I6)')'Maximum number of macro iterations',
     &                           MAXIT
       Write(LF,Fmt2//'A,T45,I6)')'Maximum number of SX iterations',
     &                           ITMAX
       Write(LF,Fmt2//'A,T45,ES10.3)')'Threshold for RASSCF energy',
     &                              THRE
       Call Put_dScalar('EThr',ThrE)
       Write(LF,Fmt2//'A,T45,ES10.3)')'Threshold for max MO rotation',
     &                              THRTE
       Write(LF,Fmt2//'A,T45,ES10.3)')'Threshold for max BLB element',
     &                              THRSX
       Write(LF,Fmt2//'A,T45,ES10.3)')'Level shift parameter',
     &                              LVSHFT
       If ( NQUNE.ne.0 ) THEN
        Write(LF,Fmt1)'Make Quasi-Newton update'
       End If
       If ( ISUPSM.ne.0 ) then
         Write(LF,Fmt1)
     &   'Supersymmetry is used to disable selected orbital rotations'
         iEnd=0
         Do iSym=1,nSym
            iStart=iEnd+1
            iEnd=iEnd+nBas(iSym)
            iTemp=0
            Do i=iStart,iEnd
               iTemp=iTemp+IXSYM(i)
            End Do
            If ( iTemp.gt.0 ) then
               Write(LF,Fmt2//'A,I3)')
     &         'Supersymmetry vector for symmetry species',iSym
               Write(LF,Fmt2//'30I3)') (IXSYM(i),i=iStart,iEnd)
            End If
         End Do
       End If
       If ( lRF ) then
         Call mma_allocate(Tmp0,nTot1+4,Label='Tmp0')
         iRc=-1
         iOpt=ibset(0,sNoOri)
         iComp=1
         iSyLbl=1
         Label='Mltpl  0'
         Call RdOne(iRc,iOpt,Label,iComp,Tmp0,iSyLbl)
         Tot_Nuc_Charge=Tmp0(nTot1+4)
         If ( iRc.ne.0 ) then
            Write(LF,*) 'InpPri: iRc from Call RdOne not 0'
            Write(LF,*) 'Label = ',Label
            Write(LF,*) 'iRc = ',iRc
            Call Abend
         Endif
         Call mma_deallocate(Tmp0)
         Tot_El_Charge=Zero
         Do iSym=1,nSym
            Tot_El_Charge=Tot_El_Charge
     &                   -2.0D0*DBLE(nFro(iSym)+nIsh(iSym))
         End Do
         Tot_El_Charge=Tot_El_Charge-DBLE(nActEl)
         Tot_Charge=Tot_Nuc_Charge+Tot_El_Charge
         iCharge=Int(Tot_Charge)
         Call PrRF(.False.,NonEq,iCharge,2)
         Write(LF,Fmt2//'A,T45,I2)')' Reaction field from state:',
     &                              IPCMROOT
       End If
       Call CollapseOutput(0,'Optimization specifications:')
       If ( RFpert ) then
         Write(LF,*)
         Write(LF,Fmt1)'Reaction field specifications:'
         Write(LF,Fmt1)'------------------------------'
         Write(LF,*)
         Write(LF,'(6X,A)')'The Reaction field is added as a '//
     &                    'perturbation and has been determined '//
     &                    'in a previous calculation'
         Write(LF,*)
       End If
       If (ICIRST.EQ.1) Then
        Write(LF,Fmt1)
     &  'Starting CI array(s) will be read from file'
       End If
       Call Put_dScalar('EThr',ThrE)
*
*---- Print out grid information in case of DFT
*
       If (KSDFT.ne.'SCF') Then
         Call Put_dScalar('DFT exch coeff',CoefX)
         Call Put_dScalar('DFT corr coeff',CoefR)
         Call Funi_Print()
         IF(IPRLEV.GE.USUAL) THEN
            Write(6,*)
            Write(6,'(6X,A)') 'DFT functional specifications'
            Write(6,'(6X,A)') '-----------------------------'
            Call libxc_version()
            Call Init_Funcs(KSDFT2)
            Call Print_Info()
            Write(6,*)
         END IF
       End If
      END IF
      Write(LF,*)

  900 CONTINUE
      Call XFlush(LF)
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
      End Subroutine InpPri
