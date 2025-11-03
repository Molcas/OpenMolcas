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
* Copyright (C) 2021, Yoshio Nishimoto                                 *
************************************************************************
      Subroutine GrdIni
C
      use caspt2_global, only: LuPT2,LuGAMMA,LuCMOPT2,LuAPT2,
     *                           do_nac,do_lindep,LUGRAD,LUSTD,iStpGrd,
     *                           idBoriMat,TraFro,
     *                           CLag,CLagFull,OLag,OLagFull,SLag,WLag,
     *                           nCLag,nOLag,nSLag,nWLag,
     *                           DPT2_tot,DPT2C_tot,DPT2_AO_tot,
     *                           DPT2C_AO_tot,DPT2Canti_tot,
     *                           FIMO_all,FIFA_all,FIFASA_all,idSDMat,
     *                           OMGDER,iTasks_grad
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp
C
C     use gugx, only: CIS
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "caspt2.fh"
#include "pt2_guga.fh"
C
      character(len=128) :: FileName
      character(len=4096) :: RealName
      Logical is_error,Exists

      real(kind=wp),allocatable :: WRK(:)

      iStpGrd = 1

      ! Define (initial) logical unit numbers for gradients files
      LUPT2    = 17 ! MCLR
      LUGAMMA  = 65 ! ERI derivatives ALASKA or 3-center RI/CD
      LUCMOPT2 = 66 ! Back-transform ALASKA
      LUSTD    = 67 ! S and T derivatives in CASPT2
      LUAPT2   = 68 ! A_PT2, 2-center derivatives for RI/CD
      LUGRAD   = 69 ! CASPT2 gradient and property

      Call PrgmTranslate('GAMMA',RealName,lRealName)
      LuGAMMA  = isFreeUnit(LuGAMMA)
      If (IfChol) Then
        LENGTH = nBas(1)
      Else
        LENGTH = nIsh(1) + nAsh(1)
      End If
      Call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),
     &                     'DIRECT','UNFORMATTED',
     &                      iost,.TRUE.,
     &                      LENGTH**2*8,'REPLACE',is_error)
      Close (LuGAMMA)

      If (.not.IfChol) Then
        Call PrgmTranslate('CMOPT2',RealName,lRealName)
        LuCMOPT2 = isFreeUnit(LuCMOPT2)
        Call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),'DIRECT',
     &                       'UNFORMATTED',iost,.FALSE.,1,'REPLACE',
     &                        is_error)
        Close (LuCMOPT2)
      End If

      CALL DANAME_MF_wa(LUSTD,'LUSTD')
      If (IfChol) CALL DANAME_MF_wa(LUAPT2,'LUAPT2')

      !! Check if this is not the first (MS-)CASPT2 call
      !! If this is the first call, compute CASPT2 energies;
      !! otherwise read many things from the PT2GRD file.
      !! PT2GRD file is always deleted when a RASSCF is finished
      !! This is written in Driver/rasscf.prgm.src, but I'm not sure
      !! this is the right way to manipulate files...
      FileName = 'PT2GRD'
      Call f_inquire(FileName,Exists)
      If (Exists) iStpGrd = 0
      ! Not sure none/mf/mf_wa
      CALL DANAME_MF_wa(LUGRAD,'LUPT2GRD')
C
      !! Allocate lagrangian terms
      ! CLag and SLag should allocate for nRoots and not nState,
      ! but for the time being we only support the case nState=nRoots
      nCLag = nconf*nState
      nOLag = NBSQT
      nSLag = nState*nState
      nWLag = NBTRI
C
      call mma_allocate(DPT2_tot    ,NBSQT,Label='DPT2_tot')
      call mma_allocate(DPT2C_tot   ,NBSQT,Label='DPT2C_tot')
      call mma_allocate(DPT2_AO_tot ,NBSQT,Label='DPT2_AO_tot')
      call mma_allocate(DPT2C_AO_tot,NBSQT,Label='DPT2C_AO_tot')
      DPT2_tot     = 0.0d+00
      DPT2C_tot    = 0.0d+00
      DPT2_AO_tot  = 0.0d+00
      DPT2C_AO_tot = 0.0d+00
C
      !! Some Lagrangians for each state are constructed in CLag or
      !! OLag. The full (sum over all states, in particular for
      !! MS-CASPT2) configuration and orbital Lagrangians are then
      !! constructed in OLagFull. For SS-CASPT2, CLag and
      !! CLagFull, for instance, will be identical.
      call mma_allocate(CLag    ,nconf,nState,Label='CLAG')
      call mma_allocate(CLagFull,nconf,nState,Label='CLAGFULL')
      call mma_allocate(OLag    ,nOLag,Label='OLAG')
      call mma_allocate(OLagFull,nOLag,Label='OLAGFULL')
      call mma_allocate(SLag    ,nState,nState,Label='SLAG')
      call mma_allocate(WLag    ,nWLag,Label='WLAG')
      CLag     = 0.0d+00
      CLagFull = 0.0d+00
      OLag     = 0.0d+00
      OLagFull = 0.0d+00
      SLag     = 0.0d+00
      WLag     = 0.0d+00
C     write(6,*) "nclag,nolag,nslag"
C     write(6,*)  nclag, nolag, nslag
C
      call mma_allocate(FIMO_all,NBSQT,Label='FIMO_all')
      call mma_allocate(FIFA_all,NBSQT,Label='FIFA_all')
      FIMO_all = 0.0d+00
      FIFA_all = 0.0d+00
C
      !! FIFASA is constructed with state-averaged density always
      !! FIFA   can be state-specific or dynamically weighted
      !! FIMO   is uniquely determined, but the basis can be
      !!        either natural or quasi-canonical
      If (IFXMS .or. IFRMS) Then
        call mma_allocate(FIFASA_all,NBSQT,Label='FIFASA_all')
        FIFASA_all = 0.0d+00
        ! norbi=norb(1)
      End If
C
      If (IFDW .and. zeta >= 0.0d0) Then
        call mma_allocate(OMGDER,nState,nState,Label='OMGDER')
        OMGDER = 0.0d+00
      End If
C
      If (do_nac) Then
        call mma_allocate(DPT2Canti_tot,NBSQT,Label='DPT2Canti_tot')
        DPT2Canti_tot = 0.0d+00
      End If
C
      MaxLen = 0
      Do iCase = 1, 11
        Do iSym = 1, nSym
          nAS = nASUP(iSym,iCase)
          MaxLen = Max(MaxLen,nAS*nAS)
        End Do
      End Do
C
      call mma_allocate(WRK,MaxLen,Label='WRK')
      WRK(:) = 0.0d+00
C
      idSD = 1
C     write (*,*) "iflindep = ", iflindeplag
      If (do_lindep) Then
        Do iCase = 1, 11
          DO iSym = 1, nSym
            idBoriMat(iSym,iCase) = idSD
            NAS=NASUP(ISYM,ICASE)
            NS=(NAS*(NAS+1))/2
            CALL DDAFILE(LuSTD,0,WRK,NS,idSD)
            idSD_ = idBoriMat(iSym,iCase)
            CALL DDAFILE(LuSTD,1,WRK,NS,idSD_)
          End Do
        End Do
      End If
C
      If (MAXIT /= 0) Then
        Do iCase = 1, 11
          Do iSym = 1, nSym
            idSDMat(iSym,iCase) = idSD
            nAS = nASUP(iSym,iCase)
            CALL DDAFILE(LuSTD,0,WRK,nAS*nAS,idSD)
            idSDer = idSDMat(iSym,iCase)
            ! idSDMat(iSym,iCase))
            CALL DDAFILE(LuSTD,1,WRK,nAS*nAS,idSDer)
          End Do
        End Do
      End If
      call mma_deallocate(WRK)
C
      if (nFroT /= 0) call mma_allocate(TraFro,nFroT**2,Label='TraFro')
      call mma_allocate(iTasks_grad,nAshT**2,Label='Tasks_grad')
      iTasks_grad(:) = 0
C
      Return

      End Subroutine GrdIni

C-----------------------------------------------------------------------

      Subroutine GrdCls(IRETURN,UEFF,U0,H0)
C
      use caspt2_global, only: iPrGlb
      use caspt2_global, only: LuPT2,LuAPT2,
     *                           do_nac,do_csf,iRoot1,iRoot2,LUGRAD,
     *                           LUSTD,TraFro,
     *                           CLag,CLagFull,OLag,OLagFull,SLag,WLag,
     *                           nOLag,nWLag,
     *                           DPT2_tot,DPT2C_tot,DPT2_AO_tot,
     *                           DPT2C_AO_tot,DPT2Canti_tot,
     *                           FIMO_all,FIFA_all,FIFASA_all,OMGDER,
     *                           iTasks_grad
      use PrintLevel, only: verbose
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "caspt2.fh"
C
      Dimension UEFF(nState,nState),U0(nState,nState),H0(nState,nState)
      Character(Len=16) mstate1
      LOGICAL DEB,Found
      logical, external :: RF_On
C
      real(kind=wp),allocatable :: HEFF1(:,:),WRK1(:,:),WRK2(:,:),
     *                             CI1(:,:)
C
      !! In case convergence of CASPT2 equation failed
      !! Call this subroutine just deallocate memory
      If (IRETURN.NE.0) GO TO 9000
C
      call mma_allocate(HEFF1,nState,nState,Label='HEFF1')
      call mma_allocate(WRK1,nState,nState,Label='WRK1')
      call mma_allocate(WRK2,nState,nState,Label='WRK2')
C
      !! Add XMS specific terms
      !! Note that CLagFull is in natural CSF basis,
      !! so everything in this subroutine has to be done in natural
      SLag(:,:) = 0.0d+00
      WRK2(:,:) = 0.0d+00
      If (IFDW .and. zeta >= 0.0d0) Then
        !! Construct Heff[1] in XMS basis
        HEFF1(:,:) = 0.0d+00
        Do ilStat = 1, nState
         HEFF1(ilStat,ilStat) = REFENE(ilStat)
        End Do
        Call DGEMM_('T','N',nState,nState,nState,
     *              1.0D+00,U0,nState,HEFF1,nState,
     *              0.0D+00,WRK1,nState)
        Call DGEMM_('N','N',nState,nState,nState,
     *              1.0D+00,WRK1,nState,U0,nState,
     *              0.0D+00,HEFF1,nState)
C
        !! Derivative of Heff[1] in XMS basis
        !! It is transformed with U0, so the contribution has to be
        !! considered when we construct the auxiliary density in the
        !! XMS-specific term
        call DWDER(OMGDER,HEFF1,SLag)
        Call DGEMM_('N','N',nState,nState,nState,
     *              1.0D+00,U0,nState,SLag,nState,
     *              0.0D+00,WRK2,nState)
        Call DGEMM_('N','T',nState,nState,nState,
     *              1.0D+00,WRK2,nState,U0,nState,
     *              0.0D+00,WRK1,nState)

        WRK2(:,:) = 0.0D+00
        Do ilStat = 1, nState
          If (DWTYPE.EQ.1) Then
            WRK2(ilStat,ilStat) = SLag(ilStat,ilStat)
          Else If (DWTYPE.EQ.2.OR.DWTYPE.EQ.3) Then
            Do jlStat = 1, nState
              WRK2(ilStat,jlStat) = SLag(ilStat,jlStat)
            End Do
          End If
          If (.not.do_nac) Then
            Do jlStat = 1, ilStat-1
             WRK1(ilStat,jlStat)=WRK1(ilStat,jlStat)+WRK1(jlStat,ilStat)
              WRK1(jlStat,ilStat) = 0.0d+00
            End Do
          End If
        End Do
        SLag(:,:) = WRK1(:,:)
      End If

      IF (IFXMS.or.IFRMS.or.(IFMSCOUP.and.do_nac.and.do_csf)) Then

        If (.not.IFXMS .and. .not.IFRMS) Then
          !! For MS-CASPT2, only the second term in eq.(68)
          U0(:,:) = 0.0D+00
          Call DCopy_(nState,[1.0D+00],0,U0,nState+1)
        End If

        CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
        CALL XMS_Grad(H0,U0,UEFF,WRK2)
        CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
        CPUT =CPTF10-CPTF0
        WALLT=TIOTF10-TIOTF0
        If (IPRGLB.ge.VERBOSE) Then
          write(6,'(a,2f10.2)')" XMS_Grad: CPU/WALL TIME=", cput,wallt
        End If
      End If

      !! Now, compute the state Lagrangian and do some projections
      !! If PCM, we have to obtain internal state rotation parameters
      !! self-consistently, so we should not call this subroutine
      If (.not.RF_On()) Call CLagFinal(CLagFull,SLag)

      !! Add MS-CASPT2 contributions
      If (IFMSCOUP) Then
        Do ilStat = 1, nState
          Do jlStat = 1, ilStat
            If (do_nac) Then
              If (.not.IFXMS .and. .not.IFRMS .and. ilstat.ne.jlstat)
     &          Cycle

              Scal = UEFF(ilStat,iRoot1)*UEFF(jlStat,iRoot2)
     *             + UEFF(jlStat,iRoot1)*UEFF(ilStat,iRoot2)
              Scal = Scal*0.5D+00
              SLag(ilStat,jlStat) = SLag(ilStat,jlStat) + Scal
C
              If (ilStat.ne.jlStat) Then
                SLag(jlStat,ilStat) = SLag(jlStat,ilStat) + Scal
              End If
            Else
              IF (IFXMS .or. IFRMS) Then
                Scal = UEFF(ilStat,iRoot1)*UEFF(jlStat,iRoot2)
                If (ilStat.ne.jlStat) Scal = Scal*2.0d+00
                SLag(ilStat,jlStat) = SLag(ilStat,jlStat) + Scal
              Else
                If (ilStat.eq.jlStat) Then
                  Scal = UEFF(ilStat,iRoot1)*UEFF(jlStat,iRoot2)
                  SLag(ilStat,jlStat) = SLag(ilStat,jlStat) + Scal
                End If
              End If
            End If
          End Do
        End Do
      End If
C
      !! Subtract the original rhs_sa.f or rhs_nac.f contribution
      !! For MS-type CASPT2, CASSCF part has to be determined by UEFF
      If (IFMSCOUP.and.iRoot1.eq.iRoot2) Then
        ilStat = MAX(iRoot1,iRoot2)
        jlStat = MIN(iRoot1,iRoot2)
        SLag(ilStat,jlStat) = SLag(ilStat,jlStat) - 1.0D+00
      End If
C
      !! Finalize the first-order transition(-like) density matrix
      !! for the CSF derivative term
      If (do_nac) Then
        If (do_csf) Then
          Call CnstAntiC(DPT2Canti_tot,UEFF,U0)
        Else
          !! Clear just in case
          DPT2Canti_tot = 0.0d+00
        End If
      End If
C
      !! Back-transform the CI Lagrangian
      !! It is in the XMS basis, so it has to be transformed to
      !! CASSCF basis to be used in Z-vector
      !! No need to do this for SLag.
      If (IFXMS .or. IFRMS) Then
        call mma_allocate(CI1,nConf,nState,Label='CI1')
        Call DGEMM_('N','T',nConf,nState,nState,
     &              1.0D+00,CLagFull,nConf,U0,nState,
     &              0.0D+00,CI1,nConf)
        CLagFull(:,:) = CI1(:,:)
        call mma_deallocate(CI1)
      End If
C
      !! Compute true unrelaxed properties for MS-CASPT2
      if ((.not.do_nac) .and. ifmscoup) CALL PRPCTL(1,UEFF,U0)
C
      LuPT2 = isFreeUnit(LuPT2)
      Call Molcas_Open(LuPT2,'PT2_Lag')

      DEB = .false.
      !! configuration Lagrangian (read in RHS_PT2)
      If (DEB) call RecPrt('CLagFull','',CLagFull,nConf,nState)
      Do j = 1, nState
        Do i = 1, nConf
          Write (LuPT2,*) CLagFull(i,j)
        End Do
      End do

      !! orbital Lagrangian (read in RHS_PT2)
#ifdef _MOLCAS_MPP_
      if (is_real_par()) then
        If (.not.King()) OLagFull(:) = 0.0d+00
        CALL GADSUM (OLagFull,nOLag)
      end if
#endif
      If (DEB) call RecPrt('OLagFull','',OLagFull,nBasT,nBasT)
      Do i = 1, nOLag
        Write (LuPT2,*) OLagFull(i)
      End Do

      !! state Lagrangian (read in RHS_PT2)
      If (DEB) call RecPrt('SLag', '', SLag, nState, nState)
      Do j = 1, nState
        Do i = 1, nState
          Write (LuPT2,*) SLag(i,j)
        End Do
      End Do

      !! renormalization contributions (read in OUT_PT2)
#ifdef _MOLCAS_MPP_
      if (is_real_par()) then
        If (.not.King()) WLag(:) = 0.0d+00
        CALL GADSUM (WLag,nWLag)
      end if
#endif
      If (DEB) call TriPrt('WLag', '', WLag, nBast)
      Do i = 1, nWLag ! = NBTRI
        Write (LuPT2,*) WLag(i)
      End Do

      !! D^PT2 in MO (read in OUT_PT2)
#ifdef _MOLCAS_MPP_
      if (is_real_par()) then
        If (.not.King()) DPT2_tot(:) = 0.0d+00
        CALL GADSUM (DPT2_tot,NBSQT)
      end if
#endif
      If (DEB) call RecPrt('DPT2', '', DPT2_tot, nBast, nBast)
      Do i = 1, NBSQT
        Write (LuPT2,*) DPT2_tot(i)
      End Do

      !! D^PT2(C) in MO (read in OUT_PT2)
#ifdef _MOLCAS_MPP_
      if (is_real_par()) then
        If (.not.King()) DPT2C_tot(:) = 0.0D+00
        CALL GADSUM (DPT2C_tot,NBSQT)
      end if
#endif
      If (DEB) call RecPrt('DPT2C', '', DPT2C_tot, nBast, nBast)
      Do i = 1, NBSQT
        Write (LuPT2,*) DPT2C_tot(i)
      End Do

      !! NAC
      If (do_nac) Then
        Do i = 1, NBSQT
          Write (LuPT2,*) DPT2Canti_tot(i)
        End Do
      End If

      !! D^PT2 in AO (not used?)
      if (RFpert.and.IFMSCOUP) then
        !! Recompute DPT2AO and DPT2AO for PCM
        Call Recompute_DPT2AO(DPT2_tot,DPT2C_tot,
     *                        DPT2_AO_tot,DPT2C_AO_tot)
      end if
      If (DEB) call TriPrt('DPT2_AO_tot', '', DPT2_AO_tot, nBast)
      Do i = 1, NBSQT
        Write (LuPT2,*) DPT2_AO_tot(i)
      End Do

      !! D^PT2(C) in AO (not used?)
      If (DEB) call TriPrt('DPT2C_AO', '', DPT2C_AO_tot, nBast)
      Do i = 1, NBSQT
        Write (LuPT2,*) DPT2C_AO_tot(i)
      End Do

      if (RFpert) then
        !! For CASPT2/PCM gradient
        !! The implicit derivative contributions have not been
        !! considered in the CASPT2 module
        if (ifmscoup) then
          ! I do not remember why this should be halved
          call daxpy_(NBTRI,0.5D+00,DPT2C_AO_tot,1,DPT2_AO_tot,1)
        else
          call daxpy_(NBTRI,1.0D+00,DPT2C_AO_tot,1,DPT2_AO_tot,1)
        end if
        Call Put_dArray('D1aoVar',DPT2_AO_tot,NBTRI)
      else
        !! not sure this is OK
        Call Put_dArray('D1aoVar',DPT2_AO_tot,0)
      end if

      ! close gradient files
      Close (LuPT2)
C
      call mma_deallocate(HEFF1)
      call mma_deallocate(WRK1)
      call mma_deallocate(WRK2)
C
 9000 CONTINUE
C
      call mma_deallocate(DPT2_tot)
      call mma_deallocate(DPT2C_tot)
      call mma_deallocate(DPT2_AO_tot)
      call mma_deallocate(DPT2C_AO_tot)
C
      call mma_deallocate(CLag)
      call mma_deallocate(CLagFull)
      call mma_deallocate(OLag)
      call mma_deallocate(OLagFull)
      call mma_deallocate(SLag)
      call mma_deallocate(WLag)
C
      call mma_deallocate(FIMO_all)
      call mma_deallocate(FIFA_all)
C
      If (IFXMS .or. IFRMS)         call mma_deallocate(FIFASA_all)
      If (IFDW .and. zeta >= 0.0d0) call mma_deallocate(OMGDER)
      If (do_nac)                   call mma_deallocate(DPT2Canti_tot)
C
      !! Prepare for MCLR
      iGo = 3
      Call Put_iScalar('SA ready',iGo)
      ! overwrites whatever was set in CASSCF with the relax
      ! root that was chosen in CASPT2
      if (do_nac) then
C       write (*,*) "NAC"
C       write (*,*) "CASSCF/Original = ", iRoot1,iRoot2
        Call Put_iScalar('Relax CASSCF root',iRoot1)
        Call Put_iScalar('Relax Original root',iRoot2)
        call Qpg_cArray('MCLR Root',Found,I)
        if (Found) then
          Call Get_cArray('MCLR Root',mstate1,16)
          if ((mstate1(8:8).eq.'0') .and. (mstate1(16:16).eq.'0')) then
            !! NAC states have not been specified
            write (mstate1,'(1X,I7,1X,I7)') iRoot1,iRoot2
            Call Put_cArray('MCLR Root',mstate1,16)
          end if
        else
          mstate1 = '****************'
          Call Put_cArray('MCLR Root',mstate1,16)
        end if
      else
C       write (*,*) "GRD"
C       write (*,*) "CASSCF/Original = ", irlxroot,irlxroot
        Call Put_iScalar('Relax CASSCF root',irlxroot)
        Call Put_iScalar('Relax Original root',irlxroot)
        mstate1 = '****************'
        Call Put_cArray('MCLR Root',mstate1,16)
      end if
C
      !! Close files
      Call DaClos(LUSTD)
      If (IfChol) Call DaClos(LUAPT2)
      Call DaClos(LUGRAD)
C
      if (nFroT /= 0) call mma_deallocate(TraFro)
      call mma_deallocate(iTasks_grad)
C
      Return
C
      End Subroutine GrdCls
C
C-----------------------------------------------------------------------
C
      Subroutine ModDip
C
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "caspt2.fh"
C
      real(kind=wp),allocatable :: DMs1(:,:),DMs2(:,:)
C
      call mma_allocate(DMs1,3,nRoots,Label='DMs1')
      call mma_allocate(DMs2,3,lRoots,Label='DMs2')
      Call Get_dArray('Last Dipole Moments',DMs2,3*LROOTS)
      Do i = 1, lRoots
        j = Root2State(i)
        If (j.eq.0) Cycle
        DMs1(:,j) = DMs2(:,i)
      End Do
      Call Put_dArray('Last Dipole Moments',DMs1,3*nROOTS)
      call mma_deallocate(DMs1)
      call mma_deallocate(DMs2)
C
      Return
C
      End Subroutine ModDip
C
C-----------------------------------------------------------------------
C
      Subroutine GradStart
C
      use caspt2_global, only:iPrGlb
      use caspt2_global, only:ipea_shift
      use caspt2_global, only: if_invar
      use PrintLevel, only: usual
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "caspt2.fh"
C
      If ((.not.if_invar) .and. (IPRGLB >= USUAL)) Then
        Write (6,*)
        Write (6,'(3X,"This is a non-invariant CASPT2 calculation")')
        If (ipea_shift /= 0.0D+00)
     *    Write (6,'(3X,"- IPEA shift is employed")')
        Write (6,'(3X,"A linear equation will be solved to obtain ",
     *                "the off-diagonal active density")')
        Write (6,*)
      End If
C
      End Subroutine GradStart
C
C-----------------------------------------------------------------------
C
      Subroutine GradPrep(UEFF,VECROT)
C
      use caspt2_global, only: iRoot1, iRoot2, jStLag
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "caspt2.fh"
C
C#include "nadc.fh"
C#include "nac.fh"
C
      Dimension UEFF(nState,nState),VECROT(nState)
C
C     H_{IJ} = <I|H^2|J>
C     we diagonalize tH_{IJ} = tilde-H_{IJ} = (H_{IJ}+H_{JI})/2
C     U^T*H*U = U_{IK}*tH_{IJ}*U_{JL}
C             = U_{IK}*(H_{IJ}+H_{JI})*U_{JL}/2 = delta_{KL}
C     Derivative of the diagonal:
C       d(U_{IK}*(H_{IJ}+H_{JI})*U_{JK})/dx/2
C       = U_{IK}*d(H_{IJ}+H_{JI})/dx*U_{JK}/2
C       = U_{IK}*dH_{IJ}/dx*U_{JK}/2 + U_{JK}*dH_{IJ}/dx*U_{IK}/2
C       = U_{IK}*dH_{IJ}/dx*U_{JK}
C     Derivative of H for off-diagonal:
C       d(U_{IK}*(H_{IJ}+H_{JI})*U_{JL})/dx/2
C       = U_{IK}*d(H_{IJ}+H_{JI})/dx*U_{JL}/2
C       = U_{IK}*dH_{IJ}/dx*U_{JL}/2 + U_{JK}*dH_{IJ}/dx*U_{IL}/2
C       = (U_{IK}*U_{JL} + U_{JK}*U_{IL}) * dH_{IJ}/dx/2
C       also
C       d(U_{IL}*(H_{IJ}+H_{JI})*U_{JK})/dx/2
C       = U_{IL}*d(H_{IJ}+H_{JI})/dx*U_{JK}/2
C       = U_{IL}*dH_{IJ}/dx*U_{JK}/2 + U_{JL}*dH_{IJ}/dx*U_{IK}/2
C
C Hij = U1i*tH11*U1j + U1i*tH12*U2j + U2i*tH21*U1j + U2i*tH22*U2j
C     = (U1i*(H11+H11)*U1j + U1i*(H12+H21)*U2j + U2i*(H21+H12)*U1j + U2i*(H22+H22)*U2j)/2
C     = (2*U1i*H1j*U12 + (U1i*U2j+U2i*U1j)H12 + (U1i*U2j+U2i*U1j)*H21 + 2*U2i*H22*U2j)/2
C     If i  = j, UIi*dHij/dx*UJj
C     If i \= j, (UIi*UJj+UJi*UIj)*dHij/dx*0.5
C
      !! Construct the rotation vector
      If (IFMSCOUP) Then
        Do iState = 1, nState
          TMP = UEFF(iState,iRoot1)*UEFF(jState,iRoot2)
     *        + UEFF(iState,iRoot2)*UEFF(jState,iRoot1)
          VECROT(iState) = TMP*0.5d+00
        End Do
        jStLag    = jState
      Else
C       write(6,*) 'jState in gradprep: ',jstate
        VECROT(jState) = 1.0D+00
        jStLag    = jState
      End If
C
      End Subroutine GradPrep
C
C-----------------------------------------------------------------------
C
      Subroutine OLagFinal(OLagLoc,Trf)
C
      use caspt2_global, only: CMOPT2
      use caspt2_global, only: OLagFull,WLag,nOLag
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp
      Implicit Real*8 (A-H,O-Z)
C
#include "caspt2.fh"
C
      Dimension OLagLoc(*),Trf(*)
      real(kind=wp),allocatable :: WRK(:),WLagLoc(:)
C
      call mma_allocate(WRK,NBSQT,Label='WRK')
      call mma_allocate(WLagLoc,NBSQT,Label='WLagLoc')
C
      WLagLoc(1:NBSQT) = 0.5D+00*OLagLoc(1:NBSQT)
C     write(6,*) "Wlag square"
C     call sqprt(wlag,nbast)
C
      !! W(MO) -> W(AO) using the quasi-canonical orbitals
      !! No need to back transform to natural orbital basis
      Call DGemm_('N','N',nBasT,nBasT,nBasT,
     *            1.0D+00,CMOPT2,nBasT,WLagLoc,nBasT,
     *            0.0D+00,WRK,nBasT)
      Call DGemm_('N','T',nBasT,nBasT,nBasT,
     *            1.0D+00,WRK,nBasT,CMOPT2,nBasT,
     *            0.0D+00,WLagLoc,nBasT)
C
      !! square -> triangle for WLag(AO)
      WRK(:) = WLagLoc(:)
      iBasTr = 1
      iBasSq = 1
      Do iSym = 1, nSym
        nBasI = nBas(iSym)
        liBasTr = iBasTr
        liBasSq = iBasSq
        Do iBasI = 1, nBasI
          Do jBasI = 1, iBasI
            liBasSq = iBasSq + iBasI-1 + nBasI*(jBasI-1)
            If (iBasI.eq.jBasI) Then
              WLagLoc(liBasTr) = WRK(liBasSq)
            Else
            liBasSq2 = iBasSq + jBasI-1 + nBasI*(iBasI-1)
            WLagLoc(liBasTr) = WRK(liBasSq)+WRK(liBasSq2)
            End If
            liBasTr = liBasTr + 1
          End Do
        End Do
        iBasTr = iBasTr + nBasI*(nBasI+1)/2
        iBasSq = iBasSq + nBasI*nBasI
      End Do
      ! accumulate W Lagrangian only for MS,XMS,XDW,RMS,
      ! but not for SS-CASPT2
      if (jState.eq.iRlxRoot .or. IFMSCOUP) then
        WLag(1:NBTRI) = WLag(1:NBTRI) + WLagLoc(1:NBTRI)
      end if
      call mma_deallocate(WLagLoc)
C
C
C
      !! Transform quasi-canonical -> natural MO basis
      !! orbital Lagrangian
      Call DGemm_('N','N',nBasT,nBasT,nBasT,
     *            1.0D+00,Trf,nBasT,OLagLoc,nBasT,
     *            0.0D+00,WRK,nBasT)
      Call DGemm_('N','T',nBasT,nBasT,nBasT,
     *            1.0D+00,WRK,nBasT,Trf,nBasT,
     *            0.0D+00,OLagLoc,nBasT)
      !! sufficient only for active
      nBasI = nBas(1)
      WRK(1:nBasI**2) = OLagLoc(1:nBasI**2)
      Call DGeSub(WRK,nBas(1),'N',
     &            WRK,nBas(1),'T',
     &            OLagLoc,nBas(1),
     &            nBas(1),nBas(1))
      ! accumulate orbital Lagrangian only for MS,XMS,XDW,RMS,
      ! but not for SS-CASPT2
      if (jState.eq.iRlxRoot .or. IFMSCOUP) then
        Call DaXpY_(nOLag,1.0D+00,OLagLoc,1,OLagFull,1)
      end if
C
      call mma_deallocate(WRK)
C
      End Subroutine OLagFinal
C
C-----------------------------------------------------------------------

      Subroutine CnstFIFAFIMO(MODE)

      use caspt2_global, only: TraFro, OLag,
     *                           FIMO_all, FIFA_all, FIFASA_all
      use caspt2_global, only: FIMO, FIFA, CMOPT2
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp

      Implicit Real*8 (A-H,O-Z)

#include "caspt2.fh"

      real(kind=wp),allocatable :: WRK1(:),WRK2(:)

      If (IfChol) Then
        !! For DF or CD, we already have FIFA and FIMO in AO,
        !! so just do AO -> MO transformation
        call mma_allocate(WRK1,NBSQT,Label='WRK1')
        call mma_allocate(WRK2,NBSQT,Label='WRK2')
        WRK1(:) = 0.0d+00
        WRK2(:) = 0.0d+00

!           write (*,*) "fifa_all,fimo_all"
!           do i = 1, 10
!           write (*,*) i,fifa_all(i),fimo_all(i)
!           end do
        iSQ = 0
        iTR = 0
        Do iSym = 1, nSym
          ! nOrbI = nOrb(iSym)
          nBasI = nBas(iSym)
          !! FIFA
          If (nFroT.eq.0) Then
            If (MODE.eq.0 .and. (IFDW.or.IFRMS)) Then
              Call SQUARE(FIFA(1+iTr),FIFASA_all(1+iSQ),1,nBasI,nBasI)
            Else If (MODE.eq.1) Then
              Call SQUARE(FIFA(1+iTr),FIFA_all(1+iSQ),1,nBasI,nBasI)
C             write (*,*) "fifa in MO"
C             call sqprt(fifa_all(1+isq),nbasi)
            End If
          Else
            Call SQUARE(FIFA_all(1+iTr),WRK1,1,nBasI,nBasI)
C             write (*,*) "fifasa in AO"
C             call sqprt(wrk1,nbasi)
            If (MODE.eq.0 .and. (IFDW.or.IFRMS)) Then
              !! with the state-average
              !! FIFASA_all will be natural basis
              Call OLagTrf(2,iSym,CMOPT2,FIFASA_all(1+iSQ),WRK1,WRK2)
C             write (*,*) "fifasa in MO"
C             call sqprt(fifasa_all(1+isq),nbasi)
            Else If (MODE.eq.1) Then
              !! with the state-specific or dynamically weighted
              !! FIFA will be quasi-canonical basis
              Call OLagTrf(2,iSym,CMOPT2,FIFA_all(1+iSQ),WRK1,WRK2)
C             write (*,*) "fifa in MO"
C             call sqprt(fifa_all(1+isq),nbasi)
              !! canonicalize frozen orbitals
              !! still under investigation, but this is something we
              !! should do to obtain "better" orbital enegies for
              !! methods using state-dependent Fock operators.
              !! We actually need to canonicalize frozen and inactive
              !! orbitals simultaneously?
              If (nFroT /= 0) Then
                CALL DCOPY_(nBasI*nBasI,WRK1,1,WRK2,1)
                CALL DIAFCK(NBAS(ISYM),FIFA_all,1,NFRO(ISYM),
     &                      TraFro,NBAS(ISYM),CMOPT2,WRK2)
                CALL DCOPY_(NBAS(ISYM)*NFRO(ISYM),WRK2,1,CMOPT2,1)
                Call OLagTrf(2,iSym,CMOPT2,FIFA_all(1+iSQ),WRK1,WRK2)
              End If
            End If
          End If
C
          !! FIMO
C         If (MODE.eq.0) Then
            If (nFroT.eq.0) Then
              Call SQUARE(FIMO(1+iTr),FIMO_all(1+iSQ),1,nBasI,nBasI)
            Else
              Call SQUARE(FIMO_all(1+iTr),WRK1,1,nBasI,nBasI)
              Call OLagTrf(2,iSym,CMOPT2,FIMO_all(1+iSQ),WRK1,OLag)
C             write (*,*) "fimo in MO"
C             call sqprt(fimo_all(1+isq),nbasi)
            End If
C         End If
          iSQ = iSQ + nBasI*nBasI
          iTR = iTR + nBasI*(nBasI+1)/2
        End Do
        call mma_deallocate(WRK1)
        call mma_deallocate(WRK2)
C
        If (IFXMS.and..not.IFDW)
     *    Call DCopy_(NBSQT,FIFA_all,1,FIFASA_all,1)
      Else
        If (nFroT.ne.0) Then
        Else
          iSQ = 0
          iTR = 0
          Do iSym = 1, nSym
            ! nOrbI = nOrb(iSym)
            nBasI = nBas(iSym)
            Call SQUARE(FIFA(1+iTr),FIFA_all(1+iSQ),1,nBasI,nBasI)
            Call SQUARE(FIMO(1+iTr),FIMO_all(1+iSQ),1,nBasI,nBasI)
            iSQ = iSQ + nBasI*nBasI
            iTR = iTR + nBasI*(nBasI+1)/2
          End Do
          If (IFXMS.and..not.IFDW)
     *      Call DCopy_(NBSQT,FIFA_all,1,FIFASA_all,1)
        End If
C
      !! XDW or RMS case: call after XDWINI
      ! If (MODE.eq.0) Then
      ! End If
C
      !! XMS case: call after GRPINI
      ! If (MODE.eq.1) Then
      ! End If
C
C     !! SS or MS case: call in dens.f
C     If (MODE.eq.2) Then
C       If (IFSADREF) Then
C       Else
C       End If
C     End If
      End If
C
      Return
C
      End Subroutine CnstFIFAFIMO
C
C-----------------------------------------------------------------------
C
      Subroutine Recompute_DPT2AO(DPT2,DPT2C,DPT2AO,DPT2CAO)
C
      use definitions, only: wp
      use stdalloc, only: mma_allocate,mma_deallocate
      use caspt2_global, only: LUONEM
C
      Implicit Real*8 (A-H,O-Z)
C
#include "caspt2.fh"
C
      real(kind=wp), intent(in) :: DPT2(NBSQT),DPT2C(NBSQT)
      real(kind=wp), intent(inout) :: DPT2AO(NBSQT),DPT2CAO(NBSQT)
      real(kind=wp), allocatable :: WRK1(:),WRK2(:),WRK3(:),WRK4(:)
C
      call mma_allocate(WRK1,NBSQT,Label='WRK1')
      call mma_allocate(WRK2,NBSQT,Label='WRK2')
      call mma_allocate(WRK3,NBSQT,Label='WRK3')
      call mma_allocate(WRK4,NBSQT,Label='WRK4')
      IDISK=IAD1M(1)
      CALL DDAFILE(LUONEM,2,WRK1,NBSQT,IDISK)

      call dcopy_(NBSQT,[0.0d+00],0,DPT2AO,1)
      call dcopy_(NBSQT,[0.0d+00],0,DPT2CAO,1)
C
      iBasTr = 1
      iBasSq = 1
      Do iSym = 1, nSym
        call OLagTrf(1,iSym,WRK1,DPT2 ,WRK3,WRK2)
        call OLagTrf(1,iSym,WRK1,DPT2C,WRK4,WRK2)
        nBasI = nBas(iSym)
        liBasTr = iBasTr
        liBasSq = iBasSq
        ljBasSq = iBasSq
        Do iBasI = 1, nBasI
          Do jBasI = 1, iBasI
            liBasSq = iBasSq + iBasI-1 + nBasI*(jBasI-1)
            ljBasSq = iBasSq + jBasI-1 + nBasI*(iBasI-1)
            If (iBasI.eq.jBasI) Then
              DPT2AO (liBasTr) = WRK3(liBasSq)
              DPT2CAO(liBasTr) = WRK4(liBasSq)
            Else
              val = WRK3(liBasSq)+WRK3(ljBasSq)
              DPT2AO (liBasTr) = val
              val = WRK4(liBasSq)+WRK4(ljBasSq)
              DPT2CAO(liBasTr) = val
            End If
            liBasTr = liBasTr + 1
          End Do
        End Do
        iBasTr = iBasTr + nBasI*(nBasI+1)/2
      End Do
C
      call mma_deallocate(WRK1)
      call mma_deallocate(WRK2)
      call mma_deallocate(WRK3)
      call mma_deallocate(WRK4)

      End Subroutine Recompute_DPT2AO
