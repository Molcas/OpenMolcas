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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************
      Subroutine GrdCls(IRETURN,nState,UEFF,U0,H0)

      use caspt2_global, only: iPrGlb
      use caspt2_global, only: LuPT2,LuAPT2,                            &
     &                           do_nac,do_csf,iRoot1,iRoot2,LUGRAD,    &
     &                           LUSTD,TraFro,                          &
     &                           CLag,CLagFull,OLag,OLagFull,SLag,WLag, &
     &                           nOLag,nWLag,                           &
     &                           DPT2_tot,DPT2C_tot,DPT2_AO_tot,        &
     &                           DPT2C_AO_tot,DPT2Canti_tot,            &
     &                           FIMO_all,FIFA_all,FIFASA_all,OMGDER,   &
     &                           iTasks_grad
      use PrintLevel, only: VERBOSE
      use stdalloc, only: mma_allocate,mma_deallocate
      use Constants, only: Zero, One, Half, Two
      use definitions, only: wp, iwp, u6
      use caspt2_module, only: REFENE, RFPERT, IfChol, IFMSCOUP, IFXMS, &
     &                         IFRMS, IFDW, NCONF, nFroT, NBAST,        &
     &                         NBTRI, NBSQT, iRlxRoot, NROOTS, DWTYPE,  &
     &                         Zeta
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif

      implicit none

      integer(kind=iwp), intent(in) :: IRETURN, nState
      real(kind=wp), intent(in) :: H0(nState,nState)
      real(kind=wp), intent(inout) :: UEFF(nState,nState),              &
     &                                U0(nState,nState)

      Character(Len=16) :: mstate1
      logical(kind=iwp) :: DEB,Found
      logical(kind=iwp), external :: RF_On

      real(kind=wp), allocatable :: HEFF1(:,:),WRK1(:,:),WRK2(:,:),     &
     &                             CI1(:,:)
      integer(kind=iwp) :: i, iGo, ilStat, j, jlStat
      integer(kind=iwp), external :: isFreeUnit
      real(kind=wp) :: Scal
      real(kind=wp) :: CPUT, WALLT, CPE, CPTF0, CPTF10, TIOE, TIOTF0,   &
     &                 TIOTF10

      !! In case convergence of CASPT2 equation failed
      !! Call this subroutine just deallocate memory
      If (IRETURN == 0) then

        call mma_allocate(HEFF1,nState,nState,Label='HEFF1')
        call mma_allocate(WRK1,nState,nState,Label='WRK1')
        call mma_allocate(WRK2,nState,nState,Label='WRK2')

        !! Add XMS specific terms
        !! Note that CLagFull is in natural CSF basis,
        !! so everything in this subroutine has to be done in natural
        SLag(:,:) = Zero
        WRK2(:,:) = Zero
        If (IFDW .and. zeta >= Zero) Then
          !! Construct Heff[1] in XMS basis
          HEFF1(:,:) = Zero
          Do ilStat = 1, nState
           HEFF1(ilStat,ilStat) = REFENE(ilStat)
          End Do
          Call DGEMM_('T','N',nState,nState,nState,                     &
     &                One,U0,nState,HEFF1,nState,                       &
     &                Zero,WRK1,nState)
          Call DGEMM_('N','N',nState,nState,nState,                     &
     &                One,WRK1,nState,U0,nState,                        &
     &                Zero,HEFF1,nState)

          !! Derivative of Heff[1] in XMS basis
          !! It is transformed with U0, so the contribution has to be
          !! considered when we construct the auxiliary density in the
          !! XMS-specific term
          call DWDER(OMGDER,HEFF1,SLag)
          Call DGEMM_('N','N',nState,nState,nState,                     &
     &                One,U0,nState,SLag,nState,                        &
     &                Zero,WRK2,nState)
          Call DGEMM_('N','T',nState,nState,nState,                     &
     &                One,WRK2,nState,U0,nState,                        &
     &                Zero,WRK1,nState)

          WRK2(:,:) = Zero
          Do ilStat = 1, nState
            If (DWTYPE == 1) Then
              WRK2(ilStat,ilStat) = SLag(ilStat,ilStat)
            Else If (DWTYPE == 2 .OR. DWTYPE == 3) Then
              Do jlStat = 1, nState
                WRK2(ilStat,jlStat) = SLag(ilStat,jlStat)
              End Do
            End If
            If (.not.do_nac) Then
              Do jlStat = 1, ilStat-1
                WRK1(ilStat,jlStat) =                                   &
     &            WRK1(ilStat,jlStat)+WRK1(jlStat,ilStat)
                WRK1(jlStat,ilStat) = Zero
              End Do
            End If
          End Do
          SLag(:,:) = WRK1(:,:)
        End If

        IF (IFXMS .or. IFRMS .or.                                       &
     &     (IFMSCOUP .and. do_nac .and. do_csf)) Then
          If (.not.IFXMS .and. .not.IFRMS) Then
            !! For MS-CASPT2, only the second term in eq.(68)
            U0(:,:) = Zero
            Call DCopy_(nState,[One],0,U0,nState+1)
          End If

          CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
          CALL XMS_Grad(H0,U0,UEFF,WRK2)
          CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
          CPUT =CPTF10-CPTF0
          WALLT=TIOTF10-TIOTF0
          If (IPRGLB >= VERBOSE) Then
           write(u6,'(a,2f10.2)')' XMS_Grad: CPU/WALL TIME=', cput,wallt
          End If
        End If

        !! Now, compute the state Lagrangian and do some projections
        !! If PCM, we have to obtain internal state rotation parameters
        !! self-consistently, so we should not call this subroutine
        If (.not.RF_On())                                               &
     &    Call CLagFinal(nConf,nRoots,nState,CLagFull,SLag)

        !! Add MS-CASPT2 contributions
        If (IFMSCOUP) Then
          Do ilStat = 1, nState
            Do jlStat = 1, ilStat
              If (do_nac) Then
                If (.not.IFXMS .and. .not.IFRMS .and. ilstat /= jlstat) &
     &            Cycle

                Scal = UEFF(ilStat,iRoot1)*UEFF(jlStat,iRoot2)          &
     &               + UEFF(jlStat,iRoot1)*UEFF(ilStat,iRoot2)
                Scal = Scal*Half
                SLag(ilStat,jlStat) = SLag(ilStat,jlStat) + Scal

                If (ilStat /= jlStat) Then
                  SLag(jlStat,ilStat) = SLag(jlStat,ilStat) + Scal
                End If
              Else
                IF (IFXMS .or. IFRMS) Then
                  Scal = UEFF(ilStat,iRoot1)*UEFF(jlStat,iRoot2)
                  If (ilStat /= jlStat) Scal = Scal*Two
                  SLag(ilStat,jlStat) = SLag(ilStat,jlStat) + Scal
                Else
                  If (ilStat == jlStat) Then
                    Scal = UEFF(ilStat,iRoot1)*UEFF(jlStat,iRoot2)
                    SLag(ilStat,jlStat) = SLag(ilStat,jlStat) + Scal
                  End If
                End If
              End If
            End Do
          End Do
        End If

        !! Subtract the original rhs_sa.f or rhs_nac.f contribution
        !! For MS-type CASPT2, CASSCF part has to be determined by UEFF
        If (IFMSCOUP .and. iRoot1 == iRoot2) Then
          ilStat = MAX(iRoot1,iRoot2)
          jlStat = MIN(iRoot1,iRoot2)
          SLag(ilStat,jlStat) = SLag(ilStat,jlStat) - One
        End If

        !! Finalize the first-order transition(-like) density matrix
        !! for the CSF derivative term
        If (do_nac) Then
          If (do_csf) Then
            Call CnstAntiC(DPT2Canti_tot,UEFF,U0)
          Else
            !! Clear just in case
            DPT2Canti_tot(:) = Zero
          End If
        End If

        !! Back-transform the CI Lagrangian
        !! It is in the XMS basis, so it has to be transformed to
        !! CASSCF basis to be used in Z-vector
        !! No need to do this for SLag.
        If (IFXMS .or. IFRMS) Then
          call mma_allocate(CI1,nConf,nState,Label='CI1')
          Call DGEMM_('N','T',nConf,nState,nState,                      &
     &                One,CLagFull,nConf,U0,nState,                     &
     &                Zero,CI1,nConf)
          CLagFull(:,:) = CI1(:,:)
          call mma_deallocate(CI1)
        End If

        !! Compute true unrelaxed properties for MS-CASPT2
        if (.not.do_nac .and. ifmscoup) CALL PRPCTL(1,UEFF,U0,nState)

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
          If (.not.King()) OLagFull(:) = Zero
          CALL GADGOP (OLagFull,nOLag,'+')
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
          If (.not.King()) WLag(:) = Zero
          CALL GADGOP (WLag,nWLag,'+')
        end if
#endif
        If (DEB) call TriPrt('WLag', '', WLag, nBast)
        Do i = 1, nWLag ! = NBTRI
          Write (LuPT2,*) WLag(i)
        End Do

        !! D^PT2 in MO (read in OUT_PT2)
#ifdef _MOLCAS_MPP_
        if (is_real_par()) then
          If (.not.King()) DPT2_tot(:) = Zero
          CALL GADGOP (DPT2_tot,NBSQT,'+')
        end if
#endif
        If (DEB) call RecPrt('DPT2', '', DPT2_tot, nBast, nBast)
        Do i = 1, NBSQT
          Write (LuPT2,*) DPT2_tot(i)
        End Do

        !! D^PT2(C) in MO (read in OUT_PT2)
#ifdef _MOLCAS_MPP_
        if (is_real_par()) then
          If (.not.King()) DPT2C_tot(:) = Zero
          CALL GADGOP (DPT2C_tot,NBSQT,'+')
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
        if (RFpert .and. IFMSCOUP) then
          !! Recompute DPT2AO and DPT2AO for PCM
          Call Recompute_DPT2AO(NBSQT,DPT2_tot,DPT2C_tot,               &
     &                          DPT2_AO_tot,DPT2C_AO_tot)
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
            call daxpy_(NBTRI,Half,DPT2C_AO_tot,1,DPT2_AO_tot,1)
          else
            call daxpy_(NBTRI,One,DPT2C_AO_tot,1,DPT2_AO_tot,1)
          end if
          Call Put_dArray('D1aoVar',DPT2_AO_tot,NBTRI)
        else
          !! not sure this is OK
          Call Put_dArray('D1aoVar',DPT2_AO_tot,0)
        end if

        ! close gradient files
        Close (LuPT2)

        call mma_deallocate(HEFF1)
        call mma_deallocate(WRK1)
        call mma_deallocate(WRK2)
      end if

      call mma_deallocate(DPT2_tot)
      call mma_deallocate(DPT2C_tot)
      call mma_deallocate(DPT2_AO_tot)
      call mma_deallocate(DPT2C_AO_tot)

      call mma_deallocate(CLag)
      call mma_deallocate(CLagFull)
      call mma_deallocate(OLag)
      call mma_deallocate(OLagFull)
      call mma_deallocate(SLag)
      call mma_deallocate(WLag)

      call mma_deallocate(FIMO_all)
      call mma_deallocate(FIFA_all)

      If (IFXMS .or. IFRMS)        call mma_deallocate(FIFASA_all)
      If (IFDW .and. zeta >= Zero) call mma_deallocate(OMGDER)
      If (do_nac)                  call mma_deallocate(DPT2Canti_tot)

      !! Prepare for MCLR
      iGo = 3
      Call Put_iScalar('SA ready',iGo)
      ! overwrites whatever was set in CASSCF with the relax
      ! root that was chosen in CASPT2
      if (do_nac) then
!       write (u6 'NAC'
!       write (u6 'CASSCF/Original = ', iRoot1,iRoot2
        Call Put_iScalar('Relax CASSCF root',iRoot1)
        Call Put_iScalar('Relax Original root',iRoot2)
        call Qpg_cArray('MCLR Root',Found,I)
        if (Found) then
          Call Get_cArray('MCLR Root',mstate1,16)
          if (mstate1(8:8) == '0' .and. mstate1(16:16) == '0') then
            !! NAC states have not been specified
            write (mstate1,'(1X,I7,1X,I7)') iRoot1,iRoot2
            Call Put_cArray('MCLR Root',mstate1,16)
          end if
        else
          mstate1 = '****************'
          Call Put_cArray('MCLR Root',mstate1,16)
        end if
      else
!       write (u6 'GRD'
!       write (u6 'CASSCF/Original = ', irlxroot,irlxroot
        Call Put_iScalar('Relax CASSCF root',irlxroot)
        Call Put_iScalar('Relax Original root',irlxroot)
        mstate1 = '****************'
        Call Put_cArray('MCLR Root',mstate1,16)
      end if

      call ModDip()

      !! Close files
      Call DaClos(LUSTD)
      If (IfChol) Call DaClos(LUAPT2)
      Call DaClos(LUGRAD)

      if (nFroT /= 0) call mma_deallocate(TraFro)
      call mma_deallocate(iTasks_grad)

      Return

      End Subroutine GrdCls
