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
      Subroutine OLagVVVO(iSym,NBSQT,lT2AO,MaxVec_PT2,DPT2AO,DPT2CAO,   &
     &                    FPT2AO,FPT2CAO,T2AO,DIA,DI,FIFA,FIMO,A_PT2)
      USE iSD_data, only: iSD
      use caspt2_global, only: LuGAMMA,LuCMOPT2,LuAPT2,OLag
      use caspt2_global, only: CMOPT2
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp, u6
      use Constants, only: Zero, One
      use caspt2_module, only: IFMSCOUP, IFXMS, IFRMS, IFDW, IFSADREF,  &
     &                         NSYM, NFRO, NISH, NASH, NSSH,            &
     &                         NBAS, NBAST, NBMX, NSTATE, JSTATE,       &
     &                         iRlxRoot
#ifdef _MOLCAS_MPP_
      use caspt2_module, only: NFROT
      use caspt2_global, only: nOLag
      USE Para_Info, ONLY: Is_Real_Par
#endif

      implicit none

      integer(kind=iwp), intent(in) :: iSym, NBSQT, lT2AO, MaxVec_PT2
      real(kind=wp), intent(in) :: DPT2AO(NBSQT), DPT2CAO(NBSQT),       &
     &  T2AO(lT2AO), DIA(NBSQT), DI(NBSQT)
      real(kind=wp), intent(inout) :: FPT2AO(NBSQT), FPT2CAO(NBSQT),    &
     &  A_PT2(MaxVec_PT2**2), FIFA(NBSQT), FIMO(NBSQT)

      character(len=4096) :: RealName
      logical(kind=iwp) :: DoRys, DoCholesky, is_error, Square
      real(kind=wp), allocatable :: T_hbf(:,:,:,:), vLag(:), WRK1(:),   &
     &                              WRK2(:)
      integer(kind=iwp), allocatable :: iOffAO(:)
      integer(kind=iwp) :: nBasX(8), KEEP(8), IRC, nSymX, id, lRealName,&
     &  iost, iSymI, iSymJ, iSymA, iSymB, nocc, i, nSSDM, nDiff, nSkal, &
     &  MaxShlAO, iSh, nBasI, jSh, nBasJ, iBas0, iBas, jBas0, jBas,     &
     &  iOcc, jOcc, loc1, loc2, j, iRec, nOrbA
      integer(kind=iwp), external :: isFreeUnit
!
!     ----- (VV|VO)
!
!     Compute L_{pq} = (pa|jb) * T_{qj}^{ab}, in particular for
!     (p,q) \in  (virt, inact+act). This operation involves
!     (VV|VO) integrals which are stored in neither disk nor memory.
!     The back-transformation is applied to occupied orbital derivatives
!     of two-electron integrals that have two virtual indices (F, G, H).
!     In principle, the algorithm is to avoid (VV|VO) integrals,
!     i.e. U_{pq} for (p,q) = (vir, inact+act), but can also be applied
!     to U_{pq} for (p,q) = (all, inact+act).
!
      call mma_allocate(vLag,nBasT*nBasT,Label='vLag')
      call mma_allocate(WRK1,nBasT*nBasT,Label='WRK1')
      Call GetOrd(IRC,Square,nSymX,nBasX,KEEP)

      vLag(:) = Zero

      Call DecideOncholesky(DoCholesky)
      If (DoCholesky) Then
        !! No need to save CMOPT2. Just save A_PT2 and B_PT2.
        !! First, save A_PT2 in LuCMOPT2
        If (IFMSCOUP .and. jState /= 1) Then
          call mma_allocate(WRK2,MaxVec_PT2**2,Label='WRK2')

          ! read A_PT2 from LUAPT2
          id = 0
          call ddafile(LUAPT2, 2, WRK2, MaxVec_PT2**2, id)

          A_PT2(1:MaxVec_PT2**2) = A_PT2(1:MaxVec_PT2**2)               &
     &      + WRK2(1:MaxVec_PT2**2)
          call mma_deallocate(WRK2)
        End If

        ! For SS-CASPT2 I should write A_PT2 on disk only
        ! for the correct iRlxRoot
        if (jState == iRlxRoot .or. IFMSCOUP) then
          ! write A_PT2 in LUAPT2
          id = 0
          call ddafile(LUAPT2,1,A_PT2,MaxVec_PT2**2,id)
        end if

        !rewind(LuGamma)
        Call PrgmTranslate('GAMMA',RealName,lRealName)
        LuGAMMA = isFreeUnit(LuGAMMA)
        Call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),            &
     &                       'DIRECT','UNFORMATTED',                    &
     &                        iost,.TRUE.,                              &
     &                        nBas(iSym)**2*8,'OLD',is_error)
      End If
      !! 2) Compute ERI (mu rho | nu sigma)
      !! 3) Quarter-transformation of ERI
      !!    (mu rho | j sigma) = sum_{j} C_{nu j} (mu rho | nu sigma)
      !! 4) Contract with AO amplitude
      !!    L_{mu i} = sum_{j,rho sigma} T_{ij}^{mu nu}*(mu rho|j sigma)
      !! the third argument is used wrongly, but should be no problem

      !! D_{mu nu} := DPT2AO
      !! FPT2AO = ((mu nu|rho sigma)-(mu rho|nu sigma)/4
      !!          -(mu sigma|nu rho)/4)*D_{rho sigma}
      isymi = 1
      isymj = 1
      isyma = 1
      isymb = 1
!     nocc = nfro(1)+nish(1)+nash(1)
      nocc = nish(1)+nash(1)
      Call VVVO_Drv(nSym,nBas,nIsh,nFro,KEEP,                           &
     &              iSym,iSymI,iSymA,iSymJ,iSymB,                       &
     &              lT2AO,T2AO,vLag,                                    &
     &              nOcc,nBasT,nBMX,                                    &
     &              CMOPT2(1+nBasT*nFro(iSymA)),                        &
     &              DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,                      &
     &              DIA,DI,FIFA,FIMO)

      !! Save the half transformed integrals on disk.
      !! It will be used in drvg1.f etc for gradient.
      If (DoCholesky) Then
        !! Do nothing
        Close (LuGAMMA)
      Else
        !! This is only for conventional calculations!!
        Call PrgmTranslate('CMOPT2',RealName,lRealName)
        LuCMOPT2 = isFreeUnit(LuCMOPT2)
        Call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),           &
     &                       'DIRECT','UNFORMATTED',                    &
     &                        iost,.FALSE.,1,'OLD',is_error)
        !! First, CMOPT2 has to be saved. The MO coefficient matrix in
        !! grvg1.f may be different from CMOPT2.
        Do i = 1, nBasT*nBasT
          Write (LuCMOPT2) CMOPT2(i)
        End Do
        write (LuCMOPT2) (nIsh(1)+nAsh(1))
        write (LuCMOPT2) (nIsh(2)+nAsh(2))
        write (LuCMOPT2) (nIsh(3)+nAsh(3))
        write (LuCMOPT2) (nIsh(4)+nAsh(4))
        write (LuCMOPT2) (nIsh(5)+nAsh(5))
        write (LuCMOPT2) (nIsh(6)+nAsh(6))
        write (LuCMOPT2) (nIsh(7)+nAsh(7))
        write (LuCMOPT2) (nIsh(8)+nAsh(8))
        write (LuCMOPT2) (nFro(1))
        write (LuCMOPT2) (nFro(2))
        write (LuCMOPT2) (nFro(3))
        write (LuCMOPT2) (nFro(4))
        write (LuCMOPT2) (nFro(5))
        write (LuCMOPT2) (nFro(6))
        write (LuCMOPT2) (nFro(7))
        write (LuCMOPT2) (nFro(8))
        !! Number of state-specific density matrix
        !! It is needed, because the separable contribution to the
        !! electron repulsion integral is different
        nSSDM = 0
        If (.not.IFSADREF) Then
          If (nState == 1) Then
            nSSDM = 0
          Else If (IFDW .or. IFRMS) Then
            !! For (X)DW, use nState density matrix
            nSSDM = nState
          Else If (IFXMS) Then
            !! For XMS, use SA density matrix
            nSSDM = 0
          Else If (IFMSCOUP) Then
            !! For MS, use nState density matrix
            nSSDM = nState
          Else
            !! Otherwise, use one SS density matrix
            nSSDM = 1
          End If
        End If
        write (LuCMOPT2) nSSDM

        Close (LuCMOPT2)

!       write(u6,*) 'mo saved'
!       call sqprt(CMOPT2,nbast)
!
        Call PrgmTranslate('GAMMA',RealName,lRealName)
        LuGAMMA = isFreeUnit(LuGAMMA)
        Call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),            &
     &                       'DIRECT','UNFORMATTED',                    &
     &                        iost,.TRUE.,                              &
     &                        nOcc*nOcc*8,'OLD',is_error)
        if (is_error) then
        write (u6,*) 'Something is wrong in opening LuGamma in olagvvvo'
          call abend()
        end if
        !  Setup for shell. Why do I have to call IniSew damn here?
        !  The number of shells should be able to be referred globally.
        nDiff=1
        DoRys=.True.
        Call IniSew(DoRys,nDiff)
        Call Nr_Shells(nSkal)
        Call Setup_iSD()
        !! see Include/info.fh
        call mma_allocate(iOffAO,nSkal+1,Label='iOffAO')
        MaxShlAO = 0
        iOffAO(1) = 0
        Do iSh = 1, nSkal
          nBasI = iSD(2,iSh)*iSD(3,iSh)
          If (nBasI > MaxShlAO) MaxShlAO = nBasI
          iOffAO(iSh+1) = iOffAO(iSh)+nBasI
        End Do
        ! nMax = MaxShlAO*MaxShlAO*nOcc*nOcc

        call mma_allocate(T_hbf,nOcc,nOcc,MaxShlAO,MaxShlAO,            &
     &                    Label='T_hbf')
        Do iSh = 1, nSkal
          !! iSD(2,iSh): number of AOs of the shell
          !! iSD(3,iSh): number of cont. func. of the shell
          nBasI = iSD(2,iSh)*iSD(3,iSh)
          Do jSh = 1, nSkal
            nBasJ = iSD(2,jSh)*iSD(3,jSh)
            Do iBas0 = 1, nBasI
              iBas = iOffAO(iSh) + iBas0
              Do jBas0 = 1, nBasJ
                jBas = iOffAO(jSh) + jBas0
                Do iOcc = 1, nOcc
                  Do jOcc = 1, nOcc
                    loc1 = jOcc-1 + (jBas-1)*nOcc                       &
     &                   + (iOcc-1)*nOcc*nBasT                          &
     &                   + (iBas-1)*nOcc*nBasT*nOcc
                    loc2 = iOcc-1 + (iBas-1)*nOcc                       &
     &                   + (jOcc-1)*nOcc*nBasT                          &
     &                   + (jBas-1)*nOcc*nBasT*nOcc
                    T_hbf(iOcc,jOcc,iBas0,jBas0)                        &
     &                = T2AO(loc1+1)+T2AO(loc2+1)
                  End Do
                End Do
                iRec = iBas+nBasT*(jBas-1)
                if (ifmscoup .and. jstate /= 1) then
                  read (lugamma,rec=irec) (wrk1(i),i=1,nocc*nocc)
                  call daxpy_(nocc*nocc,One,wrk2,1,                     &
     &                        t_hbf(1,1,ibas0,jbas0),1)
                end if
                Write (LuGamma,Rec=iRec)                                &
     &            ((T_hbf(i,j,iBas0,jBas0),i=1,nOcc),j=1,nOcc)
              End Do
            End Do
          End Do
        End Do
        call mma_deallocate(iOffAO)
        call mma_deallocate(T_hbf)
        Close (LuGAMMA)
        Call Free_iSD()
        call clssew()
      End If

      !! 5) L_{ai} = sum_{mu} C_{mu a} L_{mu i}
!     CALL DGEMM_('T','N',nSsh(iSym),nOcc,nBasT,
!    *            One,CMOPT2(1+nBasT*nOcc),nBasT,
!    *                vLag,nBasT,
!    *            One,OLAG(nOCC+1),nOrb(iSymA))
!     write(u6,*) 'olag before vvvo'
!     call sqprt(olag,nbast)
      nOrbA = nFro(iSymA)+nIsh(iSymA)+nAsh(iSymA)+nSsh(iSymA)
      If (DoCholesky) nOcc = nOrbA-nFro(iSymA)
!     write(u6,*) 'vLag'
!     call sqprt(vLag,nbast)
      CALL DGEMM_('T','N',nOrbA,nOcc,nBasT,                             &
     &            One,CMOPT2,nBasT,vLag,nBasT,                          &
     &            One,OLAG(nOrbA*nFro(iSymA)+1),nOrbA)
!     write(u6,*) 'olag after vvvo'
!     call sqprt(olag,nbast)
!
#ifdef _MOLCAS_MPP_
      If (Is_Real_Par()) Then
        if (DoCholesky) call GADGOP(OLag,nOLag,'+')
        If (nFroT == 0) Then
          CALL GADGOP (FPT2AO,nBasT**2,'+')
          CALL GADGOP (FPT2CAO,nBasT**2,'+')
        End If
      End If
#endif
      call mma_deallocate(vLag)
      call mma_deallocate(WRK1)

      END SUBROUTINE OLagVVVO
