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

      Subroutine CnstInt(Mode,INT1,INT2)

      Use CHOVEC_IO, only: NVLOC_CHOBATCH
      use caspt2_global, only: FIMO_all
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: iwp,wp
      use caspt2_module, only: IfChol, NSYM, NFRO, NISH, NASH, NASHT,   &
     &                         NBAS, NBAST, NBSQT, NBTCH, NBTCHES
!     use caspt2_module, only: NSSH
      use Constants, only: Zero, One, Half
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif

      implicit none

      integer(kind=iwp), intent(in) :: Mode
      real(kind=wp), intent(out) :: INT1(nAshT,nAshT),                  &
     &                              INT2(nAshT,nAshT,nAshT,nAshT)

      integer(kind=iwp),allocatable :: BGRP(:,:)
      real(kind=wp),allocatable :: WRK1(:),WRK2(:),KET(:)

      integer(kind=iwp), parameter :: Inactive=1, Active=2, Virtual=3
      integer(kind=iwp) ::  iSym, nFroI, nIshI, nCorI, nBasI,           &
     &  iAshI, jAshI, iSymA, iSymI, iSymB, iSymJ, JSYM, IB, IB1, IB2,   &
     &  MXBGRP, IBGRP, NBGRP, NCHOBUF, MXPIQK, NADDBUF, IBSTA, IBEND,   &
     &  NV, nKET, kAshI, lAshI, iT, iU, iTU, iV, iX, iVX, iOrb, jOrb
!     integer(kind=iwp) :: nSh(8,3)
      real(kind=wp) :: Val, SCAL

      INT1(:,:) = Zero
      Int2(:,:,:,:) = Zero

      iSym=1
      nFroI = nFro(iSym)
      nIshI = nIsh(iSym)
      nCorI = nFroI+nIshI
      nBasI = nBas(iSym)
      ! nOrbI = nOrb(iSym)

      call mma_allocate(WRK1,NBSQT,Label='WRK1')
      call mma_allocate(WRK2,NBSQT,Label='WRK2')
!
!     --- One-Electron Integral
!
      !! Read H_{\mu \nu}
      Do iAshI = 1, nAsh(iSym)
        Do jAshI = 1, nAsh(iSym)
          Val = FIMO_all(nCorI+iAshI+nBasI*(nCorI+jAshI-1))
          INT1(iAshI,jAshI) = INT1(iAshI,jAshI) + Val
        End Do
      End Do
!
!     --- Two-Electron Integral
!
      iSymA = 1
      iSymI = 1
      iSymB = 1
      iSymJ = 1

      If (IfChol) Then
!       nSh(1:nSym,Inactive) = NISH(1:nSym)
!       nSh(1:nSym,Active  ) = NASH(1:nSym)
!       nSh(1:nSym,Virtual ) = NSSH(1:nSym)
        DO JSYM=1,NSYM
          IB1=NBTCHES(JSYM)+1
          IB2=NBTCHES(JSYM)+NBTCH(JSYM)

          MXBGRP=IB2-IB1+1
          IF (MXBGRP <= 0) CYCLE
          call mma_allocate(BGRP,2,MXBGRP,Label='BGRP')
          IBGRP=1
          DO IB=IB1,IB2
           BGRP(1,IBGRP) = IB
           BGRP(2,IBGRP) = IB
           IBGRP=IBGRP+1
          END DO
          NBGRP=MXBGRP

          CALL MEMORY_ESTIMATE(JSYM,BGRP,NBGRP,                         &
     &                         NCHOBUF,MXPIQK,NADDBUF)
          call mma_allocate(KET,NCHOBUF,Label='KETBUF')
          Do IBGRP=1,NBGRP

            IBSTA=BGRP(1,IBGRP)
            IBEND=BGRP(2,IBGRP)

            NV=0
            DO IB=IBSTA,IBEND
              NV=NV+NVLOC_CHOBATCH(IB)
            END DO

            !! int2(tuvx) = (tu|vx)/2
            !! This can be computed without frozen orbitals
            Call Get_Cholesky_Vectors(Active,Active,JSYM,               &
     &                                KET,SIZE(KET),nKet,               &
     &                                IBSTA,IBEND)

            If (IBGRP == 1) SCAL = Zero
            If (IBGRP /= 1) SCAL = One
            Call DGEMM_('N','T',NASH(JSYM)**2,NASH(JSYM)**2,NV,         &
     &                  Half,KET,NASH(JSYM)**2,KET,NASH(JSYM)**2,       &
     &                  SCAL   ,INT2,NASH(JSYM)**2)
          End Do
          call mma_deallocate(KET)
          call mma_deallocate(BGRP)
        End Do
      Else
        Do iAshI = 1, nAsh(iSym)
          iOrb = nCorI+iAshI
          Do jAshI = 1, nAsh(iSym)
            jOrb = nCorI+jAshI

            Call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
            !! Put in INT2
            Do kAshI = 1, nAsh(iSym)
              Do lAshI = 1, nAsh(iSym)
                INT2(iAshI,jAshI,kAshI,lAshI)                           &
     &        = INT2(iAshI,jAshI,kAshI,lAshI)                           &
     &        + WRK1(nCorI+kAshI+nBasT*(nCorI+lAshI-1))*Half
              End Do
            End Do
          End Do
        End Do
      End If
      call mma_deallocate(WRK1)
      call mma_deallocate(WRK2)
      If (Mode == 0) Then
      Do IT = 1, nAshT
        Do iU = 1, nAshT
          iTU = iT + nAshT*(iU-1)
          Do iV = 1, nAshT
            Do iX = 1, nAshT
              iVX = iV + nAshT*(iX-1)
              If (iVX > iTU) Then
               INT2(iT,iU,iV,IX) = INT2(iT,iU,iV,iX) + INT2(iV,iX,iT,iU)
               INT2(iV,iX,iT,iU) = Zero
              End If
            End Do
          End Do
        End Do
      End Do
      End If

#ifdef _MOLCAS_MPP_
      if (is_real_par()) CALL GADGOP (INT2,nAshT**4,'+')
#endif

      Return

      End Subroutine CnstInt
