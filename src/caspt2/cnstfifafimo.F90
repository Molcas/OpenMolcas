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
      Subroutine CnstFIFAFIMO(MODE)

      use caspt2_global, only: TraFro, OLag,                            &
     &                           FIMO_all, FIFA_all, FIFASA_all
      use caspt2_global, only: FIMO, FIFA, CMOPT2
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp, iwp
      use Constants, only: Zero
      use caspt2_module, only: IfChol, IFXMS, IFRMS, IFDW, NSYM, NFRO,  &
     &                         NFROT, NBAS, NBSQT

      implicit none

      integer(kind=iwp), intent(in) :: MODE

      real(kind=wp),allocatable :: WRK1(:), WRK2(:)
      integer(kind=iwp) :: iSQ, iTr, iSym, nBasI

      If (IfChol) Then
        !! For DF or CD, we already have FIFA and FIMO in AO,
        !! so just do AO -> MO transformation
        call mma_allocate(WRK1,NBSQT,Label='WRK1')
        call mma_allocate(WRK2,NBSQT,Label='WRK2')
        WRK1(:) = Zero
        WRK2(:) = Zero

        iSQ = 0
        iTR = 0
        Do iSym = 1, nSym
          ! nOrbI = nOrb(iSym)
          nBasI = nBas(iSym)
          !! FIFA
          If (nFroT == 0) Then
            If (MODE == 0 .and. (IFDW .or. IFRMS)) Then
              Call SQUARE(FIFA(1+iTr),FIFASA_all(1+iSQ),1,nBasI,nBasI)
            Else If (MODE == 1) Then
              Call SQUARE(FIFA(1+iTr),FIFA_all(1+iSQ),1,nBasI,nBasI)
!             write (u6 'fifa in MO'
!             call sqprt(fifa_all(1+isq),nbasi)
            End If
          Else
            Call SQUARE(FIFA_all(1+iTr),WRK1,1,nBasI,nBasI)
!             write (u6 'fifasa in AO'
!             call sqprt(wrk1,nbasi)
            If (MODE == 0 .and. (IFDW .or. IFRMS)) Then
              !! with the state-average
              !! FIFASA_all will be natural basis
              Call OLagTrf(2,iSym,NBSQT,CMOPT2,FIFASA_all(1+iSQ),WRK1,  &
     &                     WRK2)
!             write (u6 'fifasa in MO'
!             call sqprt(fifasa_all(1+isq),nbasi)
            Else If (MODE == 1) Then
              !! with the state-specific or dynamically weighted
              !! FIFA will be quasi-canonical basis
              Call OLagTrf(2,iSym,NBSQT,CMOPT2,FIFA_all(1+iSQ),WRK1,    &
     &                     WRK2)
!             write (u6 'fifa in MO'
!             call sqprt(fifa_all(1+isq),nbasi)
              !! canonicalize frozen orbitals
              !! still under investigation, but this is something we
              !! should do to obtain "better" orbital enegies for
              !! methods using state-dependent Fock operators.
              !! We actually need to canonicalize frozen and inactive
              !! orbitals simultaneously?
              If (nFroT /= 0) Then
                WRK2(1:nBasI*nBasI) = WRK1(1:nBasI*nBasI)
                CALL DIAFCK(NBAS(ISYM),FIFA_all,1,NFRO(ISYM),           &
     &                      TraFro,NBAS(ISYM),CMOPT2,WRK2)
                CMOPT2(1:NBAS(ISYM)*NFRO(ISYM))                         &
     &            = WRK2(1:NBAS(ISYM)*NFRO(ISYM))
                Call OLagTrf(2,iSym,NBSQT,CMOPT2,FIFA_all(1+iSQ),WRK1,  &
     &                       WRK2)
              End If
            End If
          End If

          !! FIMO
!         If (MODE == 0) Then
            If (nFroT == 0) Then
              Call SQUARE(FIMO(1+iTr),FIMO_all(1+iSQ),1,nBasI,nBasI)
            Else
              Call SQUARE(FIMO_all(1+iTr),WRK1,1,nBasI,nBasI)
              Call OLagTrf(2,iSym,NBSQT,CMOPT2,FIMO_all(1+iSQ),WRK1,    &
     &                     OLag)
!             write (u6 'fimo in MO'
!             call sqprt(fimo_all(1+isq),nbasi)
            End If
!         End If
          iSQ = iSQ + nBasI*nBasI
          iTR = iTR + nBasI*(nBasI+1)/2
        End Do
        call mma_deallocate(WRK1)
        call mma_deallocate(WRK2)

        If (IFXMS .and. .not.IFDW)                                      &
     &    FIFASA_all(1:NBSQT) = FIFA_all(1:NBSQT)
      Else
        If (nFroT /= 0) Then
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
          If (IFXMS .and. .not.IFDW)                                    &
     &      FIFASA_all(1:NBSQT) = FIFA_all(1:NBSQT)
        End If

      !! XDW or RMS case: call after XDWINI
      ! If (MODE == 0) Then
      ! End If

      !! XMS case: call after GRPINI
      ! If (MODE == 1) Then
      ! End If

!     !! SS or MS case: call in dens.f
!     If (MODE == 2) Then
!       If (IFSADREF) Then
!       Else
!       End If
!     End If
      End If

      Return

      End Subroutine CnstFIFAFIMO
