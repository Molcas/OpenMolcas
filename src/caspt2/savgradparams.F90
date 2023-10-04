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
! Copyright (C) 2023, Yoshio Nishimoto                                 *
!***********************************************************************

Subroutine SavGradParams(Mode,IDSAVGRD)

  use caspt2_gradient, only: LUGRAD, LUSTD, do_lindep, IDBoriMat
  use definitions, only: iwp,wp,byte
  use stdalloc, only: mma_allocate, mma_deallocate

  Implicit None

#include "rasdim.fh"
#include "caspt2.fh"
#include "eqsolv.fh"
#include "pt2_guga.fh"

#include "WrkSpc.fh"

  integer(kind=iwp), intent(in) :: Mode
  integer(kind=iwp), intent(inout) :: IDSAVGRD

  integer(kind=iwp) :: IORW,ID,NIN,NAS,NIS,NNN,NMAX,ISYM,ICASE,iLUID

  real(kind=wp), allocatable :: WRK1(:)
  integer(kind=iwp), allocatable :: IWRK1(:)
  integer(kind=byte), allocatable :: idxG3(:,:)

! character(len=80) :: Label

  !! Shift the address due to SavGradParams2
  If (IDSAVGRD==0) IDSAVGRD = NSTATE + 3*NSTATE**2

  !! Decide what to do
  If (Mode.eq.1) Then
    IORW = 1 !! Write
  Else If (Mode.eq.2) Then
    IORW = 2 !! Read
  End If

  !! Save internal contractions-related quantities
  !! 1. Some integers
  !! - Number of independent vectors
  CALL IDAFILE(LUGRAD,IORW,NINDEP(1,1),8*MXCASE,IDSAVGRD)
  !! - Number of active indices
  CALL IDAFILE(LUGRAD,IORW,NASUP(1,1),8*MXCASE,IDSAVGRD)
  !! - Number of inactive + secondary indices
  CALL IDAFILE(LUGRAD,IORW,NISUP(1,1),8*MXCASE,IDSAVGRD)

  !! 2. RDMs
  !! - NG1, NG2, NG3
  Call mma_allocate(IWRK1,4,Label='IWRK1')
  If (IORW.eq.1) Then
    IWRK1(1) = NG1
    IWRK1(2) = NG2
    IWRK1(3) = NG3
    IWRK1(4) = NG3TOT
    CALL IDAFILE(LUGRAD,IORW,IWRK1,4,IDSAVGRD)
  Else If (IORW.eq.2) Then
    CALL IDAFILE(LUGRAD,IORW,IWRK1,4,IDSAVGRD)
    NG1    = IWRK1(1)
    NG2    = IWRK1(2)
    NG3    = IWRK1(3)
    NG3TOT = IWRK1(4)
  End If
  Call mma_deallocate(IWRK1)

  NMAX = 0

  Do ISYM = 1, NSYM
    Do ICASE = 1, 11
      NIN = NINDEP(ISYM,ICASE)
      NAS = NASUP(ISYM,ICASE)
      NIS = NISUP(ISYM,ICASE)
      NNN = NAS*(NAS+1)/2
      NNN = MAX(NNN,NAS*NIN)
      NNN = MAX(NNN,NIS)
      NMAX= MAX(NNN,NMAX)
    End Do
    Do ICASE = 12, 13
      NAS = NASUP(ISYM,ICASE)
      NIS = NISUP(ISYM,ICASE)
      NMAX= MAX(NAS*NIS,NMAX)
    End Do
  End Do
  NMAX = MAX(NMAX,NG3)

  Call mma_allocate(WRK1,NMAX,Label='WRK1')

  CALL mma_allocate(idxG3,6,NG3,label='idxG3')
  If (IORW.eq.1) Then
    !! NG3 index
    iLUID = 0
    CALL I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)
    CALL I1DAFILE(LUGRAD,1,idxG3,6*NG3,IDSAVGRD)

    !! D and F
    CALL PT2_GET(NG1,' GAMMA1',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG1,IDSAVGRD)
    CALL PT2_GET(NG2,' GAMMA2',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG2,IDSAVGRD)
    CALL PT2_GET(NG3,' GAMMA3',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG3,IDSAVGRD)

    CALL PT2_GET(NG1,' DELTA1',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG1,IDSAVGRD)
    CALL PT2_GET(NG2,' DELTA2',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG2,IDSAVGRD)
    CALL PT2_GET(NG3,' DELTA3',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG3,IDSAVGRD)

    !! EASUM
    WRK1(1) = EASUM
    CALL DDAFILE(LUGRAD,IORW,WRK1,1,IDSAVGRD)
  Else If (IORW.eq.2) Then
    !! NG3 index
    CALL I1DAFILE(LUGRAD,2,idxG3,6*NG3,IDSAVGRD)
    iLUID = 0
    CALL I1DAFILE(LUSOLV,1,idxG3,6*NG3,iLUID)

    !! D and F
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG1,IDSAVGRD)
    CALL PT2_PUT(NG1,' GAMMA1',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG2,IDSAVGRD)
    CALL PT2_PUT(NG2,' GAMMA2',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG3,IDSAVGRD)
    CALL PT2_PUT(NG3,' GAMMA3',WRK1)

    CALL DDAFILE(LUGRAD,IORW,WRK1,NG1,IDSAVGRD)
    CALL PT2_PUT(NG1,' DELTA1',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG2,IDSAVGRD)
    CALL PT2_PUT(NG2,' DELTA2',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG3,IDSAVGRD)
    CALL PT2_PUT(NG3,' DELTA3',WRK1)

    !! EASUM
    CALL DDAFILE(LUGRAD,IORW,WRK1,1,IDSAVGRD)
    EASUM = WRK1(1)
    CALL GETDPREF(WORK(LDREF),WORK(LPREF))
    EREF=REFENE(JSTATE)
  End If
  Call mma_deallocate(idxG3)

  !! 3. LUSBT matrices
  Do ISYM = 1, NSYM
    Do ICASE = 1, 11
      NIN = NINDEP(ISYM,ICASE)
      NAS = NASUP(ISYM,ICASE)
      NIS = NISUP(ISYM,ICASE)
      If (IORW.eq.1) Then
        !! Active overlap
        ID = IDSMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,2,WRK1,NAS*(NAS+1)/2,ID)
        CALL DDAFILE(LUGRAD,1,WRK1,NAS*(NAS+1)/2,IDSAVGRD)
        !! ST matrix
        ID = IDSTMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,2,WRK1,NAS*NIN,ID)
        CALL DDAFILE(LUGRAD,1,WRK1,NAS*NIN,IDSAVGRD)
        !! Transformation matrix (eigenvector)
        ID = IDTMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,2,WRK1,NAS*NIN,ID)
        CALL DDAFILE(LUGRAD,1,WRK1,NAS*NIN,IDSAVGRD)
        !! Eigenvalue
        ID  = IDBMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,2,WRK1,NAS,ID)
        CALL DDAFILE(LUGRAD,1,WRK1,NAS,IDSAVGRD)
        !! IS
!       CALL DDAFILE(LUSBT,2,WRK1,NIS,ID)
!       CALL DDAFILE(LUGRAD,1,WRK1,NIS,IDSAVGRD)
        if (do_lindep) then
          ID  = IDBoriMat(ISYM,ICASE)
          CALL DDAFILE(LUSTD,2,WRK1,NAS*(NAS+1)/2,ID)
          CALL DDAFILE(LUGRAD,1,WRK1,NAS*(NAS+1)/2,IDSAVGRD)
        end if
      Else If (IORW.eq.2) Then
        !! Active overlap
        CALL DDAFILE(LUGRAD,2,WRK1,NAS*(NAS+1)/2,IDSAVGRD)
        ID = IDSMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,WRK1,NAS*(NAS+1)/2,ID)
        !! ST matrix
        CALL DDAFILE(LUGRAD,2,WRK1,NAS*NIN,IDSAVGRD)
        ID = IDSTMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,WRK1,NAS*NIN,ID)
        !! Transformation matrix (eigenvector)
        CALL DDAFILE(LUGRAD,2,WRK1,NAS*NIN,IDSAVGRD)
        ID = IDTMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,WRK1,NAS*NIN,ID)
        !! Eigenvalue
        CALL DDAFILE(LUGRAD,2,WRK1,NAS,IDSAVGRD)
        ID  = IDBMAT(ISYM,ICASE)
        CALL DDAFILE(LUSBT,1,WRK1,NAS,ID)
        !! IS
!       CALL DDAFILE(LUGRAD,2,WRK1,NIS,IDSAVGRD)
!       CALL DDAFILE(LUSBT,1,WRK1,NIS,ID)
        if (do_lindep) then
          CALL DDAFILE(LUGRAD,2,WRK1,NAS*(NAS+1)/2,IDSAVGRD)
          ID  = IDBoriMat(ISYM,ICASE)
          CALL DDAFILE(LUSTD,1,WRK1,NAS*(NAS+1)/2,ID)
        end if
      End If
    End Do
  End Do

  !! 4. E2TOT
  If (IORW.eq.1) Then
    WRK1(1) = E2TOT
    CALL DDAFILE(LUGRAD,IORW,WRK1,1,IDSAVGRD)
  Else
    CALL DDAFILE(LUGRAD,IORW,WRK1,1,IDSAVGRD)
    E2TOT = WRK1(1)
  End If

  call mma_deallocate(WRK1)

  !! quasi-canonical active orbital energy
! CALL DDAFILE(LUGRAD,IORW,EPSA,NASHT,IDSAVGRD)

  !! 5. Save T-amplitude (IRHS = LURHS(1) = 51)
  Call SaveReadT1()

Contains

  Subroutine SaveReadT1()

    Implicit None

    integer(kind=iwp) :: lg_V1,NVEC

    IVECX = 2

    Do ICASE = 1, NCASES
      Do ISYM = 1, NSYM
        NIN = NINDEP(ISYM,ICASE)
        If (NIN.EQ.0) Cycle
        NIS = NISUP(ISYM,ICASE)
        NAS = NASUP(ISYM,ICASE)
        NVEC = NIN*NIS
        If (ICASE.EQ.12.OR.ICASE.EQ.13) NVEC = NAS*NIS
        If (NVEC.EQ.0) Cycle
        !! lg_V1 = T (solution; not quasi-variational)
        If (ICASE.EQ.12.OR.ICASE.EQ.13) Then
          Call RHS_ALLO(NAS,NIS,lg_V1)
          If (IORW.EQ.1) Then
            Call RHS_READ_SR(lg_V1,ICASE,ISYM,IVECX)
            CALL DDAFILE(LUGRAD,IORW,WORK(lg_V1),NAS*NIS,IDSAVGRD)
          Else If (IORW.EQ.2) Then
            CALL DDAFILE(LUGRAD,IORW,WORK(lg_V1),NAS*NIS,IDSAVGRD)
            CALL RHS_SAVE_SR(lg_V1,ICASE,ISYM,IVECX)
          End If
          CALL RHS_FREE(NAS,NIS,lg_V1)
        Else
          Call RHS_ALLO(NIN,NIS,lg_V1)
          If (IORW.EQ.1) Then
            Call RHS_READ_SR(lg_V1,ICASE,ISYM,IVECX)
            CALL DDAFILE(LUGRAD,IORW,WORK(lg_V1),NIN*NIS,IDSAVGRD)
          Else If (IORW.EQ.2) Then
            CALL DDAFILE(LUGRAD,IORW,WORK(lg_V1),NIN*NIS,IDSAVGRD)
            CALL RHS_SAVE_SR(lg_V1,ICASE,ISYM,IVECX)
          End If
          CALL RHS_FREE(NIN,NIS,lg_V1)
        End If
      End Do
    End Do

  End Subroutine SaveReadT1

End Subroutine SavGradParams

Subroutine SavGradParams2(Mode,UEFF,U0,H0)

  use caspt2_gradient, only: LUGRAD
  use definitions, only: iwp,wp

  Implicit None

#include "rasdim.fh"
#include "caspt2.fh"

  integer(kind=iwp), intent(in) :: Mode
  real(kind=wp)    , intent(inout) :: UEFF(*),U0(*),H0(*)

  integer(kind=iwp) :: IORW,ID

  !! Decide what to do
  If (Mode.eq.1) Then
    IORW = 1 !! Write
  Else If (Mode.eq.2) Then
    IORW = 2 !! Read
  End If

  ID = 0
  CALL DDAFILE(LUGRAD,IORW,ENERGY,NSTATE   ,ID)
  CALL DDAFILE(LUGRAD,IORW,UEFF  ,NSTATE**2,ID)
  CALL DDAFILE(LUGRAD,IORW,U0    ,NSTATE**2,ID)
  CALL DDAFILE(LUGRAD,IORW,H0    ,NSTATE**2,ID)

End Subroutine SavGradParams2
