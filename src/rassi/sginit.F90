!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine SGInit(nSym,nActEl,iSpin,nRasPrt,nRas,nRasEl,SGS)
      use stdalloc, only: mma_allocate, mma_deallocate
      use Struct, only: LEVEL, SGStruct
      IMPLICIT None
#include "rassi.fh"
      Integer nSym, nActEl, iSpin, nRasPrt
      Type (SGStruct), Target :: SGS
      Integer nRas(8,nRasPrt),nRasEl(nRasPrt)

      Integer, Allocatable:: Lim(:)
      Integer IAC,iErr,iLev,iRO,iSy,iSym,it,iTabs, nVert0,Lev

      Interface
      Subroutine MKDRT0(SGS)
      use struct, only: SGStruct
      IMPLICIT None
      Type(SGStruct), Target:: SGS
      End Subroutine MKDRT0
      End Interface

      Associate ( nLev => SGS%nLev, nVert => SGS%nVert, MidLev => SGS%MidLev, &
                  MVSta => SGS%MVSta, MvEnd => SGS%MVEnd,    &
                  IA0=>SGS%IA0, IB0=>SGS%IB0, IC0=>SGS%IC0, &
                  nVert0=>SGS%nVert0)

      NLEV=NASHT
! Allocate Level to Symmetry table ISm:
      Call mma_allocate(SGS%ISm,nLev,Label='SGS%ISm')
      ITABS=0
      DO ISYM=1,NSYM
        DO IT=1,NASH(ISYM)
          ITABS=ITABS+1
          ILEV=LEVEL(ITABS)
          SGS%ISM(ILEV)=ISYM
        END DO
      END DO

! Compute size of unrestricted DRT table:
      ib0=ispin-1
      ia0=(nActEl-ib0)/2
      ic0=nLev-ia0-ib0

      iErr=0
      If ((2*ia0+ib0).ne.nActEl) Then
        iErr=1
      Else If((ia0.lt.0).or.(ib0.lt.0).or.(ic0.lt.0)) Then
        iErr=1
      End If
      If(iErr.ne.0) then
        Write(6,*)' RASSI/SGINIT: Impossible input variables.'
        Write(6,*)'   nLev:',nLev
        Write(6,*)' nActEl:',nActEl
        Write(6,*)'  iSpin:',iSpin
        Write(6,*)'Program stops, sorry.'
        CALL ABEND()
      End If

! start of mkguga-like code
      IAC=MIN(IA0,IC0)
      NVERT0=((IA0+1)*(IC0+1)*(2*IB0+IAC+2))/2-(IAC*(IAC+1)*(IAC+2))/6

! Compute unrestricted DRT tables:
      CALL mma_allocate(SGS%DRT0,NVERT0,5,Label='DRT0')
      CALL mma_allocate(SGS%DOWN0,[1,NVERT0],[0,3],Label='DOWN0')
      nVert=nVert0

      SGS%DRTP => SGS%DRT0
      SGS%DOWNP => SGS%DOWN0

      CALL mkDRT0(SGS)

! Construct a restricted graph.
      Call mma_allocate(Lim,nLev,Label='Lim')
      Lim(:)=0
! Fill in the occupation limit table:
      Lev=0
      Do iRO=1,nRasPrt
        Do iSy=1,nSym
          Lev=Lev+nRas(iSy,iRO)
        End Do
        if(Lev.gt.0) Lim(Lev)=nRasEl(iRO)
      End Do

      Call mma_allocate(SGS%Ver,nVert0,Label='V11')
      Call RmVert(nLev,nVert,SGS%DRT0,SGS%Down0,Lim,SGS%Ver)
      Call mma_deallocate(Lim)

      Call mma_allocate(SGS%DRT,nVert,5,Label='SGS%DRT')
      Call mma_allocate(SGS%Down,[1,nVert],[0,3],Label='SGS%Down')
      Call mkDRT(SGS)
      Call mma_deallocate(SGS%Ver)
      Call mma_deallocate(SGS%DRT0)
      Call mma_deallocate(SGS%Down0)

! Direct Arc Weights table:
      CALL MKDAW(SGS)

! Upchain Index table:
! Reverse Arc Weights table:
      Call MKRAW(SGS)

! Level-To-Vertex table:
      CALL MKLTV(SGS)

      Call MKMID(SGS)

! Modified Arc Weights table:
      CALL MKMAW(SGS)

! The DAW, RAW tables are no longer needed:
      CALL mma_deallocate(SGS%RAW)
      CALL mma_deallocate(SGS%DAW)

      End Associate

      end Subroutine SGInit
