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
      Type (SGStruct) SGS
      Integer nRas(8,nRasPrt),nRasEl(nRasPrt)

      Integer, Allocatable:: DRT0(:), Down0(:), Tmp(:), Lim(:), NWV(:), DAW(:), RAW(:)
      Integer IA0,IB0,IC0,IAC,iErr,iLev,iRO,iSy,iSym,it,iTabs,MidLev,  &
     &        MVSta,MVEnd,NDown0,nDrt0,nLev,nTmp,nVert0,Lev,nVert

      NLEV=NASHT
! Allocate Level to Symmetry table ISm:
      Call mma_allocate(SGS%ISm,nLev,Label='SGS%ISm')
      ITABS=0
      IT=0
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

      IAC=MIN(IA0,IC0)
      NVERT0=((IA0+1)*(IC0+1)*(2*IB0+IAC+2))/2-(IAC*(IAC+1)*(IAC+2))/6
      NDRT0=5*NVERT0
      NDOWN0=4*NVERT0

! Compute unrestricted DRT tables:
      CALL mma_allocate(DRT0,NDRT0,Label='DRT0')
      CALL mma_allocate(DOWN0,NDOWN0,Label='DOWN0')
      NTMP=((NLEV+1)*(NLEV+2))/2
      CALL mma_allocate(TMP,NTMP,Label='TMP')
      CALL mkDRT0 (IA0,IB0,IC0,NVERT0,DRT0,DOWN0,NTMP,TMP)

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
      Call mma_allocate(NWV,nVert0,Label='NWV')
      nVert=nVert0
      Call RmVert(nLev,nVert,DRT0,Down0,Lim,NWV)
      Call mma_deallocate(Lim)
      Call mma_allocate(SGS%DRT,5*nVert,Label='SGS%DRT')
      Call mma_allocate(SGS%Down,4*nVert,Label='SGS%Down')
      Call mkDRT(nVert0,nVert,DRT0,Down0,NWV,SGS%DRT,SGS%Down)
      Call mma_deallocate(NWV)
      Call mma_deallocate(DRT0)
      Call mma_deallocate(Down0)

! Level-To-Vertex table:
      Call mma_allocate(SGS%LTV,nLev+2,Label='SGS%LTV')
      CALL MKLTV(nVert,nLev,SGS%DRT,SGS%LTV)

! Direct Arc Weights table:
      Call mma_allocate(DAW,5*nVert,Label='DAW')
      CALL MKDAW(nVert,SGS%DOWN,DAW)

! Upchain Index table:
      Call mma_allocate(SGS%Up,4*nVert,Label='SGS%Up')
! Reverse Arc Weights table:
      Call mma_allocate(RAW,5*nVert,Label='RAW')
! Modified Arc Weights table:
      Call mma_allocate(SGS%MAW,4*nVert,Label='SGS%MAW')
      Call MkMAW_RASSI(nLev,nVert,SGS%Down,DAW,SGS%Up,RAW,SGS%MAW,SGS%LTV,MidLev)
      MVSta=SGS%LTV(2+MidLev)
      MVEnd=SGS%LTV(1+MidLev)-1
! The DAW, RAW tables are no longer needed:
      CALL mma_deallocate(RAW)
      CALL mma_deallocate(DAW)
      CALL mma_deallocate(TMP)

! Put sizes and addresses in structure SGS:

      SGS%nLev   =nLev
      SGS%nVert  =nVert
      SGS%MidLev =MidLev
      SGS%MVSta  =MVSta
      SGS%MVEnd  =MVEnd

      end Subroutine SGInit
