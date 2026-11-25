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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

!> The split graph of a CAS/RAS CSF space and its walks (CSFs): the order of
!> the CSFs in the CI vectors of the JOBIPH/HDF5 files. Independent of the
!> CI libraries, so that SGUGA (whose SGStruct/CIStruct extend these types)
!> and Faroald (which needs the file order at its boundary) do not depend on
!> each other.
!>
!> CSFG_Graph: the distinct-row table and its arc weights (see SGStruct of the
!> SGUGA module for the components); CSFG_Walks: the half-walk case lists by
!> midvertex and symmetry and the CSF counts and offsets (see CIStruct).
!> CSFG_Make builds both from nSym, nActEl, iSpin, nLev, ISm and the RAS
!> partitioning set in the graph.
module CSF_Graph

use Molcas, only: MxGas, MxLev, MxSym
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
private

integer(kind=iwp) :: iq

! Set CSFG_nPack to the number of cases (2 bit per case) that can be packed in one integer.
#ifdef SIZE_INITIALIZATION
integer(kind=iwp), parameter :: CSFG_nPack = storage_size(1_iwp)/2-1
#elif defined (_I8_)
integer(kind=iwp), parameter:: CSFG_nPack=32-1
#else
integer(kind=iwp), parameter:: CSFG_nPack=16-1
#endif

! index for DRT
integer(kind=iwp), parameter :: CSFG_LTab = 1, CSFG_NTab = 2, CSFG_ATab = 3, CSFG_BTab = 4, CSFG_CTab = 5

type CSFG_Graph
  integer(kind=iwp) :: NSym = 0, nActEl = 0, IFRAS = 0
  integer(kind=iwp) :: IA0, IB0, IC0, iSpin, nLev, nVert, nVert0, MidLev, MVSta, MVEnd, MXUP, MXDWN
  integer(kind=iwp), allocatable :: ISm(:), DRT(:,:), DRT0(:,:), Down(:,:), Down0(:,:), Up(:,:), Ver(:), MAW(:,:), LTV(:), &
                                    DAW(:,:), RAW(:,:)
  integer(kind=iwp), pointer :: DRTP(:,:), DOWNP(:,:)
  integer(kind=iwp) :: L2ACT(MXLEV)=[(iq,iq=1,MXLEV)]
  integer(kind=iwp) :: LEVEL(MXLEV)=[(iq,iq=1,MXLEV)]
  integer(kind=iwp) :: NRAS(MxSym,MxGAS), NRASEL(MxGAS), nRsPrt=0
  ! GAS: upper bound on the electrons in the partitions 1..k (-1: none)
  integer(kind=iwp) :: NRASMX(MxGAS) = -1
end type CSFG_Graph

type CSFG_Walks
  integer(kind=iwp) :: nMidV, nIpWlk, nWalk, NUW
  integer(kind=iwp), allocatable :: NOW(:,:,:), IOW(:,:,:), NCSF(:), NOCSF(:,:,:), IOCSF(:,:,:), ICase(:)
end type CSFG_Walks

public :: CSFG_ATab, CSFG_BTab, CSFG_Count, CSFG_CTab, CSFG_Free, CSFG_Graph, CSFG_LTab, CSFG_Make, CSFG_Mk_IOW, CSFG_MkCOT, &
          CSFG_nPack, CSFG_NTab, CSFG_Walk_Next, CSFG_Walks

contains

subroutine CSFG_Make(SGS,CIS)
! PURPOSE: MAKE THE GUGA TABLES
! NOTE:    TO RETAIN THE TABLES AVAILABLE FOR LATER PURPOSES
!          THE START ADDRESSES OF OF THE ARRAYS ETC. ARE STORED IN
!          THREE USER DEFINED TYPES. Consult the sguga module for the details.

  class(CSFG_Graph), target, intent(inout) :: SGS
  class(CSFG_Walks), intent(inout) :: CIS

  ! COMPUTE TOP ROW OF THE GUGA TABLE

  call mknVert0(SGS)

  ! SET UP A FULL PALDUS DRT TABLE:
  ! (INITIALLY NO RESTRICTIONS ARE PUT UP)

  SGS%nVert = SGS%nVert0
  if (SGS%IFRAS /= 0) then
    call mma_allocate(SGS%DRT0,SGS%nVert0,5,Label='DRT0')
    call mma_allocate(SGS%DOWN0,[1,SGS%nVert0],[0,3],Label='DOWN0')
    SGS%DRTP => SGS%DRT0
    SGS%DOWNP => SGS%DOWN0
  else
    call mma_allocate(SGS%DRT,SGS%nVert,5,Label='SGS%DRT')
    call mma_allocate(SGS%DOWN,[1,SGS%nVert],[0,3],Label='SGS%DOWN')
    SGS%DRTP => SGS%DRT
    SGS%DOWNP => SGS%DOWN
  end if

  call mkDRT0(SGS)

  ! IF THIS IS A RAS CALCULATION PUT UP RESTRICTIONS BY DELETING
  ! VERTICES WHICH VIOLATE THE FORMER.

  if (SGS%IFRAS /= 0) then
    call RmVert(SGS) ! REASSEMBLE THE DRT TABLE (REMOVE DISCONNECTED VERTICES)
    call mkDRT(SGS) ! IF THIS IS A CAS CALCULATION PROCEED WITH THE UNRESTRICTED DRT TABLE
  end if

  nullify(SGS%DOWNP,SGS%DRTP)

  call MKDAW(SGS) ! COMPUTE DOWNCHAIN TABLE AND ARC WEIGHT
  call MKRAW(SGS) ! COMPUTE UPCHAIN TABLE AND REVERSE ARC WEIGHTS
  call MKLTV(SGS) ! COMPUTE LTV TABLES.
  call MKMID(SGS) ! COMPUTE MIDLEVEL AND LIMITS ON MIDVERTICE.
  call CSFG_MkCOT(SGS,CIS) ! CONSTRUCT THE CASE LIST
contains

  subroutine mknVert0(SGS)

    class(CSFG_Graph), target, intent(inout) :: SGS
    integer(kind=iwp) :: IAC

    SGS%IB0 = SGS%iSpin-1
    SGS%IA0 = (SGS%nActEl-SGS%IB0)/2
    SGS%IC0 = SGS%nLev-SGS%IA0-SGS%IB0

    if ((2*SGS%IA0+SGS%IB0 /= SGS%nActEl) .or. (SGS%IA0 < 0) .or. (SGS%IB0 < 0) .or. (SGS%IC0 < 0)) then
      write(u6,*) 'mknVert0 Error: Impossible specifications.'
      write(u6,'(1x,a,3I8)') 'NACTEL,NLEV,ISPIN:',SGS%nActEl,SGS%nLev,SGS%iSpin
      write(u6,'(1x,a,3I8)') 'IA0,IB0,IC0:      ',SGS%IA0,SGS%IB0,SGS%IC0
      write(u6,*) ' This is a severe internal error, or possibly'
      write(u6,*) ' indicates a strange input which should have been'
      write(u6,*) ' diagnosed earlier. Please submit a bug report.'
      call Abend()
    end if

!   See Eqs. (2) of 10.1002/jcc.26080

    IAC = min(SGS%IA0,SGS%IC0)
    SGS%nVert0 = ((SGS%IA0+1)*(SGS%IC0+1)*(2*SGS%IB0+IAC+2))/2-(IAC*(IAC+1)*(IAC+2))/6

  end subroutine mknVert0

subroutine mkDRT0(SGS)
  class(CSFG_Graph), target, intent(inout) :: SGS

  integer(kind=iwp) :: ADWN, BDWN, CDWN
  integer(kind=iwp) :: AUP, BUP, CUP
  integer(kind=iwp) :: LEV, ISTEP
  integer(kind=iwp) :: VUP, VSTA, VEND, VNEW, VDWN
  integer(kind=iwp) :: key, maxKey

  integer(kind=iwp), parameter :: Steps(4,0:3)=reshape( &
       [0,0,1,0, 0,1,0,1, 1,-1,1,1, 1,0,0,2],[4,4])

  ! --- bounds for encoding ---
  integer(kind=iwp) :: Amax, Bmax, Cmax
  integer(kind=iwp), allocatable :: vmap(:)

  SGS%DRTP(:,:) = 0
  SGS%DownP(:,:) = 0

  SGS%NACTEL = 2*SGS%IA0 + SGS%IB0
  SGS%nLev   = SGS%IA0 + SGS%IB0 + SGS%IC0

  ! maximum ranges (safe upper bounds)
  Amax = SGS%IA0
  Bmax = SGS%IB0 + SGS%IA0   ! loose but safe
  Cmax = SGS%IC0

  maxKey = (Amax+1)*(Bmax+1)*(Cmax+1)

  call mma_allocate(vmap,[0,maxKey-1],Label='vmap')
  vmap(:) = 0

  ! ------------------------------------------
  ! helper: inline key encoding
  ! ------------------------------------------
#   define _KEY_(A,B,C) ((A)*(Bmax+1)*(Cmax+1)+(B)*(Cmax+1)+(C))

  ! TOP vertex
  SGS%DRTP(1,CSFG_LTab) = SGS%nLev
  SGS%DRTP(1,CSFG_NTab) = SGS%NACTEL
  SGS%DRTP(1,CSFG_ATab) = SGS%IA0
  SGS%DRTP(1,CSFG_BTab) = SGS%IB0
  SGS%DRTP(1,CSFG_CTab) = SGS%IC0

  key = _KEY_(SGS%IA0,SGS%IB0,SGS%IC0)
  vmap(key) = 1

  VSTA = 1
  VEND = 1

  ! ------------------------------------------
  ! main loop
  ! ------------------------------------------
  do LEV = SGS%nLev, 1, -1
    VNEW = 0

    do VUP = VSTA, VEND
      AUP = SGS%DRTP(VUP,CSFG_ATab)
      BUP = SGS%DRTP(VUP,CSFG_BTab)
      CUP = SGS%DRTP(VUP,CSFG_CTab)

      do ISTEP = 0, 3

        ADWN = AUP - Steps(1,ISTEP)
        if (ADWN < 0) cycle

        BDWN = BUP - Steps(2,ISTEP)
        if (BDWN < 0) cycle

        CDWN = CUP - Steps(3,ISTEP)
        if (CDWN < 0) cycle

        key = _KEY_(ADWN,BDWN,CDWN)

        VDWN = vmap(key)

        if (VDWN /= 0) then
          SGS%DownP(VUP,ISTEP) = VDWN
        else
          VNEW = VNEW + 1
          VDWN = VEND + VNEW

          SGS%DownP(VUP,ISTEP) = VDWN

          SGS%DRTP(VDWN,CSFG_ATab) = ADWN
          SGS%DRTP(VDWN,CSFG_BTab) = BDWN
          SGS%DRTP(VDWN,CSFG_CTab) = CDWN

          vmap(key) = VDWN
        end if

      end do
    end do

    VSTA = VEND + 1
    VEND = VEND + VNEW

  end do

  ! ------------------------------------------
  ! finalize CSFG_LTab / CSFG_NTab
  ! ------------------------------------------
  SGS%DRTP(1:VEND,CSFG_LTab) = SGS%DRTP(1:VEND,CSFG_ATab) + &
                          SGS%DRTP(1:VEND,CSFG_BTab) + &
                          SGS%DRTP(1:VEND,CSFG_CTab)

  SGS%DRTP(1:VEND,CSFG_NTab) = 2*SGS%DRTP(1:VEND,CSFG_ATab) + &
                            SGS%DRTP(1:VEND,CSFG_BTab)

    call mma_deallocate(vmap)

#   undef _KEY_

end subroutine mkDRT0

  subroutine mkDRT(SGS)
  ! PURPOSE: USING THE UNRESTRICTED DRT TABLE GENERATED BY DRT0 AND
  !          THE MASKING ARRAY PRODUCED BY RESTR COPY ALL VALID
  !          VERTICES FROM THE OLD TO THE NEW DRT TABLE

    class(CSFG_Graph), target, intent(inout) :: SGS
    integer(kind=iwp) :: IC, ID, IDNEW, IV, IVNEW

    call mma_allocate(SGS%DRT,SGS%nVert,5,Label='DRT')
    call mma_allocate(SGS%Down,[1,SGS%nVert],[0,3],Label='SGS%DOWN')

    ! Loop over all original vertices
    do IV=1,SGS%nVert0
      IVNEW = SGS%Ver(IV) ! Pick up the new index, if iVNEW==0
      if (IVNEW == 0) cycle
      ! Move over all relevant elements from DRT0
      SGS%DRT(IVNEW,:) = SGS%DRT0(IV,1:5)
      ! Loop over step vector
      do IC=0,3
        ID = SGS%Down0(IV,IC) ! Pick up the old index of the vertex which connects to vertex IV
        IDNEW = 0
        ! If a viable arc get the new index of the vertex ID
        if (ID /= 0) IDNEW = SGS%Ver(ID)
        ! Update the downward Chain
        SGS%Down(IVNEW,IC) = IDNEW
      end do
    end do
#   ifdef _DEBUGPRINT_
    do IV=1,SGS%nVert
      write(u6,*) 'DRT(i,:)=',SGS%DRT(IV,:)
    end do
#   endif

    call mma_deallocate(SGS%Ver)
    call mma_deallocate(SGS%DRT0)
    call mma_deallocate(SGS%Down0)

  end subroutine mkDRT

  subroutine MKDAW(SGS)
  ! PURPOSE: CONSTRUCT DIRECT ARC WEIGHTS TABLE

    class(CSFG_Graph), target, intent(inout) :: SGS
    integer(kind=iwp) :: IC, IDWN, IV

    call mma_allocate(SGS%DAW,[0,SGS%nVert],[0,4],Label='SGS%DAW')

    ! Note that this is done without symmetry.
    ! BEGIN TO CONSTRUCT DOWN CHAIN TABLE

    SGS%DAW(:,:) = 0
    SGS%DAW(SGS%nVert,4) = 1

    ! Loop over all vertices
    do IV=SGS%nVert-1,1,-1
      do IC=0,3
        IDWN = SGS%Down(IV,IC)
        SGS%DAW(IV,IC+1) = SGS%DAW(IV,IC) + SGS%DAW(IDWN,4)
      end do
    end do

#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) ' DIRECT ARC WEIGHTS:'
    do IV=1,SGS%nVert
      write(u6,'(1X,I4,5X,5(1X,I6))') IV,SGS%DAW(IV,0:4)
    end do
    write(u6,*)
#   endif

  end subroutine MKDAW

  subroutine MKRAW(SGS)
  ! PURPOSE: CONSTRUCT UPCHAIN INDEX TABLE AND REVERSE ARC WEIGHTS

    class(CSFG_Graph), target, intent(inout) :: SGS
    integer(kind=iwp) :: IC, IDWN, IU, IV

    call mma_allocate(SGS%UP,[1,SGS%nVert],[0,3],Label='SGS%UP')
    call mma_allocate(SGS%RAW,[0,SGS%nVert],[0,4],Label='SGS%RAW')

    ! BEGIN BY CONSTRUCTING THE UPCHAIN TABLE IUP:

    SGS%Up(:,:) = 0
    do IU=1,SGS%nVert-1
      do IC=0,3
        IDWN = SGS%Down(IU,IC)
        if (IDWN == 0) cycle
        SGS%Up(IDWN,IC) = IU
      end do
    end do

#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) ' THE UPCHAIN TABLE IN MKRAW:'
    do IV=1,SGS%nVert
      write(u6,'(1X,I4,5X,4(1X,I6))') IV,SGS%Up(IV,0:3)
    end do
    write(u6,*)
#   endif

    ! USE UPCHAIN TABLE TO CALCULATE THE REVERSE ARC WEIGHT TABLE:

    SGS%RAW(:,:) = 0
    SGS%RAW(1,4) = 1
    do IV=2,SGS%nVert
      do IC=0,3
        IU = SGS%Up(IV,IC)
        SGS%RAW(IV,IC+1) = SGS%RAW(IV,IC) + SGS%RAW(IU,4)
      end do
    end do

#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) ' THE REVERSE ARC WEIGHT TABLE IN MKRAW:'
    do IV=1,SGS%nVert
      write(u6,'(1X,I4,5X,5(1X,I6))') IV,SGS%RAW(IV,0:4)
    end do
    write(u6,*)
#   endif

  end subroutine MKRAW

  subroutine MKLTV(SGS)

    class(CSFG_Graph), target, intent(inout) :: SGS
    integer(kind=iwp) :: IV, LEV

    !Note that the range is [-1,nLev]! An imaginary level, -1, is added.
    call mma_allocate(SGS%LTV,[-1,SGS%nLev],Label='LTV')

    ! SET UP A LEVEL-TO-VERTEX TABLE, LTV

    SGS%LTV(:) = 0

    ! Loop over all vertices
    do IV=1,SGS%nVert
      LEV = SGS%DRT(IV,CSFG_LTab) ! Pick up the level index of the vertex
      SGS%LTV(LEV) = SGS%LTV(LEV)+1 ! Increment the number of vertex for a particular level
    end do

    ! Loop over all levels, from the level corresponding to the head vertex to the level corresponding
    ! to the root vertex of the graph
    do LEV=SGS%nLev,0,-1
      SGS%LTV(LEV-1) = SGS%LTV(LEV-1)+SGS%LTV(LEV) ! Update such that  SGS%LEV(Lev) is the number of vertex at level
                                                   ! Lev and above.
    end do

    do LEV=-1,SGS%nLev-1
      SGS%LTV(LEV) = 1+SGS%LTV(LEV+1)
    end do

    ! SGS%LTV(LEV) now tabulates the index of the first vertex at level LEV if the vertices are numbered 1,..,nVert
    ! starting from the head vertex.

  end subroutine MKLTV

  subroutine MKMID(SGS)
  ! PURPOSE: FIND THE MIDLEVEL

    class(CSFG_Graph), target, intent(inout) :: SGS
    integer(kind=iwp) :: IL, IV, MINW, MV, NW

    ! USE DAW,RAW TABLES TO DETERMINE MIDLEV.
    ! THE ASSUMPTION IS THAT A BALANCED NUMBER OF UPPER/LOWER WALKS
    ! IS THE BEST CHOICE.

    !hrl 980529 fix for nLev=0 (no orbitals in any active space)
    ! Since LTV(-1:nLev)  and the statement after the loop
    ! MVSta=LTV(MidLev) we have the condition MidLev>=nLev
    ! Hence MidLev=1 is inappropriate for nLev=0
    ! MIDLEV=1

    SGS%MidLev = min(SGS%nLev,1)

    MINW = 1000000 ! Initiate to a large number
    do IL=1,SGS%nLev-1 ! Loop from level 1 (one above the root vertex) to the
                       ! level below the head vertex.

      ! SGS%RAW(IV,4) is the number of possible paths from the root vertex up to vertex IV
      ! SGS%DAW(IV,4) is the number of possible paths from the head vertex down to vertex IV
      NW = 0
      do IV=SGS%LTV(IL),SGS%LTV(IL-1)-1  ! Loop over all vertex at this level
        NW = NW+abs(SGS%RAW(IV,4)-SGS%DAW(IV,4))
      end do
      NW = abs(NW)
      if (NW >= MINW) cycle
      SGS%MidLev = IL        ! Update level index indicating the midlevel
      MINW = NW              ! Update MINW
    end do

    SGS%MVSta = SGS%LTV(SGS%MidLev)      ! Staring vertex index at level MIDLEV
    SGS%MVEnd = SGS%LTV(SGS%MidLev-1)-1  ! Ending vertex index at level MIDLEV
    CIS%nMidV = SGS%MVEnd-SGS%MVSta+1    ! Number of vertices at level MIDLEV

    ! NOW FIND THE MAX NUMBERS OF UPPER AND LOWER WALKS. RESPECTIVELY
    ! (DISREGARDING SYMMETRY)

    SGS%MxUp = 0
    SGS%MxDwn = 0
    do MV=SGS%MVSta,SGS%MVEnd ! Lopp over all vertices of level MIDLEV
      SGS%MxUp  = max(SGS%MxUp,SGS%RAW(MV,4))
      SGS%MxDwn = max(SGS%MxDwn,SGS%DAW(MV,4))
    end do

#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,'(A,I3)') ' MIDLEVEL =             ',SGS%MidLev
    write(u6,'(A,I3)') ' NUMBER OF MIDVERTICES =',CIS%NMIDV
    write(u6,'(A,I3)') ' FIRST MIDVERTEX =      ',SGS%MVSta
    write(u6,'(A,I3)') ' LAST MIDVERTEX =       ',SGS%MVEnd
    write(u6,'(A,I3)') ' MAX. NO UPPER WALKS=   ',SGS%MxUp
    write(u6,'(A,I3)') ' MAX. NO LOWER WALKS=   ',SGS%MxDwn
    write(u6,*)
#   endif

  end subroutine MKMID

subroutine RMVERT(SGS)

  class(CSFG_Graph), intent(inout) :: SGS

  integer(kind=iwp) :: IV, IC, ID
  integer(kind=iwp) :: L, N
  integer(kind=iwp) :: NCHANGES, NV, iRO, iSy, Lev
  integer(kind=iwp), allocatable :: CONN(:), Lim(:), LimMx(:)

  !-----------------------------------------
  ! Build occupation limits
  !-----------------------------------------
  call mma_allocate(Lim,SGS%nLev,Label='Lim')
  call mma_allocate(LimMx,SGS%nLev,Label='LimMx')
  Lim(:) = 0
  LimMx(:) = huge(1_iwp)

  Lev = 0
  do iRO = 1, SGS%nRsPrt
    do iSy = 1, SGS%nSym
      Lev = Lev + SGS%nRas(iSy,iRO)
    end do
    if (Lev > 0) Lim(Lev) = SGS%nRasEl(iRO)
    if ((Lev > 0) .and. (SGS%nRasMx(iRO) >= 0)) LimMx(Lev) = SGS%nRasMx(iRO)
  end do

  !-----------------------------------------
  ! Initialize vertex mask
  !-----------------------------------------
  call mma_allocate(SGS%Ver,SGS%nVert0,Label='SGS%Ver')
  SGS%Ver(:) = 1

  call mma_allocate(CONN,SGS%nVert,Label='CONN')

  !-----------------------------------------
  ! Initial pruning based on occupation
  !-----------------------------------------
  do IV = 1, SGS%nVert-1
    L = SGS%DRT0(IV,CSFG_LTab)
    N = SGS%DRT0(IV,CSFG_NTab)

    if ((N < Lim(L)) .or. (N > LimMx(L))) SGS%Ver(IV) = 0
  end do

  !-----------------------------------------
  ! Iterative pruning
  !-----------------------------------------
  do
    NCHANGES = 0

    !-------------------------------------
    ! Step 1: remove invalid downward arcs
    !-------------------------------------
    do IV = 1, SGS%nVert-1
      if (SGS%Ver(IV) == 0) then
        SGS%Down0(IV,0:3) = 0
        cycle
      end if

      do IC = 0,3
          ID = SGS%Down0(IV,IC)
          if (ID < 1) cycle

        if (SGS%Ver(ID) == 0) then
           SGS%Down0(IV,IC) = 0
           NCHANGES = NCHANGES + 1
        end if
      end do
    end do

    !-------------------------------------
    ! Step 2: remove vertices with no children
    !-------------------------------------
    do IV = 1, SGS%nVert-1
      if (SGS%Ver(IV) == 0) cycle

      if (all(SGS%Down0(IV,0:3) == 0)) then
        SGS%Ver(IV) = 0
        NCHANGES = NCHANGES + 1
      end if
    end do

    !-------------------------------------
    ! Step 3: rebuild upward connectivity
    !-------------------------------------
    CONN(:) = 0
    CONN(1) = SGS%Ver(1)

    do IV = 1, SGS%nVert-1
      if (SGS%Ver(IV) == 0) cycle

      do IC = 0,3
          ID = SGS%Down0(IV,IC)
          if (ID < 1) cycle
          if (SGS%Ver(ID) == 1) CONN(ID) = 1
      end do
    end do

    ! remove vertices not connected from above
    do IV = 1, SGS%nVert
        if ((SGS%Ver(IV) == 1) .and. (CONN(IV) == 0)) then
        SGS%Ver(IV) = 0
        NCHANGES = NCHANGES + 1
      end if
    end do

    if (NCHANGES == 0) exit

  end do

  !-----------------------------------------
  ! Safety check
  !-----------------------------------------
  if (SGS%Ver(1) == 0) then
    write(u6,*) 'RASSI/RMVERT: Too severe restrictions.'
    call Abend()
  end if

  !-----------------------------------------
  ! Renumber vertices
  !-----------------------------------------
  NV = 0
  do IV = 1, SGS%nVert
    if (SGS%Ver(IV) == 1) then
      NV = NV + 1
      SGS%Ver(IV) = NV
    end if
  end do

  SGS%nVert = NV

  call mma_deallocate(CONN)
  call mma_deallocate(Lim)
  call mma_deallocate(LimMx)

end subroutine RMVERT

end subroutine CSFG_Make

subroutine CSFG_MkCOT(SGS,CIS)
! For INIT==0
! PURPOSE: SET UP COUNTER AND OFFSET TABLES FOR WALKS AND CSFS
! NOTE:    TO GET GET VARIOUS COUNTER AND OFFSET TABLES
!          THE DOWN-CHAIN TABLE IS SCANNED TO PRODUCE ALL POSSIBLE
!          WALKS. POSSIBLY, THERE ARE MORE EFFICIENT WAYS, BUT
!          SINCE ONLY UPPER AND LOWER WALKS ARE REQUIRED
!          THEIR NUMBER IS VERY LIMITTED, EVEN FOR LARGE CASES.

! For INIT==1
! PURPOSE: CONSTRUCT THE COMPRESSED CASE-LIST, I.E.,
!          STORE THE STEP VECTOR FOR ALL POSSIBLE WALKS
!          IN THE ARRAY ICASE. GROUPS OF CSFG_nPack CASES ARE PACKED
!          INTO ONE INTEGER WORD.

class(CSFG_Graph), intent(inout) :: SGS
class(CSFG_Walks), intent(inout) :: CIS
  integer(kind=iwp) :: IC, IHALF, ILND, INIT, IPOS, ISML, ISTP, IVB, IVT, IVTEND, IVTOP, IVTSTA, IWSYM, L, LEV, LEV1, LEV2, LL, MV
# ifdef _DEBUGPRINT_
  integer(kind=iwp) :: IS, NUW
# endif
integer(kind=iwp), allocatable :: SCR(:,:)
integer(kind=iwp), parameter :: IVERT = 1, ISYM = 2, ISTEP = 3

do INIT=0,1
  if (INIT == 0) then
    ! Compute the number of integers needed to store the packed case vectors.
    CIS%nIpWlk = 1+(SGS%MidLev-1)/CSFG_nPack
    CIS%nIpWlk = max(CIS%nIpWlk,1+(SGS%nLev-SGS%MidLev-1)/CSFG_nPack)
    call mma_allocate(CIS%NOW,2,SGS%nSym,CIS%nMidV,Label='CIS%NOW')
    call mma_allocate(CIS%IOW,2,SGS%nSym,CIS%nMidV,Label='CIS%IOW')
    call mma_allocate(CIS%NOCSF,SGS%nSym,CIS%nMidV,SGS%nSym,Label='CIS%NOCSF')
    call mma_allocate(CIS%IOCSF,SGS%nSym,CIS%nMidV,SGS%nSym,Label='CIS%IOCSF')
    call mma_allocate(CIS%NCSF,SGS%nSym,Label='CIS%NCSF')

    call mma_allocate(Scr,[1,3],[0,SGS%nLev],Label='Scr')  ! First index referenced by IVERT, SYM, and ISTEP

    ! CLEAR ARRAYS IOW AND NOW

    CIS%NOW(:,:,:) = 0
    CIS%IOW(:,:,:) = 0

    ! CLEAR ARRAYS IOCSF AND NOCSF

    CIS%IOCSF(:,:,:) = 0
    CIS%NOCSF(:,:,:) = 0
  else
    !call mma_allocate(Scr,[1,3],[0,SGS%nLev],Label='Scr',safe='*')  ! First index referenced by IVERT, SYM, and ISTEP
    call mma_allocate(CIS%ICase,CIS%nWalk*CIS%nIpWlk,Label='CIS%ICase',safe='*')
    ! CLEAR ARRAY NOW. IT WILL BE RESTORED FINALLY
    CIS%NOW(:,:,:) = 0
  end if

! START MAIN LOOP OVER UPPER AND LOWER WALKS, RESPECTIVELY.

do IHALF=1,2
  ! set the loop ranges:
  ! IVSTA-IVTEND: the top vertex, or loop over midlevel vertices
  ! LEV2-LEV1: loop over levels
  if (IHALF == 1) then
    IVTSTA = 1
    IVTEND = 1
    LEV1 = SGS%nLev
    LEV2 = SGS%MidLev
  else
    IVTSTA = SGS%MVSta
    IVTEND = SGS%MVEnd
    LEV1 = SGS%MidLev
    LEV2 = 0
  end if

  ! LOOP OVER VERTICES STARTING AT TOP OF SUBGRAPH
  ! This is either the top vertex, or one of the midlevel vertices

  do IVTOP=IVTSTA,IVTEND
    ! SET CURRENT LEVEL=TOP LEVEL OF SUBGRAPH
    LEV = LEV1

    ! Store away vertex index, and initiate symmetry index, and the step vector index
    Scr(IVERT,LEV) = IVTOP
    Scr(ISYM,LEV)  =  1
    Scr(ISTEP,LEV) = -1

    do while (LEV <= LEV1)

      ! FIND FIRST POSSIBLE UNTRIED ARC DOWN FROM CURRENT VERTEX
      IVT = Scr(IVERT,LEV)  ! Pickup the current vertex index for level LEV
      ! Continue scanning the step vectors, SCR(ISTEP,LEV) is the index of the next vector to explot
      do ISTP=Scr(ISTEP,LEV)+1,3
        IVB = SGS%Down(IVT,ISTP)
        if (IVB /= 0) exit      ! Exits if step vector leads to a valid vertex below.
      end do

      ! IF NO SUCH ARC WAS POSSIBLE. GO UP ONE LEVEL AND TRY AGAIN.
      if (ISTP > 3) then
        Scr(ISTEP,LEV) = -1
        LEV = LEV+1
        cycle
      end if

      ! SUCH AN ARC WAS FOUND. WALK DOWN:
      Scr(ISTEP,LEV) = ISTP ! Store the current step vector index of level Level
      ! doubly occupied or empty orbital case are total symmetric. Singly occupied orbitals
      ! carry the symmetry of the orbital in the level.
      ISML = 1
          if (ISTP == 1 .or. ISTP == 2) ISML = SGS%ISm(LEV)

      LEV = LEV-1     ! Walk down one level

      ! Store away the accumulated symmetry, the new vertex, and initiate the step vector counter.
      Scr(ISYM,LEV) = Mul(ISML,Scr(ISYM,LEV+1))
      Scr(IVERT,LEV) = IVB
      Scr(ISTEP,LEV) = -1

      if (LEV > LEV2) cycle   ! Repeat

      ! WE HAVE NOW REACHED THE BOTTOM LEVEL. THE WALK IS COMPLETE.
      ! FIND MIDVERTEX NUMBER ORDERING NUMBER AND SYMMETRY OF THIS WALK
      MV = Scr(IVERT,SGS%MidLev)+1-SGS%MVSta ! Pick up the relative index of the midlev vertex
      IWSYM = Scr(ISYM,LEV2)                 ! Pick up the symmetry of the walk

      ILND = CIS%NOW(IHALF,IWSYM,MV) + 1
      CIS%NOW(IHALF,IWSYM,MV) = ILND   ! Increment counter for how many walks there are given the
      ! midvertex index and the total symmetry of the subwalk.
      if (INIT == 1) then
        ! CONSEQUENTLY, THE POSITION IMMEDIATELY BEFORE THIS COMPRESSED WALK:
        IPOS = CIS%IOW(IHALF,IWSYM,MV)+(ILND-1)*CIS%nIpWlk

        ! PACK THE STEPS IN GROUPS OF CSFG_nPack LEVELS PER INTEGER:
        do LL=LEV2+1,LEV1,CSFG_nPack

          IC = 0
          do L=min(LL+CSFG_nPack-1,LEV1),LL,-1
            IC = 4*IC+Scr(ISTEP,L)
          end do

          IPOS = IPOS+1
          CIS%ICase(IPOS) = IC
        end do
      end if

      ! BACK UP ONE LEVEL AND EXPLORE NEW WALKS:
      LEV = LEV+1
    end do
  end do
end do

if (INIT == 0) then
! Generate an off-set array for the CIS%NOW array
call CSFG_Mk_IOW(CIS,SGS)
call CSFG_Count(CIS,SGS)

#ifdef _DEBUGPRINT_
NUW=CIS%NUW
write(u6,*)
write(u6,*) ' TOTAL NR OF WALKS: UPPER ',NUW
write(u6,*) '                    LOWER ',CIS%nWalk-NUW
write(u6,*) '                     SUM  ',CIS%nWalk
write(u6,*)
write(u6,*) ' NR OF CONFIGURATIONS/SYMM:'
write(u6,'(8(1X,I8))') (CIS%NCSF(IS),IS=1,SGS%nSym)
write(u6,*)
write(u6,*)
write(u6,*) ' NR OF WALKS AND CONFIGURATIONS IN MkNRCOUP'
write(u6,*) ' BY MIDVERTEX AND SYMMETRY.'
do MV=1,CIS%nMidV
  write(u6,'(A,I2,A,8I6)') '  MV=',MV,'    UPPER WALKS:',(CIS%NOW(1,IS,MV),IS=1,SGS%nSym)
  write(u6,'(A,8I6)') '           LOWER WALKS:',(CIS%NOW(2,IS,MV),IS=1,SGS%nSym)
  do ISTP=1,SGS%nSym
    write(u6,'(A,I2,A,8I6)') ' ISTP=',ISTP,'  CONFIGURATIONS:',(CIS%NOCSF(IS,MV,ISTP),IS=1,SGS%nSym)
  end do
end do
#endif
end if

end do ! INIT

call mma_deallocate(Scr,safe='*')

end subroutine CSFG_MkCOT

subroutine CSFG_Count(CIS,SGS)
! CONSTRUCT COUNTER AND OFFSET TABLES FOR THE CSFS
! SEPARATED BY MIDVERTICES AND SYMMETRY.
! FORM ALSO CONTRACTED SUMS OVER MIDVERTICES.

class(CSFG_Walks), intent(inout) :: CIS
class(CSFG_Graph), intent(inout) :: SGS
integer(kind=iwp) :: ISYDWN, ISYTOT, ISYUP, MV, N

CIS%NCSF(:) = 0   ! Number of CSFs in each irrep
do ISYTOT=1,SGS%NSYM
  do MV=1,CIS%nMidV
    do ISYUP=1,SGS%NSYM
      ISYDWN = Mul(ISYTOT,ISYUP)
      N = CIS%NOW(1,ISYUP,MV)*CIS%NOW(2,ISYDWN,MV)
      CIS%NOCSF(ISYUP,MV,ISYTOT) = N
      CIS%IOCSF(ISYUP,MV,ISYTOT) = CIS%NCSF(ISYTOT)
      CIS%NCSF(ISYTOT) = CIS%NCSF(ISYTOT)+N
#     include "compiler_features.h"
#     ifdef _BUGGY_INTEL_LLVM_
      ! dummy statement to work around compiler bug, will never be executed
      if (ISYUP > 99) CIS%NCSF(ISYTOT) = -1
#     endif
    end do
  end do
end do

end subroutine CSFG_Count

subroutine CSFG_Mk_IOW(CIS,SGS)
! CONSTRUCT OFFSET TABLES FOR UPPER AND LOWER WALKS
! SEPARATED FOR EACH MIDVERTEX AND SYMMETRY
!
! CIS%IOW(1/2,ISYM,MV) is an off-set vector associated with CIS%NOW(1/2,ISYM,MV)
! CIS%NOW(1/2,ISYM,MV) is the number of upper/lower walks of symmetry ISYM that ends/starts in mid vertex MV (relative indexation).
! 1/2 indicate if the off-set is for upper or lower walks.

  class(CSFG_Walks), intent(inout) :: CIS
  class(CSFG_Graph), intent(inout) :: SGS
  integer(kind=iwp) :: ISYM, MV

CIS%NUW = 0   ! At completion the total number of upper walks.
do MV=1,CIS%nMidV
  do ISYM=1,SGS%NSYM
    CIS%IOW(1,ISYM,MV) = CIS%NUW*CIS%nIpWlk
    CIS%NUW = CIS%NUW+CIS%NOW(1,ISYM,MV)
  end do
end do
CIS%nWalk = CIS%NUW ! At completion the total number pf upper and lower walks

do MV=1,CIS%nMidV
  do ISYM=1,SGS%NSYM
    CIS%IOW(2,ISYM,MV) = CIS%nWalk*CIS%nIpWlk
    CIS%nWalk = CIS%nWalk+CIS%NOW(2,ISYM,MV)
  end do
end do
end subroutine CSFG_Mk_IOW

!> The next CSF (walk) of symmetry LSM in the split-graph order, as its step
!> vector ICS(1:nLev) (0 empty, 1 up, 2 down, 3 doubly occupied). The CSFs
!> come in blocks (MV, ISYUP), MV outermost (CSFG_Count), each block with the
!> upper walks running fastest (as GETSTEPVECTOR for symmetry 1). Start with
!> MV = ISYUP = IDWN = IUP = 1; MV = 0 after the last walk.
subroutine CSFG_Walk_Next(G,W,LSM,MV,ISYUP,IDWN,IUP,ICS)
  class(CSFG_Graph), intent(in) :: G
  class(CSFG_Walks), intent(in) :: W
  integer(kind=iwp), intent(in) :: LSM
  integer(kind=iwp), intent(inout) :: MV, ISYUP, IDWN, IUP
  integer(kind=iwp), intent(out) :: ICS(G%nLev)
  integer(kind=iwp) :: IC1, ICDPOS, ICDWN, ICUP, ICUPOS, IDW0, ISYDWN, IUW0, LEV, NDWN, NNN, NUP

  call Walk_Skip_Empty(G%nSym,W%nMidV,W%NOW,LSM,MV,ISYUP)
  if (MV == 0) then
    write(u6,'(1X,A)') 'CSFG_Walk_Next: no more walks'
    call Abend()
  end if
  ISYDWN = Mul(LSM,ISYUP)
  NUP = W%NOW(1,ISYUP,MV)
  NDWN = W%NOW(2,ISYDWN,MV)
  IUW0 = 1-W%nIpWlk+W%IOW(1,ISYUP,MV)
  IDW0 = 1-W%nIpWlk+W%IOW(2,ISYDWN,MV)
  ! lower walk
  ICDPOS = IDW0+IDWN*W%nIpWlk
  ICDWN = W%ICASE(ICDPOS)
  NNN = 0
  do LEV=1,G%MIDLEV
    NNN = NNN+1
    if (NNN == CSFG_nPack+1) then
      NNN = 1
      ICDPOS = ICDPOS+1
      ICDWN = W%ICASE(ICDPOS)
    end if
    IC1 = ICDWN/4
    ICS(LEV) = ICDWN-4*IC1
    ICDWN = IC1
  end do
  ! upper walk
  ICUPOS = IUW0+W%nIpWlk*IUP
  ICUP = W%ICASE(ICUPOS)
  NNN = 0
  do LEV=G%MIDLEV+1,G%nLev
    NNN = NNN+1
    if (NNN == CSFG_nPack+1) then
      NNN = 1
      ICUPOS = ICUPOS+1
      ICUP = W%ICASE(ICUPOS)
    end if
    IC1 = ICUP/4
    ICS(LEV) = ICUP-4*IC1
    ICUP = IC1
  end do
  ! next indices: upper walks fastest, then lower walks, then the blocks
  if (IUP == NUP) then
    IUP = 1
    if (IDWN == NDWN) then
      IDWN = 1
      call Walk_Next_Block(G%nSym,W%nMidV,MV,ISYUP)
    else
      IDWN = IDWN+1
    end if
  else
    IUP = IUP+1
  end if
  call Walk_Skip_Empty(G%nSym,W%nMidV,W%NOW,LSM,MV,ISYUP)

end subroutine CSFG_Walk_Next

!> The next (MV, ISYUP) block of CSFG_Walk_Next: ISYUP fastest; MV = 0 after
!> the last one. (Module procedures with explicit arguments: internal
!> procedures that changed the host's MV/ISYUP made nvfortran loop forever.)
subroutine Walk_Next_Block(nSym,nMidV,MV,ISYUP)
  integer(kind=iwp), intent(in) :: nSym, nMidV
  integer(kind=iwp), intent(inout) :: MV, ISYUP
  if (ISYUP == nSym) then
    ISYUP = 1
    if (MV == nMidV) then
      MV = 0
    else
      MV = MV+1
    end if
  else
    ISYUP = ISYUP+1
  end if
end subroutine Walk_Next_Block

!> Advance (MV, ISYUP) to the next block with CSFs of symmetry LSM.
subroutine Walk_Skip_Empty(nSym,nMidV,NOW,LSM,MV,ISYUP)
  integer(kind=iwp), intent(in) :: nSym, nMidV, NOW(:,:,:), LSM
  integer(kind=iwp), intent(inout) :: MV, ISYUP
  do while (MV /= 0)
    if (NOW(1,ISYUP,MV)*NOW(2,Mul(LSM,ISYUP),MV) > 0) exit
    call Walk_Next_Block(nSym,nMidV,MV,ISYUP)
  end do
end subroutine Walk_Skip_Empty

!> Free the arrays of a graph and its walks built by CSFG_Make.
subroutine CSFG_Free(G,W)
  class(CSFG_Graph), intent(inout) :: G
  class(CSFG_Walks), intent(inout) :: W
  call mma_deallocate(G%ISm,safe='*')
  call mma_deallocate(G%DRT,safe='*')
  call mma_deallocate(G%DRT0,safe='*')
  call mma_deallocate(G%Down,safe='*')
  call mma_deallocate(G%Down0,safe='*')
  call mma_deallocate(G%Up,safe='*')
  call mma_deallocate(G%Ver,safe='*')
  call mma_deallocate(G%MAW,safe='*')
  call mma_deallocate(G%LTV,safe='*')
  call mma_deallocate(G%DAW,safe='*')
  call mma_deallocate(G%RAW,safe='*')
  call mma_deallocate(W%NOW,safe='*')
  call mma_deallocate(W%IOW,safe='*')
  call mma_deallocate(W%NCSF,safe='*')
  call mma_deallocate(W%NOCSF,safe='*')
  call mma_deallocate(W%IOCSF,safe='*')
  call mma_deallocate(W%ICase,safe='*')
end subroutine CSFG_Free

end module CSF_Graph
