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
! Copyright (C) Valera Veryazov                                        *
!               Luca De Vico                                           *
!***********************************************************************

subroutine CHARGE(NSYM,NBAS,BNAME,CMO,OCCN,SMAT,iCase,FullMlk,lSave)

use SpinAV, only: Do_SpinAV, DSc
use UnixInfo, only: ProgName
use define_af, only: AngTp, iTabMx
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: NSYM, NBAS(NSYM), iCase
character(len=LenIn8), intent(in) :: BNAME(*)
real(kind=wp), intent(in) :: CMO(*), OCCN(*), SMAT(*)
logical(kind=iwp), intent(in) :: FullMlk, lSave
integer(kind=iwp) :: AtomA, AtomB, i0, iAB, iAng, iB, iBlo, i, iEnd, iix, iixx, ik, ikk, iM, iMN, IMO, IO, iPair, iPL, IS, ISMO, &
                     IST, iStart, iSum, iSwap, iSyLbl, ISYM, IT, ix, J, jAng, jEnd, jM, jx, k, l, lqSwap, MY, MYNUC, MYTYP, NB, &
                     nBas2, NBAST, nNuc, NPBonds, nScr, NXTYP, NY, NYNUC, NYTYP, tNUC
real(kind=wp) :: BO, BOThrs, DET, DMN, TCh, TERM, xsg
#ifdef _DEBUGPRINT_
real(kind=wp) :: E
#endif
logical(kind=iwp) :: DMN_SpinAV, DoBond
character(len=len(ProgName)) :: PName
character(len=8) :: TMP
integer(kind=iwp), allocatable :: center(:), ICNT(:), ITYP(:), nStab(:)
real(kind=wp), allocatable, save :: Bonds(:), Chrg(:), D(:,:), D_blo(:), D_tmp(:,:), DS(:,:), DSswap(:,:), Fac(:), P(:,:), &
                                    PInv(:,:), Q2(:), QQ(:,:), QSUM(:), QSUM_TOT(:), qSwap(:), S(:,:), S_blo(:), S_tmp(:,:), Scr(:)
character(len=LenIn4), allocatable :: LblCnt4(:)
character(len=LenIn), allocatable :: CNAME(:)
character(len=8), allocatable :: TNAME(:), TSwap(:)
character(len=*), parameter :: AufBau(19) = ['01s',                   &
                                             '02s',            '02p', &
                                             '03s',            '03p', &
                                             '04s',      '03d','04p', &
                                             '05s',      '04d','05p', &
                                             '06s','04f','05d','06p', &
                                             '07s','05f','06d','07p']
integer(kind=iwp), external :: iPrintLevel
real(kind=wp), external :: DDot_
logical(kind=iwp), external :: Reduce_Prt
character(len=LenIn8), external :: Clean_BName
!character(len=4), allocatable :: TLbl(:)
!character(len=LenIn), allocatable :: LblCnt(:)
#include "WrkSpc.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
iPL = iPrintLevel(-1)
if (Reduce_Prt() .and. (iPL < 3)) iPL = 0

if (iCase /= 0) then
  if (iPL >= 2) then
    write(u6,*)
    call CollapseOutput(1,'   Molecular charges:')
    write(u6,'(3X,A)') '   ------------------'
    write(u6,*)
  end if
end if

xsg = -One
DMN_SpinAV = .false.
if ((iCase == 0) .or. (iCase == 1)) then
  if (Do_SpinAV) DMN_SpinAV = .true.
  if (iCase == 1) xsg = One
end if
!                                                                      *
!***********************************************************************
!                                                                      *

!----------------------------------------------------------------------*
! Get the name of the calling module.                                  *
! If CPFMCPF no bond analysis is done.                                 *
!----------------------------------------------------------------------*

PName = ProgName
call Upcase(PName)
PName = adjustl(PName)
iEnd = 1
do
  if (PName(iEnd:iEnd) == ' ') exit
  iEnd = iEnd+1
end do

DoBond = .true.
if (PName(1:iEnd) == 'CPF') DoBond = .false.

!----------------------------------------------------------------------*
! Set the Mulliken Bond Order threshold for printout                   *
!----------------------------------------------------------------------*

BOThrs = Half

!----------------------------------------------------------------------*
! GET THE TOTAL NUMBER OF BASIS FUNCTIONS                              *
!----------------------------------------------------------------------*

NBAST = 0
do I=1,NSYM
  NBAST = NBAST+NBAS(I)
end do

!----------------------------------------------------------------------*
! Find the list of unique center labels                                *
!----------------------------------------------------------------------*

call Get_iScalar('Unique atoms',nNuc)
call mma_allocate(CNAME,nNuc,label='CNAME')
call mma_allocate(nStab,nNuc,label='nStab')
call Get_cArray('Unique Atom Names',CNAME,LenIn*nNuc)
call Get_iArray('nStab',nStab,nNuc)

!----------------------------------------------------------------------*
! Find the center label for each basis function                        *
!----------------------------------------------------------------------*

call mma_allocate(ICNT,NBAST,label='ICNT')
ICNT(:) = -1
do I=1,NBAST
  do J=1,NNUC
    if (BNAME(I)(1:LenIn) == CNAME(J)) then
      ICNT(I) = J
      exit
    end if
  end do
end do

!----------------------------------------------------------------------*
! Find the type label for each basis function                          *
!----------------------------------------------------------------------*

call mma_allocate(ITYP,nBAST,label='ITYP')
ITYP(:) = 0
NXTYP = 0
outer: do I=1,NBAST
  if (ICNT(I) < 0) cycle outer  ! skip pseudo center
  do J=1,I-1
    if (BNAME(I)(LenIn1:LenIn8) == BNAME(J)(LenIn1:LenIn8)) then
      ITYP(I) = ITYP(J)
      cycle outer
    end if
  end do
  NXTYP = NXTYP+1
  ITYP(I) = NXTYP
end do outer
call mma_allocate(TNAME,NXTYP,label='TNAME')
TNAME(:) = ''
do I=1,NBAST
  J = ITYP(I)
  if (J == 0) cycle
  if (TNAME(J) /= '') cycle
  TNAME(J) = BNAME(I)(LenIn1:LenIn8)
end do

lqSwap = NNUC+NNUC*NXTYP

if (iCase == 0) then
  ! instead of printing charges we dump everything into a memory
  ! same with DS matrix

  call mma_allocate(qSwap,lqSwap,label='CHRG_SWP')

  if (DoBond) then
    call mma_allocate(DSswap,NBAST,NBAST,label='DSswap')
  end if

end if

!----------------------------------------------------------------------*
! Do some trivial sorting of the type labels                           *
!----------------------------------------------------------------------*
!
! Sort with respect to radial index

ix = 0
jx = 0
do i=1,NxTyp-1
  ix = ichar(TNAME(i)(1:1))-ichar('1')+1
  ix = 10*ix+ichar(TNAME(i)(2:2))-ichar('1')+1
  ! Put polarization and diffuse functions last
  if (tName(i)(1:1) == '*') ix = 100
  do j=i+1,NxTyp
    jx = ichar(TNAME(j)(1:1))-ichar('1')+1
    jx = 10*jx+ichar(TNAME(j)(2:2))-ichar('1')+1
    if (tName(j)(1:1) == '*') jx = 100
    if (ix > jx) then
      iSwap = ix
      ix = jx
      jx = iSwap
      TMP = TNAME(i)
      TNAME(i) = TNAME(j)
      TNAME(j) = TMP
    end if
  end do
end do

! Sort with respect to angular index

iAng = 0
jAng = 0
ix = 1
do
  iix = ichar(tName(ix)(1:1))
  iixx = ichar(tName(ix)(2:2))
  jx = ix
  do i=min(ix+1,NxTyp),NxTyp
    if ((ichar(tName(i)(1:1)) == iix) .and. (ichar(tName(i)(2:2)) == iixx)) jx = i
  end do

  do i=ix,jx-1
    do k=0,iTabMx
      if (AngTp(k) == tName(i)(3:3)) iAng = k
    end do
    do j=i+1,jx
      do l=0,iTabMx
        if (AngTp(l) == tName(j)(3:3)) jAng = l
      end do
      if (iAng > jAng) then
        iSwap = iAng
        iAng = jAng
        jAng = iSwap
        TMP = TNAME(i)
        TNAME(i) = TNAME(j)
        TNAME(j) = TMP
      end if
    end do
  end do
  !write(u6,*) ' Sorted n subrange'
  !do i=ix,jx
  !  write(u6,*) TName(i)
  !end do

  ! Now sort with respect to the magnetic index

  iEnd = jx
  iStart = ix
  do
    do k=0,iTabMx
      if (AngTp(k) == tName(iStart)(3:3)) iAng = k
    end do
    jEnd = iStart
    do i=min(iStart+1,iEnd),iEnd
      if (tName(i)(3:3) == AngTp(iAng)) jEnd = i
    end do

    i0 = ichar('1')-1

    iM = 0
    jM = 0
    if (iAng == 1) then
      do i=iStart,jEnd-1
        if (tName(i)(4:4) == 'x') iM = 1
        if (tName(i)(4:4) == 'z') iM = 0
        if (tName(i)(4:4) == 'y') iM = -1
        do j=i+1,jEnd
          if (tName(j)(4:4) == 'x') jM = 1
          if (tName(j)(4:4) == 'z') jM = 0
          if (tName(j)(4:4) == 'y') jM = -1
          if (jM > iM) then
            iSwap = iM
            iM = jM
            jM = iSwap
            TMP = TNAME(i)
            TNAME(i) = TNAME(j)
            TNAME(j) = TMP
          end if
        end do
      end do
    else if (iAng >= 2) then
      do i=iStart,jEnd-1
        iM = ichar(tName(i)(4:4))-i0
        iM = 10*iM+ichar(tName(i)(5:5))-i0
        if (tName(i)(6:6) == '-') iM = -iM
        do j=i+1,jEnd
          jM = ichar(tName(j)(4:4))-i0
          jM = 10*jM+ichar(tName(j)(5:5))-i0
          if (tName(j)(6:6) == '-') jM = -jM
          if (jM > iM) then
            iSwap = iM
            iM = jM
            jM = iSwap
            TMP = TNAME(i)
            TNAME(i) = TNAME(j)
            TNAME(j) = TMP
          end if
        end do
      end do
    end if

    if (jEnd == iEnd) exit
    iStart = jEnd+1
  end do

  if (jx == NxTyp) exit
  ix = jx+1
end do

! Sort according to AufBau

call mma_allocate(TSwap,NxTyp,label='TSwap')
iStart = 1
do iAB=1,19
  do i=1,NxTyp
    if (TName(i)(1:3) == AufBau(iAB)) then
      TSwap(iStart) = TName(i)
      TName(i) = '        '
      iStart = iStart+1
    end if
  end do
end do
do i=1,NxTyp
  if (TName(i) /= '        ') then
    TSwap(iStart) = TName(i)
    TName(i) = '        '
    iStart = iStart+1
  end if
end do
TName(:) = TSwap(:)
call mma_deallocate(TSwap)

do I=1,NBAST
  if (ICNT(I) < 0) cycle  ! skip pseudo center
  do J=1,NXTYP
    if (BNAME(I)(LenIn1:LenIn8) == TNAME(J)) then
      ITYP(I) = J
      exit
    end if
  end do
end do

!----------------------------------------------------------------------*
! Get the total number of atoms tNUC, regardless of symmetry           *
!----------------------------------------------------------------------*

call Get_iScalar('LP_nCenter',tNUC)

!----------------------------------------------------------------------*
! Bond analysis initialization                                         *
!----------------------------------------------------------------------*

if (DoBond) then

!----------------------------------------------------------------------*
! In case of symmetry we need the desymmetrization matrix,             *
! for the bond order calculation only.                                 *
!----------------------------------------------------------------------*

  if (nSym > 1) then

    call mma_allocate(P,NBAST,NBAST,label='P')
    call mma_allocate(PInv,NBAST,NBAST,label='PInv')
    call Get_dArray('SM',P,NBAST**2)
#ifdef _DEBUGPRINT_
    call RecPrt('SM',' ',P,NBAST,NBAST)
#endif
    call MINV(P,PInv,DET,NBAST)
#ifdef _DEBUGPRINT_
    call RecPrt('SMInv',' ',PInv,NBAST,NBAST)
#endif
    call DGeTMi(PInv,NBAST,NBAST)
  end if

  ! Pick up index array of which center a basis function belongs to.
  ! If no symmetry, it is the same as ICNT(I).

  call mma_allocate(center,NBAST,label='center')
  call Get_iArray('Center Index',center,NBAST)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !--------------------------------------------------------------------*
  ! Initialize symmetric density D_tmp and overlap S_tmp matrices,     *
  ! block D_blo and S_blo matrices (if symmetry)                       *
  ! plus asymmetric D, S and DS matrices                               *
  !--------------------------------------------------------------------*

  call mma_allocate(D_tmp,NBAST,NBAST,label='D_tmp')
  call mma_allocate(S_tmp,NBAST,NBAST,label='S_tmp')
  call mma_allocate(D,NBAST,NBAST,label='D')
  call mma_allocate(S,NBAST,NBAST,label='S')
  call mma_allocate(DS,NBAST,NBAST,label='DS')
  D_tmp(:,:) = Zero
  S_tmp(:,:) = Zero
  D(:,:) = Zero
  S(:,:) = Zero
  DS(:,:) = Zero

  if (nSym > 1) then
    nBas2 = 0
    do I=1,nsym
      nBas2 = nBas2+nBas(i)*nBas(i)
    end do
    call mma_allocate(D_blo,nBas2,label='D_blo')
    call mma_allocate(S_blo,nBas2,label='S_blo')
    D_blo(:) = Zero
    S_blo(:) = Zero
  end if

  !--------------------------------------------------------------------*
  !   Find the center label for each atom, regardless of symmetry      *
  !--------------------------------------------------------------------*

  ! Just atom label. It's a double of the next one,
  ! but someone could find it usefull in future

  !call Get_Name_All(TLbl)

  ! Atom labels plus symmetry generator

  call mma_allocate(LblCnt4,tNUC,label='LblCnt4')
  call Get_cArray('LP_L',LblCnt4,LenIn4*tNUC)
  !do i=1,tNUC
  !  LblCnt(i)(1:LenIn) = LblCnt4(i)(1:LenIn)
  !end do

  !--------------------------------------------------------------------*
  ! Initialize bond order vector                                       *
  !--------------------------------------------------------------------*

  NPBonds = tNUC*(tNUC-1)/2
  call mma_allocate(Bonds,NPBonds,label='Bonds')
  Bonds(:) = Zero

  !--------------------------------------------------------------------*
  ! End of Bond analysis initialization                                *
  !--------------------------------------------------------------------*

end if

!----------------------------------------------------------------------*
! Compute Mulliken atomic charges for each center and basis            *
! function type                                                        *
!----------------------------------------------------------------------*

call mma_allocate(QQ,NXTYP,NNUC,label='QQ')
QQ(:,:) = Zero
IB = 0
IS = 0
IMO = 0
do ISYM=1,NSYM
  NB = NBAS(ISYM)
  if (NB /= 0) then
    IMN = 0
    do MY=1,NB
      do NY=1,MY
        IMN = IMN+1
        DMN = Zero
        ISMO = IMO
        do IO=1,NB
          DMN = DMN+OCCN(IO+IB)*CMO(ISMO+MY)*CMO(ISMO+NY)
          ISMO = ISMO+NB
        end do

        if (DMN_SpinAV) then
          DMN = DMN+xsg*DSc((NY-1)*NB+MY+IMO)
        end if

        if (DoBond) then
          ! Save the Density matrix element (my,ny) and (ny,my) in D_tmp
          ! Save the Overlap matrix element (my,ny) and (ny,my) in S_tmp
          D_tmp(MY+IB,NY+IB) = DMN
          D_tmp(NY+IB,MY+IB) = DMN
          S_tmp(MY+IB,NY+IB) = SMAT(IMN+IS)
          S_tmp(NY+IB,MY+IB) = SMAT(IMN+IS)
        end if

        if (MY /= NY) DMN = Two*DMN
        MYNUC = ICNT(MY+IB)
        NYNUC = ICNT(NY+IB)
        MYTYP = ITYP(MY+IB)
        NYTYP = ITYP(NY+IB)
        if (MY == NY) then
          TERM = SMAT(IMN+IS)*DMN
          if (MYNUC > 0) QQ(MYTYP,MYNUC) = QQ(MYTYP,MYNUC)+TERM
        else
          TERM = Half*SMAT(IMN+IS)*DMN
          if (MYNUC > 0) QQ(MYTYP,MYNUC) = QQ(MYTYP,MYNUC)+TERM
          if (NYNUC > 0) QQ(NYTYP,NYNUC) = QQ(NYTYP,NYNUC)+TERM
        end if
      end do
    end do
    IB = IB+NB
    IS = IS+(NB+NB**2)/2
    IMO = IMO+NB**2
  end if
end do

!----------------------------------------------------------------------*
! Density and overlap matrix handling for bond order                   *
!----------------------------------------------------------------------*

if (DoBond) then

# ifdef _DEBUGPRINT_
  call RecPrt('Density Matrix = ',' ',D_tmp,NBAST,NBAST)
  call RecPrt('Overlap Matrix = ',' ',S_tmp,NBAST,NBAST)
  E = Zero
  do I=1,NBAST
    do J=1,NBAST
      E = E+D_tmp(I,J)*S_tmp(I,J)
    end do
  end do
  write(u6,*)
  write(u6,*) 'Number of electrons as sum of D and S elements = ',E
# endif

  ! In case of symmetry, we desymmetrize D and S through D_blo and S_blo

  if (nSym > 1) then
    iBlo = 1
    iSum = 0
    do i=1,NSYM
      if (nbas(i) /= 0) then
        do j=1,nbas(i)
          do k=1,nbas(i)
            D_blo(iBlo) = D_tmp(iSum+k,iSum+j)
            S_blo(iBlo) = S_tmp(iSum+k,iSum+j)
            iBlo = iBlo+1
          end do
        end do
        iSum = iSum+nbas(i)
      end if
    end do

#   ifdef _DEBUGPRINT_
    write(u6,*) 'D_blo = '
    do i=1,nBas2
      write(u6,*) D_blo(I)
    end do
    write(u6,*) 'S_blo = '
    do i=1,nBas2
      write(u6,*) S_blo(I)
    end do
#   endif

    nScr = NBAST**2
    iSyLbl = 1
    call mma_allocate(Scr,nScr,label='Scr')
    call Desymmetrize(D_blo,nBas2,Scr,nScr,D,nBas,NBAST,P,nSym,iSyLbl)
    call Desymmetrize(S_blo,nBas2,Scr,nScr,S,nBas,NBAST,PInv,nSym,iSyLbl)
    call mma_deallocate(Scr)

    ! Otherwise we simply copy D and S tmp into D and S

  else
    D(:,:) = D_tmp(:,:)
    S(:,:) = S_tmp(:,:)
  end if

# ifdef _DEBUGPRINT_
  write(u6,*) 'After Desymmetrization'
  !call RecPrt('Density Matrix = ',' ',D,NBAST,NBAST)
  !call RecPrt('Overlap Matrix = ',' ',S,NBAST,NBAST)
  write(u6,*) 'Dens=',DDot_(nBast**2,D,1,D,1),DDot_(nBast**2,D,1,[One],0)
  write(u6,*) 'Ovrl=',DDot_(nBast**2,S,1,S,1),DDot_(nBast**2,S,1,[One],0)
  write(u6,*) 'DO  =',DDot_(nBast**2,S,1,D,1)
  E = Zero
  do I=1,NBAST
    do J=1,NBAST
      E = E+D(I,J)*S(I,J)
    end do
  end do
  write(u6,*)
  write(u6,*) 'Number of electrons as sum of D by S elements = ',E
# endif

  ! Finally, we compute the DS matrix as product of D and S

  call DGEMM_('N','N',NBAST,NBAST,NBAST,One,D,NBAST,S,NBAST,Zero,DS,NBAST)

# ifdef _DEBUGPRINT_
  call RecPrt('DS Matrix = ',' ',DS,NBAST,NBAST)
  E = Zero
  do I=1,NBAST
    E = E+DS(I,I)
  end do
  write(u6,*)
  write(u6,*) 'Number of electrons as sum of the DS diagonal = ',E
# endif

  ! in case of first call for UHF we dump everything only

  if (iCase == 0) then
    DSswap(:,:) = DS(:,:)
  end if

  ! in case of second call for UHF we add what dumped before
  ! and release swap memory

  if (iCase == 1) then
    DS(:,:) = DS(:,:)+DSswap(:,:)
    call mma_deallocate(DSswap)

#   ifdef _DEBUGPRINT_
    call RecPrt('DS Matrix = ',' ',DS,NBAST,NBAST)
    E = Zero
    do I=1,NBAST
      E = E+DS(I,I)
    end do
    write(u6,*)
    write(u6,*) 'Number of electrons as sum of the DS diagonal = ',E
#   endif

  end if

end if

!----------------------------------------------------------------------*
! Compute gross atomic charges                                         *
!----------------------------------------------------------------------*

call mma_allocate(QSUM,nNuc,label='QSUM')
call mma_allocate(QSUM_TOT,nNuc,label='QSUM_TOT')
QSUM(:) = Zero
do I=1,NNUC
  do J=1,NXTYP
    QSUM(I) = QSUM(I)+QQ(J,I)
  end do
end do
! if iCase=0, or 1 we need to put/get QSUM
if (iCase == 0) then
  qSwap(1:NNUC) = QSUM(:)
else if (iCase == 1) then
  QSUM_TOT(:) = QSUM(:)+qSwap(1:NNUC)
else if (iCase >= 2) then
  QSUM_TOT(:) = QSUM(:)
end if

!----------------------------------------------------------------------*
! Pick up the nuclear charge                                           *
!----------------------------------------------------------------------*

if (iCase /= 0) then
  call mma_allocate(Chrg,nNuc,label='Chrg')
  call Get_dArray('Effective nuclear charge',Chrg,nNuc)
  Chrg(:) = Chrg(:)*real(nSym/nStab(:),kind=wp)-QSUM_TOT(:)
  if (lSave) call Put_dArray('Mulliken Charge',Chrg,nNuc)
end if

!----------------------------------------------------------------------*
! Compute the 'Mulliken' Bond Order                                    *
!----------------------------------------------------------------------*

if (DoBond .and. (tNUC > 1) .and. (iCase >= 1)) then

# ifdef _DEBUGPRINT_
  write(u6,*) 'nPBonds,tNuc=',nPBonds,tNuc
  do MY=1,NBAST
    AtomA = center(MY)
    write(u6,*) 'AtomA,My=',AtomA,My
  end do
# endif
  do MY=1,NBAST
    AtomA = center(MY)
    if (ICNT(MY) <= 0) cycle    ! skip pseudo center
    do NY=1,MY
      AtomB = center(NY)
      if (ICNT(NY) <= 0) cycle  ! skip pseudo center
      if (AtomA == AtomB) cycle ! same atom

      iPair = (max(AtomA,AtomB)-1)*(max(AtomA,AtomB)-2)/2+min(AtomA,AtomB)

      Bonds(iPair) = Bonds(iPair)+DS(MY,NY)*DS(NY,MY)

#     ifdef _DEBUGPRINT_
      write(u6,*) 'Bond Number=',iPair
      write(u6,*) 'Atom numbers = ',AtomA,AtomB
      write(u6,*) 'Bond number = ',iPair,'bond order = ',Bonds(iPair)
      write(u6,*) 'DS(MY,NY) =',DS(MY,NY)
      write(u6,*) 'DS(NY,MY) =',DS(NY,MY)
#     endif
    end do
  end do

  ! distant atoms could have negative bond order, set to zero

  do I=1,NPBonds
    if (Bonds(I) < Zero) Bonds(I) = Zero
  end do

# ifdef _DEBUGPRINT_
  write(u6,*) 'Bond order vector'
  call TriPrt('Bonds','(10F10.5)',Bonds,tNUC-1)
# endif

end if

!----------------------------------------------------------------------*
! Printout section                                                     *
!----------------------------------------------------------------------*

call mma_allocate(Fac,nNuc,label='Fac')
Fac(:) = real(nStab(:),kind=wp)/real(nSym,kind=wp)

if (iCase == 0) then
  ! first call for UHF, so just dump numbers to swap
  IEND = 0
  ik = 1
  do IST=1,NNUC,6
    IEND = min(IEND+6,NNUC)
    do IT=1,NXTYP
      do j=IST,IEND
        qSwap(NNUC+ik) = QQ(IT,J)
        ik = ik+1
      end do
    end do
  end do
end if

if ((iCase == 1) .and. (iPL >= 2)) then
  ! second call, make a real print out
  if (FullMlk) then
    write(u6,'(6X,A)') 'Mulliken charges per centre and basis function type'
    write(u6,'(6X,A)') '---------------------------------------------------'
  else
    write(u6,'(6X,A)') 'Mulliken charges per centre'
    write(u6,'(6X,A)') '---------------------------'
  end if
  IEND = 0
  ik = 1
  ikk = 1
  call mma_allocate(Q2,nNuc,label='Q2')
  do IST=1,nNuc,6
    IEND = min(IEND+6,nNuc)
    write(u6,*)
    write(u6,'(14X,6(14X,A,4X))') (CNAME(I),I=IST,IEND)
    write(u6,'(14X,6(A12,A12))') (' alpha','  beta',I=IST,IEND)
    do IT=1,NXTYP
      do J=IST,IEND
        Q2(J) = qSwap(NNUC+ik)
        ik = ik+1
      end do
      if (FullMlk) then
        write(u6,'(5X,A8,12F12.4)') Clean_BName(TNAME(IT),0),(Fac(j)*Q2(J),Fac(j)*QQ(IT,J),J=IST,IEND)
      end if
    end do

    do J=IST,IEND
      Q2(J) = qSwap(ikk)
      ikk = ikk+1
    end do

    write(u6,'(6X,A,12F12.4)') 'Total  ',(Fac(i)*Q2(I),Fac(i)*QSUM(I),I=IST,IEND)
    write(u6,'(6X,A,6(6X,F12.4,6X))') 'Total  ',(Fac(i)*(Q2(I)+QSUM(I)),I=IST,IEND)
    write(u6,*)
    write(u6,'(6X,A,6(5X,F12.4,7X))') 'Charge ',(Fac(i)*Chrg(I),I=IST,IEND)
  end do
  call mma_deallocate(Q2)
  write(u6,*)
  write(u6,'(6X,A,F12.6)') 'Total electronic charge=',DDot_(nNuc,[One],0,QSum_TOT,1)
  write(u6,*)
  write(u6,'(6X,A,F12.6)') 'Total            charge=',DDot_(nNuc,[One],0,Chrg,1)

end if
if (iCase == 1) then
  call mma_deallocate(Chrg)
  call mma_deallocate(qSwap)
end if

if (((iCase == 2) .and. (iPL >= 2)) .or. ((iCase == 3) .and. (iPL >= 2))) then
  ! icase=2 for usual mulliken, =2 for spin population.

  if (FullMlk) then
    if (iCase == 2) then
      write(u6,'(6X,A)') 'Mulliken charges per centre and basis function type'
    else
      write(u6,'(6X,A)') 'Mulliken spin population per centre and basis function type'
    end if
    write(u6,'(6X,A)') '---------------------------------------------------'
  else
    if (iCase == 2) then
      write(u6,'(6X,A)') 'Mulliken charges per centre'
    else
      write(u6,'(6X,A)') 'Mulliken spin population per centre'
    end if
    write(u6,'(6X,A)') '---------------------------'
  end if

  IEND = 0
  do IST=1,nNuc,12
    IEND = min(IEND+12,nNuc)
    write(u6,*)
    write(u6,'(14X,12(2X,A))') (CNAME(I),I=IST,IEND)
    if (FullMlk) then
      do IT=1,NXTYP
        write(u6,'(5X,A8,12F8.4)') Clean_BName(TNAME(IT),0),(Fac(j)*QQ(IT,J),J=IST,IEND)
      end do
    end if
    write(u6,'(6X,A,12F8.4)') 'Total  ',(Fac(i)*QSUM(I),I=IST,IEND)
    if (iCase /= 3) then
      write(u6,*)
      write(u6,'(6X,A,12F8.4)') 'N-E    ',(Fac(i)*Chrg(I),I=IST,IEND)
    end if
  end do
  if (iCase == 3) then
    write(u6,*)
    write(u6,'(6X,A,F12.6)') 'Total electronic spin=',DDot_(nNuc,[One],0,QSum,1)
  else
    write(u6,*)
    write(u6,'(6X,A,F12.6)') 'Total electronic charge=',DDot_(nNuc,[One],0,QSum,1)
    write(u6,*)
    TCh = DDot_(nNuc,[One],0,Chrg,1)
    write(u6,'(6X,A,F12.6)') 'Total            charge=',DDot_(nNuc,[One],0,Chrg,1)
    call xml_dDump('FormalCharge','Total charge','a.u',0,[TCh],1,1)
  end if
end if
if (iCase >= 2) then
  call mma_deallocate(Chrg)
end if

! Mulliken bond order print

if (iPL > 2) then
  if ((iCase >= 1) .and. (iCase <= 2) .and. (tNUC > 1) .and. DoBond) then
    write(u6,*)
    write(u6,'(6X,A)') 'Mulliken Bond Order analysis'
    write(u6,'(6X,A)') '----------------------------'
    write(u6,'(6X,A,F5.3,A)') 'Only bonds with order larger than ',BOThrs,' are printed'
    write(u6,*)
    if (nSym > 1) then
      write(u6,'(8X,A)') 'Atom A:Gen.   Atom B:Gen.   Bond Order'
    else
      write(u6,'(8X,A)') 'Atom A        Atom B        Bond Order'
    end if
    do I=1,tNUC-1
      do J=I+1,tNUC
        iPair = (J-1)*(J-2)/2+I
        BO = Bonds(iPair)
        if (BO >= BOThrs) then
          write(u6,'(8X,2(A,4X),F7.3)') LblCnt4(I),LblCnt4(J),BO
        end if
      end do
    end do
    write(u6,*)
  end if
end if

if (iCase /= 0) then
  if (iPL >= 2) then
    call CollapseOutput(0,'   Molecular charges:')
    write(u6,*)
  end if
end if

if (DoBond) then
  if (nSym > 1) then
    call mma_deallocate(P)
    call mma_deallocate(PInv)
    call mma_deallocate(D_blo)
    call mma_deallocate(S_blo)
  end if
  call mma_deallocate(LblCnt4)
  call mma_deallocate(center)
  call mma_deallocate(D_tmp)
  call mma_deallocate(S_tmp)
  call mma_deallocate(D)
  call mma_deallocate(S)
  call mma_deallocate(DS)
  call mma_deallocate(bonds)
end if
call mma_deallocate(CNAME)
call mma_deallocate(nStab)
call mma_deallocate(ICNT)
call mma_deallocate(ITYP)
call mma_deallocate(TNAME)
call mma_deallocate(QQ)
call mma_deallocate(QSUM)
call mma_deallocate(QSUM_TOT)
call mma_deallocate(Fac)

return

end subroutine CHARGE
