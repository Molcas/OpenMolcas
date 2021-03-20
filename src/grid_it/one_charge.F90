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
! Copyright (C) 2011, Francesco Aquilante                              *
!***********************************************************************

subroutine One_CHARGE(NSYM,NBAS,UBNAME,CMO,OCCN,SMAT,iCase,FullMlk,MXTYP,QQ,nNuc)

use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: NSYM, NBAS(NSYM), iCase, MXTYP, nNuc
character(len=LenIn8), intent(in) :: UBNAME(*)
real(kind=wp), intent(in) :: CMO(*), OCCN(*), SMAT(*)
logical(kind=iwp), intent(in) :: FullMlk
real(kind=wp), intent(out) :: QQ(MXTYP,nNuc)
#include "angtp.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: AtomA, AtomB, i, i0, iAB, iAng, IB, iBlo, ICNT(MXBAS), iEnd, iix, iixx, ik, ikk, iM, IMN, IMO, iNuc, IO, &
                     ip_center, ip_Charge, iPair, ipBonds, ipD, ipD_blo, ipD_tmp, ipDS, iPL, ipP, ipPInv, ipS, ipS_blo, ipS_tmp, &
                     ipScr, IS, ISING, ISMO, IST, iStart, iSum, iSwap, iSyLbl, ISYM, IT, ITYP(MXBAS), ix, J, jAng, jEnd, jM, &
                     jPair, jx, k, l, lqSwap, MY, MYNUC, MYTYP, NB, nBas2, NBAST, NDIM, NPBonds, nScr, nStab(MxAtom), NXTYP, NY, &
                     NYNUC, NYTYP, tNUC
real(kind=wp) :: BO, BOThrs, Det, DMN, Q2(MXATOM), QSUM(MXATOM), QSUM_TOT(MXATOM), QSUMI, TERM
logical(kind=iwp) :: DoBond
character(len=100) :: ProgName
character(len=8) TNAME(MXTYP), TMP, TSwap(MXTYP)
!character(len=4) TLbl(MXATOM)
character(len=LenIn) :: CNAME(MXATOM)
character(len=LenIn4) :: LblCnt4(MxAtom)
integer(kind=iwp), save :: ipqswap, ipDSswap
integer(kind=iwp), external :: iPrintLevel
logical(kind=iwp), external :: Reduce_Prt
character(len=100), external :: Get_ProgName
character(len=LenIn8), external :: Clean_Bname
character(len=3), parameter :: AufBau(19) = ['01s',                   &
                                             '02s',            '02p', &
                                             '03s',            '03p', &
                                             '04s',      '03d','04p', &
                                             '05s',      '04d','05p', &
                                             '06s','04f','05d','06p', &
                                             '07s','05f','06d','07p']

!                                                                      *
!***********************************************************************
!                                                                      *
!---- Statement function

real(kind=wp) :: Fac
Fac(i) = real(nStab(i),kind=wp)/real(nSym,kind=wp)
!                                                                      *
!***********************************************************************
!                                                                      *
iPL = iPrintLevel(-1)
if (Reduce_Prt() .and. iPL < 3) iPL = 0
!                                                                      *
!***********************************************************************
!                                                                      *
do i=1,mxTyp
  TName(i) = '        '
end do

!----------------------------------------------------------------------*
!     Get the name of the calling module.                              *
!     If CPFMCPF no bond analysis is done.                             *
!----------------------------------------------------------------------*

ProgName = Get_ProgName()
call Upcase(ProgName)
call LeftAd(ProgName)
iEnd = 1
93 if (ProgName(iEnd:iEnd) /= ' ') then
  iEnd = iEnd+1
  Go To 93
end if

DoBond = .false.

!----------------------------------------------------------------------*
!     Set the Mulliken Bond Order threshold for printout               *
!----------------------------------------------------------------------*

BOThrs = Half

!----------------------------------------------------------------------*
!     GET THE TOTAL NUMBER OF BASIS FUNCTIONS AND CHECK LIMITS         *
!----------------------------------------------------------------------*

NBAST = 0
do I=1,NSYM
  NBAST = NBAST+NBAS(I)
end do
if (NBAST > MXBAS) goto 991

!----------------------------------------------------------------------*
!     Find the list of unique center labels                            *
!----------------------------------------------------------------------*

call Get_cArray('Unique Atom Names',CNAME,LenIn*nNuc)
call Get_iArray('nStab',nStab,nNuc)

!----------------------------------------------------------------------*
!     Find the center label for each basis function                    *
!----------------------------------------------------------------------*

do I=1,NBAST
  ICNT(I) = -1
  do J=1,NNUC
    if (UBNAME(I)(1:LenIn) == CNAME(J)) ICNT(I) = J
  end do
end do

!----------------------------------------------------------------------*
!     Find the type label for each basis function                      *
!----------------------------------------------------------------------*

NXTYP = 0
call ICopy(nBAST,[0],0,ITYP,1)
do I=1,NBAST
  if (ICNT(I) < 0) Go To 99  ! skip pseudo center
  do J=1,NXTYP
    if (J > MxTyp) then
      write(u6,*) 'Charge: J.gt.MxTyp'
      write(u6,*) 'J=',J
      write(u6,*) 'MxTyp=',MxTyp
      write(u6,*) 'Increase MxType and recompile!'
      call Abend()
    end if
    if (UBNAME(I)(LenIn1:LenIn8) == TNAME(J)) then
      ITYP(I) = J
      Go To 99
    end if
  end do
  NXTYP = NXTYP+1
  TNAME(NXTYP) = UBNAME(I)(LenIn1:LenIn8)

  ITYP(I) = NXTYP
99 continue
end do

lqSwap = NNUC+NNUC*NXTYP

if (iCase == 0) then
  ! instead of printing charges we dump everything into a memory same with DS matrix

  call GetMem('CHRG_SWP','ALLO','REAL',ipqSwap,lqSwap)

  if (DoBond) then
    call Allocate_Work(ipDSswap,(NBAST*NBAST))
  end if

end if

!----------------------------------------------------------------------*
!     Do some trivial sorting of the type labels                       *
!----------------------------------------------------------------------*

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
666 iix = ichar(tName(ix)(1:1))
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
!  Write(u6,*) TName(i)
!end do
!
! Now sort with respect to the magnetic index
!
iEnd = jx
iStart = ix
777 do k=0,iTabMx
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

if (jEnd /= iEnd) then
  iStart = jEnd+1
  Go To 777
end if

if (jx /= NxTyp) then
  ix = jx+1
  Go To 666
end if

! Sort according to AufBau

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
do i=1,NxTyp
  TName(i) = TSwap(i)
end do

do I=1,NBAST
  if (ICNT(I) < 0) Go To 98  ! skip pseudo center
  do J=1,NXTYP
    if (UBNAME(I)(LenIn1:LenIn8) == TNAME(J)) then
      ITYP(I) = J
      Go To 98
    end if
  end do
98 continue
end do

!----------------------------------------------------------------------*
!     Get the total number of atoms tNUC, regardless of symmetry       *
!----------------------------------------------------------------------*

call Get_iScalar('LP_nCenter',tNUC)

!----------------------------------------------------------------------*
!     Bond analysis initialization                                     *
!----------------------------------------------------------------------*

if (DoBond) then

  !--------------------------------------------------------------------*
  !   In case of symmetry we need the desymmetrization matrix,         *
  !   for the bond order calculation only.                             *
  !--------------------------------------------------------------------*

  if (nSym > 1) then

    call Allocate_Work(ipP,NBAST**2)
    call Allocate_Work(ipPInv,NBAST**2)
    call Get_dArray('SM',Work(ipP),NBAST**2)
#   ifdef _DEBUGPRINT_
    call RecPrt('SM',' ',Work(ipP),NBAST,NBAST)
#   endif
    call MINV(Work(ipP),Work(ipPInv),ISING,DET,NBAST)
#   ifdef _DEBUGPRINT_
    call RecPrt('SMInv',' ',Work(ipPInv),NBAST,NBAST)
#   endif
    call DGeTMi(Work(ipPInv),NBAST,NBAST)
  end if

  ! Pick up index array of which center a basis function belongs to.
  ! If no symmetry, it is the same as ICNT(I).

  call Allocate_iWork(ip_center,NBAST)
  call Get_iArray('Center Index',iWork(ip_center),NBAST)

  !*********************************************************************

  !--------------------------------------------------------------------*
  !   Initialize symmetric density D_tmp and overlap S_tmp matrices,   *
  !   block D_blo and S_blo matrices (if symmetry)                     *
  !   plus asymmetric D, S and DS matrices                             *
  !--------------------------------------------------------------------*

  call Allocate_Work(ipD_tmp,(NBAST*NBAST))
  call Allocate_Work(ipS_tmp,(NBAST*NBAST))
  call Allocate_Work(ipD,(NBAST*NBAST))
  call Allocate_Work(ipS,(NBAST*NBAST))
  call Allocate_Work(ipDS,(NBAST*NBAST))
  do I=1,(NBAST*NBAST)
    Work(ipD_tmp+I-1) = Zero
    Work(ipS_tmp+I-1) = Zero
    Work(ipD+I-1) = Zero
    Work(ipS+I-1) = Zero
    Work(ipDS+I-1) = Zero
  end do

  if (nSym > 1) then
    nBas2 = 0
    do I=1,nsym
      nBas2 = nBas2+nBas(i)*nBas(i)
    end do
    call Allocate_Work(ipD_blo,nBas2)
    call Allocate_Work(ipS_blo,nBas2)
    do I=1,nBas2
      Work(ipD_blo+I-1) = Zero
      Work(ipS_blo+I-1) = Zero
    end do
  end if

  !--------------------------------------------------------------------*
  !   Find the center label for each atom, regardless of symmetry      *
  !--------------------------------------------------------------------*

  ! Just atom label. It's a double of the next one,
  ! but someone could find it usefull in future

  !call Get_LblCnt_All(TLbl)

  ! Atom labels plus symmetry generator

  call Get_cArray('LP_L',LblCnt4,(LenIn4)*tNUC)
  !do i=1,tNUC
  !  LblCnt(i)(1:LenIn) = LblCnt4(i)(1:LenIn)
  !end do

  !--------------------------------------------------------------------*
  !   Initialize bond order vector                                     *
  !--------------------------------------------------------------------*

  NPBonds = tNUC*(tNUC-1)/2
  call Allocate_Work(ipBonds,NPBonds)
  call FZero(Work(ipBonds),nPBonds)

  !--------------------------------------------------------------------*
  !   End of Bond analysis initialization                              *
  !--------------------------------------------------------------------*

end if

!----------------------------------------------------------------------*
!     Compute Mulliken atomic charges for each center and basis        *
!     function type                                                    *
!----------------------------------------------------------------------*

NDIM = NXTYP*NNUC
call FZero(QQ,nDim)
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

        if (DoBond) then
          !  Save the Density matrix element (my.ny) and (ny,my) in work(ipD_tmp)
          !  Save the Overlap matrix element (my.ny) and (ny,my) in work(ipS_tmp)
          Work(ipD_tmp+(NY+IB-1)*NBAST+MY+IB-1) = DMN
          Work(ipD_tmp+(MY+IB-1)*NBAST+NY+IB-1) = DMN
          Work(ipS_tmp+(NY+IB-1)*NBAST+MY+IB-1) = SMAT(IMN+IS)
          Work(ipS_tmp+(MY+IB-1)*NBAST+NY+IB-1) = SMAT(IMN+IS)
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
!     Density and overlap matrix handling for bond order               *
!----------------------------------------------------------------------*
if (DoBond) then

# ifdef _DEBUGPRINT_
  call RecPrt('Density Matrix = ',' ',Work(ipD_tmp),NBAST,NBAST)
  call RecPrt('Overlap Matrix = ',' ',Work(ipS_tmp),NBAST,NBAST)
  E = Zero
  do I=1,NBAST
    do J=1,NBAST
      E = E+Work(ipD_tmp+(J-1)*NBAST+I-1)*Work(ipS_tmp+(J-1)*NBAST+I-1)
    end do
  end do
  write(u6,*)
  write(u6,*) 'Number of electrons as sum of D and S elements = ',E
# endif

  ! In case of symmetry, we desymmetrize D and S through D_blo and S_blo
  if (nSym > 1) then
    iBlo = 0
    iSum = 0
    do i=1,NSYM
      if (nbas(i) /= 0) then
        do j=0,nbas(i)-1
          do k=0,nbas(i)-1
            Work(ipD_blo+iBlo) = Work(ipD_tmp+(j+iSum)*NBAST+iSum+k)
            Work(ipS_blo+iBlo) = Work(ipS_tmp+(j+iSum)*NBAST+iSum+k)
            iBlo = iBlo+1
          end do
        end do
        iSum = iSum+nbas(i)
      end if
    end do

#   ifdef _DEBUGPRINT_
    write(u6,*) 'D_blo = '
    do i=1,nBas2
      write(u6,*) (Work(ipD_blo+I-1))
    end do
    write(u6,*) 'S_blo = '
    do i=1,nBas2
      write(u6,*) (Work(ipS_blo+I-1))
    end do
#   endif

    nScr = MXBAS*NBAST
    iSyLbl = 1
    call Allocate_Work(ipScr,nScr)
    call Desymmetrize(Work(ipD_blo),nBas2,Work(ipScr),nScr,Work(ipD),nBas,NBAST,Work(ipP),nSym,iSyLbl)
    call Free_Work(ipScr)

    call Allocate_Work(ipScr,nScr)
    call Desymmetrize(Work(ipS_blo),nBas2,Work(ipScr),nScr,Work(ipS),nBas,NBAST,Work(ipPInv),nSym,iSyLbl)
    call Free_Work(ipScr)

  ! Otherwise we simply copy D and S tmp into D and S

  else
    call dcopy_(nBasT**2,Work(ipD_tmp),1,Work(ipD),1)
    call dcopy_(nBasT**2,Work(ipS_tmp),1,Work(ipS),1)
    !do I=1,NBAST*NBAST
    !   Work(ipD+I-1) = Work(ipD_tmp+I-1)
    !   Work(ipS+I-1) = Work(ipS_tmp+I-1)
    !end do
  end if

# ifdef _DEBUGPRINT_
  write(u6,*) 'After Desymmetrization'
  !call RecPrt('Density Matrix = ', ' ', Work(ipD), NBAST, NBAST)
  !call RecPrt('Overlap Matrix = ', ' ', Work(ipS), NBAST, NBAST)
  write(u6,*) 'Dens=',DDot_(nBast**2,Work(ipD),1,Work(ipD),1),DDot_(nBast**2,Work(ipD),1,[One],0)
  write(u6,*) 'Ovrl=',DDot_(nBast**2,Work(ipS),1,Work(ipS),1),DDot_(nBast**2,Work(ipS),1,[One],0)
  write(u6,*) 'DO  =',DDot_(nBast**2,Work(ipS),1,Work(ipD),1)
  E = Zero
  do I=1,NBAST
    do J=1,NBAST
      E = E+Work(ipD+(J-1)*NBAST+I-1)*Work(ipS+(J-1)*NBAST+I-1)
    end do
  end do
  write(u6,*)
  write(u6,*) 'Number of electrons as sum of D by S elements = ',E
# endif

  ! Finally, we compute the DS matrix as product of D and S

  call DGEMM_('N','N',NBAST,NBAST,NBAST,One,Work(ipD),NBAST,Work(ipS),NBAST,Zero,Work(ipDS),NBAST)

# ifdef _DEBUGPRINT_
  call RecPrt('DS Matrix = ',' ',Work(ipDS),NBAST,NBAST)
  E = Zero
  do I=1,NBAST
    E = E+Work(ipDS+(I-1)*NBAST+I-1)
  end do
  write(u6,*)
  write(u6,*) 'Number of electrons as sum of the DS diagonal = ',E
# endif

  ! in case of first call for UHF we dump everything only

  if (iCase == 0) then
    do I=1,NBAST
      Work(ipDSswap+I-1) = Work(ipDS+I-1)
    end do
  end if

  ! in case of second call for UHF we add what dumped before and release swap memory

  if (iCase == 1) then
    do I=1,NBAST
      Work(ipDS+I-1) = Work(ipDS+I-1)+Work(ipDSswap+I-1)
    end do
    call Free_Work(ipDSswap)

#   ifdef _DEBUGPRINT_
    call RecPrt('DS Matrix = ',' ',Work(ipDS),NBAST,NBAST)
    E = Zero
    do I=1,NBAST
      E = E+Work(ipDS+(I-1)*NBAST+I-1)
    end do
    write(u6,*)
    write(u6,*) 'Number of electrons as sum of the DS diagonal = ',E
#   endif

  end if

end if

!----------------------------------------------------------------------*
!     Compute gross atomic charges                                     *
!----------------------------------------------------------------------*

do I=1,NNUC
  QSUMI = Zero
  do J=1,NXTYP
    QSUMI = QSUMI+QQ(J,I)
  end do
  QSUM(I) = QSUMI
end do
! if iCase=0, or 1 we need to put/get QSUM
do i=1,NNUC
  if (iCase == 0) Work(ipqSwap+i-1) = QSUM(I)
  if (iCase == 1) QSUM_TOT(I) = QSUM(I)+Work(ipqSwap+i-1)
  if (iCase >= 2) QSUM_TOT(I) = QSUM(I)
end do

!----------------------------------------------------------------------*
!     Pick up the nuclear charge                                       *
!----------------------------------------------------------------------*

if (iCase /= 0) then
  call Allocate_Work(ip_Charge,nNuc)
  call Get_dArray('Effective nuclear charge',Work(ip_Charge),nNuc)
  do iNuc=0,nNuc-1
    Work(ip_Charge+iNuc) = Work(ip_Charge+iNuc)*(nSym/nStab(iNuc+1))
  end do
  call DaXpY_(nNuc,-One,QSUM_TOT,1,Work(ip_Charge),1)
end if

!----------------------------------------------------------------------*
!     Compute the 'Mulliken' Bond Order                                *
!----------------------------------------------------------------------*

if (DoBond .and. (tNUC > 1) .and. (iCase >= 1)) then

# ifdef _DEBUGPRINT_
  write(u6,*) 'nPBonds,tNuc=',nPBonds,tNuc
  do MY=1,NBAST
    AtomA = iWork(ip_center+MY-1)
    write(u6,*) 'AtomA,My=',AtomA,My
  end do
# endif
  do MY=1,NBAST
    AtomA = iWork(ip_center+MY-1)
    if (ICNT(MY) <= 0) Go To 95    ! skip pseudo center
    do NY=1,MY
      AtomB = iWork(ip_center+NY-1)
      if (ICNT(NY) <= 0) Go To 94  ! skip pseudo center
      if (AtomA == AtomB) Go To 94  ! same atom

      iPair = (max(AtomA,AtomB)-1)*(max(AtomA,AtomB)-2)/2+min(AtomA,AtomB)
      jPair = ipBonds-1+iPair

      Work(jPair) = Work(jPair)+Work(ipDS+(NY-1)*NBAST+MY-1)*Work(ipDS+(MY-1)*NBAST+NY-1)

#     ifdef _DEBUGPRINT_
      write(u6,*) 'Bond Number=',iPair
      write(u6,*) 'Atom numbers = ',AtomA,AtomB
      write(u6,*) 'Bond number = ',iPair,'bond order = ',Work(jPair)
      write(u6,*) 'Work(ipDS+ (NY-1) * NBAST + MY-1) =',Work(ipDS+(NY-1)*NBAST+MY-1)
      write(u6,*) 'Work(ipDS+ (MY-1) * NBAST + NY-1) =',Work(ipDS+(MY-1)*NBAST+NY-1)
#     endif
94    continue
    end do
95  continue
  end do

  ! distant atoms could have negative bond order, set to zero

  do I=1,NPBonds
    if (Work(ipBonds+I-1) < Zero) Work(ipBonds+I-1) = Zero
  end do

# ifdef _DEBUGPRINT_
  write(u6,*) 'Bond order vector'
  call TriPrt('Bonds','(10F10.5)',Work(ipBonds),tNUC-1)
# endif

end if

!----------------------------------------------------------------------*
!     Printout section                                                 *
!----------------------------------------------------------------------*

if (iCase == 0) then
  ! first call for UHF, so just dump numbers to swap
  IEND = 0
  ik = 0
  do IST=1,NNUC,6
    IEND = min(IEND+6,NNUC)
    do IT=1,NXTYP
      do j=IST,IEND
        Work(ipqSwap+NNUC+ik) = QQ(IT,J)
        ik = ik+1
      end do
    end do
  end do
end if

if (iCase == 1 .and. iPL >= 2) then
  ! second call, make a real print out
  if (FullMlk) then
    write(u6,'(6X,A)') 'Mulliken charges per centre and basis function type'
    write(u6,'(6X,A)') '---------------------------------------------------'
  else
    write(u6,'(6X,A)') 'Mulliken charges per centre'
    write(u6,'(6X,A)') '---------------------------'
  end if
  IEND = 0
  ik = 0
  ikk = 0
  do IST=1,nNuc,6
    IEND = min(IEND+6,nNuc)
    write(u6,*)
    write(u6,'(14X,6(14X,A,4X))') (CNAME(I),I=IST,IEND)
    write(u6,'(14X,6(A12,A12))') (' alpha','  beta',I=IST,IEND)
    do IT=1,NXTYP
      do J=IST,IEND
        Q2(J) = Work(ipqSwap+NNUC+ik)
        ik = ik+1
      end do
      if (FullMlk) then
        write(u6,'(5X,A8,12F12.4)') Clean_BName(TNAME(IT),0),(Fac(j)*Q2(J),Fac(j)*QQ(IT,J),J=IST,IEND)
      end if
    end do

    do J=IST,IEND
      Q2(J) = Work(ipqSwap+ikk)
      ikk = ikk+1
    end do

    write(u6,'(6X,A,12F12.4)') 'Total  ',(Fac(i)*Q2(I),Fac(i)*QSUM(I),I=IST,IEND)
    write(u6,'(6X,A,6(6X,F12.4,6X))') 'Total  ',(Fac(i)*(Q2(I)+QSUM(I)),I=IST,IEND)
    write(u6,*)
    write(u6,'(6X,A,6(5X,F12.4,7X))') 'Charge ',(Fac(i)*Work(ip_Charge+I-1),I=IST,IEND)
  end do
  write(u6,*)
  !Write(u6,'(6X,A,F12.6)') 'Total electronic charge=',DDot_(nNuc,[One],0,QSum_TOT,1)
  write(u6,*)
  !Write(u6,'(6X,A,F12.6)') 'Total            charge=',DDot_(nNuc,[One],0,Work(ip_Charge),1)

end if
if (iCase == 1) then
  call Free_Work(ip_Charge)
  call GetMem('CHRG_SWP','FREE','REAL',ipqSwap,lqSwap)
end if

if ((iCase == 2 .and. iPL >= 2) .or. (iCase == 3 .and. iPL >= 2)) then
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
      !write(u6,'(6X,A,12F8.4)')'N-E    ',(Fac(i)*Work(ip_Charge+I-1),I=IST,IEND)
    end if
  end do
  if (iCase == 3) then
    write(u6,*)
    !write(u6,'(6X,A,F12.6)') 'Total electronic spin=',DDot_(nNuc,[One],0,QSum,1)
  else
    write(u6,*)
    !write(u6,'(6X,A,F12.6)') 'Total electronic charge=',DDot_(nNuc,[One],0,QSum,1)
    write(u6,*)
    !TCh = DDot_(nNuc,[One],0,Work(ip_Charge),1)
    !write(u6,'(6X,A,F12.6)') 'Total            charge=',DDot_(nNuc,[One],0,Work(ip_Charge),1)
    !call xml_dDump('FormalCharge','Total charge','a.u',0,TCh,1,1)
  end if
end if
if (iCase >= 2) then
  call Free_Work(ip_Charge)
end if

!  Mulliken bond order print

if (iPL <= 2) Go To 9999
if ((iCase >= 1) .and. (iCase <= 2) .and. (tNUC > 1) .and. (DoBond)) then
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
      BO = Work(ipBonds-1+iPair)
      if (BO >= BOThrs) then
        write(u6,'(8X,2(A,4X),F7.3)') LblCnt4(I),LblCnt4(J),BO
      end if
    end do
  end do
  write(u6,*)
end if

9999 continue
if (DoBond) then
  if (nSym > 1) then
    call Free_Work(ipP)
    call Free_Work(ipPInv)
    call Free_Work(ipD_blo)
    call Free_Work(ipS_blo)
  end if
  call Free_iWork(ip_center)
  call Free_Work(ipD_tmp)
  call Free_Work(ipS_tmp)
  call Free_Work(ipD)
  call Free_Work(ipS)
  call Free_Work(ipDS)
  call Free_Work(ipBonds)
end if

return

!----------------------------------------------------------------------*
!     Error Exits                                                      *
!----------------------------------------------------------------------*

991 write(u6,'(/6X,A)') 'The number of basis functions exceeds the present limit'
call Abend()
!992 write(u6,'(/6X,A)') 'The number of basis functions exceeds the present limit'
!call Abend()
!993 write(u6,'(/6X,A)') 'Warning: Total charge is not equal to number of electrons'
!call Abend()

return

end subroutine One_CHARGE
