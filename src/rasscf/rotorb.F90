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

subroutine rotorb(cmoo,cmon,c,x,x2,y,thmax,FA)
! RASSCF program: version IBM-3090: SX section
!
! PURPOSE: Calculation of the rotation matrix X from the super-CI
!          coefficient matrix C, formation of exp(x), and rotation
!          of the orbitals with this matrix.
!          These orbitals are used in the next
!          RASSCF iteration as starting orbitals.
!          Called from SXCTL
!
!      ********** IBM-3090 MOLCAS Release: 90 02 22 **********

use Index_Functions, only: nTri_Elem
use gas_data, only: iDoGAS, NGAS, NGSSH
use rasscf_global, only: CMAX, iXSym, PURIFY, ROTMAX
use PrintLevel, only: DEBUG, TERSE, VERBOSE
use output_ras, only: IPRLOC
use general_data, only: NASH, NBAS, NDEL, NFRO, NISH, NORB, NRS1, NRS2, NSSH, NSYM, NTOT2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: cmoo(*), cmon(*), c(*), x(*), x2(*), y(*), THMAX, FA(*)
integer(kind=iwp) :: I, IB, ICORE, ICOREGAS(0:4,0:4), IDAMP, IDAMPGAS(0:4,0:4), IGAS, II, IJ, IO, iOff, iOrb, iPrLev, iSpace, IST, &
                     ISTBM, ISTMO, ISTMO1, ISUM, ISYM, jPr, jSPace, MOType, NACI, NACJ, NAE, NAO, NB, NBO, ND, NDB, NEO, NF, NFB, &
                     NI, NII, NIO, NIO1, NJ, NO, NOC, NOC1, NP, NR
real(kind=wp) :: COREGAS(10), DAMPGAS(10), TERM, THM, Xn, XX
logical(kind=iwp) :: iFrzAct
real(kind=wp), allocatable :: SqFA(:), Unt(:)
real(kind=wp), parameter :: ACC = 1.0e-13_wp, Thrs = 1.0e-14_wp

IPRLEV = IPRLOC(4)
if (IPRLEV >= DEBUG) then
  write(u6,*) ' Entering ROTORB'

  write(u6,*)
  write(u6,*) 'FI+FA in RotOrb by Unitary transform'
  write(u6,*) ' --------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iOrb = nOrb(iSym)
    call TriPrt(' ',' ',FA(iOff),iOrb)
    iOff = iOff+nTri_Elem(iOrb)
  end do

end if

istbm = 1
istmo1 = 0
ib = 0
rotmax = Zero
cmax = Zero
thmax = Zero
iOff = 1
! A long loop over symmetry
do isym=1,nsym
  no = norb(isym)
  nb = nbas(isym)
  nf = nfro(isym)
  io = ib+nf
  nio = nish(isym)
  nao = nash(isym)
  noc = nio+nao
  neo = nssh(isym)
  nae = nao+neo
  if ((noc == 0) .or. (nae == 0) .or. (nio+neo == 0)) then

    if (IPRLEV >= VERBOSE) write(u6,*) 'No rotations active for symmetry=',isym
    !No rotations active in this symmetry
    !Move orbitals from old to new set

    cmon(istmo1+1:istmo1+nb**2) = cmoo(istmo1+1:istmo1+nb**2)
  else
    istmo = istmo1+nf*nb

    ! Form the quadratic matrix X from the coefficients C

    x(1:no**2) = Zero
    do nr=1,noc
      do np=max(nr+1,nio+1),no
        jpr = no*(nr-1)+np
        !jrp = no*(np-1)+nr
        xx = c(istbm+nr+noc*(np-nio-1))

        ! Any numerical information smaller than Acc is considered
        ! numerical noise and is ignored.

        if (abs(xx) < Acc) cycle

        x(jpr) = xx
        !x(jrp) = -xx
      end do
    end do

    ! Set zero to matrix elements corresponding to orbitals
    ! not allowed to rotate into each other.

    ij = 0
    do ni=1,no
      do nj=1,no
        ij = ij+1
        if (ixsym(ni+io) /= ixsym(nj+io)) x(ij) = Zero
      end do
    end do

    ! Freeze Active orbitals for difficult orbital optimization
    ! This is activated by FRAC keyword in RASSCF module

    iFrzAct = .false.
    if (iFrzAct) then
      ij = 0
      do ni=1,no
        do nj=1,no
          ij = ij+1
          !io = ib+nf ! offset counting all nbas of previous sym and nfro of current sym.
          if ((ni > nio) .and. (ni < nio+nao) .or. (nj > nio) .and. (nj < nio+nao)) x(ij) = Zero
        end do
      end do
    end if

    ! For optimization of the core-hole we want to damp/eliminate
    ! certain rotations to find the local minimum.

    ! Damp/eliminate certain rotations inside RAS/GAS.

    ! Only loop over the active spaces

    IDAMP = 0
    IDAMPGAS = 0
    DAMPGAS = Zero
    ICORE = 0
    ICOREGAS = 0
    COREGAS = 0.1_wp
    !COREGAS(1) = Zero
    !COREGAS(2) = Zero
    ICOREGAS(0,1) = 1
    ICOREGAS(1,0) = 2
    !ICOREGAS(0,2) = 1
    !ICOREGAS(2,0) = 2
    ICOREGAS(2,1) = 3
    ICOREGAS(1,2) = 4
    ICOREGAS(3,1) = 5
    ICOREGAS(1,3) = 6
    ICOREGAS(1,4) = 1
    ICOREGAS(4,1) = 2
    !write(u6,*) 'NOC,NO',NOC,NO
    if ((ICORE == 1) .or. (IDAMP == 1)) then ! New keywords

      IJ = 0
      do NACI=1,NO ! Notice only loop over inactive and active
        if (IDOGAS) then ! Find which GAS this index belong to
          ISPACE = 0
          if (NACI <= NIO) then
            ISPACE = 0
          else if (NACI > NOC) then
            ISPACE = NGAS+1
          else
            ISUM = 0
            do IGAS=1,NGAS
              if ((NGSSH(IGAS,ISYM)+NIO+ISUM) >= NACI) then
                ISPACE = IGAS
                exit
              else
                ISUM = ISUM+NGSSH(IGAS,ISYM)
              end if
            end do
          end if
        else ! Find which RAS this index belong to
          if (NACI <= NIO) then
            ISPACE = 0
          else if (NACI > NOC) then
            ISPACE = 4
          else if ((NRS1(ISYM)+NIO) >= NACI) then
            ISPACE = 1
          else if ((NRS1(ISYM)+NRS2(ISYM)+NIO) >= NACI) then
            ISPACE = 2
          else
            ISPACE = 3
          end if
        end if

        do NACJ=1,NO ! Over all orbitals due to counting
          !write(u6,*) 'NACI,NACJ',NACI,NACJ
          if (IDOGAS) then ! Find which GAS this index belong to
            JSPACE = 0
            if (NACJ <= NIO) then
              JSPACE = 0
            else if (NACJ > NOC) then
              JSPACE = NGAS+1
              !IJ = IJ+1
              !cycle
            else
              ISUM = 0
              do IGAS=1,NGAS
                if ((NGSSH(IGAS,ISYM)+NIO+ISUM) >= NACJ) then
                  JSPACE = IGAS
                  exit
                else
                  ISUM = ISUM+NGSSH(IGAS,ISYM)
                end if
              end do
            end if
          else ! Find which RAS this index belong to
            JSPACE = 0
            if (NACJ <= NIO) then
              JSPACE = 0
            else if (NACJ > NOC) then
              JSPACE = 4
              !IJ = IJ+1
              !write(u6,*) 'IJ',IJ
              !cycle
            else if ((NRS1(ISYM)+NIO) >= NACJ) then
              JSPACE = 1
            else if ((NRS1(ISYM)+NRS2(ISYM)+NIO) >= NACJ) then
              JSPACE = 2
            else
              JSPACE = 3
            end if
          end if
          IJ = IJ+1
          write(u6,*) 'NACI,NACJ',NACI,NACJ
          write(u6,*) 'ISPACE,JSPACE',ISPACE,JSPACE
          write(u6,*) 'IJ',IJ
          if (IDAMP == 1) then
            if (IDAMPGAS(ISPACE,JSPACE) > 0) X(IJ) = X(IJ)*DAMPGAS(IDAMPGAS(ISPACE,JSPACE))
          else if (ICORE == 1) then
            if (ICOREGAS(ISPACE,JSPACE) > 0) then
              write(u6,*) ' damping'
              write(u6,*) 'ICOREGAS(ISPACE,JSPACE)',ICOREGAS(ISPACE,JSPACE)
              write(u6,*) 'fac',COREGAS(ICOREGAS(ISPACE,JSPACE))
              X(IJ) = X(IJ)*COREGAS(ICOREGAS(ISPACE,JSPACE))
            end if
          end if
        end do

      end do
      ! temp print
      ij = 0
      do ni=1,no
        do nj=1,no
          ij = ij+1
          write(u6,*) 'ij,x(ij)',ij,x(ij)
        end do
      end do

    end if

    ! Now form the unitary matrix exp(X)

    call exp_eig(no,x,thm)
    thmax = max(thmax,thm)

    ! Check for largest non diagonal element

    ij = 0
    do ni=1,no
      do nj=1,no
        ij = ij+1
        if (abs(x(ij)) < Thrs) then
          x(ij) = Zero
          cycle
        end if
        if (ni == nj) cycle
        if (abs(x(ij)) > abs(rotmax)) rotmax = x(ij)
      end do
    end do

    ! Check for large rotations and phase of new orbital

    ii = 1
    ist = 0
    do ni=1,no
      if (ni <= nio) then
        motype = 1
        xn = Zero
        do nii=1,nio
          xn = xn+x(ist+nii)**2
        end do
        if (IPRLEV >= TERSE) then
          if (xn < Half) then
            call WarningMessage(1,'Large orbital rotation.')
            write(u6,1010) ni,isym,motype,xn
          end if
        end if
      end if
      if ((ni > nio) .and. (ni <= noc)) then
        motype = 2
        xn = Zero
        nio1 = nio+1
        do nii=nio1,noc
          xn = xn+x(ist+nii)**2
        end do
        if (IPRLEV >= TERSE) then
          if (xn < Half) then
            call WarningMessage(1,'Large orbital rotation.')
            write(u6,1010) ni,isym,motype,xn
          end if
        end if
      end if
      if (ni > noc) then
        motype = 3
        noc1 = noc+1
        xn = Zero
        do nii=noc1,no
          xn = xn+x(ist+nii)**2
        end do
        if (IPRLEV >= TERSE) then
          if (xn < Half) then
            call WarningMessage(1,'Large orbital rotation.')
            write(u6,1010) ni,isym,motype,xn
          end if
        end if
      end if
      ii = ii+no+1
      ist = ist+no
    end do

    ! Transformation of the Fock matrix according to the orbital rotation matrix.
    ! This step is not required other than for FCIQMC as orbital energies could be
    ! used for choosing the reference determinant.

    if (iprlev >= debug) then
      call recprt('X in RotOrb',' ',x,no,no)
      write(u6,*) 'FA for sym = ',iSym
      write(u6,*) 'iOff is set to = ',iOff
      call TriPrt(' ',' ',FA(iOff),no)
    end if

    call mma_allocate(Unt,no*no,Label='Unt')
    call mma_allocate(SqFA,no*no,Label='SqFA')
    Unt(:) = Zero
    SqFA(:) = Zero
    call Square(FA(iOff),SqFA,1,no,no)
    if (iprlev >= debug) call recprt('Square FA in RotOrb',' ',SqFA,no,no)
    call DGEMM_('N','N',no,no,no,One,x,no,SqFA,no,Zero,Unt,no)

    call DGEMM_('N','T',no,no,no,One,Unt,no,x,no,Zero,SqFA,no)

    call Fold_Mat(1,[no],SqFA,FA(iOff))

    iOff = iOff+nTri_Elem(no)

    call mma_deallocate(Unt)
    call mma_deallocate(SqFA)

    ! Print output in MO-basis

    ! Transform new orbitals to AO-basis

    call DGEMM_('N','N',nb,no,no,One,cmoo(istmo+1),nb,x,no,Zero,cmon(istmo+1),nb)

    ! Calculate max change in occupied molecular orbital coefficient

    nbo = noc*nb
    do i=1,nbo
      term = cmon(istmo+i)-cmoo(istmo+i)
      if (abs(term) > abs(cmax)) cmax = term
    end do

  end if

  ! Move frozen orbitals from old to new set of orbitals

  if (nf /= 0) then
    nfb = nf*nb
    cmon(istmo1+1:istmo1+nfb) = cmoo(istmo1+1:istmo1+nfb)
  end if

  ! MOVE DELETED ORBITALS.

  nd = ndel(isym)
  if (nd /= 0) then
    ist = istmo1+nb*(nf+no)+1
    ndb = nd*nb
    cmon(ist:ist+ndb-1) = cmoo(ist:ist+ndb-1)
  end if

  istbm = istbm+noc*nae
  istmo1 = istmo1+nb*nb
  ib = ib+nb
end do
!^ End of the long loop over symmetry

if (PURIFY(1:4) == 'ATOM') call SPHPUR(CMON)
if (PURIFY(1:6) == 'LINEAR') call LINPUR(CMON)

! Orthogonalize new MO's and move them back to CMOO

call supsch(x2,cmoo,cmon)
! PAM07: In ortho, cmoo is scratch and will be destroyed.
call ortho_rasscf(x2,cmoo,cmon,y)
cmoo(1:ntot2) = cmon(1:ntot2)

if (iprlev >= debug) then
  write(u6,*)
  write(u6,*) ' >>> Exit RotOrb <<< '
  write(u6,*)
end if

1010 format(6x,'Molecular orbital',i4,' of symmetry',i2,' MO space',i2,'  weight is',f12.6)

end subroutine rotorb
