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
! Copyright (C) 1994,1996,2014, Per Ake Malmqvist                      *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine NEWFOCK(FIFA,NFIFA,CMO,NCMO,DREF,nDREF)
! Purpose: Modify the standard fock matrix for experimental
! purposes. The string variable FOCKTYPE (character*8) has a
! keyword value given as input. The experimental user modifies
! this routine to suit his purposes, and links with the rest
! of the program. The Fock matrix is given as call parameter
! and returned after modification.
! To define the modified Fock matrix, a number of arrays on
!  LUONE may be useful. In addition, the active 1- and 2-
! electron density matrices, and the inactive Fock matrix
! FIMO, are available in workspace at DREF,
! PREF,FIMO, and FIFA.
! Two-electron integrals involving non-frozen, non-deleted
! orbitals, at most two secondary, are available from
! subroutines COUL and EXCH (See).
! Meaningful modifications must define a Fock operator
! which is invariant to transformations among inactives,
! among actives, and among virtuals.
! Coded 94-01-31 by Malmqvist, for CASPT2 MOLCAS-3.
! Modif 96-10-06 by Malmqvist, restructured, options added.
! Modif 14-03-19 by Malmqvist, restructured, options removed.

use caspt2_global, only: iPrGlb
use PrintLevel, only: USUAL
use caspt2_module, only: FockType, IfChol, nAES, nAMx, nAsh, nAshT, nIMx, nIsh, nOMx, nOrb, nOSqT, nSMx, nSsh, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NFIFA, NCMO, nDREF
real(kind=wp), intent(inout) :: FIFA(NFIFA)
real(kind=wp), intent(in) :: CMO(NCMO), DREF(nDREF)
integer(kind=iwp) :: I, IA, IATOT, IDDVX, IDREF, IDTT, IDTU, IDUT, IFGFOCK, II, IP, IQ, IR, IS, ISC, ISYM, ISYMPQ, ISYMRS, IT, &
                     ITABS, ITTOT, ITU, IU, IUABS, IUTOT, IV, IX, J, KFIFA, LEIG, LEV1, LEV2, LINT, LSC, LSC1, LSC2, LSCR, LXAI, &
                     LXPQ, LXQP, MA, MI, MTRES, N3, NA, NA2, NA3, NAS, NASQES, NASQT, NATR, NATRES, NI, NIA, NINTBUF, NO, NOSQES, &
                     NOTRES, NS, NSCR, NSCR1, NSCR2, NSCR3, NSQES, NTRES
real(kind=wp) :: D, DDVX, EIGVAL, VAL, VALTU, VALUT, X
real(kind=wp), allocatable :: INTBUF(:), DSQ(:), DD(:), DDTR(:), TWOMDSQ(:), XMAT(:), SC(:)

! I never meant to cause you any sorrow
if (FOCKTYPE == 'STANDARD') return

! Options MC and MC2 removed, PAM March 2014.
!CPAM96 The option FOCKTYPE='MC' added 961006. This option will
!C replace the active/active block with the MCSCF Fock matrix
!C while zeroing any non-diagonal blocks. This option is usually
!C quite ridiculous, but it can be used in very particular cases
!C when all active orbitals are singly occupied.
!CPAM96 The option FOCKTYPE='MC2' added 961006. Similar to the
!C above, but using as  active/active block the matrix
!C    D**(-1/2) FMC D**(-1/2)

IFGFOCK = 0
if ((FOCKTYPE == 'G1') .or. (FOCKTYPE == 'G2') .or. (FOCKTYPE == 'G3')) IFGFOCK = 1

if ((IFGFOCK == 0) .or. (NASHT <= 0)) return
if (IPRGLB >= USUAL) write(u6,*) ' THE FOCK MATRIX IS MODIFIED BY KEYWORD FOCKTYPE=',FOCKTYPE

! Determine sizes of areas for memory allocation
NASQT = 0
NATR = 0
do ISYM=1,NSYM
  NI = NISH(ISYM)
  NA = NASH(ISYM)
  NS = NSSH(ISYM)
  NO = NORB(ISYM)
  NASQT = NASQT+NA**2
  NATR = NATR+(NA*(NA+1))/2
end do
NINTBUF = NOMX**2
NSCR1 = NAMX*max(2*NAMX,NIMX,NSMX)
NSCR2 = 3*NAMX*(NAMX+1)
NSCR3 = NAMX*(3*NAMX+1)
NSCR = NSCR1
if (FOCKTYPE == 'G2') NSCR = NSCR2
if (FOCKTYPE == 'G3') NSCR = NSCR3

! Allocate memory: Integral buffer and scratch array:
call mma_allocate(INTBUF,2*NINTBUF,LABEL='INTBUF')
LINT = 1
LSCR = LINT+NINTBUF
call mma_allocate(SC,NSCR,Label='SC')
LSC = 1
! Form symmetry-packed squares of density matrix DSQ, and
! similarly (2I-DSQ)
call mma_allocate(DSQ,NASQT,Label='DSQ')
call mma_allocate(TWOMDSQ,NASQT,LABEL='TWOMDSQ')
! Symmetry-packed triangles of D*(2I-D)
call mma_allocate(DDTR,NATR,Label='DDTR')
! Temporary use of single square symmetry-block:
call mma_allocate(DD,NAMX**2,Label='DD')
! The exchange matrix, A(pq)=sum over rs of (ps,rq)*DD(rs)
call mma_allocate(XMAT,NOSQT,Label='XMAT')
NSQES = 0
do ISYM=1,NSYM
  NA = NASH(ISYM)
  do IT=1,NA
    ITABS = IT+NAES(ISYM)
    do IU=1,NA
      IUABS = IU+NAES(ISYM)
      IDREF = (ITABS*(ITABS-1))/2+IUABS
      D = DREF(IDREF)
      IDTU = NSQES+IT+NA*(IU-1)
      IDUT = NSQES+IU+NA*(IT-1)
      DSQ(IDTU) = D
      DSQ(IDUT) = D
    end do
  end do
  do I=1,NA*NA
    IDTU = NSQES+I
    TWOMDSQ(IDTU) = -DSQ(IDTU)
  end do
  do I=1,NA*NA,(NA+1)
    IDTT = NSQES+I
    TWOMDSQ(IDTT) = Two-DSQ(IDTT)
  end do
  NSQES = NSQES+NA**2
end do

! Create the matrix DDTR =D(2I-D) (triangular symmetry blocks)
! Use also temporary DD, single symmetry blocks of D*(2I-D):
NTRES = 0
NSQES = 1
do ISYM=1,NSYM
  NA = NASH(ISYM)
  if (NA > 0) then
    N3 = (NA*(NA+1))/2
    call DGEMM_('N','N',NA,NA,NA,One,DSQ(NSQES),NA,TWOMDSQ(NSQES),NA,Zero,DD,NA)
    call TRIANG(NA,DD)
    DDTR(NTRES+1:NTRES+N3) = DD(1:N3)
    NTRES = NTRES+N3
    NSQES = NSQES+NA**2
  end if
end do

! Calculation of the exchange matrix, A(pq)=sum over rs of (ps,rq)*DD(rs)
XMAT(:) = Zero
if (IfChol) then
  call Cho_Amatrix(XMAT,NOSQT,CMO,NCMO,DDTR,NATR)
else
  NOSQES = 0
  do ISYMPQ=1,NSYM
    NI = NISH(ISYMPQ)
    NA = NASH(ISYMPQ)
    NS = NSSH(ISYMPQ)
    NO = NORB(ISYMPQ)
    if (NO > 0) then
      MTRES = 0
      do ISYMRS=1,NSYM
        MI = NISH(ISYMRS)
        MA = NASH(ISYMRS)
        do IV=1,MA
          IR = IV+MI
          do IX=1,IV
            IS = IX+MI
            call EXCH(ISYMPQ,ISYMRS,ISYMPQ,ISYMRS,IR,IS,INTBUF,INTBUF(LSCR))
            IDDVX = MTRES+(IV*(IV-1))/2+IX
            DDVX = DDTR(IDDVX)
            if (IR == IS) DDVX = Half*DDVX
            XMAT(NOSQES+1:NOSQES+NO**2) = XMAT(NOSQES+1:NOSQES+NO**2)+DDVX*INTBUF(1:NO**2)
          end do
        end do
        MTRES = MTRES+(MA*(MA+1))/2
      end do
      do IP=2,NO
        do IQ=1,IP-1
          LXPQ = NOSQES+IP+NO*(IQ-1)
          LXQP = NOSQES+IQ+NO*(IP-1)
          VAL = Half*(XMAT(LXPQ)+XMAT(LXQP))
          XMAT(LXPQ) = VAL
          XMAT(LXQP) = VAL
        end do
      end do
      NOSQES = NOSQES+NO**2
    end if
  end do
end if

! Determine the correction to the Fock matrix
select case (FOCKTYPE)

  case ('G1')
    ! Focktype=g1 case. A very long IF block.
    NOSQES = 0
    NOTRES = 0
    NASQES = 0
    do ISYMPQ=1,NSYM
      NI = NISH(ISYMPQ)
      NA = NASH(ISYMPQ)
      NS = NSSH(ISYMPQ)
      NO = NORB(ISYMPQ)
      NIA = NI+NA
      NAS = NA+NS
      if (NIA*NAS > 0) then

        ! the active-inactive block
        if (NA*NI > 0) then
          call DGEMM_('N','N',NA,NI,NA,One,TWOMDSQ(1+NASQES),NA,XMAT(1+NOSQES+NI),NO,Zero,SC,NA)
          do IT=1,NA
            ITTOT = NI+IT
            do II=1,NI
              KFIFA = NOTRES+(ITTOT*(ITTOT-1))/2+II
              ISC = IT+NA*(II-1)
              FIFA(KFIFA) = FIFA(KFIFA)-SC(ISC)
            end do
          end do
        end if

        ! the secondary-inactive block
        if (NS*NI > 0) then
          do IA=1,NS
            IATOT = NI+NA+IA
            do II=1,NI
              KFIFA = NOTRES+(IATOT*(IATOT-1))/2+II
              LXAI = NOSQES+IATOT+NO*(II-1)
              FIFA(KFIFA) = FIFA(KFIFA)-Two*XMAT(LXAI)
            end do
          end do
        end if

        ! the active-active block
        if (NA > 0) then
          IX = NOSQES+NI+1+NO*NI
          call DGEMM_('N','N',NA,NA,NA,One,DSQ(1+NASQES),NA,XMAT(IX),NO,Zero,SC(LSC+NA*NA),NA)
          call DGEMM_('N','N',NA,NA,NA,One,SC(LSC+NA*NA),NA,TWOMDSQ(1+NASQES),NA,Zero,SC,NA)
          do IT=1,NA
            ITTOT = NI+IT
            do IU=1,IT
              IUTOT = NI+IU
              KFIFA = NOTRES+(ITTOT*(ITTOT-1))/2+IUTOT
              VALTU = SC(IT+NA*(IU-1))
              VALUT = SC(IU+NA*(IT-1))
              FIFA(KFIFA) = FIFA(KFIFA)-Half*(VALTU+VALUT)
            end do
          end do
        end if

        ! the secondary-active block
        if (NS*NA > 0) then
          IX = NOSQES+NI+NA+1+NO*NI
          call DGEMM_('N','N',NS,NA,NA,One,XMAT(IX),NO,DSQ(1+NASQES),NA,Zero,SC,NS)
          do IA=1,NS
            IATOT = NI+NA+IA
            do IT=1,NA
              ITTOT = NI+IT
              KFIFA = NOTRES+(IATOT*(IATOT-1))/2+ITTOT
              FIFA(KFIFA) = FIFA(KFIFA)-SC(IA+NS*(IT-1))
            end do
          end do
        end if

      end if
      NOSQES = NOSQES+NO**2
      NOTRES = NOTRES+(NO*(NO+1))/2
      NASQES = NASQES+NA**2
    end do
    ! Focktype=g1 case ends.
  case ('G2','G3')
    ! Focktype=g2 or g3
    NOSQES = 0
    NOTRES = 0
    NASQES = 0
    NATRES = 0
    do ISYMPQ=1,NSYM
      NI = NISH(ISYMPQ)
      NA = NASH(ISYMPQ)
      NO = NORB(ISYMPQ)
      NIA = NI+NA
      NA2 = NA*NA
      NA3 = (NA+NA2)/2
      if (NA > 0) then

        ! Determine the matrix blocks of the correction to
        ! the Fock matrix and add them to the Fock matrix

        ! Form the selection matrix as a temporary square matrix.
        ! Compute it by spectral resolution.
        ! First, form a copy of the triangular D(2I-D) matrix block,
        ! and diagonalize it. The DDTR copy at LSC:
        SC(1:NA3) = DDTR(NATRES+1:NATRES+NA3)
        ! A unit matrix at LEV1, to become eigenvectors:
        LEV1 = LSC+NA2
        call unitmat(SC(LEV1),NA)
        ! A call to NIDiag diagonalizes the triangular matrix:
        call NIDiag(SC,SC(LEV1),NA,NA)
        call JACORD(SC,SC(LEV1),NA,NA)
        ! Make a copy of the eigenvector matrix:
        LEV2 = LEV1+NA2
        SC(LEV2:LEV2+NA2-1) = SC(LEV1:LEV1+NA2-1)
        ! Put eigenvalues at LEIG:
        LEIG = LEV2+NA2
        call VEIG(NA,SC,SC(LEIG))
        ! Now scale the second array of eigenvectors with any required
        ! function of the eigenvalues:
        do J=1,NA
          EIGVAL = SC(LEIG-1+J)
          if (FOCKTYPE == 'G2') then
            X = sqrt(max(Zero,EIGVAL))
          else
            X = EIGVAL
          end if
          do I=1,NA
            SC(LEV2-1+I+NA*(J-1)) = X*SC(LEV2-1+I+NA*(J-1))
          end do
        end do
        ! Now the selection matrix can be formed, at LSC:
        call DGEMM_('N','T',NA,NA,NA,One,SC(LEV1),NA,SC(LEV2),NA,Zero,SC(LSC),NA)
        ! Obviously, the FOCKTYPE=G3 case can be obtained by just
        ! squaring the DDTR block into SC.

        ! Focktype=g2 or g3
        IX = NOSQES+NI+1+NO*NI
        LSC1 = LSC+NA2
        LSC2 = LSC1+NA2
        call DGEMM_('N','N',NA,NA,NA,One,SC(LSC),NA,XMAT(IX),NO,Zero,SC(LSC1),NA)
        call DGEMM_('N','N',NA,NA,NA,One,SC(LSC1),NA,SC(LSC),NA,Zero,SC(LSC2),NA)
        do IT=1,NA
          ITTOT = NI+IT
          do IU=1,IT
            IUTOT = NI+IU
            KFIFA = NOTRES+(ITTOT*(ITTOT-1))/2+IUTOT
            ITU = IT+NA*(IU-1)
            FIFA(KFIFA) = FIFA(KFIFA)-SC(LSC2-1+ITU)
          end do
        end do

      end if
      NOSQES = NOSQES+NO**2
      NOTRES = NOTRES+(NO*(NO+1))/2
      NASQES = NASQES+NA**2
      NATRES = NATRES+(NA*(NA+1))/2
    end do
end select

call mma_deallocate(SC)
call mma_deallocate(INTBUF)
call mma_deallocate(DSQ)
call mma_deallocate(TWOMDSQ)
call mma_deallocate(DD)
call mma_deallocate(DDTR)
call mma_deallocate(XMAT)

end subroutine NEWFOCK
