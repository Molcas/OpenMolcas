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

subroutine TR2CTL(CMO)
! Two-electron integral transformation program: control section
!
! Purpose: Set up of memory locations (decide if out of core is
!          needed. Loop over symmetry blocks of AO integrals.
!          The transformation routine TRAMO is called for each
!          symmetry block of integrals.

use motra_global, only: Debug, FnTwoAO, FnTwoMO, IAD13, iPrint, ISP, ISQ, ISR, ISS, LMOP, LMOQ, LMOR, LMOS, LTUVX, LuTwoAO, &
                        LuTwoMO, NBP, NBPQ, NBQ, NBR, NBRS, NBS, NOP, NOQ, NOR, NOS, NOVX, nBas, nFro, nOrb, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6, RtoB

implicit none
real(kind=wp), intent(in) :: CMO(*)
#include "tratoc.fh"
integer(kind=iwp) :: IBATCH, INTBUF, IOPT, IRC, ISTBS, ISTSQ(8), ISYM, KEEP(8), KEEPP, KEEPQ, KEEPR, KEEPS, KEEPT, MEMX, NB1, NB2, &
                     NBSX(8), NORBP, NSP, NSPQ, NSPQR, NSPQRS, NSYM2, NSQ, NSR, NSS, NSSM, NW1, NW2, NW3, NW4
real(kind=wp) :: CPE, CPT, TIOE, TIOT
logical(kind=iwp) :: FoundTwoEls, DoDirect, DoCholesky, ISQUAR
integer(kind=iwp), allocatable :: iDsk(:,:)
real(kind=wp), allocatable :: W1(:), W2(:), W3(:), W4(:), W5(:)
integer(kind=iwp), external :: mma_avmem

! Set time at start of transformation

call SETTIM()

! Initiate unit LUTWOMO (two-electron integrals in MO basis)

call DANAME_MF(LUTWOMO,FNTWOMO)
IAD13 = 0
iTraToc(:) = 0
call iDAFILE(LUTWOMO,1,iTraToc,nTraToc,IAD13)

! Initiate unit LUTWOAO (two-electron integrals in AO basis)

call f_Inquire(FnTwoAo,FoundTwoEls)
call DecideOnDirect(.false.,FoundTwoEls,DoDirect,DoCholesky)
if (.not. DoCholesky) then
  IOPT = 0
  call OPNORD(IRC,IOPT,FNTWOAO,LUTWOAO)
end if

call GETORD(IRC,ISQUAR,NSYM2,NBSX,KEEP)

! Compare content of 1el and 2el integral file

if (NSYM2 /= NSYM) then
  write(u6,*) 'Tr2Ctl: NSYM2 /= NSYM'
  write(u6,*) 'NSYM2=',NSYM2
  write(u6,*) 'NSYM=',NSYM
  call Abend()
end if
do ISYM=1,NSYM
  NB1 = NBAS(ISYM)
  NB2 = NBSX(ISYM)
  if (NB1 /= NB2) then
    write(u6,*) 'Tr2Ctl: NB1 /= NB2'
    write(u6,*) 'NB1=',NB1
    write(u6,*) 'NB2=',NB2
    call Abend()
  end if
end do

! Precompute start points

ISTBS = 1
do ISYM=1,NSYM
  ISTSQ(ISYM) = ISTBS+NBAS(ISYM)*NFRO(ISYM)
  ISTBS = ISTBS+NBAS(ISYM)*NBAS(ISYM)
end do

! Loop over quadruples of symmetries (nsp,nsq,nsr,nss)
! Note that the integrals on LUTWOAO have to be sorted in the
! same order as the loop structure below.

if (iPrint >= 0) write(u6,2000)
IBATCH = 0
do NSP=1,NSYM
  NBP = NBAS(NSP)
  NOP = NORB(NSP)
  KEEPP = KEEP(NSP)
  LMOP = ISTSQ(NSP)
  ISP = NSP
  do NSQ=1,NSP
    NBQ = NBAS(NSQ)
    NOQ = NORB(NSQ)
    KEEPQ = KEEP(NSQ)
    LMOQ = ISTSQ(NSQ)
    NSPQ = ieor(NSP-1,NSQ-1)+1
    ISQ = NSQ
    do NSR=1,NSP
      NBR = NBAS(NSR)
      NOR = NORB(NSR)
      KEEPR = KEEP(NSR)
      LMOR = ISTSQ(NSR)
      NSPQR = ieor(NSPQ-1,NSR-1)+1
      ISR = NSR
      NSSM = NSR
      if (NSR == NSP) NSSM = NSQ
      do NSS=1,NSSM
        NBS = NBAS(NSS)
        NOS = NORB(NSS)
        KEEPS = KEEP(NSS)
        LMOS = ISTSQ(NSS)
        NSPQRS = ieor(NSPQR-1,NSS-1)+1
        if (NSPQRS /= 1) cycle
        IBATCH = IBATCH+1
        ISS = NSS

        ! Check the loop conditions and skip transformation step if possible

        KEEPT = KEEPP+KEEPQ+KEEPR+KEEPS
        NORBP = NOP*NOQ*NOR*NOS
        if (NORBP == 0) cycle
        if (KEEPT /= 0) then
          write(u6,*) 'Tr2Ctl: NORBP /= 0 .AND. KEEPT /= 0'
          write(u6,*) 'NORBP=',NORBP
          write(u6,*) 'KEEPT=',KEEPT
          call Abend()
        end if

        ! Define matrix sizes

        if (ISR == ISS) then
          NBPQ = (NBP+NBP**2)/2
          NBRS = (NBR+NBR**2)/2
          NOVX = (NOR+NOR**2)/2
        else
          NBPQ = NBP*NBQ
          NBRS = NBR*NBS
          NOVX = NOR*NOS
        end if

        ! Set up dynamic memory

        ! INTBUF = Size of output buffer.
        !INTBUF = 256*256
        INTBUF = 16*256*256
        NW1 = 2*nTraBuf
        NW2 = max(INTBUF,NBP*NOQ,NBQ*NOP)
        ! NW2 is size of 'X1' in TRAMO, used as LBUF in call to RDORD and RDORD_.
        ! LBUF-1 must be at least 'klB':
        NW2 = max(NW2,NBRS+1,NBPQ+1)
        NW3 = max(NBR**2,NBP**2,NOQ*NBP,NOVX)
        NW4 = max(NBR*NOS,NBQ*NOP)
        call mma_allocate(W1,NW1,label='OUTBUF')
        call mma_allocate(W2,NW2,label='X1')
        call mma_allocate(W3,NW3,label='X2')
        call mma_allocate(W4,NW4,label='X3')
        call mma_allocate(iDsk,3,nOVX,label='iDsk')
        MEMX = int(mma_avmem()*0.9_wp,kind=iwp)/RtoB

        if (DoCholesky) then ! save some space for GenInt
          MEMX = max(MEMX-MEMX/10,0)
        end if

        call mma_allocate(W5,MEMX,label='VXPQ')

        if (MEMX < NOVX) then
          write(u6,*) 'Tr2Ctl: MEMX < NOVX'
          write(u6,*) 'MEMX=',MEMX
          write(u6,*) 'NOVX=',NOVX
          call Abend()
        end if

        ! transform the symmetry block (ISP,ISQ|ISR,ISS)

        iTraToc(IBATCH) = IAD13
        ! NW2 is size of 'X1' in TRAMO, used as LBUF in call to RDORD and RDORD_.
        call TRAMO(NW2,W1,nW1,W2,nW2,W3,nW3,W4,nW4,W5,MEMX,CMO,iDsk,nOVX)
        call TIMING(CPT,CPE,TIOT,TIOE)
        if (iPrint >= 0) write(u6,2100) ISP,ISQ,ISR,ISS,NBP,NBQ,NBR,NBS,NOP,NOQ,NOR,NOS,LTUVX,CPE,TIOE
        call Xflush(u6)

        ! Deallocate work space

        call mma_deallocate(W1)
        call mma_deallocate(W2)
        call mma_deallocate(W3)
        call mma_deallocate(W4)
        call mma_deallocate(iDsk)
        call mma_deallocate(W5)

        ! End of loop over quadruples of symmetries

      end do
    end do
  end do
end do
call TIMING(CPT,CPE,TIOT,TIOE)
if (iPrint >= 0) write(u6,2200) CPT,TIOT

! Close LUTWOAO

if (.not. DoCholesky) then
  call CLSORD(IRC)
end if

! Write address block to LUTWOMO

IAD13 = 0
call iDAFILE(LUTWOMO,1,iTraToc,nTraToc,IAD13)
if ((IPRINT >= 5) .or. (DEBUG /= 0)) then
  write(u6,'(6X,A)') 'DISK ADRESSES FOR SYMMETRY BLOCKS OF TRANSFORMED TWO-ELECTRON INTEGRALS'
  write(u6,'(6X,10I8)') iTraToc(:)
end if

call DACLOS(LUTWOMO)

return

2000 format(/7X,'SYMMETRY',2X,'BASIS FUNCTIONS',6X,' ORBITALS',6X,'INTEGRALS   CPU(SEC)  I/O(SEC)')
2100 format(7X,4I2,1X,4I4,2X,4I4,3X,I9,F11.2,F10.2)
2200 format(/6X,' TOTAL CPU TIME(SEC)',F8.2,'TOTAL I/O TIME(SEC)',F8.2)

end subroutine TR2CTL
