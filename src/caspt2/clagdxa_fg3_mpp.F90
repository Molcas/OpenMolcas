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

#ifdef _MOLCAS_MPP_
subroutine CLagDXA_FG3_MPP(ISYM,NASHT,NG3,lg_BDER,lg_SDER,DG1,DG2,DG3,DF1,DF2,DF3,DEPSA,G2,idxG3)

use Symmetry_Info, only: Mul
use SUPERINDEX, only: KTUV
use definitions, only: iwp, RtoB, wp, byte
use stdalloc, only: mma_allocate, mma_deallocate, mma_MaxDBLE
use Para_Info, only: Is_Real_Par, nProcs
use caspt2_module, only: IASYM, EPSA, NTUVES
use Constants, only: Zero

implicit none
#include "global.fh"
#include "mafdecls.fh"
integer(kind=iwp), intent(in) :: ISYM, NASHT, NG3, lg_BDER, lg_SDER
real(kind=wp), intent(inout) :: DG1(NASHT,NASHT), DG2(NASHT,NASHT,NASHT,NASHT), DG3(NG3), DF1(NASHT,NASHT), &
                                DF2(NASHT,NASHT,NASHT,NASHT), DF3(NG3), DEPSA(NASHT,NASHT)
real(kind=wp), intent(in) :: G2(NASHT,NASHT,NASHT,NASHT)
integer(kind=byte), intent(in) :: idxG3(6,NG3)
integer(kind=iwp), allocatable :: INDI(:), INDJ(:), NELBsav(:), NELSsav(:)
real(kind=wp), allocatable :: BUFFB(:), BUFFS(:)
integer(kind=iwp) :: NG3MAX, MAXMEM, iscal, MAXBUF, NG3B, NBUF, NBLOCKS, IBLOCK, IG3STA, IG3END, NtotELB, NtotELS, iG3, NELB, &
                     NELS, iT, iU, iV, iX, iY, iZ, iST, iSU, iSV, iSX, iSY, iSZ, ituvs, ixyzs, iTU, iVX, iYZ, jSym, IROW, ICOL, i, &
                     IW
real(kind=wp) :: G3VAL, F3VAL

! Since we are stuck with collective calls to MPI_Alltoallv in
! order to gather the elements, each process needs to loop over
! the same number of blocks.
NG3MAX = NG3
call GAIGOP_SCAL(NG3MAX,'max')
if (NG3MAX == 0) return

! The global SC matrix has already been allocated, so we need to
! find out how much memory is left for buffering (4 equally sized
! buffers for sending and receiving values and indices)
call mma_MaxDBLE(MAXMEM)
! we need two real and two integer values per element
iscal = (iwp*2+wp*2)/RtoB
!MAXBUF=MIN(NINT(0.95D0*MAXMEM)/4,2000000000/8)
MAXBUF = min(nint(0.95_wp*MAXMEM,kind=iwp)/iscal,2000000000/8)
MAXBUF = MAXBUF-2*NG3 !! for NELBsav and NELSsav

! Loop over blocks NG3B of NG3, so that 12*NG3B < MAXBUF/NPROCS.
! This guarantees that e.g. if all processes send all their data
! to one other, that process receives NPROCS*NG3B*12 elements
! in the receive buffer.
NG3B = MAXBUF/(NPROCS*16) !! originally 12, but not sure why
NG3B = min(NG3B,NG3MAX)
call GAIGOP_SCAL(NG3B,'min')
NBUF = 16*NG3B

call mma_allocate(BUFFB,NBUF,Label='BUFFB')
call mma_allocate(BUFFS,NBUF,Label='BUFFS')
call mma_allocate(INDI,NBUF,Label='INDI')
call mma_allocate(INDJ,NBUF,Label='INDJ')

NBLOCKS = (NG3MAX-1)/NG3B+1
do IBLOCK=1,NBLOCKS
  IG3STA = 1+(IBLOCK-1)*NG3B
  IG3END = min(IG3STA+NG3B-1,NG3)
  call mma_allocate(NELBsav,IG3END-IG3STA+1,Label='NELBsav')
  call mma_allocate(NELSsav,IG3END-IG3STA+1,Label='NELSsav')

  ! Second pass fills the buffers with values and indices
  NtotELB = 0 ! Number of total elements of BDER
  NtotELS = 0 ! Number of total elements of SDER
  do iG3=IG3STA,IG3END
    NELB = 0 ! Number of elements of BDER for iG3
    NELS = 0 ! Number of elements of SDER for iG3
    iT = idxG3(1,iG3)
    iU = idxG3(2,iG3)
    iV = idxG3(3,iG3)
    iX = idxG3(4,iG3)
    iY = idxG3(5,iG3)
    iZ = idxG3(6,iG3)
    iST = IASYM(iT)
    iSU = IASYM(iU)
    iSV = IASYM(iV)
    iSX = IASYM(iX)
    iSY = IASYM(iY)
    iSZ = IASYM(iZ)
    ituvs = Mul(IST,Mul(ISU,ISV))
    ixyzs = Mul(ISX,Mul(ISY,ISZ))
    if (ituvs == ixyzs) then
      iTU = iT+NASHT*(iU-1)
      iVX = iV+NASHT*(iX-1)
      iYZ = iY+NASHT*(iZ-1)
      !F3VAL = F3(iG3)
      !-SVC20100829: 12 equivalent cases, of which the second
      !  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
      !  - G(tuvxyz) -> SA(xut,vyz)
      jSYM = Mul(IASYM(iX),Mul(IASYM(iU),IASYM(iT)))
      if (jSYM == iSYM) then
        IROW = KTUV(iX,iU,iT)-nTUVES(jSYM)
        ICOL = KTUV(iV,iY,iZ)-nTUVES(jSYM)
        NELB = NELB+1
        NtotELB = NtotELB+1
        INDI(NtotELB) = IROW
        INDJ(NtotELB) = ICOL
      end if
      if ((iTU /= iVX) .or. (iVX /= iYZ)) then
        if ((iTU /= iVX) .and. (iTU /= iYZ) .and. (iVX /= iYZ)) then
          !  - G(vxtuyz) -> SA(uxv,tyz)
          jSYM = Mul(IASYM(iU),Mul(IASYM(iX),IASYM(iV)))
          if (jSYM == iSYM) then
            IROW = KTUV(iU,iX,iV)-nTUVES(jSYM)
            ICOL = KTUV(iT,iY,iZ)-nTUVES(jSYM)
            NELB = NELB+1
            NtotELB = NtotELB+1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          end if
          !  - G(yzvxtu) -> SA(xzy,vtu)
          jSYM = Mul(IASYM(iX),Mul(IASYM(iZ),IASYM(iY)))
          if (jSYM == iSYM) then
            IROW = KTUV(iX,iZ,iY)-nTUVES(jSYM)
            ICOL = KTUV(iV,iT,iU)-nTUVES(jSYM)
            NELB = NELB+1
            NtotELB = NtotELB+1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          end if
          !  - G(tuyzvx) -> SA(zut,yvx)
          jSYM = Mul(IASYM(iZ),Mul(IASYM(iU),IASYM(iT)))
          if (jSYM == iSYM) then
            IROW = KTUV(iZ,iU,iT)-nTUVES(jSYM)
            ICOL = KTUV(iY,iV,iX)-nTUVES(jSYM)
            NELB = NELB+1
            NtotELB = NtotELB+1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          end if
        end if
        !  - G(yztuvx) -> SA(uzy,tvx)
        jSYM = Mul(IASYM(iU),Mul(IASYM(iZ),IASYM(iY)))
        if (jSYM == iSYM) then
          IROW = KTUV(iU,iZ,iY)-nTUVES(jSYM)
          ICOL = KTUV(iT,iV,iX)-nTUVES(jSYM)
          NELB = NELB+1
          NtotELB = NtotELB+1
          INDI(NtotELB) = IROW
          INDJ(NtotELB) = ICOL
        end if
        !  - G(vxyztu) -> SA(zxv,ytu)
        jSYM = Mul(IASYM(iZ),Mul(IASYM(iX),IASYM(iV)))
        if (jSYM == iSYM) then
          IROW = KTUV(iZ,iX,iV)-nTUVES(jSYM)
          ICOL = KTUV(iY,iT,iU)-nTUVES(jSYM)
          NELB = NELB+1
          NtotELB = NtotELB+1
          INDI(NtotELB) = IROW
          INDJ(NtotELB) = ICOL
        end if
      end if
      if (((iT /= iU) .or. (iV /= iX) .or. (iY /= iZ)) .and. ((iT /= iU) .or. (iV /= iZ) .or. (iX /= iY)) .and. &
          ((iX /= iV) .or. (iT /= iZ) .or. (iU /= iY)) .and. ((iZ /= iY) .or. (iV /= iU) .or. (iX /= iT))) then
        !  - G(utxvzy) -> SA(vtu,xzy)
        jSYM = Mul(IASYM(iV),Mul(IASYM(iT),IASYM(iU)))
        if (jSYM == iSYM) then
          IROW = KTUV(iV,iT,iU)-nTUVES(jSYM)
          ICOL = KTUV(iX,iZ,iY)-nTUVES(jSYM)
          NELB = NELB+1
          NtotELB = NtotELB+1
          INDI(NtotELB) = IROW
          INDJ(NtotELB) = ICOL
        end if
        if ((iTU /= iVX) .or. (iVX /= iYZ)) then
          if ((iTU /= iVX) .and. (iTU /= iYZ) .and. (iVX /= iYZ)) then
            !  - G(xvutzy) -> SA(tvx,uzy)
            jSYM = Mul(IASYM(iT),Mul(IASYM(iV),IASYM(iX)))
            if (jSYM == iSYM) then
              IROW = KTUV(iT,iV,iX)-nTUVES(jSYM)
              ICOL = KTUV(iU,iZ,iY)-nTUVES(jSYM)
              NELB = NELB+1
              NtotELB = NtotELB+1
              INDI(NtotELB) = IROW
              INDJ(NtotELB) = ICOL
            end if
            !  - G(zyxvut) -> SA(vyz,xut)
            jSYM = Mul(IASYM(iV),Mul(IASYM(iY),IASYM(iZ)))
            if (jSYM == iSYM) then
              IROW = KTUV(iV,iY,iZ)-nTUVES(jSYM)
              ICOL = KTUV(iX,iU,iT)-nTUVES(jSYM)
              NELB = NELB+1
              NtotELB = NtotELB+1
              INDI(NtotELB) = IROW
              INDJ(NtotELB) = ICOL
            end if
            !  - G(utzyxv) -> SA(ytu,zxv)
            jSYM = Mul(IASYM(iY),Mul(IASYM(iT),IASYM(iU)))
            if (jSYM == iSYM) then
              IROW = KTUV(iY,iT,iU)-nTUVES(jSYM)
              ICOL = KTUV(iZ,iX,iV)-nTUVES(jSYM)
              NELB = NELB+1
              NtotELB = NtotELB+1
              INDI(NtotELB) = IROW
              INDJ(NtotELB) = ICOL
            end if
          end if
          !  - G(zyutxv) -> SA(tyz,uxv)
          jSYM = Mul(IASYM(iT),Mul(IASYM(iY),IASYM(iZ)))
          if (jSYM == iSYM) then
            IROW = KTUV(iT,iY,iZ)-nTUVES(jSYM)
            ICOL = KTUV(iU,iX,iV)-nTUVES(jSYM)
            NELB = NELB+1
            NtotELB = NtotELB+1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          end if
          !  - G(xvzyut) -> SA(yvx,zut)
          jSYM = Mul(IASYM(iY),Mul(IASYM(iV),IASYM(iX)))
          if (jSYM == iSYM) then
            IROW = KTUV(iY,iV,iX)-nTUVES(jSYM)
            ICOL = KTUV(iZ,iU,iT)-nTUVES(jSYM)
            NELB = NELB+1
            NtotELB = NtotELB+1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          end if
        end if
      end if
    end if
    nels = nelb
    NELBsav(iG3-IG3STA+1) = NELB
    NELSsav(iG3-IG3STA+1) = NELS
  end do
  ntotels = ntotelb

  call GA_GATHER(lg_BDER,BUFFB,INDI,INDJ,NtotELB)
  call GA_GATHER(lg_SDER,BUFFS,INDI,INDJ,NtotELS)

  ! Finally, fill the local chunk of the SC matrix (block of rows)
  ! with the received values at their appropriate place.
  NtotELB = 0
  NtotELS = 0
  do iG3=IG3STA,IG3END
    iT = idxG3(1,iG3)
    iU = idxG3(2,iG3)
    iV = idxG3(3,iG3)
    iX = idxG3(4,iG3)
    iY = idxG3(5,iG3)
    iZ = idxG3(6,iG3)
    NELB = NELBsav(iG3-IG3STA+1)
    NELS = NELSsav(iG3-IG3STA+1)

    F3VAL = Zero
    do i=1,NELB
      NtotELB = NtotELB+1
      F3VAL = F3VAL-BUFFB(NtotELB)
    end do

    G3VAL = Zero
    do i=1,NELS
      NtotELS = NtotELS+1
      G3VAL = G3VAL-BUFFS(NtotELS)
    end do
    G3VAL = G3VAL-(EPSA(iU)+EPSA(iY))*F3VAL

    DF3(iG3) = DF3(iG3)+F3VAL
    DG3(iG3) = DG3(iG3)+G3VAL

    !! DEPSA is done in DF3_DEPSA_MPP

    !! remaining F3 and G3 transformation in mkfg3
    if (iY == iX) then
      DF2(iT,iU,iV,iZ) = DF2(iT,iU,iV,iZ)-F3VAL
      DG2(iT,iU,iV,iZ) = DG2(iT,iU,iV,iZ)-EPSA(iU)*F3VAL
      do iW=1,nAshT
        DEPSA(iU,iW) = DEPSA(iU,iW)-F3VAL*G2(iT,iW,iV,iZ)
      end do
      DG2(iT,iU,iV,iZ) = DG2(iT,iU,iV,iZ)-G3VAL
    end if
    if (iV == iU) then
      DF2(iT,iX,iY,iZ) = DF2(iT,iX,iY,iZ)-F3VAL
      DG2(iT,iX,iY,iZ) = DG2(iT,iX,iY,iZ)-EPSA(iY)*F3VAL
      do iW=1,nAshT
        DEPSA(iW,iY) = DEPSA(iW,iY)-F3VAL*G2(iT,iX,iW,iZ)
      end do
      DG2(iT,iX,iY,iZ) = DG2(iT,iX,iY,iZ)-G3VAL
    end if
    if (iY == iU) then
      DF2(iV,iX,iT,iZ) = DF2(iV,iX,iT,iZ)-F3VAL
      DG2(iV,iX,iT,iZ) = DG2(iV,iX,iT,iZ)-EPSA(iU)*F3VAL
      DG2(iV,iX,iT,iZ) = DG2(iV,iX,iT,iZ)-G3VAL
    end if
    DEPSA(iY,iU) = DEPSA(iY,iU)-F3VAL*G2(iV,iX,iT,iZ)
    if ((iY == iX) .and. (iV == iU)) then
      DF1(iT,iZ) = DF1(iT,iZ)-F3VAL
      DG1(iT,iZ) = DG1(iT,iZ)-G3VAL
    end if
  end do

  call mma_deallocate(NELBsav)
  call mma_deallocate(NELSsav)
end do ! end loop over blocks of G3 values

call mma_deallocate(BUFFB)
call mma_deallocate(BUFFS)
call mma_deallocate(INDI)
call mma_deallocate(INDJ)

return

end subroutine CLagDXA_FG3_MPP

#elif defined (NAGFOR)
! Some compilers do not like empty files
subroutine empty_CLagDXA_FG3_MPP()
end subroutine empty_CLagDXA_FG3_MPP
#endif
