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

subroutine HAM3(OP0,OP1,NOP2,OP2,NOP3,OP3,ISYCI,CI,SGM,NCI)
! Purpose: Compute and add a contribution to SGM which is
! obtained from a sum of zero- one- two- and three-electron
! operators acting on wave function CI.

! Note that the coefficients in OP1 and OP2 must have been
! modified by adding elements from OP2 and OP3, as done in
! subroutine MODOP.

! Presently symmetry blocking is disregarded for OP2, OP3, but
! index pair C permutation symmetry is used.
! NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
! NOP3=(NASHT**2+2 over 3)  (Binomial coefficient)

use sguga, only: sg_epq_psi
use Index_Functions, only: iTri, nTri3_Elem
use Symmetry_Info, only: Mul
use caspt2_global, only: CIS, EXS, SGS
use Molcas, only: MxLev
use caspt2_module, only: IASYM, ISCF, MxCI, NACTEL, NASHT, NCONF, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NOP2, NOP3, ISYCI, NCI
real(kind=wp), intent(in) :: OP0, OP1(NASHT,NASHT), OP2(NOP2), OP3(NOP3), CI(NCI)
real(kind=wp), intent(inout) :: SGM(NCI)
integer(kind=iwp) :: I, IATOG(MXLEV), ISTU, ISVX, ISVXYZ, ISYM, ISYM1, ISYM2, ISYZ, IT, ITABS, ITMIN, ITU, ITUVXYZ, IU, IV, IVMIN, &
                     IVX, IVXYZ, IX, IY, IYZ, IZ, LEVT, LEVU, LEVV, LEVX, LEVY, LEVZ, nLev, NSGM1, NSGM2
real(kind=wp) :: OCCNO, X
real(kind=wp), allocatable :: SGM1(:), SGM2(:)

nLev = SGS%nLev

if (NCONF == 0) return
if (abs(OP0) > 1.0e-15_wp) SGM(1:NCONF) = SGM(1:NCONF)+OP0*CI(1:NCONF)
if (NACTEL == 0) return

! Unless this is a special-case wave function, reserve space
! for intermediate results of elementary excitations.
if (ISCF == 0) then
  call mma_allocate(SGM1,MXCI,Label='SGM1')
  if (NACTEL >= 2) call mma_allocate(SGM2,MXCI,Label='SGM2')
end if
! Special cases:
OCCNO = Zero
if (ISCF == 1) OCCNO = Two
if (ISCF == 2) OCCNO = One

! Create reorder table giving the GUGA level, i.e. CI-coupling
! ordinal number of each active orbital.
ITABS = 0
do ISYM=1,NSYM
  do I=1,NLEV
    if (SGS%ISM(I) == ISYM) then
      ITABS = ITABS+1
      IATOG(ITABS) = I
    end if
  end do
end do

do IZ=1,NASHT
  do IY=1,NASHT
    IYZ = IY+(IZ-1)*NASHT
    ISYZ = Mul(IASYM(IY),IASYM(IZ))
    ISYM1 = Mul(ISYZ,ISYCI)
    NSGM1 = CIS%NCSF(ISYM1)
    if (NSGM1 == 0) cycle
    if (ISCF == 0) then
      ! The general case:
      ! Compute SGM1:=E(IY,IZ) PSI
      SGM1(1:nSGM1) = Zero
      LEVY = IATOG(IY)
      LEVZ = IATOG(IZ)
      call SG_Epq_Psi(SGS,CIS,EXS,LEVY,LEVZ,One,ISYCI,CI,SGM1)
      ! Add non-zero 1-el contribution to SGM:
      if (ISYZ == 1) then
        X = OP1(IY,IZ)
        if (abs(X) > 1.0e-15_wp) SGM(1:nConf) = SGM(1:nConf)+X*SGM1(1:nConf)
      end if
    else
      ! Closed-shell or hi-spin case:
      if (IY /= IZ) cycle
      X = OCCNO*OP1(IY,IZ)
      SGM(1) = SGM(1)+X*CI(1)
    end if
    if (NACTEL == 1) cycle
    do IX=IZ,NASHT
      IVMIN = 1
      if (IX == IZ) IVMIN = IY
      do IV=IVMIN,NASHT
        IVX = IV+(IX-1)*NASHT
        ISVX = Mul(IASYM(IV),IASYM(IX))
        ISVXYZ = Mul(ISVX,ISYZ)
        IVXYZ = iTri(IVX,IYZ)
        ISYM2 = Mul(ISVX,ISYM1)
        NSGM2 = CIS%NCSF(ISYM2)
        if (NSGM2 == 0) cycle
        if (ISCF == 0) then
          ! The general case:
          ! Compute SGM2:=E(IV,IX) SGM1
          SGM2(1:nSGM2) = Zero
          LEVV = IATOG(IV)
          LEVX = IATOG(IX)
          call SG_Epq_Psi(SGS,CIS,EXS,LEVV,LEVX,One,ISYM1,SGM1,SGM2)
          ! Add non-zero 2-el contribution to SGM:
          if (ISVXYZ == 1) then
            X = OP2(IVXYZ)
            if (abs(X) > 1.0e-15_wp) SGM(1:nConf) = SGM(1:nConf)+X*SGM2(1:nConf)
          end if
        else
          ! Closed-shell or hi-spin case:
          if (IY /= IZ) cycle
          if (IV /= IX) cycle
          X = (OCCNO**2)*OP2(IVXYZ)
          SGM(1) = SGM(1)+X*CI(1)
        end if
        if (NACTEL == 2) cycle
        do IU=IX,NASHT
          ITMIN = 1
          if (IU == IX) ITMIN = IV
          do IT=ITMIN,NASHT
            ITU = IT+(IU-1)*NASHT
            ISTU = Mul(IASYM(IT),IASYM(IU))
            if (ISTU /= ISVXYZ) cycle
            ITUVXYZ = nTri3_Elem(ITU-1)+IVXYZ
            X = OP3(ITUVXYZ)
            if (abs(X) < 1.0e-15_wp) cycle
            ! Add non-zero 3-el contribution to SGM:
            if (ISCF == 0) then
              LEVT = IATOG(IT)
              LEVU = IATOG(IU)
              call SG_Epq_Psi(SGS,CIS,EXS,LEVT,LEVU,X,ISYM2,SGM2,SGM)
            else
              ! Closed-shell or hi-spin case:
              if (IT /= IU) cycle
              if (IV /= IX) cycle
              if (IY /= IZ) cycle
              X = (OCCNO**3)*OP3(ITUVXYZ)
              SGM(1) = SGM(1)+X*CI(1)
            end if
          end do
        end do
      end do
    end do
  end do
end do

! Deallocate temporary arrays, if any:
if (ISCF == 0) then
  call mma_deallocate(SGM1)
  if (NACTEL >= 2) call mma_deallocate(SGM2)
end if

end subroutine HAM3
