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

subroutine MODRHS(IVEC,FIMO,NFIMO)

use definitions, only: iwp, wp
use SUPERINDEX, only: KTUV, KTU
use stdalloc, only: mma_allocate, mma_deallocate
use caspt2_module, only: NSYM, NINDEP, NTUV, NISH, NASH, NAES, NACTEL, NASHT, NTUVES, NORB, NSSH, NISUP, NASUP, NTUES

implicit none
integer(kind=iwp), intent(in) :: IVEC, NFIMO
real(kind=wp), intent(in) :: FIMO(NFIMO)
real(kind=wp), allocatable :: WA(:), WC(:), WD(:)
integer(kind=iwp) ICASE, IFOFF, ISYM, NAS, NIS, NWA, ISYJ, ISYT, NAT, NIT, IT, ITTOT, ITABS, IJ, IVABS, IW1, IW2, IA, IAJ, IATOT, &
                  ISYU, IU, IUABS, IUU, IWD, IX, IXABS, IXTOT, IYABS, IYYW, lg_A, lg_C, lg_D, NAJ, NAX, NIJ, NIX, NO, NSJ, NSX, &
                  NWC, NWD, IYYWA
real(kind=wp) ONEADD, SUM, value

!**************************************************************
! Case A:
ICASE = 1
IFOFF = 0
do ISYM=1,NSYM
  NAS = NTUV(ISYM)
  NIS = NISH(ISYM)
  NWA = NAS*NIS
  if (NINDEP(ISYM,1)*NWA /= 0) then
    call mma_allocate(WA,NWA,Label='WA')
    call RHS_ALLO(NAS,NIS,lg_A)
    ! Read W from disk:
    call RHS_READ(NAS,NIS,lg_A,ICASE,ISYM,IVEC)
    call RHS_GET(NAS,NIS,lg_A,WA)
    ! Insert one-electron contribution to coupling <A|0>:
    ! WA(tvv,j)=FIMO(t,j)/NACTEL (+two-electron part)
    ISYJ = ISYM
    ISYT = ISYM
    NAT = NASH(ISYT)
    NIT = NISH(ISYT)
    do IT=1,NAT
      ITTOT = NIT+IT
      ITABS = NAES(ISYT)+IT
      do IJ=1,NIT
        value = FIMO(IFOFF+(ITTOT*(ITTOT-1))/2+IJ)/dble(max(1,NACTEL))
        do IVABS=1,NASHT
          IW1 = KTUV(ITABS,IVABS,IVABS)-NTUVES(ISYM)
          IW2 = IJ
          WA(IW1+NAS*(IW2-1)) = WA(IW1+NAS*(IW2-1))+value
        end do
      end do
    end do
    call RHS_PUT(NAS,NIS,lg_A,WA)
    ! Put W on disk:
    call RHS_SAVE(NAS,NIS,lg_A,ICASE,ISYM,IVEC)
    call RHS_FREE(lg_A)
    call mma_deallocate(WA)

  end if
  ! End of loop over ISYM.
  NO = NORB(ISYM)
  IFOFF = IFOFF+(NO*(NO+1))/2
end do

!**************************************************************
! Case C:
ICASE = 4
IFOFF = 0
do ISYM=1,NSYM
  NAS = NTUV(ISYM)
  NIS = NSSH(ISYM)
  NWC = NAS*NIS
  if (NINDEP(ISYM,4)*NWC /= 0) then
    call mma_allocate(WC,NWC,LABEL='WC')
    call RHS_ALLO(NAS,NIS,lg_C)
    ! Read W from disk:
    call RHS_READ(NAS,NIS,lg_C,ICASE,ISYM,IVEC)
    call RHS_GET(NAS,NIS,lg_C,WC)
    ! Insert one-electron contribution to coupling <C|0>:
    ! WC(xuu,a)=(FIMO(a,x)-sum((ay,yx), y=1,NASHT) )/NACTEL (+ two-el part)
    NIX = NISH(ISYM)
    NAX = NASH(ISYM)
    NSX = NSSH(ISYM)
    do IX=1,NAX
      IXTOT = NIX+IX
      IXABS = NAES(ISYM)+IX
      do IA=1,NSX
        IATOT = NIX+NAX+IA
        SUM = FIMO(IFOFF+(IATOT*(IATOT-1))/2+IXTOT)
        do IYABS=1,NASHT
          IYYW = KTUV(IYABS,IYABS,IXABS)-NTUVES(ISYM)
          IYYWA = IYYW+NAS*(IA-1)
          SUM = SUM-WC(IYYWA)
        end do
        ONEADD = SUM/dble(max(1,NACTEL))
        do IUABS=1,NASHT
          IW1 = KTUV(IXABS,IUABS,IUABS)-NTUVES(ISYM)
          IW2 = IA
          WC(IW1+NAS*(IW2-1)) = WC(IW1+NAS*(IW2-1))+ONEADD
        end do
      end do
    end do
    call RHS_PUT(NAS,NIS,lg_C,WC)
    ! Put W on disk:
    call RHS_SAVE(NAS,NIS,lg_C,ICASE,ISYM,IVEC)
    call RHS_FREE(lg_C)
    call mma_deallocate(WC)

  end if
  ! End of loop over ISYM.
  NO = NORB(ISYM)
  IFOFF = IFOFF+(NO*(NO+1))/2
end do

!**************************************************************
! Case D1:
ICASE = 5
ISYM = 1

NAS = NASUP(ISYM,5)
NIS = NISUP(ISYM,5)
NWD = NAS*NIS
if (NINDEP(ISYM,5)*NWD /= 0) then
  call mma_allocate(WD,NWD,LABEL='WD')
  call RHS_ALLO(NAS,NIS,lg_D)
  ! Read W from disk:
  call RHS_READ(NAS,NIS,lg_D,ICASE,ISYM,IVEC)
  call RHS_GET(NAS,NIS,lg_D,WD)

  ! Insert one-electron contribution to coupling <D1|0>:
  ! Compute WD1(vv,aj)=FIMO(a,j)/NACTEL (+ two-el part)
  IFOFF = 0
  IAJ = 0
  do ISYJ=1,NSYM
    NIJ = NISH(ISYJ)
    NAJ = NASH(ISYJ)
    NSJ = NSSH(ISYJ)
    do IA=1,NSJ
      IATOT = NIJ+NAJ+IA
      do IJ=1,NIJ
        ONEADD = FIMO(IFOFF+(IATOT*(IATOT-1))/2+IJ)/dble(max(1,NACTEL))
        IAJ = IAJ+1
        do ISYU=1,NSYM
          do IU=1,NASH(ISYU)
            IUABS = NAES(ISYU)+IU
            IUU = KTU(IUABS,IUABS)-NTUES(ISYM)
            IWD = IUU+NAS*(IAJ-1)
            WD(IWD) = WD(IWD)+ONEADD
          end do
        end do
      end do
    end do
    NO = NORB(ISYJ)
    IFOFF = IFOFF+(NO*(NO+1))/2
  end do
  call RHS_PUT(NAS,NIS,lg_D,WD)

  ! Put W on disk:
  call RHS_SAVE(NAS,NIS,lg_D,ICASE,ISYM,IVEC)
  call RHS_FREE(lg_D)
  call mma_deallocate(WD)
end if

end subroutine MODRHS
