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

use Index_Functions, only: iTri, nTri_Elem
use SUPERINDEX, only: KTU, KTUV
use general_data, only: NACTEL, NASH
use caspt2_module, only: NAES, NASHT, NASUP, NINDEP, NISH, NISUP, NORB, NSSH, NSYM, NTUES, NTUV, NTUVES
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IVEC, NFIMO
real(kind=wp), intent(in) :: FIMO(NFIMO)
integer(kind=iwp) :: IA, IAJ, IATOT, ICASE, IFOFF, IJ, ISYJ, ISYM, ISYT, ISYU, IT, ITABS, ITTOT, IU, IUABS, IUU, IVABS, IW1, IW2, &
                     IWD, IX, IXABS, IXTOT, IYABS, IYYW, IYYWA, lg_A, lg_C, lg_D, NAJ, NAS, NAT, NAX, NIJ, NIS, NIT, NIX, NO, NSJ, &
                     NSX, NWA, NWC, NWD
real(kind=wp) :: ONEADD, rSUM, Val
real(kind=wp), allocatable :: WA(:), WC(:), WD(:)

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
        Val = FIMO(IFOFF+iTri(ITTOT,IJ))/real(max(1,NACTEL),kind=wp)
        do IVABS=1,NASHT
          IW1 = KTUV(ITABS,IVABS,IVABS)-NTUVES(ISYM)
          IW2 = IJ
          WA(IW1+NAS*(IW2-1)) = WA(IW1+NAS*(IW2-1))+Val
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
  IFOFF = IFOFF+nTri_Elem(NO)
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
        rSUM = FIMO(IFOFF+iTri(IATOT,IXTOT))
        do IYABS=1,NASHT
          IYYW = KTUV(IYABS,IYABS,IXABS)-NTUVES(ISYM)
          IYYWA = IYYW+NAS*(IA-1)
          rSUM = rSUM-WC(IYYWA)
        end do
        ONEADD = rSUM/real(max(1,NACTEL),kind=wp)
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
  IFOFF = IFOFF+nTri_Elem(NO)
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
        ONEADD = FIMO(IFOFF+iTri(IATOT,IJ))/real(max(1,NACTEL),kind=wp)
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
    IFOFF = IFOFF+nTri_Elem(NO)
  end do
  call RHS_PUT(NAS,NIS,lg_D,WD)

  ! Put W on disk:
  call RHS_SAVE(NAS,NIS,lg_D,ICASE,ISYM,IVEC)
  call RHS_FREE(lg_D)
  call mma_deallocate(WD)
end if

end subroutine MODRHS
