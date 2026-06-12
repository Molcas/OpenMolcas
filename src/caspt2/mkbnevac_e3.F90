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
! Copyright (C) 2026, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine MKBNEVAC_E3(nAshT,NG3,Hact,Htilde,Gact,G1,G2,G3,idxG3)

use Index_Functions, only: nTri_Elem
use caspt2_global, only: LUSBT
use caspt2_module, only: NINDEP, NSYM, NTUV
use fake_GA, only: Allocate_GA_Array, Deallocate_GA_Array, GA_Arrays
use SC_NEVPT2, only: IDBMAT_NEVPT2
#ifdef _MOLCAS_MPP_
use Index_Functions, only: iTri
use Para_Info, only: Is_Real_Par
use GA_Wrapper, only: DBL_MB, GA_NodeId
use Definitions, only: u6
#endif
use Definitions, only: wp, iwp, byte

implicit none
integer(kind=iwp), intent(in) :: nAshT, NG3
real(kind=wp), intent(in) :: Hact(nAshT,nAshT), Htilde(nAshT,nAshT), Gact(nAshT,nAshT,nAshT,nAshT), G1(nAshT,nAshT), &
                             G2(nAshT,nAshT,nAshT,nAshT), G3(NG3)
integer(kind=byte), intent(in) :: idxG3(6,NG3)
integer(kind=iwp) :: ICASE1, ICASE4, idisk, ISYM, lg_BA, lg_BC, NASA, NASC, NBA, NBC, NINA, NINC
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: IHIA, IHIC, ii, ILOA, ILOC, JHIA, JHIC, jj, JLOA, JLOC, LDA, LDC, lg_BASQ, lg_BCSQ, MA, MC, MYRANK
#endif

ICASE1 = 1
ICASE4 = 4
do ISYM=1,NSYM
  NINA = NINDEP(ISYM,ICASE1)
  NINC = NINDEP(ISYM,ICASE4)
  if ((NINA == 0) .and. (NINC == 0)) cycle
  NASA = NTUV(ISYM)
  NASC = NTUV(ISYM)
  NBA = nTri_Elem(NASA)
  NBC = nTri_Elem(NASC)
  if ((NBA <= 0) .and. (NBC <= 0)) cycle

    !! Here uses triangular arrays
  lg_BA = Allocate_GA_Array(NBA,'BA')
  lg_BC = Allocate_GA_Array(NBC,'BC')

  ! fill in the 3-el parts
  call MKBNEVAC_E3x(ISYM,nAshT,NG3,NBA,NBC,GA_Arrays(lg_BA)%A(:),GA_Arrays(lg_BC)%A(:),G1,G2,G3,Hact,Htilde,Gact,idxG3)

  ! this is done intentionally before GADGOP
  ! summation will be done in MKBNEVAC_E4
  idisk = IDBMAT_NEVPT2(iSym,iCase1,1)
  call DDAFILE(LUSBT,1,GA_Arrays(lg_BA)%A(1),NBA,IDISK)
  idisk = IDBMAT_NEVPT2(iSym,iCase4,1)
  call DDAFILE(LUSBT,1,GA_Arrays(lg_BC)%A(1),NBC,IDISK)

# ifdef _MOLCAS_MPP_
  if (is_real_par()) then
    call GADGOP(GA_Arrays(lg_BA)%A(:),NBA,'+')
    call GADGOP(GA_Arrays(lg_BC)%A(:),NBC,'+')

    MYRANK = GA_NODEID()
    ! Put the partial contribution to the squared matrix
    call PSBMAT_GETMEM('BA',lg_BASQ,NASA)
    call GA_DISTRIBUTION(LG_BASQ,MYRANK,ILOA,IHIA,JLOA,JHIA)
    if ((JLOA /= 0) .and. (JHIA-JLOA+1 /= NASA)) then
      write(u6,*) 'MKBA: MISMATCH IN RANGE OF THE SUPERINDICES'
      call ABEND()
    end if
    if ((ILOA > 0) .and. (JLOA > 0)) then
      call GA_ACCESS(LG_BASQ,ILOA,IHIA,JLOA,JHIA,MA,LDA)
      do ii=iloa,ihia
        do jj=jloa,jhia
          DBL_MB(MA+ii-iloa+LDA*(jj-jloa)) = GA_Arrays(lg_BA)%A(iTri(ii,jj))
        end do
      end do
      call GA_RELEASE_UPDATE(LG_BASQ,ILOA,IHIA,JLOA,JHIA)
    end if
    call PSBMAT_WRITE('B',iCase1,iSYM,lg_BASQ,NASA)
    call PSBMAT_FREEMEM(lg_BASQ)

    call PSBMAT_GETMEM('BC',lg_BCSQ,NASC)
    call GA_DISTRIBUTION(LG_BCSQ,MYRANK,ILOC,IHIC,JLOC,JHIC)
    if ((JLOC /= 0) .and. (JHIC-JLOC+1 /= NASC)) then
      write(u6,*) 'MKBC: MISMATCH IN RANGE OF THE SUPERINDICES'
      call ABEND()
    end if
    if ((ILOC > 0) .and. (JLOC > 0)) then
      call GA_ACCESS(LG_BCSQ,ILOC,IHIC,JLOC,JHIC,MC,LDC)
      do ii=iloc,ihic
        do jj=jloc,jhic
          DBL_MB(MC+ii-iloc+LDC*(jj-jloc)) = GA_Arrays(lg_BC)%A(iTri(ii,jj))
        end do
      end do
      call GA_RELEASE_UPDATE(LG_BCSQ,ILOC,IHIC,JLOC,JHIC)
    end if
    call PSBMAT_WRITE('B',iCase4,iSYM,lg_BCSQ,NASC)
    call PSBMAT_FREEMEM(lg_BCSQ)
  else
# endif
    call PSBMAT_WRITE('B',iCase1,iSYM,lg_BA,NASA)
    call PSBMAT_WRITE('B',iCase4,iSYM,lg_BC,NASC)
# ifdef _MOLCAS_MPP_
  end if
# endif
  call Deallocate_GA_Array(lg_BA)
  call Deallocate_GA_Array(lg_BC)
end do

end subroutine MKBNEVAC_E3
