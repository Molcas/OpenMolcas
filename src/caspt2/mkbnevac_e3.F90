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
  use caspt2_global, only: LUSBT
  use caspt2_module, only: NINDEP, NSYM, NTUV
  use Constants, only: Zero
  use stdalloc, only: mma_allocate, mma_deallocate
  use definitions, only: iwp,wp,byte
  use fake_GA, only: GA_Arrays, Allocate_GA_Array, Deallocate_GA_Array
  use SC_NEVPT2, only: IDBMAT_NEVPT2
#ifdef _MOLCAS_MPP_
  USE Para_Info, ONLY: Is_Real_Par
  use GA_Wrapper, only: DBL_MB, GA_NodeId
  use definitions, only: u6
#endif

  implicit none

  integer(kind=iwp), intent(in) :: nAshT, NG3
  real(kind=wp), intent(in) :: Hact(nAshT,nAshT), Htilde(nAshT,nAshT), Gact(nAshT,nAshT,nAshT,nAshT), &
                               G1(nAshT,nAshT), G2(nAshT,nAshT,nAshT,nAshT), G3(NG3)
  integer(kind=byte), intent(in) :: idxG3(6,NG3)

  integer(kind=iwp) :: ICASE1, ICASE4, idisk, NASA, NASC, NBA, NBC, NINA, NINC, ISYM, lg_BA, lg_BC
#ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: IHIA, IHIC, ILOA, ILOC, JHIA, JHIC, JLOA, JLOC, LDA, LDC, lg_BASQ, lg_BCSQ, MA, MC, MYRANK
  integer(kind=iwp) :: ii, jj
#endif

  ICASE1=1
  ICASE4=4
  DO ISYM=1,NSYM
    NINA=NINDEP(ISYM,ICASE1)
    NINC=NINDEP(ISYM,ICASE4)
    if (NINA == 0 .and. NINC == 0) cycle
    NASA=NTUV(ISYM)
    NASC=NTUV(ISYM)
    NBA=NASA*(NASA+1)/2
    NBC=NASC*(NASC+1)/2
    if (NBA <= 0 .and. NBC <= 0) cycle

    !! Here uses triangular arrays
    lg_BA = Allocate_GA_Array(NBA,'BA')
    lg_BC = Allocate_GA_Array(NBC,'BC')

    ! fill in the 3-el parts
    CALL MKBNEVAC_E3x(ISYM,nAshT,NG3,NBA,NBC,GA_Arrays(lg_BA)%A(:),GA_Arrays(lg_BC)%A(:),G1,G2,G3,Hact,Htilde,Gact,idxG3)

    ! this is done intentionally before GADGOP
    ! summation will be done in MKBNEVAC_E4
    idisk = IDBMAT_NEVPT2(iSym,iCase1,1)
    CALL DDAFILE(LUSBT,1,GA_Arrays(lg_BA)%A(1),NBA,IDISK)
    idisk = IDBMAT_NEVPT2(iSym,iCase4,1)
    CALL DDAFILE(LUSBT,1,GA_Arrays(lg_BC)%A(1),NBC,IDISK)

#ifdef _MOLCAS_MPP_
    if (is_real_par()) then
      call GADGOP(GA_Arrays(lg_BA)%A(:),NBA,'+')
      call GADGOP(GA_Arrays(lg_BC)%A(:),NBC,'+')
    end if
#endif

#ifdef _MOLCAS_MPP_
    IF (IS_REAL_PAR()) THEN
      MYRANK = GA_NODEID()
      ! Put the partial contribution to the squared matrix
      CALL PSBMAT_GETMEM('BA',lg_BASQ,NASA)
      CALL GA_DISTRIBUTION (LG_BASQ,MYRANK,ILOA,IHIA,JLOA,JHIA)
      IF (JLOA /= 0 .AND. (JHIA-JLOA+1) /= NASA) THEN
        WRITE(u6,*) 'MKBA: MISMATCH IN RANGE OF THE SUPERINDICES'
        CALL ABEND()
      END IF
      IF (ILOA > 0 .AND. JLOA > 0) THEN
        CALL GA_ACCESS (LG_BASQ,ILOA,IHIA,JLOA,JHIA,MA,LDA)
        do ii = iloa, ihia
          do jj = jloa, jhia
            if (ii >= jj) then
              DBL_MB(MA+ii-iloa+LDA*(jj-jloa)) = GA_Arrays(lg_BA)%A(ii*(ii-1)/2+jj)
            else
              DBL_MB(MA+ii-iloa+LDA*(jj-jloa)) = GA_Arrays(lg_BA)%A(jj*(jj-1)/2+ii)
            end if
          end do
        end do
        CALL GA_RELEASE_UPDATE (LG_BASQ,ILOA,IHIA,JLOA,JHIA)
      END IF
      CALL PSBMAT_WRITE('B',iCase1,iSYM,lg_BASQ,NASA)
      CALL PSBMAT_FREEMEM(lg_BASQ)

      CALL PSBMAT_GETMEM('BC',lg_BCSQ,NASC)
      CALL GA_DISTRIBUTION (LG_BCSQ,MYRANK,ILOC,IHIC,JLOC,JHIC)
      IF (JLOC /= 0 .AND. (JHIC-JLOC+1) /= NASC) THEN
        WRITE(u6,*) 'MKBC: MISMATCH IN RANGE OF THE SUPERINDICES'
        CALL ABEND()
      END IF
      IF (ILOC > 0 .AND. JLOC > 0) THEN
        CALL GA_ACCESS (LG_BCSQ,ILOC,IHIC,JLOC,JHIC,MC,LDC)
        do ii = iloc, ihic
          do jj = jloc, jhic
            if (ii >= jj) then
              DBL_MB(MC+ii-iloc+LDC*(jj-jloc)) = GA_Arrays(lg_BC)%A(ii*(ii-1)/2+jj)
            else
              DBL_MB(MC+ii-iloc+LDC*(jj-jloc)) = GA_Arrays(lg_BC)%A(jj*(jj-1)/2+ii)
            end if
          end do
        end do
        CALL GA_RELEASE_UPDATE (LG_BCSQ,ILOC,IHIC,JLOC,JHIC)
      END IF
      CALL PSBMAT_WRITE('B',iCase4,iSYM,lg_BCSQ,NASC)
      CALL PSBMAT_FREEMEM(lg_BCSQ)
    ELSE
#endif
      CALL PSBMAT_WRITE('B',iCase1,iSYM,lg_BA,NASA)
      CALL PSBMAT_WRITE('B',iCase4,iSYM,lg_BC,NASC)
#ifdef _MOLCAS_MPP_
    END IF
#endif
    Call Deallocate_GA_Array(lg_BA)
    Call Deallocate_GA_Array(lg_BC)
  END DO

end subroutine MKBNEVAC_E3
