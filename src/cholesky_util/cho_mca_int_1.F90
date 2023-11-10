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

subroutine CHO_MCA_INT_1(IJ,KL,XINT,LINT,LOCPRT)
!
! Purpose: calculate shell-quadruple (IJ|KL) and return
!          them in XINT.
!
! Notes:
!    LOCPRT: flag for printing the shell quadruple to output;
!            output format differs depending on IFCSEW.

use Index_Functions, only: iTri, nTri_Elem
use Integral_Interfaces, only: Int_PostProcess, Integral_WrOut_Cho
use Cholesky, only: IFCSew, iSP2F, LuPri, nBstSh, ShA, ShAB, ShB, ShC, ShCD, ShD
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IJ, KL, LINT
real(kind=wp), intent(inout) :: XINT(LINT)
logical(kind=iwp), intent(in) :: LOCPRT
integer(kind=iwp) :: I, II, IIJJ, J, JJ, K, KK, KKLL, KOFF, L, LL, NUMI, NUMIJ, NUMJ, NUMK, NUML
character(len=*), parameter :: SECNAM = 'CHO_MCA_INT_1'

! Initializations.
! ----------------

call CHO_INVPCK(ISP2F(IJ),I,J,.true.)
call CHO_INVPCK(ISP2F(KL),K,L,.true.)

SHCD = IJ
SHAB = KL
SHC = I
SHD = J
SHA = K
SHB = L

! Calculate integrals.
! --------------------

Int_PostProcess => Integral_WrOut_Cho
call EVAL_IJKL(I,J,K,L,XINT,LINT)
Int_PostProcess => null()

! Print integrals.
! ----------------

if (LOCPRT) then

  if (IFCSEW == 1) then

    write(LUPRI,'(//,5X,A,A,4I5,A)') SECNAM,': shell quadruple ',I,J,K,L,':'

    NUMI = NBSTSH(I)
    NUMJ = NBSTSH(J)
    NUMK = NBSTSH(K)
    NUML = NBSTSH(L)
    if (I == J) then
      NUMIJ = nTri_Elem(NUMI)
    else
      NUMIJ = NUMI*NUMJ
    end if

    if (K == L) then
      do LL=1,NUML
        do KK=1,LL
          KKLL = ITRI(KK,LL)
          if (I == J) then
            do JJ=1,NUMJ
              do II=1,JJ
                IIJJ = ITRI(II,JJ)
                KOFF = NUMIJ*(KKLL-1)+IIJJ
                write(LUPRI,*) '(',I,J,K,L,') [',II,JJ,KK,LL,'] = ',XINT(KOFF)
              end do
            end do
          else
            do JJ=1,NUMJ
              do II=1,NUMI
                IIJJ = NUMI*(JJ-1)+II
                KOFF = NUMIJ*(KKLL-1)+IIJJ
                write(LUPRI,*) '(',I,J,K,L,') [',II,JJ,KK,LL,'] = ',XINT(KOFF)
              end do
            end do
          end if
        end do
      end do
    else
      do LL=1,NUML
        do KK=1,NUMK
          KKLL = NUMK*(LL-1)+KK
          if (I == J) then
            do JJ=1,NUMJ
              do II=1,JJ
                IIJJ = ITRI(II,JJ)
                KOFF = NUMIJ*(KKLL-1)+IIJJ
                write(LUPRI,*) '(',I,J,K,L,') [',II,JJ,KK,LL,'] = ',XINT(KOFF)
              end do
            end do
          else
            do JJ=1,NUMJ
              do II=1,NUMI
                IIJJ = NUMI*(JJ-1)+II
                KOFF = NUMIJ*(KKLL-1)+IIJJ
                write(LUPRI,*) '(',I,J,K,L,') [',II,JJ,KK,LL,'] = ',XINT(KOFF)
              end do
            end do
          end if
        end do
      end do
    end if

  else if ((IFCSEW == 2) .or. (IFCSEW == 3)) then

    call CHO_PRTINT(IJ,KL,XINT,LINT)

  else

    write(LUPRI,*) SECNAM,': IFCSEW=',IFCSEW
    call CHO_QUIT(SECNAM//': IFCSEW out of bounds!',103)

  end if

end if

end subroutine CHO_MCA_INT_1
