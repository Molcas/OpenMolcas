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

subroutine CHO_MCA_INT_1_DBG2()
!
! Purpose: test symmetry of integral matrix.

use Symmetry_Info, only: Mul
use Index_Functions, only: iTri, nTri_Elem
use Cholesky, only: IFCSew, iSP2F, LuPri, Mx2Sh, nBstSh, nnShl
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: IA, IAB, IABCD, IABMN, IABMX, IB, IC, ICD, ICDMN, ICDMX, ID, IERR, ISHLA, ISHLAB, ISHLB, ISHLC, ISHLCD, &
                     ISHLD, ISYMA, ISYMAB, ISYMB, ISYMC, ISYMCD, ISYMD, ITST, LINT, LINTT, LSEW, NERR, NTST, NUMAB, NUMCD
real(kind=wp) :: ERRMAX, ERRMIN, TST
real(kind=wp), allocatable, target :: INT1(:)
real(kind=wp), pointer :: pINT1(:), pINT2(:)
real(kind=wp), parameter :: THR = 1.0e-14_wp
logical(kind=iwp), parameter :: PRTINT = .false.
character(len=*), parameter :: SECNAM = 'CHO_MCA_INT_1_DBG2'
integer(kind=iwp), external :: CHO_ISAOSH

write(LUPRI,*)
write(LUPRI,*)
write(LUPRI,*) SECNAM,': testing integral matrix symmetry'
write(LUPRI,*)

! Force computation of full shell quadruple.
! ------------------------------------------

if (IFCSEW /= 1) then
  write(LUPRI,*) SECNAM,': WARNING: resetting IFCSEW from ',IFCSEW,' to 1.'
  IFCSEW = 1
end if

LINTT = 2*MX2SH*MX2SH
call mma_allocate(INT1,LINTT,Label='INT1')
call mma_maxDBLE(LSEW)
call XSETMEM_INTS(LSEW)

NTST = 0
NERR = 0
do ISHLAB=1,NNSHL

  call CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.true.)
  if (ISHLB == ISHLA) then
    NUMAB = nTri_Elem(NBSTSH(ISHLA))
  else
    NUMAB = NBSTSH(ISHLA)*NBSTSH(ISHLB)
  end if

  do ISHLCD=1,ISHLAB

    call CHO_INVPCK(ISP2F(ISHLCD),ISHLC,ISHLD,.true.)
    if (ISHLD == ISHLC) then
      NUMCD = nTri_Elem(NBSTSH(ISHLC))
    else
      NUMCD = NBSTSH(ISHLC)*NBSTSH(ISHLD)
    end if
    LINT = NUMCD*NUMAB
    pINT1(1:LINT) => INT1(1:LINT)
    pINT2(1:LINT) => INT1(LINT+1:LINT+LINT)

    ! Calculate integrals (CD|AB).
    ! ----------------------------

    pINT1(:) = Zero
    call CHO_MCA_INT_1(ISHLCD,ISHLAB,pINT1,LINT,PRTINT)

    ! Calculate integrals (AB|CD).
    ! ----------------------------

    pINT2(:) = Zero
    call CHO_MCA_INT_1(ISHLAB,ISHLCD,pINT2,LINT,PRTINT)

    ! Compare.
    ! --------

    ITST = 0
    IERR = 0
    call CHO_MCA_INT1_1_DBG2_CMP(pINT1,pINT2,NUMCD,NUMAB,ERRMIN,ICDMN,IABMN,ERRMAX,ICDMX,IABMX,ITST,IERR,THR,.false.,LUPRI)
    NTST = NTST+ITST
    NERR = NERR+IERR

    write(LUPRI,*) '#sym. errors for (',ISHLC,ISHLD,'|',ISHLA,ISHLB,'): ',IERR,' #tests: ',ITST
    if (IERR /= 0) then
      !write(LUPRI,*) '    Here is the shell quadruple in INT1:'
      !call CHO_OUTPUT(pINT1,1,NUMCD,1,NUMAB,NUMCD,NUMAB,1,LUPRI)
      !write(LUPRI,*) '    And the shell quadruple in INT2:'
      !call CHO_OUTPUT(pINT2,1,NUMAB,1,NUMCD,NUMAB,NUMCD,1,LUPRI)
      do IB=1,NBSTSH(ISHLB)
        ISYMB = CHO_ISAOSH(IB,ISHLB)
        do IA=1,NBSTSH(ISHLA)
          ISYMA = CHO_ISAOSH(IA,ISHLA)
          ISYMAB = MUL(ISYMA,ISYMB)
          if (ISHLB == ISHLA) then
            IAB = ITRI(IA,IB)
          else
            IAB = NBSTSH(ISHLA)*(IB-1)+IA
          end if
          do ID=1,NBSTSH(ISHLD)
            ISYMD = CHO_ISAOSH(ID,ISHLD)
            do IC=1,NBSTSH(ISHLC)
              ISYMC = CHO_ISAOSH(IC,ISHLC)
              ISYMCD = MUL(ISYMC,ISYMD)
              if (ISHLC == ISHLD) then
                ICD = ITRI(IC,ID)
              else
                ICD = NBSTSH(ISHLC)*(ID-1)+IC
              end if
              IABCD = NUMCD*(IAB-1)+ICD
              TST = abs(pINT1(IABCD))
              if ((TST > Zero) .and. (ISYMCD /= ISYMAB)) then
                write(LUPRI,*) 'Symmetry break!!'
                write(LUPRI,*) 'element ',ICD,IAB,' is non-zero: ',pINT1(IABCD)
                write(LUPRI,*) 'Symmetry is: ',MUL(ISYMCD,ISYMAB)
              end if
            end do
          end do
        end do
      end do
    end if

  end do
end do

call XRLSMEM_INTS()
call mma_deallocate(INT1)
nullify(pINT1)
nullify(pINT2)

write(LUPRI,*) '***END OF ',SECNAM,': #tests: ',NTST,' #errors: ',NERR

end subroutine CHO_MCA_INT_1_DBG2
