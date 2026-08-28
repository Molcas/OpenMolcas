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

subroutine Process_RHS_Block(ITI,ITP,ITK,ITQ,nCase,Cho_Bra,nBra,Cho_Ket,nKet,nSh,JSYM,IVEC,NV)

use Symmetry_Info, only: Mul
use PrintLevel, only: DEBUG
use AddRHS, only: ADDRHSA, ADDRHSB, ADDRHSC, ADDRHSD1, ADDRHSD2, ADDRHSE, ADDRHSF, ADDRHSG, ADDRHSH
#ifdef _MOLCAS_MPP_
use ADDRHS_STRIPED, only: ADDRHSA_STRIPED, ADDRHSB_STRIPED, ADDRHSC_STRIPED, ADDRHSD1_STRIPED, ADDRHSD2_STRIPED, &
                          ADDRHSE_STRIPED, ADDRHSF_STRIPED, ADDRHSG_STRIPED, ADDRHSH_STRIPED
#endif
use caspt2_global, only: BUFF, idxb, iParRHS, iPrGlb, PIQK
use caspt2_module, only: NSYM
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ITI, ITP, ITK, ITQ, nBra, nKet, nSh(8,3), JSYM, iVec, nV
character(len=2), intent(in) :: nCase
real(kind=wp), intent(in) :: Cho_Bra(nBra), Cho_Ket(nKet)
integer(kind=iwp) :: ISYI, ISYK, ISYP, ISYQ, KPI, KQK, LBRASM, LKETSM, mxPIQK, NBRASM, nBuff, NI, NK, NKETSM, NP, NPI, NPIQK, NQ, &
                     NQK

mxPIQK = size(PIQK)
nBuff = size(BUFF)

if (iPrGlb >= DEBUG) write(u6,*) 'Processing RHS block '//nCase

LBRASM = 1
do ISYI=1,NSYM
  NI = NSH(ISYI,ITI)
  if (NI == 0) cycle
  ISYP = Mul(ISYI,JSYM)
  NP = NSH(ISYP,ITP)
  if (NP == 0) cycle
  NPI = NP*NI
  NBRASM = NPI*NV

  LKETSM = 1
  do ISYK=1,NSYM
    NK = NSH(ISYK,ITK)
    if (NK == 0) cycle
    ISYQ = Mul(ISYK,JSYM)
    NQ = NSH(ISYQ,ITQ)
    if (NQ == 0) cycle
    NQK = NQ*NK
    NKETSM = NQK*NV

    ! SVC: we need an NPI*NQK to store the 2-electron integrals, and 2
    ! buffers (values+indices) for sorting them.  Later, we can try to get
    ! rid of the buffer that stores the values and only use an index buffer
    ! and the two-electron integrals for the scatter operation.  For the
    ! buffer, any size can be taken, but assuming there is enough memory
    ! available, it's set to the size of the two-electron integrals unless
    ! larger than some predefined maximum buffer size.
    ! The striped algorithm builds only the local columns and chunks internally, so
    ! it takes the whole buffer and never needs the full (pi,qk) block.
    if (iParRHS == 4) then
      NPIQK = MXPIQK
    else
      NPIQK = NPI*NQK
      if (NPIQK > MXPIQK) then
        if (nCase == 'H') then
          KPI = MXPIQK/NQK
          NPIQK = KPI*NQK
        else if (nCase == 'G') then
          KQK = MXPIQK/NPI
          NPIQK = NPI*KQK
        else
          write(u6,*) ' NPIQK > MXPIQK and case != G or H'
          write(u6,'(A,A2)') ' CASE =   ',nCase
          write(u6,'(A,I12)') ' NPIQK =  ',NPIQK
          write(u6,'(A,I12)') ' MXPIQK = ',MXPIQK
          write(u6,*) ' This should not happen, please report.'
          call AbEnd()
        end if
      end if
    end if

    !-SVC: sanity check
    if (NPIQK <= 0) then
      write(u6,'(1X,A)') ' ADDRHS: zero-sized NPIQK'
      call AbEnd()
    end if

#   ifdef _MOLCAS_MPP_
    if (iParRHS == 4) then
      ! striped RHS construction
      select case (nCase)
        case ('A ')
          call ADDRHSA_STRIPED(JSYM,ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,Cho_Bra(LBRASM),Cho_Ket(LKETSM),NV)
        case ('B ')
          call ADDRHSB_STRIPED(JSYM,ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,Cho_Bra(LBRASM),Cho_Ket(LKETSM),NV)
        case ('D1')
          call ADDRHSD1_STRIPED(JSYM,ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,Cho_Bra(LBRASM),Cho_Ket(LKETSM),NV)
        case ('H ')
          call ADDRHSH_STRIPED(JSYM,ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,Cho_Bra(LBRASM),Cho_Ket(LKETSM),NV)
        case ('C ')
          call ADDRHSC_STRIPED(JSYM,ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,Cho_Bra(LBRASM),Cho_Ket(LKETSM),NV)
        case ('F ')
          call ADDRHSF_STRIPED(JSYM,ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,Cho_Bra(LBRASM),Cho_Ket(LKETSM),NV)
        case ('D2')
          call ADDRHSD2_STRIPED(JSYM,ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,Cho_Bra(LBRASM),Cho_Ket(LKETSM),NV)
        case ('G ')
          call ADDRHSG_STRIPED(JSYM,ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,Cho_Bra(LBRASM),Cho_Ket(LKETSM),NV)
        case ('E ')
          call ADDRHSE_STRIPED(JSYM,ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,Cho_Bra(LBRASM),Cho_Ket(LKETSM),NV)
        case default
          call Abend()
      end select
    else
#   endif
      select case (nCase)
        case ('A ')
          call ADDRHSA(IVEC,JSYM,ISYI,ISYK,NP,NI,NQ,NK,PIQK,nBuff,Buff,idxb,Cho_Bra(LBRASM),Cho_Ket(LKETSM),NV)
        case ('B ')
          call ADDRHSB(IVEC,JSYM,ISYI,ISYK,NP,NI,NQ,NK,PIQK,nBuff,Buff,idxb,Cho_Bra(LBRASM),Cho_Ket(LKETSM),NV)
        case ('D1')
          call ADDRHSD1(IVEC,JSYM,ISYI,ISYK,NP,NI,NQ,NK,PIQK,nBuff,Buff,idxb,Cho_Bra(LBRASM),Cho_Ket(LKETSM),NV)
        case ('H ')
          call ADDRHSH(IVEC,JSYM,ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,nBuff,Buff,idxb,Cho_Bra(LBRASM),Cho_Ket(LKETSM),NV)
        case ('C ')
          call ADDRHSC(IVEC,JSYM,ISYI,ISYK,NP,NI,NQ,NK,PIQK,nBuff,Buff,idxb,Cho_Bra(LBRASM),Cho_Ket(LKETSM),NV)
        case ('F ')
          call ADDRHSF(IVEC,JSYM,ISYI,ISYK,NP,NI,NQ,NK,PIQK,nBuff,Buff,idxb,Cho_Bra(LBRASM),Cho_Ket(LKETSM),NV)
        case ('D2')
          call ADDRHSD2(IVEC,JSYM,ISYI,ISYK,NP,NI,NQ,NK,PIQK,nBuff,Buff,idxb,Cho_Bra(LBRASM),Cho_Ket(LKETSM),NV)
        case ('G ')
          call ADDRHSG(IVEC,JSYM,ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,nBuff,Buff,idxb,Cho_Bra(LBRASM),Cho_Ket(LKETSM),NV)
        case ('E ')
          call ADDRHSE(IVEC,JSYM,ISYI,ISYK,NP,NI,NQ,NK,PIQK,nBuff,Buff,idxb,Cho_Bra(LBRASM),Cho_Ket(LKETSM),NV)
        case default
          call Abend()
      end select
#   ifdef _MOLCAS_MPP_
    end if
#   endif

    LKETSM = LKETSM+NKETSM
  end do
  LBRASM = LBRASM+NBRASM
end do

end subroutine Process_RHS_Block
