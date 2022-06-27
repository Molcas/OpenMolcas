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
! Copyright (C) 2004,2005, Giovanni Ghigo                              *
!               2021, Roland Lindh                                     *
!***********************************************************************
!  ChoMP2_TraCtl
!
!> @brief
!>   Driver for the generation of the two-electrons integrals file (``MOLINT``)
!>   from the AO-based Cholesky Full Vectors for MBPT2 program
!> @author Giovanni Ghigo
!>
!> @details
!> This routine is similar to ::Cho_TraCtl routine but only the
!> transformed vectors type ``C`` (TCVC) are generated. See ::Cho_TraCtl
!> and related routines for more details.
!>
!> @note
!> The number of frozen and deleted MO used in the post-SCF
!> must be written in the RunFile in \c nFroPT and \c nDelPT arrays.
!>
!> @param[in] LUINTM Unit number of two-electrons integrals file (``MOLINT``)
!> @param[in] CMO    MO coefficients
!> @param[in] NCMO   Total number of MO coefficients
!***********************************************************************

subroutine ChoMP2_TraCtl(LUINTM,CMO,NCMO)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden                                    *
!           October 2004                                               *
! Modified for Cholesky-MP2 May 2005                                   *
!----------------------------------------------------------------------*
! This is the routine for the transformation from AO basis to MO basis *
! of the Cholesky Full Vectors & Two-electrons integral generation.    *
! Reads the Cholesky Full Vectors (CHFV), generate new Transformed     *
! Cholesky Full Vectors (TCVx) then generate the Two-electrons inte-   *
! gral file (MOLINT).                                                  *
!  <p,q|k,l>  where p,q: All MO; k,l: Occupied                         *
! THIS CODE IS ONLY FOR MBPT2 AND IT IS SPLIT FROM THE GENERAL CODE    *
! BUT IT IS STILL INTEGRATED AND DEPENDENT ON THE GENERAL CODE         *
!***********************************************************************

use Cho_Tra, only: DoCoul, DoFull, DoTCVA, IAD2M, IfTest, nAsh, nBas, nDel, nFro, nIsh, nOrb, nOsh, NumCho, nSsh, nSym, TCVX, &
                   TCVXist
use Symmetry_Info, only: Mul
use stdalloc, only: mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LUINTM, NCMO
real(kind=wp), intent(in) :: CMO(NCMO)
integer(kind=iwp) :: i, iAddrIAD2M, iBatch, iStrtVec_AB, iSym, iSymA, iSymAI, iSymB, iSymBJ, iSymI, iSymJ, iSymL, jSym, LenIAD2M, &
                     lUCHFV, nBasT, nBatch, nData, nFVec, NumV, nVec
real(kind=wp) :: CPE, CPU0, CPU1, CPU2, CPU3, CPU4, CPU_Gen, CPU_Tot, CPU_Tra, TIO0, TIO1, TIO2, TIO3, TIO4, TIO_Gen, TIO_Tot, &
                 TIO_Tra, TIOE
logical(kind=iwp) :: Found
character(len=6) :: CHName
character(len=4), parameter :: CHNm = 'CHFV'

IfTest = .false.

call Timing(CPU0,CPE,TIO0,TIOE)

!**** INITIALIZATION ***************************************************

! Define what has to be calculated.
!  DoTCVA flag for the generation of TCVA
!  DoFull flag for the generation of TCVF
!  DoCoul flag for the generation of coulomb integrals

! MBPT2
DoTCVA = .false.
DoFull = .false.
DoCoul = .false.

! Get Informations from RunFile.
!  The following informations must be passed to the Cholesky
!  transformation routine through RunFile. COMMON blocks could
!  not be used due to several conflicts.
call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)
call Get_iArray('nFroPT',nFro,nSym)
call Get_iArray('nDelPT',nDel,nSym)
call Get_iArray('nIsh',nIsh,nSym)
!PAM07 call Get_iArray('nAsh',nAsh,nSym)
! Replaced by:
do i=1,nSym
  nAsh(i) = 0
end do
call qpg_iArray('nAsh',Found,nData)
if (Found .and. (nData == nSym)) then
  call Get_iArray('nAsh',nAsh,nSym)
end if
! End of replacement
call Get_iArray('NumCho',NumCho,nSym)
nBasT = 0
do i=1,nSym
  nBasT = nBasT+nBas(i)
  nOrb(i) = nBas(i)-nFro(i)-nDel(i)
  nIsh(i) = nIsh(i)-nFro(i)
  nOsh(i) = nIsh(i)+nAsh(i)
  nSsh(i) = nOrb(i)-nOsh(i)
end do

! Initialize information arrays.

TCVXist(:,:,:) = .false. ! TCVx existing flag.

! The Address Field for MOLINT:
LenIAD2M = 3*36*36
do i=1,36*36
  IAD2M(1,i) = 0
  IAD2M(2,i) = 0
  IAD2M(3,i) = 0
end do
iAddrIAD2M = 0
call iDaFile(LUINTM,1,IAD2M,LenIAD2M,iAddrIAD2M)

! The Timing:
CPU_Tra = Zero
TIO_Tra = Zero
CPU_Gen = Zero
TIO_Gen = Zero

!**** START LOOP iSymL on TOTAL SYMMETRY of L (iSym * jSym) ************
do iSymL=1,nSym

  ! Re-Initialize the TCVx
  TCVXist(:,:,:) = .false. ! TCVx existing flag.

  call Mem_Est(iSymL,nVec,nFVec)
  if ((nVec <= 0) .or. (nFVec <= 0)) then
    write(u6,*)
    write(u6,*) ' ************************************'
    write(u6,*) ' *  Insufficient memory for batch ! *'
    write(u6,*) ' ************************************'
    write(u6,*)
    call Abend()
  end if
  nBatch = (NumCho(iSymL)-1)/nVec+1

  ! START LOOP iBatch
  do iBatch=1,nBatch
    if (iBatch == nBatch) then
      NumV = NumCho(iSymL)-nVec*(nBatch-1)
    else
      NumV = nVec
    end if

    ! Start Transformation of Cholesy Vectors  CHFV -> TCVx
    call Timing(CPU1,CPE,TIO1,TIOE)

    ! Start Loop on CHFV-iSym-jSym
    do iSym=1,nSym
      lUCHFV = -1
      if (nBas(iSym) > 0) then
        do jSym=1,iSym
          lUCHFV = -1
          if ((nBas(jSym) > 0) .and. (Mul(iSym,jSym) == iSymL)) then
            lUCHFV = 7
            iStrtVec_AB = nVec*(iBatch-1)+1
            write(CHName,'(A4,I1,I1)') CHNm,iSym,jSym
            call dAName_MF_WA(lUCHFV,CHName)
            if (iSym == jSym) then
              call ChoMP2_TraS(iSym,jSym,NumV,CMO,NCMO,lUCHFV,iStrtVec_AB,nFVec)
            else
              call ChoMP2_TraA(iSym,jSym,NumV,CMO,NCMO,lUCHFV,iStrtVec_AB,nFVec)
            end if
            call dAClos(lUCHFV)
          end if
        end do
      end if
    end do
    ! End Loop on CHFV-iSym-jSym

    call Timing(CPU2,CPE,TIO2,TIOE)
    CPU_Tra = CPU_Tra+CPU2-CPU1
    TIO_Tra = TIO_Tra+TIO2-TIO1
    ! End Transformation of Cholesky Vectors  CHFV -> TCVx

    ! Start Generation of Integrals files  TCVx -> MOLINT

    ! Start Loop on I, J, A, B Symmetries
    do iSymI=1,nSym
      do iSymJ=1,iSymI
        do iSymA=1,nSym
          do iSymB=1,nSym
            iSymAI = Mul(iSymA,iSymI)
            iSymBJ = Mul(iSymB,iSymJ)

            if ((iSymAI == iSymL) .and. (iSymBJ == iSymL)) call ChoMP2_TwoEl(iBatch,NumV,LUINTM,iAddrIAD2M,iSymI,iSymJ,iSymA,iSymB)

          end do
        end do
      end do
    end do
    ! End Loop on I, J, A, B Symmetries

    do iSym=1,size(TCVX,2)
      do jSym=1,size(TCVX,3)

        if (allocated(TCVX(3,iSym,jSym)%A)) call mma_deallocate(TCVX(3,iSym,jSym)%A)

      end do
    end do
    call Timing(CPU3,CPE,TIO3,TIOE)
    CPU_Gen = CPU_Gen+CPU3-CPU2
    TIO_Gen = TIO_Gen+TIO3-TIO2
    ! End Generation of Two-Electrons Integrals files

  end do
  ! END LOOP iBatch

end do
!**** END MAIN LOOP iSymL on TOTAL SYMMETRY of L (iSym * jSym) *********

iAddrIAD2M = 0
call iDaFile(LUINTM,1,IAD2M,LenIAD2M,iAddrIAD2M)

write(u6,*) 'TIMING INFORMATION:   CPU(s)   %CPU   Elapsed(s)'
write(u6,'(A,F9.2,1X,F6.1,1X,F12.2)') ' Transformation     ',CPU_Tra,1.0e2_wp*(CPU_Tra+5.0e-5_wp)/(TIO_Tra+5.0e-5_wp),TIO_Tra
write(u6,'(A,F9.2,1X,F6.1,1X,F12.2)') ' Generation         ',CPU_Gen,1.0e2_wp*(CPU_Gen+5.0e-5_wp)/(TIO_Gen+5.0e-5_wp),TIO_Gen
call Timing(CPU4,CPE,TIO4,TIOE)
CPU_Tot = CPU4-CPU0
TIO_Tot = TIO4-TIO0
write(u6,'(A,F9.2,1X,F6.1,1X,F12.2)') ' TOTAL              ',CPU_Tot,1.0e2_wp*(CPU_Tot+5.0e-5_wp)/(TIO_Tot+5.0e-5_wp),TIO_Tot
call XFlush(u6)

return

end subroutine ChoMP2_TraCtl
