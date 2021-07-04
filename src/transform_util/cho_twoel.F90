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
! Copyright (C) 2004, Giovanni Ghigo                                   *
!***********************************************************************
!  Cho_TwoEl
!
!> @brief
!>   Driver for the generation of the two-electrons integrals file (``MOLINT``) from the Transformed Cholesky vectors TCVx
!> @author Giovanni Ghigo
!>
!> @details
!> The routine generates the symmetry
!> block \p iSymI, \p iSymJ, \p iSymA, \p iSymB of two-electrons integrals. The
!> number of integrals to generate for each  \p iSymI, \p iSymJ couple is
!> defined by \c LenInt. The exch-1 and exch-2 integrals are then
!> generated for the occupied MO calling ::Cho_GenE.
!>
!> @param[in]     iBatch                  Main batch current number
!> @param[in]     NumV                    Number of Cholesky vectors to transform in the current batch
!> @param[in]     LUINTM                  Unit number of two-electrons integrals file (``MOLINT``)
!> @param[in,out] iAddrIAD2M              Current disk address in ``MOLINT``
!> @param[in]     iSymI,iSymJ,iSymA,iSymB Symmetry block of the two-electrons integrals
!***********************************************************************

subroutine Cho_TwoEl(iBatch,numV,LUINTM,iAddrIAD2M,iSymI,iSymJ,iSymA,iSymB)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden                                    *
!           October-November 2004                                      *
!----------------------------------------------------------------------*
! Routine for the generation of the new Two-electrons integrals files  *
! (MOLINT) from the Transformed Cholesky Full Vectors (TCVX).          *
!                                                                      *
! A,B are MO indices, counting only non-frozen and non-deleted.        *
! I,J are occupied MO indices, only non-frozen and non-deleted.        *
! (AB/IJ) ARE ALWAYS GENERATED                                         *
! (AI/BJ) IF ISP >= ISR                                                *
! (AI/JB) IF ISP > ISS AND ISP /= ISQ                                  *
! (IA/BJ) IF ISQ > ISR AND ISP /= ISQ                                  *
! (IA/JB) IF ISQ >= ISS AND ISP /= ISQ                                 *
!                                                                      *
!   IAD2M CONTAINS START ADRESS FOR EACH TYPE OF INTEGRALS:            *
!    IAD2M(1,ISPQRS)   COULOMB INTEGRALS <AB|IJ>                       *
!    IAD2M(2,ISPQRS)   EXCHANGE INTEGRALS <AB|IJ> FOR iSymI > iSymJ    *
!    IAD2M(3,ISPQRS)   EXCHANGE INTEGRALS <AB|IJ> FOR iSymI < iSymJ    *
!   THE LAST ADRESS IS ZERO IF iSymI = iSymJ                           *
!                                                                      *
!***********************************************************************

use Cho_Tra, only: DoCoul, DoExc2, DoTCVA, IAD2M, IfTest, nOrb, nOsh, nSsh, nSym, SubBlocks, TCVXist
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iBatch, numV, LUINTM, iSymI, iSymJ, iSymA, iSymB
integer(kind=iwp), intent(inout) :: iAddrIAD2M
integer(kind=iwp) :: i, iAddrIAD2Mij, iEndJ, iI, iIJAB, iJ, iType, j, nA, nB, nN_AB, nN_Ex1, nN_Ex2, nN_IJ, nSymP
real(kind=wp), allocatable :: AddCou(:), AddEx1(:), AddEx2(:), AddEx2t(:)

nSymP = (nSym**2+nSym)/2
call LenInt(iSymI,iSymJ,iSymA,iSymB,nN_IJ,nN_AB,nN_Ex1,nN_Ex2)
!-----------------------------------------------------------------------
if (IfTest) then
  write(u6,*)
  write(u6,'(A,4I3,A,I8,A,I9,A,I9)') '    * [CGG:Cho_TwoEl]: SYMMETRY BLOCK < A B | I J >',iSymA,iSymB,iSymI,iSymJ,': nN_AB=', &
                                     nN_AB,', nN_Ex1=',nN_Ex1,', nN_Ex2=',nN_Ex2
  if (nN_IJ*(nN_AB+nN_Ex1+nN_Ex2) == 0) write(u6,*) '                      Nothing to do!'
  call XFlush(u6)
end if
!-----------------------------------------------------------------------
if (nN_IJ*(nN_AB+nN_Ex1+nN_Ex2) > 0) then
  iIJAB = ((iSymI**2-iSymI)/2+iSymJ-1)*nSymP+(iSymA**2-iSymA)/2+iSymB

  !**** START GENERATION of COULOMB INTEGRALS **************************
  if (DoCoul .and. (nN_AB > 0)) then
    !call Def_SubBlockC(iSymA,iSymB)
    call Def_SubBlockE(iSymA,iSymB)
    !-------------------------------------------------------------------
    if (IfTest) then
      write(u6,*)
      write(u6,*) '    Generation of Coulomb Integrals'
      write(u6,*) '       Available TCVx for Cou: '
      do iType=1,6
        if (TCVXist(iType,iSymA,iSymI)) write(u6,*) '       -TCV x=',iType,' Sym=',iSymA,iSymI
        if (TCVXist(iType,iSymB,iSymJ) .and. (iSymA /= iSymB)) write(u6,*) '       -TCV x=',iType,' Sym=',iSymB,iSymJ
      end do
      write(u6,*)
      write(u6,*) '       SubBlocks to create for Cou: '
      do i=1,3
        do j=1,3
          if (SubBlocks(i,j)) write(u6,*) '       -SB(',i,',',j,')'
        end do
      end do
      call XFlush(u6)
    end if
    !-------------------------------------------------------------------
    if (iBatch == 1) then
      IAD2M(1,iIJAB) = iAddrIAD2M
    else
      iAddrIAD2M = IAD2M(1,iIJAB)
    end if
    ! Start Loop on i, j
    iAddrIAD2Mij = iAddrIAD2M
    do iI=1,nOsh(iSymI)
      if (iSymI == iSymJ) then
        iEndJ = iI
      else
        iEndJ = nOsh(iSymJ)
      end if
      do iJ=1,iEndJ
        !---------------------------------------------------------------
        if (IfTest) then
          write(u6,*)
          write(u6,*) '   Coulomb Integrals for |ij> pair',iI,iJ,'  iAddrIAD2Mij=',iAddrIAD2Mij
          call XFlush(u6)
        end if
        !---------------------------------------------------------------
        call mma_allocate(AddCou,nN_AB,Label='AddCou')
        if (iBatch > 1) then
          call dDaFile(LUINTM,2,AddCou,nN_AB,iAddrIAD2Mij)
          ! Reload Int
          iAddrIAD2Mij = iAddrIAD2Mij-nN_AB
        else
          AddCou(:) = Zero
        end if
        call Cho_GenC(iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV,AddCou,nN_AB,nN_Ex1)
        call GAdSum(AddCou,nN_AB)
        call dDaFile(LUINTM,1,AddCou,nN_AB,iAddrIAD2Mij)
        call mma_deallocate(AddCou)
      end do
    end do
    ! End Loop on i, j
    iAddrIAD2M = iAddrIAD2Mij ! Last written+1 in iAddrIAD2M
  end if
  !**** END GENERATION of COULOMB INTEGRALS ****************************

  !**** START GENERATION of EXCHANGE-1 INTEGRALS ***********************
  if (nN_Ex1 > 0) then
    call Def_SubBlockE(iSymA,iSymB)
    !-------------------------------------------------------------------
    !If(IfTest) then
    !  write(u6,*)
    !  write(u6,*) '    Generation of Exchange-1 Integrals'
    !  write(u6,*) '       Available TCVx for Ex-1: '
    !  do iType=1,MxTCVx
    !    if (TCVXist(iType,iSymA,iSymI)) write(u6,*) '       -TCV x=',iType,' Sym=',iSymA,iSymI
    !    if (TCVXist(iType,iSymB,iSymJ) .and. (iSymA /= iSymB)) write(u6,*) '       -TCV x=',iType,' Sym=',iSymB,iSymJ
    !  end do
    !  write(u6,*)
    !  write(u6,*) '       SubBlocks to create for Ex-1: '
    !  do i=1,3
    !    do j=1,3
    !      if (SubBlocks(i,j)) write(u6,*) '       -SB(',i,',',j,')'
    !    end do
    !  end do
    !  call XFlush(u6)
    !end if
    !-------------------------------------------------------------------
    if (iBatch == 1) then
      IAD2M(2,iIJAB) = iAddrIAD2M
    else
      iAddrIAD2M = IAD2M(2,iIJAB)
    end if
    ! Start Loop on i, j
    iAddrIAD2Mij = iAddrIAD2M
    do iI=1,nOsh(iSymI)
      if (iSymI == iSymJ) then
        iEndJ = iI
      else
        iEndJ = nOsh(iSymJ)
      end if
      do iJ=1,iEndJ
        !---------------------------------------------------------------
        !if (IfTest) then
        !  write(u6,*)
        !  write(u6,*) '   Excha-1 Integrals for |ij> pair',iI,iJ,'  iAddrIAD2Mij=',iAddrIAD2Mij,'   #'
        !  call XFlush(u6)
        !end if
        !---------------------------------------------------------------
        call mma_allocate(AddEx1,nN_Ex1,Label='AddEx1')
        if (iBatch > 1) then
          ! Reload Int
          call dDaFile(LUINTM,2,AddEx1,nN_Ex1,iAddrIAD2Mij)
          iAddrIAD2Mij = iAddrIAD2Mij-nN_Ex1
        else
          AddEx1(:) = Zero
        end if
        call Cho_GenE(iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV,AddEx1,nN_Ex1)
        call GAdSum(AddEx1,nN_Ex1)
        call dDaFile(LUINTM,1,AddEx1,nN_Ex1,iAddrIAD2Mij)
        call mma_deallocate(AddEx1)
      end do
    end do
    ! End Loop on i, j
    iAddrIAD2M = iAddrIAD2Mij
  end if
  !**** END GENERATION of EXCHANGE-1 INTEGRALS *************************

  !**** START GENERATION of EXCHANGE-2 INTEGRALS ***********************
  if ((nN_Ex2 > 0) .and. DoExc2) then
    iIJAB = ((iSymI**2-iSymI)/2+iSymJ-1)*nSymP+(iSymB**2-iSymB)/2+iSymA
    call Def_SubBlockE(iSymA,iSymB)
    ! ------------------------------------------------------------------
    !if (IfTest) then
    !  write(u6,*)
    !  write(u6,*) '    Generation of Exchange-2 Integrals'
    !  write(u6,*) '       Available TCVx for Ex-2: '
    !  do iType=1,MxTCVx
    !    if (TCVXist(iType,iSymA,iSymI)) write(u6,*) '       -TCV x=',iType,' Sym=',iSymA,iSymI
    !    if (TCVXist(iType,iSymB,iSymJ) .and. (iSymA /= iSymB)) write(u6,*) '       -TCV x=',iType,' Sym=',iSymB,iSymJ
    !  end do
    !  write(u6,*)
    !  write(u6,*) '       SubBlocks to create for Ex-2: '
    !  do i = 1, 3
    !    do j = 1, 3
    !      if (SubBlocks(i,j)) write(u6,*) '       -SB(',i,',',j,')'
    !    end do
    !  end do
    !  call XFlush(u6)
    !end if
    ! ------------------------------------------------------------------
    if (iBatch == 1) then
      IAD2M(3,iIJAB) = iAddrIAD2M
    else
      iAddrIAD2M = IAD2M(3,iIJAB)
    end if
    ! Start Loop on i, j
    iAddrIAD2Mij = iAddrIAD2M
    do iI=1,nOsh(iSymI)
      if (iSymI == iSymJ) then
        iEndJ = iI
      else
        iEndJ = nOsh(iSymJ)
      end if
      do iJ=1,iEndJ
        !---------------------------------------------------------------
        !if (IfTest) then
        !  write(u6,*)
        !  write(u6,*) '   Excha-2 Integrals for |ij> pair',iI,iJ,'  iAddrIAD2Mij=',iAddrIAD2Mij,'   #'
        !  call XFlush(u6)
        !end if
        !---------------------------------------------------------------
        if (DoTCVA) then
          nA = nOrb(iSymA)
          nB = nOrb(iSymB)
        else
          nA = nSsh(iSymA)
          nB = nSsh(iSymB)
        end if
        call mma_allocate(AddEx2,nN_Ex2,Label='AddEx2')
        call mma_allocate(AddEx2t,nN_Ex2,Label='AddEx2t')
        if (iBatch > 1) then
          ! Reload Int
          call dDaFile(LUINTM,2,AddEx2,nN_Ex2,iAddrIAD2Mij)
          iAddrIAD2Mij = iAddrIAD2Mij-nN_Ex2
          call Trnsps(nA,nB,AddEx2,AddEx2t)
        else
          AddEx2t(:) = Zero
        end if
        call Cho_GenE(iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV,AddEx2t,nN_Ex2)
        call Trnsps(nB,nA,AddEx2t,AddEx2)
        !---------------------------------------------------------------
        !if (IfTest) then
        !  write(u6,*) '    Integrals from Production code:'
        !  write(u6,'(8F10.6)') (AddEx2(i),i=1,nN_Ex2)
        !  call XFlush(u6)
        !end if
        !---------------------------------------------------------------
        call GAdSum(AddEx2,nN_Ex2)
        call dDaFile(LUINTM,1,AddEx2,nN_Ex2,iAddrIAD2Mij)
        call mma_deallocate(AddEx2t)
        call mma_deallocate(AddEx2)
      end do
    end do
    ! End Loop on i, j
    iAddrIAD2M = iAddrIAD2Mij
  end if
  !**** END GENERATION of EXCHANGE-2 INTEGRALS *************************
end if

return

end subroutine Cho_TwoEl
