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
! Copyright (C) 2005, Giovanni Ghigo                                   *
!***********************************************************************
!  Mem_Est
!
!> @brief
!>   The routine makes an estimation of the optimal memory allocation and defines the
!>   maximum number of Cholesky vectors to transform both for the main (\p nVec) and
!>   the inner (\p nFVec) batch procedures
!> @author Giovanni Ghigo
!>
!> @param[in]  iSymL Symmetry of the Cholesky vector
!> @param[out] nVec  Number of Cholesky vectors to transform in the batch procedure
!> @param[out] nFVec Number of Cholesky vectors to transform in the inner batch procedure
!***********************************************************************

subroutine Mem_Est(iSymL,nVec,nFVec)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Torino University, Italy                                   *
!           January-February, July 2005                                *
!----------------------------------------------------------------------*
! Routine for the estimation of the memory usage and the number (nVec) *
! of Cholesky Vectors transformable.                                   *
! The routine also define (calling Def_TCVx) which TCVx create.        *
!***********************************************************************

use Cho_Tra, only: IfTest, nAsh, nBas, nIsh, nOrb, nSsh, nSym, NumCho, TCVXist
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iSymL
integer(kind=iwp), intent(out) :: nVec, nFVec
integer(kind=iwp) :: iSym, iSymA, iSymAI, iSymB, iSymBJ, iSymI, iSymJ, jSym, Len_ABSq, Len_XAb, Len_XAj, Len_XAu, Len_XBi, &
                     Len_XBt, Len_YAb, Len_YAj, Len_YAu, LenCHFV, LenTCVx, LenTmpTra, MaxInt, MaxNum, MaxSize, MaxSlice, MemAlloc, &
                     MemFree, MemFree0, MemMin, MemPerVec2, MEMX, Nab, Naj, Nau, Nbi, Nbt, Nij, Niu, Nji, Njt, nN_AB, nN_Ex1, &
                     nN_Ex2, nN_IJ, Ntj, Ntu, Nui, Nut
character(len=16) :: Frmt1, Frmt3

nVec = 0  ! Nr. of transformable vectors for batch.
nFVec = 0 ! Nr. of full vectors for inner transformation batch

! Memory for Reading (will change when reading reduced set)
LenCHFV = 0   ! Mem for Full Vectors (p.v. = per vector)
! Memory for Transformation
LenTCVx = 0   ! Mem for TCVx (p.v.)
LenTmpTra = 0 ! Mem for First-half Transformation
! Memory for Generation
MaxInt = 0    ! Mem for integrals
MaxSlice = 0  ! Mem for Lx (p.v.)

! Memory for Transformation
Nij = 0  ! A
Nji = 0  ! A"
Ntj = 0  ! B
Nui = 0  ! B"
Naj = 0  ! C
Nbi = 0  ! C"
Ntu = 0  ! D
Nut = 0  ! D"
Nau = 0  ! E
Nbt = 0  ! E"
Nab = 0  ! F
Njt = 0  ! G
Niu = 0  ! G"

do iSym=1,nSym
  if (nBas(iSym) > 0) then
    do jSym=1,iSym
      if ((nBas(jSym) > 0) .and. (Mul(iSym,jSym) == iSymL)) then

        Len_YAj = 0  ! iSym=jSym: A, B, C
        Len_YAu = 0  ! iSym=jSym: D, E
        Len_YAb = 0  ! iSym=jSym: F
        Len_XAj = 0  ! iSym=/=jSym: A/A", B, C
        Len_XAu = 0  ! iSym=/=jSym: D/D", E
        Len_XAb = 0  ! iSym=/=jSym: F
        Len_XBi = 0  ! iSym=/=jSym: B",C"
        Len_XBt = 0  ! iSym=/=jSym: E"
        Len_ABSq = 0 ! Squared CHFV

        call Def_TCVx(iSym,jSym)

        if (iSym == jSym) then

          ! TCV-A :
          if (TCVXist(1,iSym,jSym)) then
            Len_YAj = nBas(iSym)*nIsh(jSym)
            Nij = Nij+nIsh(iSym)*nIsh(jSym)
          end if
          ! TCV-B :
          if (TCVXist(2,iSym,jSym)) then
            Len_YAj = nBas(iSym)*nIsh(jSym)
            Ntj = Ntj+(nAsh(iSym)*nIsh(jSym))
            ! TCV-G :
            Njt = Njt+(nIsh(jSym)*nAsh(iSym))
          end if
          ! TCV-C :
          if (TCVXist(3,iSym,jSym)) then
            Len_YAj = nBas(iSym)*nIsh(jSym)
            Naj = Naj+nSsh(iSym)*nIsh(jSym)
          end if
          ! TCV-D :
          if (TCVXist(4,iSym,jSym)) then
            Len_YAu = nBas(iSym)*nAsh(jSym)
            Ntu = Ntu+nAsh(iSym)*nAsh(jSym)
          end if
          ! TCV-E :
          if (TCVXist(5,iSym,jSym)) then
            Len_YAu = nBas(iSym)*nAsh(jSym)
            Nau = Nau+nSsh(iSym)*nAsh(jSym)
          end if
          ! TCV-F :
          if (TCVXist(6,iSym,jSym)) then
            Len_YAb = nBas(iSym)*nSsh(jSym)
            Nab = Nab+nSsh(iSym)*nssh(jSym)
          end if

          LenCHFV = LenCHFV+(nBas(iSym)*(nBas(jSym)+1)/2)
          Len_ABSq = nBas(iSym)*nBas(jSym)
          LenTmpTra = max(LenTmpTra,(Len_ABSq+Len_YAj+Len_YAu+Len_YAb))

        else

          ! TCV-A :
          if (TCVXist(1,iSym,jSym)) then
            Len_XAj = nBas(iSym)*nIsh(jSym)
            Nij = Nij+(nIsh(iSym)*nIsh(jSym))
            Nji = Nji+(nIsh(jSym)*nIsh(iSym))
          end if
          ! TCV-B :
          if (TCVXist(2,iSym,jSym)) then
            Len_XAj = nBas(iSym)*nIsh(jSym)
            Ntj = Ntj+(nAsh(iSym)*nIsh(jSym))
            ! TCV-G :
            Njt = Njt+(nIsh(jSym)*nAsh(iSym))
          end if
          if (TCVXist(2,jSym,iSym)) then
            Len_XBi = nBas(jSym)*nIsh(iSym)
            Nui = Nui+(nAsh(jSym)*nIsh(iSym))
            ! TCV-G :
            Niu = Niu+(nIsh(iSym)*nAsh(jSym))
          end if
          ! TCV-C :
          if (TCVXist(3,iSym,jSym)) then
            Len_XAj = nBas(iSym)*nIsh(jSym)
            Naj = Naj+nSsh(iSym)*nIsh(jSym)
          end if
          if (TCVXist(3,jSym,iSym)) then
            Len_XBi = nBas(jSym)*nIsh(iSym)
            Nbi = Nbi+nSsh(jSym)*nIsh(iSym)
          end if
          ! TCV-D :
          if (TCVXist(4,iSym,jSym)) then
            Len_XAu = nBas(iSym)*nAsh(jSym)
            Ntu = Ntu+(nAsh(iSym)*nAsh(jSym))
            Nut = Nut+(nAsh(jSym)*nAsh(iSym))
          end if
          ! TCV-E :
          if (TCVXist(5,iSym,jSym)) then
            Len_XAu = nBas(iSym)*nAsh(jSym)
            Nau = Nau+nSsh(iSym)*nAsh(jSym)
          end if
          if (TCVXist(5,jSym,iSym)) then
            Len_XBt = nBas(jSym)*nAsh(iSym)
            Nbt = Nbt+nSsh(jSym)*nAsh(iSym)
          end if
          ! TCV-F :
          if (TCVXist(6,iSym,jSym)) then
            Len_XAb = nBas(iSym)*nSsh(jSym)
            Nab = Nab+nSsh(iSym)*nSsh(jSym)
          end if

          LenCHFV = LenCHFV+(nBas(iSym)*nBas(jSym))
          LenTmpTra = max(LenTmpTra,(Len_XAj+Len_XAu+Len_XAb+Len_XBi+Len_XBt))

        end if

      end if
    end do ! jSym
  end if
end do ! iSym

LenTCVx = Nij+Nji+Ntj+Nui+Naj+Nbi+Ntu+Nut+Nau+Nbt+Nab+Niu+Njt

! Memory for Generation
do iSymI=1,nSym
  do iSymJ=1,iSymI
    do iSymA=1,nSym
      do iSymB=1,iSymA
        iSymAI = Mul(iSymA,iSymI)
        iSymBJ = Mul(iSymB,iSymJ)
        call LenInt(iSymI,iSymJ,iSymA,iSymB,nN_IJ,nN_AB,nN_Ex1,nN_Ex2)
        if ((iSymAI == iSymL) .and. (iSymBJ == iSymL) .and. (nN_IJ*nN_AB > 0)) then
          MaxInt = max(MaxInt,max(nN_AB,max(nN_Ex1,2*nN_Ex2)))
          MaxSlice = max(MaxSlice,nOrb(iSymA)+nOrb(iSymB))
        end if
      end do
    end do
  end do
end do

call mma_maxDBLE(MEMX)
MemFree0 = max(MEMX-MEMX/10,0)

MemPerVec2 = LenTCVx+MaxSlice

MemMin = LenCHFV+LenTmpTra+LenTCVx+MaxSlice+2*MaxInt
MemFree = MemFree0-MemMin

nVec = min((MemFree/MemPerVec2+1),NumCho(iSymL))
MemFree = MemFree-(nVec-1)*MemPerVec2
nFVec = min((MemFree/LenCHFV+1),NumCho(iSymL))
MemAlloc = MemMin+(nVec-1)*MemPerVec2+(nFVec-1)*LenCHFV

if ((.not. IfTest) .and. (nVec == NumCho(iSymL)) .and. (nFVec > 0)) return

MaxNum = max(MemFree0,MemAlloc)
MaxSize = max(9,int(log10(real(MaxNum,kind=wp)))+1)
write(Frmt1,'(I16)') MaxSize
Frmt1 = '(A,1X,I'//trim(adjustl(Frmt1))//')'
write(Frmt3,'(I16)') MaxSize-5
Frmt3 = '(A,1X,I'//trim(adjustl(Frmt3))//',A)'

write(u6,*)
write(u6,'(A)') repeat('-',60)
write(u6,105) ' MEM for TRANSF/GENER of VECTORS in Sym:',iSymL
write(u6,*)
write(u6,Frmt1) ' Mem (p.v.) for CHFV            :',LenCHFV
write(u6,Frmt1) ' Tmp Mem for transformation     :',LenTmpTra
write(u6,Frmt1) ' Mem (p.v.) for TCVx            :',LenTCVx
write(u6,Frmt1) ' Max Tmp Mem (p.v.) for generat.:',MaxSlice
write(u6,Frmt1) ' Max Mem for Integrals (twice)  :',2*MaxInt
write(u6,Frmt1)
write(u6,Frmt1) ' TOTAL AVAILABLE MEMORY         :',MemFree0
write(u6,Frmt1) ' Minimal Memory required        :',MemMin
write(u6,Frmt1) ' Memory Required by generation  :',(nVec-1)*MemPerVec2
write(u6,Frmt1) ' Memory Available for transform.:',MemFree
write(u6,Frmt1) ' Memory Required by transform.  :',(nFVec-1)*LenCHFV
write(u6,Frmt1) ' TOTAL ALLOCATED MEMORY         :',MemAlloc
write(u6,Frmt1) ' Unemployed memory              :',MemFree0-MemAlloc
write(u6,*)
write(u6,120) ' Max nr. of Transformed vectors  :',nVec
write(u6,120) ' Max nr. of vectors in sub-batch :',nFVec
write(u6,120) ' Total Number of Cholesky vectors:',NumCho(iSymL)
write(u6,*)
write(u6,*) 'ESTIMATED MEMORY REQUIREMENTS'
write(u6,Frmt3) '  Minimum:',2+MemMin/119000,' MB'
write(u6,Frmt3) '  Normal :',2+(MemMin+MemPerVec2*(NumCho(iSymL)-1))/119000,' MB'
write(u6,Frmt3) '  Maximum:',2+(MemMin+(MemPerVec2+LenCHFV)*(NumCho(iSymL)-1))/119000,' MB'
write(u6,'(A)') repeat('-',20)
write(u6,*)
call XFlush(u6)
if (nVec < NumCho(iSymL)) then
  write(u6,*) ' Batch procedure used. Increase memory if possible!'
  write(u6,*)
  write(u6,'(A)') repeat('-',20)
  call XFlush(u6)
end if

return

105 format(A,1X,I1)
120 format(A,1X,I6)

end subroutine Mem_Est
