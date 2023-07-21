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
! Copyright (C) 2004,2008, Thomas Bondo Pedersen                       *
!***********************************************************************

subroutine CD_Decomposer(CD_Col,CD_Vec,MxNumCho,Thr,Span,MxQual,ThrNeg,ThrFail,Diag,Qual,Buf,iPivot,iQual,nDim,lBuf,NumCho,irc)
!
! Thomas Bondo Pedersen, October 2004.
! Modified to compute at most MxNumCho vectors,
!    Thomas Bondo Pedersen, January 2008.
!
! Purpose: Cholesky decompose a matrix.
!          Stop decomposition when either
!          1) max. diag <= Thr
!          2) NumCho = MxNumCho
!
! To use criterion 1) only (standard procedure),
! simply set MxNumCho = nDim.
! To use criterion 2) only,
! simply set Thr=1.0d-20 (i.e. zero)
!
! Note: do *not* call this routine directly;
!       use ChoDec(...) or ChoDec_MxVec instead
!       (see those routines for documentation).
!       This routine contains implicit assumptions
!       that are checked by ChoDec and ChoDec_MxVec!!!
!
! Error codes, irc:
!    0 : all OK
!  301 : too few qualified (probably a bug)
!  302 : insufficient buffer size, lBuf
!  303 : too negative diagonal encountered
!        (matrix non-positive definite!)

implicit real*8(a-h,o-z)
external CD_Col    ! external routine for matrix columns
external CD_Vec    ! external routine for Cholesky vectors
real*8 Diag(nDim), Qual(nDim,MxQual), Buf(lBuf)
integer iPivot(nDim), iQual(MxQual)
character*13 SecNam
parameter(SecNam='CD_Decomposer')
logical Last

irc = 0

iPass = 0
mPass = MxNumCho
do while (iPass < mPass)

  ! Update counter.
  ! ---------------

  iPass = iPass+1

  ! Find max. diagonal.
  ! -------------------

  Dmax = Diag(1)
  do i=2,nDim
    Dmax = max(Dmax,Diag(i))
  end do

  ! Check for convergence.
  ! ----------------------

  if ((Dmax > Thr) .and. (NumCho < MxNumCho)) then

    ! Find largest diagonal elements > DiaMin.
    ! I.e., qualify columns.
    ! ========================================

    nQual = min(MxQual,MxNumCho-NumCho)
    DiaMin = max(Dmax*Span,Thr)
    call CD_DiaMax(Diag,nDim,iPivot,iQual,nQual,DiaMin)

    if (nQual < 1) then ! this would be a bug...
      irc = 301
      Go To 1  ! exit
    end if

    ! Get qualified columns from external routine.
    ! ============================================

    call CD_Col(Qual,nDim,iQual,nQual,Buf,lBuf)

    ! Subtract previous vectors (if any).
    ! ===================================

    if (NumCho > 0) then

      MinBuf = nDim+nQual
      nVec = min(NumCho,lBuf/MinBuf)
      if (nVec < 1) then  ! insufficient buffer size
        irc = 302
        Go To 1 ! exit
      else
        nBatch = (NumCho-1)/nVec+1
      end if

      do iBatch=1,nBatch

        if (iBatch == nBatch) then
          NumV = NumCho-nVec*(nBatch-1)
        else
          NumV = nVec
        end if

        iVec1 = nVec*(iBatch-1)+1
        lVec = nDim*NumV

        kOffV = 1
        kOffQ = kOffV+lVec

        iOpt = 2
        call CD_Vec(iVec1,NumV,Buf(kOffV),lBuf-kOffV+1,nDim,iOpt)

        do jVec=1,NumV
          do i=1,nQual
            kOff1 = kOffQ+NumV*(i-1)+jVec-1
            kOff2 = kOffV+nDim*(jVec-1)+iQual(i)-1
            Buf(kOff1) = Buf(kOff2)
          end do
        end do

        call DGEMM_('N','N',nDim,nQual,NumV,-1.0d0,Buf(kOffV),nDim,Buf(kOffQ),NumV,1.0d0,Qual(1,1),nDim)

      end do

    end if

    ! Decompose.
    ! ==========

    MxVec = min(nQual,lBuf/nDim)
    iDump = 0
    iChoMx = nQual
    iCho = 0
    do while (iCho < iChoMx)

      ! Find max. among qualified.
      ! --------------------------

      Dx = Diag(iQual(1))
      ix = 1
      do i=2,nQual
        if (Diag(iQual(i)) > Dx) then
          Dx = Diag(iQual(i))
          ix = i
        end if
      end do

      Last = (Dx < DiaMin) .or. (Dx <= Thr)
      if (.not. Last) then

        ! Calculate new vector.
        ! ---------------------

        Factor = 1.0d0/sqrt(Dx)
        do i=1,nDim
          if (Diag(i) == 0.0d0) then
            Qual(i,ix) = 0.0d0
          else
            Qual(i,ix) = Factor*Qual(i,ix)
          end if
        end do

        ! Update diagonal and find new max.
        ! ---------------------------------

        Diag(1) = Diag(1)-Qual(1,ix)*Qual(1,ix)
        xm = Diag(1)
        do i=2,nDim
          Diag(i) = Diag(i)-Qual(i,ix)*Qual(i,ix)
          xm = max(xm,Diag(i))
        end do

        ! Zero treated diagonal and find new DiaMin.
        ! ------------------------------------------

        Diag(iQual(ix)) = 0.0d0
        DiaMin = max(xm*Span,Thr)

        ! Zero negative diagonals (quit if too negative).
        ! -----------------------------------------------

        do i=1,nDim
          if (Diag(i) < ThrNeg) then
            if (Diag(i) < ThrFail) then
              irc = 303
              Go To 1 ! exit (too negative diagonal)
            else
              Diag(i) = 0.0d0
            end if
          end if
        end do

        ! Subtract this vector from qualified columns.
        ! --------------------------------------------

        do i=1,nQual
          if (Diag(iQual(i)) /= 0.0d0) then
            Factor = -Qual(iQual(i),ix)
            call dAXPY_(nDim,Factor,Qual(1,ix),1,Qual(1,i),1)
          end if
        end do

        ! Store vector in buffer.
        ! -----------------------

        kOff = nDim*iDump+1
        call dCopy_(nDim,Qual(1,ix),1,Buf(kOff),1)

        ! Update counter.
        ! ---------------

        iDump = iDump+1

      end if

      ! Dump vectors to external routine CD_Vec.
      ! ----------------------------------------

      if (Last .or. (iDump == MxVec)) then
        if (iDump > 0) then
          iVec1 = NumCho+1
          iOpt = 1
          call CD_Vec(iVec1,iDump,Buf,lBuf,nDim,iOpt)
          NumCho = NumCho+iDump
        end if
        if (Last) then
          iCho = iChoMx+1 ! break iCho loop (next pass)
        else
          iDump = 0
          iCho = iCho+1
        end if
      end if

    end do

  else ! converged; break while loop

    iPass = mPass+1

  end if

end do

1 continue

end subroutine CD_Decomposer
