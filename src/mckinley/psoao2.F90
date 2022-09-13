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
! Copyright (C) 1990,1992, Roland Lindh                                *
!               1990, IBM                                              *
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine PSOAO2(nSO,MemPrm,MemM,iAnga,iCmpa,iAO,iFnc,iBas,iBsInc,jBas,jBsInc,kBas,kBsInc,lBas,lBsInc,iPrim,iPrInc,jPrim,jPrInc, &
                  kPrim,kPrInc,lPrim,lPrInc,nAco,Mem1,Mem2,Mem3,Mem4,MemX,MemPSO,MemFck,nFT,nCMO,MemFin,MemBuffer,iMemB)
!***********************************************************************
!                                                                      *
!  Object: to partion the SO and AO block. It will go to some length   *
!          before it will start and break up the SO block. This will   *
!          reduce the total flop count. However, as we are breaking up *
!          the AO block this will affect the vectorization. Hence, at  *
!          some point it will actually be better to recompute the      *
!          primitives.                                                 *
!          Current stratergy:                                          *
!          1. Reduce the size of the density matrix and buffer so that *
!             it fits into memory.                                     *
!                                                                      *
!          2. Start reducing the length of the primitives in the order *
!             lPrim,jPrim.                                             *
!                                                                      *
!          3. Reduce the size of the SO block by reducing the number of*
!             basis functions in the order lBas, jBas.                 *
!                                                                      *
!          4. Reduce the size of the Buffer.                           *
!                                                                      *
!          5. Reduce kBas,iBas                                         *
!                                                                      *
!          6. Terminate run telling job max and min of additional      *
!             memory needed to perform the calculation.                *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!             Modified to first order derivatives. January '92         *
!             Anders Bernhardsson Theoretical chemistry, Lund 1995     *
!***********************************************************************

! Memory map for mckinley
!
!---------------------------------------------------------------------------
!|      |               |       |               |               |          |
!|REAL  |  P TRANSF     | RYSG2 |   Transf      | FCK GENERAT   |MO Transf |
!|      |               |       |               |               |          |
!---------------------------------------------------------------------------
!|      |               |       |               |               |Scratch   |
!|  MX  |               |       | 9*abcd*ijkl   |    SS         |space     |
!|      |               |       |               |               |          |
!---------------------------------------------------------------------------
!|      |               |       |               |               |          |
!|  M3  |Scratch space  |Memrys |Scratch space  |    SS         | Scratch  |
!|      |               |       |               |               | space    |
!---------------------------------------------------------------------------
!|      |MEM4 (half tr) |*******|***************|               |          |
!|  M2  |               |       |               |               |  SS      |
!|      |PSO transf     |       |               |Scratch space  |          |
!---------------------------------------------------------------------------
!|      |               |       |               |               |          |
!|  M1  |      P        |   *   |      *        |     *         |    *     |
!|      |               |       |               |               |          |
!---------------------------------------------------------------------------
!|      |      ?        |       |               |               |          |
!|Buffer|***************|*******|Transformed    |***************|**********|
!|      |               |       |integrals      |               |          |
!---------------------------------------------------------------------------

use McKinley_global, only: nMethod, RASSCF
use Index_Functions, only: nTri_Elem1
use Gateway_global, only: force_part_p !, force_part_c
use SOAO_Info, only: iAOtSO
use pso_stuff, only: lPSO
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSO, MemPrm, MemM, iAnga(4), iCmpa(4), iAO(4), iBas, jBas, kBas, lBas, iPrim, jPrim, kPrim, &
                                 lPrim, nAco, iMemB
integer(kind=iwp), intent(out) :: iFnc(4), iBsInc, jBsInc, kBsInc, lBsInc, iPrInc, jPrInc, kPrInc, lPrInc, Mem1, Mem2, Mem3, Mem4, &
                                  MemX, MemPSO, MemFck, nFT, nCMO, MemFin, MemBuffer
#include "pstat.fh"
integer(kind=iwp) :: i1, iiBas(4), iCmp, iFac, iTmp1, j, jCmp, jPam, kCmp, kSOInt, la, lb, lc, lCmp, ld, mabcd, Mem0, MemAux, &
                     MemCntrct, MemDep, MemF, MemMax, MemMO, MemRys, MemScr, MemSph, MemTrn, nabcd, nFac, nijkl, nMax, nMaxC, &
                     nPam(4,0:7), nTmp1, nTmp2
logical(kind=iwp) :: Fail, QiBas, QjBas, QjPrim, QkBas, QlBas, QlPrim
integer(kind=iwp), external :: MemTra

!iRout = 10
!iPrint = nPrint(iRout)
la = iAnga(1)
lb = iAnga(2)
lc = iAnga(3)
ld = iAnga(4)
iCmp = iCmpa(1)
jCmp = iCmpa(2)
kCmp = iCmpa(3)
lCmp = iCmpa(4)
iTotal = iTotal+1
mabcd = nTri_Elem1(la)*nTri_Elem1(lb)*nTri_Elem1(lc)*nTri_Elem1(ld)
nabcd = iCmp*jCmp*kCmp*lCmp

!if (force_part_c) then
!  iBsInc = (iBas+1)/2
!  jBsInc = (jBas+1)/2
!  kBsInc = (kBas+1)/2
!  lBsInc = (lBas+1)/2
!else
iBsInc = iBas
jBsInc = jBas
kBsInc = kBas
lBsInc = lBas
!end if
if (force_part_p) then
  jPrInc = (jPrim+1)/2
  !lPrInc = (lPrim+1)/2
else
  jPrInc = jPrim
  !lPrInc = lPrim
end if
iPrInc = iPrim
kPrInc = kPrim
lPrInc = lPrim
MemBuffer = iMemB
MemMax = MemM-MemBuffer

do
  nijkl = iBsInc*jBsInc*kBsInc*lBsInc
  QjPrim = .false.
  QlPrim = .true.
  QiBas = .false.
  QjBas = .false.
  QkBas = .false.
  QlBas = .false.
  Mem0 = MemMax

  ! Picked MO coeff

  if (nMethod == RASSCF) then
    nCMO = nACO*kCmp*kBas+nACO*lCmp*lBas
  else
    nCMO = 0
  end if

  ! Area for integral storage before transforming them to FM/MO
  ! and place for the CMOs

  MemFin = 9*nijkl*nabcd
  if (MemFin+ncmo+1 > Mem0) then
    MaxReq = max(MaxReq,nCMO+MemFin+1-Mem0)
    QlPrim = .false.
    call Change(iBas,iBsInc,QiBas,kBas,kBsInc,QkBas,jBas,jBsInc,QjBas,lBas,lBsInc,QlBas,jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim, &
                MaxReq,Fail)
    if (Fail) then
      write(u6,*) 'PSOAO2: memory partitioning failed!'
      write(u6,*) '        Restart with more memory!'
      call Abend()
    end if
    cycle
  end if
  ! Subtract one additional word for getmem's internal error check
  Mem0 = Mem0-MemFin-nCMO-1

  !---------------------------------------------------------------------

  ! *** Work1 ***

  ! Memory for 2nd order density matrix in SO basis.

  kSOInt = nSO*iBsInc*jBsInc*kBsInc*lBsInc
  Mem1 = kSOInt

  ! Allocate memory for MO to SO/AO transformation
  ! of the 2nd order density matrix for this shell quadruplet.
  ! and area for AO/SO transformation of Fock matrix.

  if (lPSO) then
    iiBas(1) = iBsInc
    iiBas(2) = jBsInc
    iiBas(3) = kBsInc
    iiBas(4) = lBsInc
    nPam(:,:) = 0
    MemPSO = 1
    nTmp2 = 0
    !call IecPrt('iiBas',iiBas,1,4)

    do jPam=1,4
      iTmp1 = 0
      nTmp1 = 0
      do j=0,nIrrep-1
        do i1=1,iCmpa(jPam)
          if (iAOtSO(iAO(jPam)+i1,j) > 0) then
            nPam(jPam,j) = nPam(jPam,j)+iiBas(jPam)
            nTmp1 = nTmp1+iiBas(jPam)
            iTmp1 = iTmp1+1
          end if
        end do
      end do
      MemPSO = MemPSO*nTmp1
      nTmp2 = nTmp2+nTmp1
      iFnc(jPam) = iTmp1
    end do
    MemScr = MemTra(nPam)
    nFac = 4
    nTmp2 = nTmp2+4
  else
    MemScr = 0
    MemPSO = 0
    nFac = 0
    nTmp2 = 0
  end if
  MemAux = MemPSO+MemScr+nFac*S%nDim+nTmp2+4
  if (Mem1+1+MemAux > Mem0) then
    MaxReq = max(MaxReq,Mem1+1+MemAux-Mem0)
    QjPrim = .false.
    QlPrim = .false.
    QiBas = .false.
    QjBas = .false.
    QkBas = .false.
    QlBas = .true.
    call Change(iBas,iBsInc,QiBas,kBas,kBsInc,QkBas,jBas,jBsInc,QjBas,lBas,lBsInc,QlBas,jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim, &
                MaxReq,Fail)
    if (Fail) then
      write(u6,*) 'PSOAO2: memory partitioning failed!'
      write(u6,*) '        Restart with more memory!'
      call Abend()
    end if
    cycle
  end if
  Mem0 = Mem0-Mem1-1
  !---------------------------------------------------------------------
  !
  ! MemFck: Target for generating the symmetrized Fock Matrix.
  !          Distributed localy.
  ! MemMo : Target for generating the  MO integrals
  !          Distributed localy.
  ! Whole work area is used, if work area is big enough work3 is increased
  !
  !----------------------------------------------------------------------

  MemDep = nijkl*nabcd
  ! Temp+S1+S2
  MemFck = 2*MemDep+max(MemDep,nijkl+max(iBsInc*lBsInc,jBsInc*lBsInc,iBsInc*kBsInc,jBsInc*kBsInc))
  nFT = iBsInc*jBsInc*iCmp*jCmp+kBsInc*lBsInc*kCmp*lCmp+iBsInc*kBsInc*iCmp*kCmp+jBsInc*lBsInc*jCmp*lCmp+iBsInc*lBsInc*iCmp*lCmp+ &
        jBsInc*kBsInc*jCmp*kCmp
  MemFck = MemFck+nFT
  if (nmethod == RASSCF) then

    !3 scratch spaces, sorted integrals and translation

    nMaxC = nACO
    MemFck = MemFck+2*nMaxC
    nMax = max(iCmp*iBsInc,jCmp*jBsInc,kCmp*kBsInc,lcmp*lBsInc)
    nMax = max(nMax,nMaxC)
    MemMO = 3*nMax**4+10*nabcd*nijkl
  else
    MemMo = 0
  end if

  !---------------------------------------------------------------------

  ! *** Work2 and Work4 ***

  ! Memory for 2nd order density matrix in contracted basis
  !   (both cartesian and spherical harmonic) and in primitive basis.
  ! MemDeP: Target for desymmetrization
  ! MemTrn: Scratch and target for decontraction
  ! MemAux: Contracted 2nd order density matrix (if partial decon.)
  ! MemSph: transformation spherical harmonics to cartesian, source and target.

  MemDeP = nabcd*nijkl
  MemTrn = mabcd*max(iBsInc*jBsInc*kBsInc*lBsInc,iPrInc*jPrInc*kBsInc*lBsInc,iPrInc*jPrInc*kPrInc*lPrInc)
  MemTrn = MemTrn+1

  ! If partial decontraction we need to keep the contracted 2nd
  ! order density matrix. (Work4)
  if ((jPrInc /= jPrim) .or. (lPrInc /= lPrim)) then
    MemAux = mabcd*iBsInc*jBsInc*kBsInc*lBsInc
  else
    MemAux = 0
  end if
  MemSph = mabcd*iBsInc*jBsInc*kBsInc*lBsInc
  Mem2 = max(MemTrn+MemAux,MemDeP,MemSph)
  MemFck = MemFck-Mem2
  MemMO = MemMo-Mem2
  if (Mem2+1 > Mem0) then
    MaxReq = max(MaxReq,Mem2+1-Mem0)
    call Change(iBas,iBsInc,QiBas,kBas,kBsInc,QkBas,jBas,jBsInc,QjBas,lBas,lBsInc,QlBas,jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim, &
                MaxReq,Fail)
    if (Fail) then
      write(u6,*) 'PSOAO2: memory partitioning failed!'
      write(u6,*) '        Restart with more memory!'
      call Abend()
    end if
    cycle
  end if
  ! Subtract one additional word for getmem's internal error check
  Mem0 = Mem0-Mem2-1
  MemX = 9*mabcd*iBsInc*jBsInc*kBsInc*lBsInc
  MemFck = MemFck-MemX
  MemMO = MemMo-MemX
  if (MemX+1 > Mem0) then
    MaxReq = max(MaxReq,MemX+1-Mem0)
    call Change(iBas,iBsInc,QiBas,kBas,kBsInc,QkBas,jBas,jBsInc,QjBas,lBas,lBsInc,QlBas,jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim, &
                MaxReq,Fail)
    if (Fail) then
      write(u6,*) 'PSOAO2: memory partitioning failed!'
      write(u6,*) '        Restart with more memory!'
      call Abend()
    end if
    cycle
  end if
  Mem0 = Mem0-MemX

  ! *** Work3 and Work5 ***

  ! Scratch for decontraction and transformation to spherical gaussian.
  ! Working array for Rysg2.
  ! Scratch area for resolving degeneracies due to the double coset
  ! treatement of the symmetry.
  ! MemTrn: Scratch for decontraction
  ! MemRys: Scratch for calualation of primitive integral gradients.

  iFac = 1
  if (mabcd /= 1) iFac = 2
  MemF = 9*nabcd*nijkl
  MemTrn = mabcd*max(iPrInc*jBsInc*kBsInc*lBsInc,iPrInc*jPrInc*kPrInc*lBsInc,iPrInc*jPrInc*kPrInc*lPrInc*iFac)
  MemRys = MemPrm*iPrInc*jPrInc*kPrInc*lPrInc+80

  ! Scratch space for contraction of the integrals

  MemCntrct = 9*mabcd*(max(iBsInc*jBsInc*kBsInc*lPrInc,iBsInc*jPrInc*kPrInc*lPrInc)+iBsInc*jBsInc*kPrInc*lPrInc)

  MemFck = max(0,MemFck)
  MemMo = max(0,MemMo)
  Mem3 = max(MemMO,MemFck,MemTrn,MemRys,2*MemF,MemF+MemCntrct)
  if (Mem3+1 <= Mem0) exit
  MaxReq = max(MaxReq,Mem3+1-Mem0)
  call Change(iBas,iBsInc,QiBas,kBas,kBsInc,QkBas,jBas,jBsInc,QjBas,lBas,lBsInc,QlBas,jPrim,jPrInc,QjPrim,lPrim,lPrInc,QlPrim, &
              MaxReq,Fail)
  if (Fail) then
    write(u6,*) 'PSOAO2: memory partitioning failed!'
    write(u6,*) '        Restart with more memory!'
    call Abend()
  end if
end do
! Subtract one additional word for getmem's internal error check
Mem0 = Mem0-Mem3-1

! Work4, if used, is placed at the end of Work2
if ((jPrInc /= jPrim) .or. (lPrInc /= lPrim)) then
  Mem4 = MemAux
else
  Mem4 = Mem2
end if

return

end subroutine PSOAO2
