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
! Copyright (C) 2005, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 2005  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine EQCTL1()
! On return, the following data sets will be defined and stored
! on LUSOLV.
! At position IVEC=IRHS, the RHS array, in SR representation.
! At position IVEC=IVECX, the solution array, in SR representation.
! At position IVEC=IVECR, the residual array, in SR representation.
! At position IVEC=IVECC, the solution array, in contravariant rep.
! At position IVEC=IVECC2, the solution array, in covariant repr.
! At position IVEC=IVECW, the RHS array, in contravariant repr.

use EQSOLV, only: iRHS, iVecc, iVecc2, iVecR, iVecW, iVecX, MxSct, ModVec, IDSMat, IDBMat, IDTMat, IDSTMat
use caspt2_global, only: do_grad, IDSCT, LUSBT, LUSOLV
use caspt2_module, only: HZERO, MxCase, nASup, nCases, nG2, nG3Tot, nInDep, nISup, nSym
use SC_NEVPT2, only: IDBMAT_NEVPT2
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp, u6, ItoB, RtoI

implicit none
integer(kind=iwp) :: ICASE, IDS, IDS1, IDS2, IDUM(1), IDV, iPad, iPARDIV, ISCT, ISYM, IVEC, LADDR, LENGTH, LSTA, MXWRT, NAS, NB, &
                     NBD, NCOEF, NG3MAX, NID, NIDSCT, NIN, NIS, NISCT, NS, NT
real(kind=wp) :: DUMMY(1)
integer(kind=iwp), parameter :: MXBLK = 40*256*256, MXVEC = 6

IRHS = 1
IVECX = 2
IVECR = 3
IVECC = 4
IVECC2 = 5
IVECW = 6

!SVC: MODVEC isn't used in the Cholesky version, and neither in the
! sigma routines any more. Probably only in MKRHS...
MXSCT = 1
do ICASE=1,NCASES
  do ISYM=1,NSYM
    MODVEC(ISYM,ICASE) = 0
    NAS = NASUP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    NCOEF = NAS*NIS
    if (NCOEF > 0) then
      ! Module lengths for reading/writing:
      MODVEC(ISYM,ICASE) = max(1,min(MXBLK/NAS,NIS))
      MXSCT = max(MXSCT,1+(NIS-1)/MODVEC(ISYM,ICASE))
    end if
  end do
end do
NIDSCT = MXSCT*8*MXCASE*MXVEC
call mma_allocate(IDSCT,NIDSCT,Label='IDSCT')

#ifdef _DEBUG
!SVC: when using Cholesky decomposition, the actual use of the RHS
! vector sizes is automatically controlled in RHSALL2. Furthermore, the
! sigma routines now use the full RHS size.
if (.not. IFCHOL) then
  write(u6,*)
  write(u6,*) ' Size of vector buffers for coefficient arrays.'
  write(u6,*) ' ICASE ISYM    NROW     NCOL     NBLK       Size'
  NVCMX = 0
  do ICASE=1,NCASES
    do ISYM=1,NSYM
      NROW = NASUP(ISYM,ICASE)
      NCOL = MODVEC(ISYM,ICASE)
      NVC = NROW*NCOL
      NVCMX = max(NVCMX,NVC)
      NBLK = 0
      if (NVC > 0) NBLK = 1+(NISUP(ISYM,ICASE)-1)/NCOL
      write(u6,'(2x,I2,3x,I2,3x,I16,3X,I16,3X,I16,3x,I16)') ICASE,ISYM,NROW,NCOL,NBLK,NVC
    end do
  end do
  write(u6,*)
  write(u6,*) ' Largest vector buffer size:',NVCMX
  write(u6,*)
else
  write(u6,*)
  write(u6,*) ' Sizes of the coefficient arrays.'
  write(u6,'(2X,A4,2X,A4,2X,5X,A4,5X,2X,5X,A4,5X,2X,5X,A4,5X)') 'CASE','SYM','NROW','NCOL','SIZE'
  NVCMX = 0
  do ICASE=1,NCASES
    do ISYM=1,NSYM
      NROW = NASUP(ISYM,ICASE)
      NCOL = NISUP(ISYM,ICASE)
      NVC = NROW*NCOL
      NVCMX = max(NVCMX,NVC)
      write(u6,'(2x,I4,2x,I4,2x,I16,2X,I16,2X,I16)') ICASE,ISYM,NROW,NCOL,NVC
    end do
  end do
  write(u6,*)
  write(u6,'(A,I14)') ' Largest vector size: ',NVCMX
  write(u6,*)
end if
#endif

IDV = 0
if (do_grad) then
  !! idxG3 matrix is needed for computing Lagrangian. Here, the
  !! shift avoids the matrix overwritten in PCOLLVEC -> SOLV2DRA
  NG3MAX = iPARDIV(NG3TOT,NG2)
  iPad = ItoB-mod(6*NG3MAX,ItoB)
  IDV = 6*NG3MAX+iPad
end if
do IVEC=1,MXVEC
  do ICASE=1,NCASES
    do ISYM=1,NSYM
      NAS = NASUP(ISYM,ICASE)
      NIS = NISUP(ISYM,ICASE)
      ! NINDEP() IS HERE SET PROVISIONALLY. IT WILL
      ! BE ADJUSTED FOR LINEAR DEPENDENCE LATER.
      NCOEF = NAS*NIS
      MXWRT = max(1,NAS*MODVEC(ISYM,ICASE))
      NISCT = NCOEF/MXWRT+min(1,mod(NCOEF,MXWRT))
      LADDR = 1+MXSCT*(ISYM-1+8*(ICASE-1+MXCASE*(IVEC-1)))
      if (NISCT == 0) IDSCT(LADDR) = IDV
      if (NISCT > MXSCT) then
        write(u6,*) 'EQCTL1 : NISCT= ',NISCT,' > MXSCT= ',MXSCT
        write(u6,*) 'Please, increase MXSCT in eqsolv'
        write(u6,*) 'Do not forget to recompile Molcas afterwards.'
        call Abend()
      end if
      do ISCT=1,NISCT
        LSTA = MXWRT*(ISCT-1)
        LENGTH = RtoI*min(NCOEF-LSTA,MXWRT)
        LADDR = ISCT+MXSCT*(ISYM-1+8*(ICASE-1+MXCASE*(IVEC-1)))
        IDSCT(LADDR) = IDV
        call IDAFILE(LUSOLV,0,IDUM,LENGTH,IDV)
      end do
    end do
  end do
end do

! IDSCT(ISCT,ISYM,ICASE,IVEC) now gives the disk address, on LUSOLV,
! to the ISCT section of the (ISYM,ICASE) block of a vector
! containing expansion C coefficients of excitation operators.
! IVEC=1..MXVEC enumerates the different vectors needed to solve
! the equations.
IDSMAT(:,:) = -1
IDS = 0
do ICASE=1,NCASES
  do ISYM=1,NSYM
    NIN = NINDEP(ISYM,ICASE)
    if (NIN > 0) then
      IDSMAT(ISYM,ICASE) = IDS
      NAS = NASUP(ISYM,ICASE)
      !NS = nTri_Elem(NAS)
      NS = NAS**2
      if (ICASE == 12) NS = 1
      if (ICASE == 13) NS = 1
      call DDAFILE(LUSBT,0,DUMMY,NS,IDS)
    end if
  end do
end do

! IDSMAT(ISYM,ICASE) gives the disk address on LUSBT to the
! (ISYM,ICASE) section of the overlap matrix for the basis
! of excitation wave function terms: those that result from
! the basic excitation operators acting on Psi0.
do ICASE=1,NCASES
  do ISYM=1,NSYM
    NAS = NASUP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    !NB = nTri_Elem(NAS)
    NB = NAS**2
    if (ICASE == 12) NB = 1
    if (ICASE == 13) NB = 1
    IDBMAT(ISYM,ICASE) = IDS
    NBD = NAS
    NID = NIS
    NT = NAS**2
    if (ICASE == 12) NT = 1
    if (ICASE == 13) NT = 1
    IDS1 = IDS
    if (NB > 0) call DDAFILE(LUSBT,0,DUMMY,NB,IDS1)
    IDS2 = IDS
    if (NBD > 0) call DDAFILE(LUSBT,0,DUMMY,NBD,IDS2)
    if (NID > 0) call DDAFILE(LUSBT,0,DUMMY,NID,IDS2)
    if (NBD > 0) call DDAFILE(LUSBT,0,DUMMY,NBD,IDS2)
    if (NID > 0) call DDAFILE(LUSBT,0,DUMMY,NID,IDS2)
    IDTMAT(ISYM,ICASE) = IDS2
    if (NT > 0) call DDAFILE(LUSBT,0,DUMMY,NT,IDS2)
    IDSTMAT(ISYM,ICASE) = IDS2
    if (NT > 0) call DDAFILE(LUSBT,0,DUMMY,NT,IDS2)
    IDS = max(IDS1,IDS2)
    if (HZERO == 'DYALL') then
      IDBMAT_NEVPT2(ISYM,ICASE,1) = IDS
      if (NT > 0) call DDAFILE(LUSBT,0,DUMMY,NT,IDS)
      IDBMAT_NEVPT2(ISYM,ICASE,2) = IDS
      if (NID > 0) call DDAFILE(LUSBT,0,DUMMY,NID,IDS)
    end if
  end do
end do

! IDBMAT(ISYM,ICASE) is similar to IDSMAT() but gives disk address
! to the (ISYM,ICASE) diagonal block of matrix elements of H0
! over the excitation wave function terms.
! IDTMAT() similarly addresses transformation matrices that
! orthonormalize the S matrix blocks.

end subroutine EQCTL1
