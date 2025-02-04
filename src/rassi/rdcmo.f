************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE RDCMO_RASSI(JOB,CMO)
#ifdef _HDF5_
      USE mh5, ONLY: mh5_is_hdf5, mh5_open_file_r, mh5_fetch_dset,
     &               mh5_close_file
#endif
      use rassi_aux, only: ipglob
      use stdalloc, only: mma_allocate, mma_deallocate
      use Cntrl, only: NJOB, PRORB, JBNAME
      use cntrl, only: IDCMO, iTOC15, LuIph
      use Symmetry_Info, only: nSym=>nIrrep
      use rassi_data, only: NCMO,NBASF,NOSH

      IMPLICIT NONE
#ifdef _HDF5_
      integer :: refwfn_id
#endif

      INTEGER JOB
      REAL*8 CMO(NCMO)

      INTEGER I, IAD, IDISK, ISY, L1, L2, LEN, NB, NBUF
      Real*8, Allocatable:: BUF(:)

      CMO(:)=0.0D0
      IF(JOB.LT.1 .OR. JOB.GT.NJOB) THEN
        WRITE(6,*)' RDCMO_RASSI: Invalid JOB parameter.'
        WRITE(6,*)' JOB, NJOB:',JOB,NJOB
        CALL ABEND()
      END IF
      IF(IPGLOB.GE.4) THEN
        WRITE(6,*)' RDCMO_RASSI called for file '//TRIM(JBNAME(JOB))
      END IF
C READ ORBITAL COEFFICIENTS FROM INTERFACE. ORIGINALLY ALL
C CMO COEFFS, INCLUDING VIRTUALS, WERE WRITTEN CONTIGUOUSLY.
      NBUF=0
      DO I=1,NSYM
        NBUF=NBUF+NBASF(I)**2
      END DO
      CALL mma_allocate(BUF,NBUF,Label='BUF')

#ifdef _HDF5_
************************************************************************
*
* For HDF5 formatted job files
*
************************************************************************
      If (mh5_is_hdf5(jbname(job))) Then
        refwfn_id = mh5_open_file_r(jbname(job))
        call mh5_fetch_dset(refwfn_id,'MO_VECTORS',BUF)
        call mh5_close_file(refwfn_id)
      Else
#endif
************************************************************************
*
* For JOBIPH/JOBMIX formatted job files
*
************************************************************************
        CALL DANAME(LUIPH,JBNAME(JOB))
        IAD=0
        CALL IDAFILE(LUIPH,2,ITOC15,30,IAD)
        IDISK=IDCMO(JOB)
        CALL DDAFILE(LUIPH,2,BUF,NBUF,IDISK)
        CALL DACLOS(LUIPH)
************************************************************************
#ifdef _HDF5_
      End If
#endif

      IF(IPGLOB.GE.5) THEN
        write(6,*)' Reading CMO'
        write(6,*)' NBUF=',NBUF
        write(6,*)' Array read in:'
        write(6,'(1x,5f16.8)')(buf(i),i=1,nbuf)
      END IF
      L1=1
      L2=1
      DO ISY=1,NSYM
        NB=NBASF(ISY)
        LEN=NOSH(ISY)*NB
        IF(LEN.GT.0) CALL DCOPY_(LEN,BUF(L1),1,CMO(L2),1)
        L2=L2+LEN
        L1=L1+NB**2
      END DO
      IF(IPGLOB.GE.5) THEN
        write(6,*)' Gathered CMO from array.'
        write(6,*)' NCMO=',NCMO
        write(6,'(1x,5f16.8)')(CMO(i),i=1,ncmo)
      END IF
      CALL mma_deallocate(BUF)

      IF (IPGLOB.gt.0 .and. PRORB) THEN
        WRITE(6,*)
        CALL WRMAT('MO ORBITAL COEFFICIENTS:',
     *               1,NBASF,NOSH,NCMO,CMO)
      END IF

      END SUBROUTINE RDCMO_RASSI
