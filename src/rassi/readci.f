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
      SUBROUTINE READCI(ISTATE,ISGS,ICIS,NCI,CI)
      use rassi_global_arrays, only: JBNUM, LROOT
#ifdef _HDF5_
      USE mh5, ONLY: mh5_is_hdf5, mh5_open_file_r, mh5_exists_attr,
     &               mh5_fetch_attr, mh5_fetch_dset_array_real
#endif
      IMPLICIT NONE
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='READCI')
#include "rasdim.fh"
#include "cntrl.fh"
#include "Files.fh"
#include "WrkSpc.fh"
#include "rassi.fh"
#include "jobin.fh"
#include "Struct.fh"
#include "SysDef.fh"
#ifdef _HDF5_
      integer :: refwfn_id
      integer :: root2state(mxroot), IDXCI
#endif

      INTEGER ISTATE
      INTEGER ISGS(NSGSIZE), ICIS(NCISIZE)
      INTEGER NCI
      REAL*8 CI(NCI)

      INTEGER I, IAD, IDISK, JOB, LROOT1, LSYM


      IF(ISTATE.LT.1 .OR. ISTATE.GT.NSTATE) THEN
        WRITE(6,*)'RASSI/READCI: Invalid ISTATE parameter.'
        WRITE(6,*)' ISTATE, NSTATE:',ISTATE,NSTATE
        CALL ABEND()
      END IF
      JOB=JBNUM(ISTATE)
      LROOT1=LROOT(ISTATE)

#ifdef _HDF5_
************************************************************************
*
* For HDF5 formatted job files
*
************************************************************************
      If (mh5_is_hdf5(jbname(job))) Then
        refwfn_id = mh5_open_file_r(jbname(job))
        if (mh5_exists_attr(refwfn_id, 'ROOT2STATE')) then
          call mh5_fetch_attr(refwfn_id,'ROOT2STATE',root2state)
          IDXCI = root2state(lroot1)
        else
          IDXCI = lroot1
        end if
        IF (IDXCI.LE.0.OR.IDXCI.GT.NROOTS(JOB)) THEN
          call WarningMessage(2,'Invalid CI array index, abort!')
          call AbEnd
        END IF
        call mh5_fetch_dset_array_real(refwfn_id,
     &         'CI_VECTORS',CI,[NCI,1],[0,IDXCI-1])
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
        IDISK=ITOC15(4)
        DO I=1,LROOT1-1
          CALL DDAFILE(LUIPH,0,CI,NCI,IDISK)
        END DO
        CALL DDAFILE(LUIPH,2,CI,NCI,IDISK)
        CALL DACLOS(LUIPH)
************************************************************************
#ifdef _HDF5_
      End If
#endif

      IF(IPGLOB.gt.SILENT .and. PRCI) THEN
        WRITE(6,*)' READCI called for state ',ISTATE
        WRITE(6,*)' This is on JobIph nr.',JOB
        WRITE(6,*)' JobIph file name:',JBNAME(JOB)
        WRITE(6,*)' It is root nr.',LROOT(ISTATE)
        WRITE(6,*)' Its length NCI=',NCI
        WRITE(6,*)' Its symmetry  =',IRREP(JOB)
        WRITE(6,*)' Spin multiplic=',MLTPLT(JOB)
        LSYM=IRREP(JOB)
        CALL PRWF(ISGS,ICIS,LSYM,CI,CITHR)
      END IF

      RETURN
      END
