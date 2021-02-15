************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2005, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 2005  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE EQCTL1
      IMPLICIT REAL*8 (A-H,O-Z)
C On return, the following data sets will be defined and stored
C on LUSOLV.
C At position IVEC=IRHS, the RHS array, in SR representation.
C At position IVEC=IVECX, the solution array, in SR representation.
C At position IVEC=IVECR, the residual array, in SR representation.
C At position IVEC=IVECC, the solution array, in contravariant rep.
C At position IVEC=IVECC2, the solution array, in covariant repr.
C At position IVEC=IVECW, the RHS array, in contravariant repr.
#include "rasdim.fh"
#include "caspt2.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
      DIMENSION DUMMY(1),IDUM(1)

      IRHS  =1
      IVECX =2
      IVECR =3
      IVECC =4
      IVECC2=5
      IVECW =6

CSVC: MODVEC isn't used in the Cholesky version, and neither in the
C sigma routines any more. Probably only in MKRHS...
      MXSCT=1
      DO ICASE=1,NCASES
        DO ISYM=1,NSYM
          MODVEC(ISYM,ICASE)=0
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          NCOEF=NAS*NIS
          IF(NCOEF.GT.0) THEN
C Module lengths for reading/writing:
            MODVEC(ISYM,ICASE)=MAX(1,MIN(MXBLK/NAS,NIS))
            MXSCT=MAX(MXSCT,1+(NIS-1)/MODVEC(ISYM,ICASE))
          ENDIF
        END DO
      END DO
      NIDSCT=MXSCT*8*MXCASE*MXVEC
      CALL GETMEM('IDSCT','ALLO','INTE',LIDSCT,NIDSCT)

#ifdef _DEBUG
CSVC: when using Cholesky decomposition, the actual use of the RHS
C vector sizes is automatically controlled in RHSALL2. Furthermore, the
C sigma routines now use the full RHS size.
      IF (.NOT.IFCHOL) THEN
        WRITE(6,*)
        WRITE(6,*)' Size of vector buffers for coefficient arrays.'
        WRITE(6,*)' ICASE ISYM    NROW     NCOL     NBLK       Size'
        NVCMX=0
        DO ICASE=1,NCASES
          DO ISYM=1,NSYM
            NROW=NASUP(ISYM,ICASE)
            NCOL=MODVEC(ISYM,ICASE)
            NVC=NROW*NCOL
            NVCMX=MAX(NVCMX,NVC)
            NBLK=0
            IF(NVC.GT.0) NBLK=1+(NISUP(ISYM,ICASE)-1)/NCOL
            WRITE(6,'(2x,I2,3x,I2,3x,I16,3X,I16,3X,I16,3x,I16)')
     &                     ICASE,ISYM,NROW,NCOL,NBLK,NVC
          END DO
        END DO
        WRITE(6,*)
        WRITE(6,*)' Largest vector buffer size:',NVCMX
        WRITE(6,*)
      ELSE
        WRITE(6,*)
        WRITE(6,*)' Sizes of the coefficient arrays.'
        WRITE(6,'(2X,A4,2X,A4,2X,5X,A4,5X,2X,5X,A4,5X,2X,5X,A4,5X)')
     &          'CASE','SYM','NROW','NCOL','SIZE'
        NVCMX=0
        DO ICASE=1,NCASES
          DO ISYM=1,NSYM
            NROW=NASUP(ISYM,ICASE)
            NCOL=NISUP(ISYM,ICASE)
            NVC=NROW*NCOL
            NVCMX=MAX(NVCMX,NVC)
            WRITE(6,'(2x,I4,2x,I4,2x,I16,2X,I16,2X,I16)')
     &                     ICASE,ISYM,NROW,NCOL,NVC
          END DO
        END DO
        WRITE(6,*)
        WRITE(6,'(A,I14)')' Largest vector size: ',NVCMX
        WRITE(6,*)
      ENDIF
#endif

      IDV=0
      DO IVEC=1,MXVEC
        DO ICASE=1,NCASES
          DO ISYM=1,NSYM
            NAS=NASUP(ISYM,ICASE)
            NIS=NISUP(ISYM,ICASE)
C NINDEP() IS HERE SET PROVISIONALLY. IT WILL
C BE ADJUSTED FOR LINEAR DEPENDENCE LATER.
            NCOEF=NAS*NIS
            MXWRT=Max(1,NAS*MODVEC(ISYM,ICASE))
            NISCT=NCOEF/MXWRT+Min(1,Mod(NCOEF,MXWRT))
            LADDR=LIDSCT+MXSCT*(ISYM-1+8*(ICASE-1+MXCASE*(IVEC-1)))
            IF (NISCT.eq.0) IWORK(LADDR)=IDV
            If (NISCT.gt.MXSCT) Then
              write(6,*) 'EQCTL1 : NISCT= ',NISCT,' > MXSCT= ',MXSCT
              write(6,*) 'Please, increase MXSCT in eqsolv.fh'
              write(6,*) 'Do not forget to recompile Molcas afterwards.'
              Call Abend()
            EndIf
            DO ISCT=1,NISCT
              LSTA=MXWRT*(ISCT-1)
              LENGTH=RtoI*MIN(NCOEF-LSTA,MXWRT)
              LADDR=LIDSCT-1+ISCT+MXSCT*(ISYM-1+8*
     &                         (ICASE-1+MXCASE*(IVEC-1)))
              IWORK(LADDR)=IDV
              CALL IDAFILE(LUSOLV,0,IDUM,LENGTH,IDV)
            END DO
          END DO
        END DO
      END DO

C IDSCT(ISCT,ISYM,ICASE,IVEC) now gives the disk address, on LUSOLV,
C to the ISCT section of the (ISYM,ICASE) block of a vector
C containing expansion C coefficients of excitation operators.
C IVEC=1..MXVEC enumerates the different vectors needed to solve
C the equations.
      IDS=0
      DO ICASE=1,NCASES
        DO ISYM=1,NSYM
          IDSMAT(ISYM,ICASE)=IDS
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.GT.0) THEN
            NAS=NASUP(ISYM,ICASE)
*           NS=(NAS*(NAS+1))/2
            NS=NAS**2
            IF(ICASE.EQ.12) NS=1
            IF(ICASE.EQ.13) NS=1
            CALL DDAFILE(LUSBT ,0,DUMMY,NS,IDS)
          END IF
        END DO
      END DO

C IDSMAT(ISYM,ICASE) gives the disk address on LUSBT to the
C (ISYM,ICASE) section of the overlap matrix for the basis
C of excitation wave function terms: those that result from
C the basic excitation operators acting on Psi0.
      DO ICASE=1,NCASES
        DO ISYM=1,NSYM
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
*         NB=(NAS*(NAS+1))/2
          NB=NAS**2
          IF(ICASE.EQ.12) NB=1
          IF(ICASE.EQ.13) NB=1
          IDBMAT(ISYM,ICASE)=IDS
          NBD=NAS
          NID=NIS
          NT=NAS**2
          IF(ICASE.EQ.12) NT=1
          IF(ICASE.EQ.13) NT=1
          IDS1=IDS
         IF(NB.GT.0)then
         end if
          IF(NB.GT.0)CALL DDAFILE(LUSBT ,0,DUMMY,NB ,IDS1)
          IDS2=IDS
          IF(NBD.GT.0)CALL DDAFILE(LUSBT ,0,DUMMY,NBD,IDS2)
          IF(NID.GT.0)CALL DDAFILE(LUSBT ,0,DUMMY,NID,IDS2)
          IF(NBD.GT.0)CALL DDAFILE(LUSBT ,0,DUMMY,NBD,IDS2)
          IF(NID.GT.0)CALL DDAFILE(LUSBT ,0,DUMMY,NID,IDS2)
          IDTMAT(ISYM,ICASE)=IDS2
          IF(NT.GT.0)CALL DDAFILE(LUSBT ,0,DUMMY,NT,IDS2)
          IDSTMAT(ISYM,ICASE)=IDS2
          IF(NT.GT.0)CALL DDAFILE(LUSBT ,0,DUMMY,NT,IDS2)
          IDS=MAX(IDS1,IDS2)
        END DO
      END DO

C IDBMAT(ISYM,ICASE) is similar to IDSMAT() but gives disk address
C to the (ISYM,ICASE) diagonal block of matrix elements of H0
C over the excitation wave function terms.
C IDTMAT() similarly addresses transformation matrices that
C orthonormalize the S matrix blocks.

      RETURN
      END
