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
      SUBROUTINE get_Umn(PHP,EnIn,DHAM,
     &                   IPCSF,IPCNF,MXPDIM,
     &                   DTOC,IPRODT,ICONF,
     &                   IREFSM,ONEBOD,ECORE,NACTOB,
     &                   NCONF,NEL,NAEL,NBEL,
     &                   NPCSF,NPCNF,DIAG,TUVX,
     &                   iterSplit,ITER,
     &                   NTEST,ExFac,IREOTS)
C
C     ARGUMENTS :
C     ===========
C     PHP    : AA Block Hamiltonian un-dressed                  (Output)
C     EnIn   : energy value of the root selected                (Input)
C     DHAM   : Dressed AA block Hamiltonian                     (Output)
C     IPCSF  : CSF's order - Index Array -                      (Input)
C     IPCNF  : CNF's order - Index Array -                      (Input)
C     MXPDIM : Total number of CSFs                             (Input)
C     DTOC   : Transformation matrix between CSF's and DET's    (input)
C     IPRODT : Prototype determinants                           (input)
C     ICONF  : List of configurations                           (input)
C     IREFSM : symmetry of considered CI space                  (input)
C     Onebod : one body hamilton matrix in rectangular form     (input)
C     ECORE  : Core energy                                      (input)
C     NACTOB : Number of active orbitals                        (input)
C     NCONF  : Number of CNFs of symmetry IREFSM                (Input)
C     NEL    : total number of active electrons                 (input)
C     NAEL   : number of alpha active electron                  (input)
C     NBEL   : number of beta active electron                   (input)
C     NPCSF  : Number of CSFs in AA block                       (Input)
C     NPCNF  : Number of CNFs in AA block                       (Input)
C     DIAG   : Hamilton diagonal over CSF's                     (Input)
C     TUVX   : Two-electron integrals (MO space)                (Input)
C     NTEST  :
C     ExFac  :
C     IREOTS : Type => symmetry reordering array
C
      IMPLICIT REAL*8 (A-H,O-Z)

#include "spinfo.fh"
#include "WrkSpc.fh"

      DIMENSION PHP(NPCSF*(NPCSF+1)/2),DHAM(NPCSF*(NPCSF+1)/2)
      DIMENSION IPCSF(*),IPCNF(*),DIAG(*)
      DIMENSION DTOC(*),IPRODT(*),ICONF(*)
      DIMENSION ONEBOD(*)
      DIMENSION TUVX(*), IREOTS(*)

      IF( NTEST.GE.30 ) THEN
        WRITE(6,*) ' Input in get_Umn '
        WRITE(6,*) ' ================== '
        WRITE(6,*)
     &  ' Number of CNFs ',NCONF
        WRITE(6,*)
     &  ' Number of CSFs ',MXPDIM
        WRITE(6,*) ' Configurations included : '
        CALL IWRTMA(IPCNF,1,NCONF,1,NCONF)
        WRITE(6,*) ' CSFs included : '
        CALL IWRTMA(IPCSF,1,MXPDIM,1,MXPDIM)
        write(6,*) ' Number of CNFs in AA block:',NPCNF
        write(6,*) ' Number of CSFs in AA block:',NPCSF
      END IF

      Call FZero(DHAM,NPCSF*(NPCSF+1)/2)
* construct the Dressed Hamiltonian matrix
      MXCSFC = 0
      DO ITYP = 1, NTYP
        MXCSFC = MAX(MXCSFC,NCSFTP(ITYP))
      END DO
*      write(6,*) 'MXCSFC = ', MXCSFC
      call getmem('AuxDia','ALLO','REAL',ipAuxD,MXCSFC)
      call getmem('AuxVer', 'ALLO','REAL',ipAuxV, MXCSFC*NPCSF)
      call getmem('AuxCopy','ALLO','REAL',ipAuxC, MXCSFC*NPCSF)
      CALL GETMEM('EXHSCR','MAX','REAL',LW2,MXXWS)
      MXXWS = MXXWS/2
      CALL GETMEM('EXHSCR','ALLO','REAL',LW2,MXXWS)

      KACONF = LW2
      KLCONF = KACONF+NEL
      KRCONF = KLCONF+NEL
      KLAUXD = KRCONF+NEL
      KLPHPS = KLAUXD+MXCSFC*MXCSFC
      KLFREE = KLPHPS+MXCSFC*MXCSFC

C Shape of PHP matrix:
C
C           CNFL (n)
C
C      1  2  4  7 11 16 22 ...
C  C   -  3  5  8 12 17 23 ...
C  N   -  -  6  9 13 18 24 ...
C  F   -  -  - 10 14 19 25 ...
C  R   -  -  -  - 15 20 26 ...
C (m)  -  -  -  -  - 21 27 ...
C      -  -  -  -  -  - 28 ...
C      -  -  -  -  -  -  - ...

************************************************************************
* Now a loop over alpha will start to calculate:                       *
*   1) BB-Block DIAGONAL element 1/(En-H(alpha,alpha))                 *
*   2) AB-Block array H(m,alpha)
************************************************************************
      if (ITER.eq.1.and.iterSplit.eq.1) goto 29555
      IIAB = 1
      Do iAlpha = NPCNF+1,NCONF
*     ^ Loop over alpha
        Call FZero(Work(KLAUXD),MXCSFC*MXCSFC)
*        write(6,*)'iAlpha = ', iAlpha
        iKACONF=ip_of_iWork_d(Work(KACONF))
        CALL GETCNF_LUCIA(iWork(iKACONF),IATYP,IPCNF(iAlpha),
     &                    ICONF,IREFSM,NEL)
        NCSFA = NCSFTP(IATYP)
*        write(6,*)'NCSFA = ', NCSFA
        iKACONF=ip_of_iWork_d(Work(KACONF))
        CALL CNHCN(iWork(iKACONF),IATYP,iWork(iKACONF),
     &             IATYP,Work(KLAUXD),Work(KLFREE),
     &             NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB,TUVX,
     &             NTEST,ExFac,IREOTS)
        DO IIA = 1, NCSFA
          ILAI = IIA*IIA
          IF( NTEST.GE.30 ) write(6,*)'ILAI =', ILAI
*          Work(ipAuxD+IIA-1) = Work(KLAUXD+ILAI-1)
*          write(6,*) 'Work(ipAuxD+IIA-1)',Work(ipAuxD+IIA-1)
          Work(ipAuxD+IIA-1) = 1.0d0/
     &                           (EnIn-Work(KLAUXD+ILAI-1))
          IF( NTEST.GE.30 )
     &         write(6,*) 'Work(ipAuxD+IIA-1)',Work(ipAuxD+IIA-1)
        END DO
****************** 2) AB-Block Array (alpha Column) ********************
        IILB = 1
        do Mindex = 1, NPCNF
*       ^ Loop over AB-Block
          Call FZero(Work(KLAUXD),MXCSFC*MXCSFC)
*          write(6,*)'Mindex',Mindex
          iKLCONF=ip_of_iWork_d(Work(KLCONF))
          CALL GETCNF_LUCIA(iWork(iKLCONF),ILTYP,IPCNF(Mindex),
     &                      ICONF,IREFSM,NEL)
          NCSFL = NCSFTP(ILTYP)
*          write(6,*)'NCSFL = ', NCSFL
          iKACONF=ip_of_iWork_d(Work(KACONF))
          iKLCONF=ip_of_iWork_d(Work(KLCONF))
          CALL CNHCN(iWork(iKACONF),IATYP,iWork(iKLCONF),
     &               ILTYP,Work(KLAUXD),Work(KLFREE),
     &               NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB,TUVX,
     &               NTEST,ExFac,IREOTS)
          IF( NTEST.GE.30 ) then
            write(6,*)'M_Alpha elements'
            call wrtmat(Work(KLAUXD),MXCSFC,MXCSFC,MXCSFC,MXCSFC)
          end if
          DO IIL = 1, NCSFL
            DO IIA = 1, NCSFA
              IILACT = IILB-1+IIL
              IIAACT = NPCSF*(IIA-1)
              ILAI = (IIL-1)*NCSFA+IIA
*             ILAI = (IIA-1)*MXCSFC+IIL
*             ILAI = (IIL-1)*MXCSFC+IIA
              ILAOV = IILACT+IIAACT
              Work(ipAuxV+ILAOV-1) = Work(KLAUXD+ILAI-1)
*              write(6,*)'ILAI, ILAOV =', ILAI,ILAOV
              IF( NTEST.GE.30 )
     &         write(6,*) 'Work(ipAuxV+ILAOV-1)',Work(ipAuxV+ILAOV-1)
              Work(ipAuxC+ILAOV-1) = Work(ipAuxV+ILAOV-1)*
     &                               Work(ipAuxD+IIA-1)
              IF( NTEST.GE.30 )
     &            write(6,*) 'Work(ipAuxC+ILAOV-1)',Work(ipAuxC+ILAOV-1)
            END DO
          END DO
          IILB = IILB + NCSFL
        end do
*       ^ End loop over AB-Block
        IF( NTEST.GE.30 ) then
          write(6,*)'AB-Block Vertical Vector'
          call wrtmat(Work(ipAuxV),NPCSF,NCSFA,NPCSF,NCSFA)
          write(6,*)'AB-Block Vertical Vector times Daa'
          call wrtmat(Work(ipAuxC),NPCSF,NCSFA,NPCSF,NCSFA)
        end if
************************************************************************
        call dGeMM_Tri('N','T',NPCSF,NPCSF,NCSFA,
     &                   1.0d0,Work(ipAuxC),NPCSF,
     &                         Work(ipAuxV),NPCSF,
     &                   1.0d0,DHAM,NPCSF)
        IF( NTEST.GE.30 )
     &        call TRIPRT('correction to the AA block',' ',DHAM,NPCSF)
        IIAB = IIAB + NCSFA
      End Do
*     ^ End of the loop over iAlpha
************************************************************************
*    1. A-Block matrix element                                         *
************************************************************************
29555 continue
      IILB = 1
      do Nindex = 1, NPCNF
*     ^ Loop over the AA-block (vertical index)
        IF( NTEST.GE.30 ) write(6,*)'Nindex',Nindex
*       write(6,*)'IILB',IILB
        iKLCONF=ip_of_iWork_d(Work(KLCONF))
        CALL GETCNF_LUCIA(iWork(iKLCONF),ILTYP,IPCNF(Nindex),ICONF,
     &                    IREFSM,NEL)
        NCSFL = NCSFTP(ILTYP)
*        write(6,*)'NCSFL = ', NCSFL

        IIRB = 1
        do Mindex = 1, Nindex
*       ^ Loop over the AA-block (horizontal index)
          Call FZero(Work(KLPHPS),MXCSFC*MXCSFC)
*          write(6,*)'Nindex,Mindex', Nindex,Mindex
          iKRCONF=ip_of_iWork_d(Work(KRCONF))
          CALL GETCNF_LUCIA(iWork(iKRCONF),IRTYP,IPCNF(Mindex),
     &                      ICONF,IREFSM,NEL)
          NCSFR = NCSFTP(IRTYP)
*          write(6,*)'NCSFR = ', NCSFR
          iKLCONF=ip_of_iWork_d(Work(KLCONF))
          iKRCONF=ip_of_iWork_d(Work(KRCONF))
          CALL CNHCN(iWork(iKLCONF),ILTYP,iWork(iKRCONF),IRTYP,
     &               Work(KLPHPS),Work(KLFREE),NAEL,NBEL,
     &               ECORE,ONEBOD,IPRODT,DTOC,NACTOB,TUVX,NTEST,
     &               ExFac,IREOTS)
          IF( NTEST.GE.30 ) then
            write(6,*)'AA block elements'
            call wrtmat(Work(KLPHPS),MXCSFC,MXCSFC,MXCSFC,MXCSFC)
          End if
          DO IIL = 1, NCSFL
            IF(IILB.EQ.IIRB) THEN
              IIRMAX = IIL
            ELSE
              IIRMAX = NCSFR
            END IF
            DO IIR = 1, IIRMAX
              IIRACT = IIRB-1+IIR
              IILACT = IILB-1+IIL
*             ILRI = (IIR-1)*MXCSFC+IIL
              ILRI = (IIR-1)*NCSFL+IIL
*             ^ Forse questo e' quello giusto; la precedente formula
*               potrebbe essere fonte di BUGS! Vedremo!
              ILRO = ((IILACT*IILACT-IILACT)/2)+IIRACT
              PHP(ILRO) = Work(KLPHPS+ILRI-1)
*              write(6,*)'ILRI, ILRO =', ILRI,ILRO
*              write(6,*) 'PHP(ILRO)',PHP(ILRO)
            END DO
          END DO
          IIRB = IIRB + NCSFR
        end do
*       ^End loop over the AA-block (horizontal index)
        IILB = IILB + NCSFL
      end do
*     ^End loop over the AA-block (vertical index)

************************************************************************
* Let's add Hmn (PHP) to the correction (DHAM)
          call daxpy_(NPCSF*(NPCSF+1)/2,1.0d0,PHP,1,DHAM,1)
************************************************************************
      IF( NTEST.GE.30 ) then
        write(6,*)'AA-Block matrix un-dressed'
        call wrtmat(PHP,NPCSF*(NPCSF+1)/2,1,NPCSF*(NPCSF+1)/2,1)
        write(6,*)'AA-Block matrix dressed'
        call wrtmat(DHAM,NPCSF*(NPCSF+1)/2,1,NPCSF*(NPCSF+1)/2,1)
      end if

      IF( NTEST.GE.30 ) then
        call TRIPRT('AA block Hamiltonian Matrix un-dressed',' ',
     &              PHP,NPCSF)
        call TRIPRT('Dressed AA block Hamiltonian Matrix',' ',
     &              DHAM,NPCSF)
      END IF
      CALL GETMEM('EXHSCR','FREE','REAL',LW2,MXXWS)
      call getmem('AuxCopy','FREE','REAL',ipAuxC,MXCSFC*NPCSF)
      call getmem('AuxVer', 'FREE','REAL',ipAuxV, MXCSFC*NPCSF)
      call getmem('AuxDia', 'FREE','REAL',ipAuxD,MXCSFC)

      RETURN
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(DIAG)
      END
