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
      SUBROUTINE get_Cm(IPCSF,IPCNF,
     &                  MXPDIM,NCONF,
     &                  NPCSF,NPCNF,
     &                  Cn,lrootSplit,EnFin,
     &                  DTOC,IPRODT,ICONF,
     &                  IREFSM,ONEBOD,ECORE,
     &                  NACTOB,NEL,NAEL,NBEL,
     &                  DIAG,TUVX,
     &                  NTEST,ExFac,IREOTS,
     &                  FordSplit,
     &                  Ctot)

************** Author : GLMJ *****************
C
C     Obtain Cm coefficients out of the AA Block for root = lrootSplit
C
C     ARGUMENTS :
C     ===========
C     IPCSF      : CSF's order - Index Array -                    (Input)
C     IPCNF      : CNF's order - Index Array -                    (Input)
C     MXPDIM     : Total number of CSFs                           (Input)
C     NCONF      : Total Number of CNFs of symmetry IREFSM        (Input)
C     NPCSF      : Number of CSFs in AA block                     (Input)
C     NPCNF      : Number of CNFs in AA block                     (Input)
C     Cn         : AA Block CI-Coefficients for root selected     (Input)
C     lRootSplit : computed root                                  (Input)
C     EnFin      : Final Energy for the root selected             (Input)
C     DTOC       : Transformation matrix between CSF's and DET's  (input)
C     IPRODT     : Prototype determinants                         (input)
C     ICONF      : List of configurations                         (input)
C     IREFSM     : symmetry of considered CI space                (input)
C     Onebod     : one body hamilton matrix in rectangular form   (input)
C     ECORE      : Core energy                                    (input)
C     NACTOB     : Number of active orbitals                      (input)
C     NEL        : total number of active electrons               (input)
C     NAEL       : number of alpha active electron                (input)
C     NBEL       : number of beta active electron                 (input)
C     DIAG       : Hamilton diagonal over CSFs                    (Input)
C     TUVX       : Two-electron integrals (MO space)
C     IREOTS     : Type => symmetry reordering array
C     Ctot       : Vector of all nConf CI-coeff for a single root (Output)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "spinfo.fh"
#include "WrkSpc.fh"

C
      DIMENSION Ctot(MXPDIM),Cn(NPCSF),IPCSF(MXPDIM),IPCNF(NCONF)
      DIMENSION DIAG(*)
      DIMENSION DTOC(*),IPRODT(*),ICONF(*)
      DIMENSION ONEBOD(*)
      DIMENSION TUVX(*), IREOTS(*)
      LOGICAL FordSplit


      IF( NTEST.GE.30 ) THEN
        WRITE(6,*) ' Input in get_Cm '
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' Total Number of CNFs ',NCONF
        WRITE(6,*) ' Total Number of CSFs ',MXPDIM
        WRITE(6,*) ' CNFs included : '
        CALL IWRTMA(IPCNF,1,NCONF,1,NCONF)
        WRITE(6,*) ' CSFs included : '
        CALL IWRTMA(IPCSF,1,MXPDIM,1,MXPDIM)
        write(6,*) ' Number of CNFs in AA block:',NPCNF
        write(6,*) ' Number of CSFs in AA block:',NPCSF
        write(6,*)'Cn Coefficients'
        call wrtmat(Cn,NPCSF,1,NPCSF,1)
      END IF
************************************************************************
* get Cm coefficients corrected to the firts order                     *
* according to Lowdin's equations                                      *
************************************************************************
      If(FOrdSplit) then
        call get_Cm_(IPCSF,IPCNF,
     &               MXPDIM,NCONF,
     &               NPCSF,NPCNF,
     &               Cn,lrootSplit,EnFin,
     &               DTOC,IPRODT,ICONF,
     &               IREFSM,ONEBOD,ECORE,
     &               NACTOB,NEL,NAEL,NBEL,
     &               DIAG,TUVX,
     &               NTEST,ExFac,IREOTS,
     &               Ctot)
        Return
      End If

************************************************************************
* get Cm coefficients corrected to the zeroth order                    *
* according to Lowdin's equations                                      *
************************************************************************
      call Fzero(Ctot,MXPDIM)

      MXCSFC = 0
      DO ITYP = 1, NTYP
        MXCSFC = MAX(MXCSFC,NCSFTP(ITYP))
      END DO

      IBblockV = MXPDIM-NPCSF
      call getmem('AuxDia' ,'ALLO','REAL', ipAuxD,  MXCSFC)
      call getmem('AuxGa'  ,'ALLO','REAL',ipAuxGa,  MXCSFC)
      call getmem('AuxGaTi','ALLO','REAL',ipAuxGaTi,MXCSFC)
      call getmem('AuxVer' ,'ALLO','REAL',ipAuxV,   MXCSFC*NPCSF)
      CALL GETMEM('EXHSCR','MAX','REAL',LW2,MXXWS)
      CALL GETMEM('EXHSCR','ALLO','REAL',LW2,MXXWS)

      KACONF = LW2
      KLCONF = KACONF+NEL
      KLAUXD = KLCONF+NEL
      KLFREE = KLAUXD+MXCSFC*MXCSFC

C Shape of matrix:
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
C
C Splitting:
C
C          |
C     AA   |   AB
C  --------|--------
C     BA   |   BB
C          |

      call cwtime(C_AlphaLoop1,W_AlphaLoop1)
      C_ComputeH_AB = 0.0d0
      W_ComputeH_AB = 0.0d0
      C_ComputeH_BB = 0.0d0
      W_ComputeH_BB = 0.0d0
      C_Oper = 0.0d0
      W_Oper = 0.0d0

      IIAB = 1
      Do iAlpha = NPCNF+1,NCONF
C       Call FZero(Work(KLAUXD),MXCSFC*MXCSFC)
        call cwtime(C_computeH_AB1,W_computeH_AB1)
        IF( NTEST.GE.30 )  write(6,*)'iAlpha = ', iAlpha
        CALL GETCNF_LUCIA(Work(KACONF),IATYP,IPCNF(iAlpha),
     &                    ICONF,IREFSM,NEL)
        NCSFA = NCSFTP(IATYP)
        IF( NTEST.GE.30 ) write(6,*)'NCSFA = ', NCSFA
************************************************************************
*                        BB-Block DIAGONAL Elements                    *
************************************************************************
        CALL CNHCN(Work(KACONF),IATYP,Work(KACONF),
     &             IATYP,Work(KLAUXD),Work(KLFREE),
     &             NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB,TUVX,
     &             NTEST,ExFac,IREOTS)
        IF( NTEST.GE.30 ) then
          write(6,*)'Alpha_Alpha elements in BB-block'
          call wrtmat(Work(KLAUXD),MXCSFC,MXCSFC,MXCSFC,MXCSFC)
        end if
        DO IIA = 1, NCSFA
          ILAI = IIA*IIA
          Work(ipAuxD+IIA-1) = Work(KLAUXD+ILAI-1)
*         Write(6,*)'ILAI =', ILAI
          IF( NTEST.GE.30 ) then
            write(6,*) 'Work(ipAuxD+IIA-1)',Work(ipAuxD+IIA-1)
          end if
        END DO

        IILB = 1
        do Mindex = 1, NPCNF
*       ^ Loop over AB-Block
C         Call FZero(Work(KLAUXD),MXCSFC*MXCSFC)
          IF( NTEST.GE.30 ) then
            write(6,*)'Mindex in AB-Block',Mindex
          End if
          CALL GETCNF_LUCIA(Work(KLCONF),ILTYP,IPCNF(Mindex),
     &                      ICONF,IREFSM,NEL)
          NCSFL = NCSFTP(ILTYP)
          IF( NTEST.GE.30 ) write(6,*)'NCSFL = ', NCSFL
          CALL CNHCN(Work(KACONF),IATYP,Work(KLCONF),
     &               ILTYP,Work(KLAUXD),Work(KLFREE),
     &               NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB,TUVX,
     &               NTEST,ExFac,IREOTS)
          IF( NTEST.GE.30 ) then
            write(6,*)'M_Alpha elements'
            call wrtmat(Work(KLAUXD),MXCSFC,MXCSFC,MXCSFC,MXCSFC)
          End if
          DO IIA = 1, NCSFA
            DO IIL = 1, NCSFL
              IILACT = IILB-1+IIL
              IIAACT = NPCSF*(IIA-1)
              ILAI = (IIL-1)*NCSFA+IIA
*              ILAI = (IIA-1)*MXCSFC+IIL
              ILAOV = IILACT+IIAACT
              Work(ipAuxV+ILAOV-1) = Work(KLAUXD+ILAI-1)
              IF( NTEST.GE.30 ) then
                Write(6,*)'ILAI, ILAOV =', ILAI,ILAOV
                write(6,*) 'Work(ipAuxV+ILAOV-1)',Work(ipAuxV+ILAOV-1)
              END IF
            END DO
          END DO
          IILB = IILB + NCSFL
        end do
*       ^ End loop over AB-Block
        call cwtime(C_computeH_AB2,W_computeH_AB2)
        C_ComputeH_AB = C_ComputeH_AB + C_computeH_AB2 - C_computeH_AB1
        W_ComputeH_AB = W_ComputeH_AB + W_computeH_AB2 - W_computeH_AB1

        IF( NTEST.GE.30 )  then
          write(6,*)'AB-Block Vertical Vector'
          call wrtmat(Work(ipAuxV),NPCSF,NCSFA,NPCSF,NCSFA)
        end if

************************************************************************
*    Compute:                                                          *
*    1.  G_Alpha_tilde                                                 *
*    2.  G_Alpha                                                       *
************************************************************************
        call cwtime(C_oper1,W_oper1)
        do IIA = 1,NCSFA
          Work(ipAuxGaTi+IIA-1) =
     &      ddot_(NPCSF,Work(ipAuxV+(IIA-1)*NPCSF),1,Cn,1)
          Work(ipAuxGa+IIA-1) =
     &      Work(ipAuxGaTi+IIA-1)/(EnFin-Work(ipAuxD+IIA-1))
          IF( NTEST.GE.30 ) then
            write(6,*) 'Work(ipAuxGaTi+IIA-1)',Work(ipAuxGaTi+IIA-1)
            write(6,*) 'Work(ipAuxGa  +IIA-1)',Work(ipAuxGa  +IIA-1)
          end if
        end do
        call cwtime(C_oper2,W_oper2)
        C_oper = C_oper + C_oper2 - C_oper1
        W_oper = W_oper + W_oper2 - W_oper1

        do IIA = 1, NCSFA
          Ctot(NPCSF+IIAB+IIA-1)=
     &           Ctot(NPCSF+IIAB+IIA-1)+ Work(ipAuxGa+IIA-1)
          IF( NTEST.GE.30 ) then
            write(6,*)'Ctot '
            call wrtmat(Ctot,MXPDIM,1,MXPDIM,1)
          end if
        end do
        call cwtime(C_oper2,W_oper2)
        C_oper = C_oper + C_oper2 - C_oper1
        W_oper = W_oper + W_oper2 - W_oper1

        IIAB = IIAB + NCSFA
      End Do
*     ^ End of the loop over iAlpha
      IF( NTEST.GE.30 ) then
       call cwtime(C_AlphaLoop2,W_AlphaLoop2)
       write(6,*) 'Total time needed to get_Cm in Alpha Loop'
       write(6,*) 'CPU timing : ', C_AlphaLoop2 - C_AlphaLoop1
       write(6,*) 'W. timing  : ', W_AlphaLoop2 - W_AlphaLoop1

       write(6,*) 'Total time to read H_AB :'
       write(6,*) 'CPU timing : ', C_ComputeH_AB
       write(6,*) 'W. timing  : ', W_ComputeH_AB

       write(6,*) 'Total time to calculate (ddot+dscal+daxpy) :'
       write(6,*) 'CPU timing : ', C_Oper
       write(6,*) 'W. timing  : ', W_Oper
      end if
      call dcopy_(NPCSF,Cn,1,Ctot,1)
      IF( NTEST.GE.30 ) then
        write(6,*)'final Ctot vector'
        call wrtmat(Ctot,MXPDIM,1,MXPDIM,1)
      end if
      CALL GETMEM('EXHSCR','FREE','REAL',LW2,MXXWS)
      call getmem('AuxVer' ,'FREE','REAL',ipAuxV,   MXCSFC*NPCSF)
      call getmem('AuxGaTi','FREE','REAL',ipAuxGaTi,MXCSFC)
      call getmem('AuxGa'  ,'FREE','REAL',ipAuxGa,  MXCSFC)
      call getmem('AuxDia' ,'FREE','REAL', ipAuxD,  MXCSFC)

      RETURN
      END

      SUBROUTINE get_Cm_(IPCSF,IPCNF,
     &                   MXPDIM,NCONF,
     &                   NPCSF,NPCNF,
     &                   Cn,lrootSplit,EnFin,
     &                   DTOC,IPRODT,ICONF,
     &                   IREFSM,ONEBOD,ECORE,
     &                   NACTOB,NEL,NAEL,NBEL,
     &                   DIAG,TUVX,
     &                   NTEST,ExFac,IREOTS,
     &                   Ctot)
C
C     Obtain Cm coefficients out of the AA Block for root = lrootSplit
C     to the first order in Lowdin equation
C
C     ARGUMENTS :
C     ===========
C     IPCSF      : CSF's order - Index Array -                    (Input)
C     IPCNF      : CNF's order - Index Array -                    (Input)
C     MXPDIM     : Total number of CSFs                           (Input)
C     NCONF      : Total Number of CNFs of symmetry IREFSM        (Input)
C     NPCSF      : Number of CSFs in AA block                     (Input)
C     NPCNF      : Number of CNFs in AA block                     (Input)
C     Cn         : AA Block CI-Coefficients for root selected     (Input)
C     lRootSplit : computed root                                  (Input)
C     EnFin      : Final Energy for the root selected             (Input)
C     DTOC       : Transformation matrix between CSF's and DET's  (input)
C     IPRODT     : Prototype determinants                         (input)
C     ICONF      : List of configurations                         (input)
C     IREFSM     : symmetry of considered CI space                (input)
C     Onebod     : one body hamilton matrix in rectangular form   (input)
C     ECORE      : Core energy                                    (input)
C     NACTOB     : Number of active orbitals                      (input)
C     NEL        : total number of active electrons               (input)
C     NAEL       : number of alpha active electron                (input)
C     NBEL       : number of beta active electron                 (input)
C     DIAG       : Hamilton diagonal over CSFs                    (Input)
C     TUVX       : Two-electron integrals (MO space)
C     IREOTS     : Type => symmetry reordering array
C     Ctot       : Vector of all nConf CI-coeff for a single root (Output)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "spinfo.fh"
#include "WrkSpc.fh"

C
      DIMENSION Ctot(MXPDIM),Cn(NPCSF),IPCSF(MXPDIM),IPCNF(NCONF)
      DIMENSION DIAG(*)
      DIMENSION DTOC(*),IPRODT(*),ICONF(*)
      DIMENSION ONEBOD(*)
      DIMENSION TUVX(*), IREOTS(*)


      IF( NTEST.GE.30 ) THEN
        WRITE(6,*) ' Input in get_Cm '
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' Total Number of CNFs ',NCONF
        WRITE(6,*) ' Total Number of CSFs ',MXPDIM
        WRITE(6,*) ' CNFs included : '
        CALL IWRTMA(IPCNF,1,NCONF,1,NCONF)
        WRITE(6,*) ' CSFs included : '
        CALL IWRTMA(IPCSF,1,MXPDIM,1,MXPDIM)
        write(6,*) ' Number of CNFs in AA block:',NPCNF
        write(6,*) ' Number of CSFs in AA block:',NPCSF
        write(6,*)'Cn Coefficients'
        call wrtmat(Cn,NPCSF,1,NPCSF,1)
      END IF

      call Fzero(Ctot,MXPDIM)

      MXCSFC = 0
      DO ITYP = 1, NTYP
        MXCSFC = MAX(MXCSFC,NCSFTP(ITYP))
      END DO

      IBblockV = MXPDIM-NPCSF
      call getmem('AuxDia' ,'ALLO','REAL', ipAuxD,  MXCSFC)
      call getmem('AuxGa'  ,'ALLO','REAL',ipAuxGa,  MXCSFC)
      call getmem('AuxGaTi','ALLO','REAL',ipAuxGaTi,MXCSFC)
      call getmem('AuxVer' ,'ALLO','REAL',ipAuxV,   MXCSFC*NPCSF)
      call getmem('AuxBB'  ,'ALLO','REAL',ipAuxBB,  MXCSFC*IBblockV)
      CALL GETMEM('EXHSCR','MAX','REAL',LW2,MXXWS)
      CALL GETMEM('EXHSCR','ALLO','REAL',LW2,MXXWS)

      KACONF = LW2
      KLCONF = KACONF+NEL
      KLAUXD = KLCONF+NEL
      KLFREE = KLAUXD+MXCSFC*MXCSFC

C Shape of matrix:
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
C
C Splitting:
C
C          |
C     AA   |   AB
C  --------|--------
C     BA   |   BB
C          |

      call cwtime(C_AlphaLoop1,W_AlphaLoop1)
      C_ComputeH_AB = 0.0d0
      W_ComputeH_AB = 0.0d0
      C_ComputeH_BB = 0.0d0
      W_ComputeH_BB = 0.0d0
      C_Oper = 0.0d0
      W_Oper = 0.0d0

      IIAB = 1
      Do iAlpha = NPCNF+1,NCONF
C       Call FZero(Work(KLAUXD),MXCSFC*MXCSFC)
        call cwtime(C_computeH_AB1,W_computeH_AB1)
        IF( NTEST.GE.30 )  write(6,*)'iAlpha = ', iAlpha
        CALL GETCNF_LUCIA(Work(KACONF),IATYP,IPCNF(iAlpha),
     &                    ICONF,IREFSM,NEL)
        NCSFA = NCSFTP(IATYP)
        IF( NTEST.GE.30 ) write(6,*)'NCSFA = ', NCSFA
************************************************************************
*                        BB-Block DIAGONAL Elements                    *
************************************************************************
        CALL CNHCN(Work(KACONF),IATYP,Work(KACONF),
     &             IATYP,Work(KLAUXD),Work(KLFREE),
     &             NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB,TUVX,
     &             NTEST,ExFac,IREOTS)
        IF( NTEST.GE.30 ) then
          write(6,*)'Alpha_Alpha elements in BB-block'
          call wrtmat(Work(KLAUXD),MXCSFC,MXCSFC,MXCSFC,MXCSFC)
        end if
        DO IIA = 1, NCSFA
          ILAI = IIA*IIA
          Work(ipAuxD+IIA-1) = Work(KLAUXD+ILAI-1)
*         Write(6,*)'ILAI =', ILAI
          IF( NTEST.GE.30 ) then
            write(6,*) 'Work(ipAuxD+IIA-1)',Work(ipAuxD+IIA-1)
          end if
        END DO

        IILB = 1
        do Mindex = 1, NPCNF
*       ^ Loop over AB-Block
C         Call FZero(Work(KLAUXD),MXCSFC*MXCSFC)
          IF( NTEST.GE.30 ) then
            write(6,*)'Mindex in AB-Block',Mindex
          End if
          CALL GETCNF_LUCIA(Work(KLCONF),ILTYP,IPCNF(Mindex),
     &                      ICONF,IREFSM,NEL)
          NCSFL = NCSFTP(ILTYP)
          IF( NTEST.GE.30 ) write(6,*)'NCSFL = ', NCSFL
          CALL CNHCN(Work(KACONF),IATYP,Work(KLCONF),
     &               ILTYP,Work(KLAUXD),Work(KLFREE),
     &               NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB,TUVX,
     &               NTEST,ExFac,IREOTS)
          IF( NTEST.GE.30 ) then
            write(6,*)'M_Alpha elements'
            call wrtmat(Work(KLAUXD),MXCSFC,MXCSFC,MXCSFC,MXCSFC)
          End if
          DO IIA = 1, NCSFA
            DO IIL = 1, NCSFL
              IILACT = IILB-1+IIL
              IIAACT = NPCSF*(IIA-1)
              ILAI = (IIL-1)*NCSFA+IIA
*              ILAI = (IIA-1)*MXCSFC+IIL
              ILAOV = IILACT+IIAACT
              Work(ipAuxV+ILAOV-1) = Work(KLAUXD+ILAI-1)
              IF( NTEST.GE.30 ) then
                Write(6,*)'ILAI, ILAOV =', ILAI,ILAOV
                write(6,*) 'Work(ipAuxV+ILAOV-1)',Work(ipAuxV+ILAOV-1)
              END IF
            END DO
          END DO
          IILB = IILB + NCSFL
        end do
*       ^ End loop over AB-Block
        call cwtime(C_computeH_AB2,W_computeH_AB2)
        C_ComputeH_AB = C_ComputeH_AB + C_computeH_AB2 - C_computeH_AB1
        W_ComputeH_AB = W_ComputeH_AB + W_computeH_AB2 - W_computeH_AB1

        IF( NTEST.GE.30 )  then
          write(6,*)'AB-Block Vertical Vector'
          call wrtmat(Work(ipAuxV),NPCSF,NCSFA,NPCSF,NCSFA)
        end if

************************************************************************
*    Compute:                                                          *
*    1.  G_Alpha_tilde                                                 *
*    2.  G_Alpha                                                       *
************************************************************************
        call cwtime(C_oper1,W_oper1)
        do IIA = 1,NCSFA
          Work(ipAuxGaTi+IIA-1) =
     &      ddot_(NPCSF,Work(ipAuxV+(IIA-1)*NPCSF),1,Cn,1)
          Work(ipAuxGa+IIA-1) =
     &      Work(ipAuxGaTi+IIA-1)/(EnFin-Work(ipAuxD+IIA-1))
          IF( NTEST.GE.30 ) then
            write(6,*) 'Work(ipAuxGaTi+IIA-1)',Work(ipAuxGaTi+IIA-1)
            write(6,*) 'Work(ipAuxGa  +IIA-1)',Work(ipAuxGa  +IIA-1)
          end if
        end do
        call cwtime(C_oper2,W_oper2)
        C_oper = C_oper + C_oper2 - C_oper1
        W_oper = W_oper + W_oper2 - W_oper1

************************************************************************
*  Loop over BB block to get the first order correction to Cm.         *
*  If not required by the USER (by using a keyword) it will be skipped!*
************************************************************************

        xmaxGaTi = 0.0d0
        do IIA = 1, NCSFA
          xmaxGaTi = MAX(xmaxGaTi,abs(Work(ipAuxGaTi+IIA-1)))
        end do
        if(xmaxGaTi.ge.1.0d-12) then
*       ^ let's try this trick to make SplitCAS faster!
          call cwtime(C_computeH_BB1,W_computeH_BB1)
          call Fzero(Work(ipAuxBB),MXCSFC*IBblockV)
          IILB = 1
          do Mindex = NPCNF+1, NCONF
*         ^ Loop over BB-Block
C           Call FZero(Work(KLAUXD),MXCSFC*MXCSFC)
*           IF( NTEST.GE.30 ) write(6,*)'Mindex in BB-Block',Mindex
            CALL GETCNF_LUCIA(Work(KLCONF),ILTYP,IPCNF(Mindex),
     &                        ICONF,IREFSM,NEL)
            NCSFL = NCSFTP(ILTYP)
*           IF( NTEST.GE.30 ) write(6,*)'NCSFL = ', NCSFL
            CALL CNHCN(Work(KACONF),IATYP,Work(KLCONF),
     &                 ILTYP,Work(KLAUXD),Work(KLFREE),
     &                 NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB,TUVX,
     &                 NTEST,ExFac,IREOTS)
*           IF( NTEST.GE.30 ) then
*             write(6,*)'M_Alpha elements in BB-block'
*             call wrtmat(Work(KLAUXD),MXCSFC,MXCSFC,MXCSFC,MXCSFC)
*           end if
            DO IIA = 1, NCSFA
              DO IIL = 1, NCSFL
                IILACT = IILB-1+IIL
*               IIAACT = NPCSF*(IIA-1)
                IIAACT = iBblockV*(IIA-1)
                ILAI = (IIL-1)*NCSFA+IIA
*               ILAI = (IIA-1)*MXCSFC+IIL
                ILAOV = IILACT+IIAACT
                Work(ipAuxBB+ILAOV-1) = Work(KLAUXD+ILAI-1)
                if(iAlpha.eq.Mindex.and.IIA.eq.IIL) then
*                 write(6,*)'iAlpha, IIAB', iAlpha, IIAB
                  Work(ipAuxBB+ILAOV-1)=0.0d0
                end if
*               Write(6,*)'ILAI, ILAOV =', ILAI,ILAOV
*               IF( NTEST.GE.30 ) then
*                 write(6,*) 'Work(ipAuxBB+ILAOV-1)',Work(ipAuxBB+ILAOV-1)
*               END IF
              END DO
            END DO
            IILB = IILB + NCSFL
          end do
*         ^ End Loop over BB-Block
*         IF( NTEST.GE.30 ) then
*           write(6,*)'BB-Block Vertical Vector'
*           call wrtmat(Work(ipAuxBB),iBblockV,NCSFA,iBblockV,NCSFA)
*         end if
          call cwtime(C_computeH_BB2,W_computeH_BB2)
          C_ComputeH_BB= C_ComputeH_BB + C_computeH_BB2 - C_computeH_BB1
          W_ComputeH_BB= W_ComputeH_BB + W_computeH_BB2 - W_computeH_BB1

          call cwtime(C_oper1,W_oper1)
          do IIA = 1, NCSFA
            call dscal_(iBblockV,Work(ipAuxGa + IIA-1),
     &                          Work(ipAuxBB +(IIA-1)*iBblockV),1)
            IF( NTEST.GE.30 ) then
              write(6,*)'BB-Block Vertical Vector times Ga'
              call wrtmat(Work(ipAuxBB),iBblockV,NCSFA,iBblockV,NCSFA)
            end if
           call daxpy_(iBblockV,1.0d0,Work(ipAuxBB +(IIA-1)*iBblockV),1,
     &                                Ctot(NPCSF+1),1)
            IF( NTEST.GE.30 ) then
              write(6,*)'Ctot correction'
              call wrtmat(Ctot,MXPDIM,1,MXPDIM,1)
            end if
          end do
        end if
*       ^ End of trick to make SplitCAS faster

        do IIA = 1, NCSFA
          Ctot(NPCSF+IIAB+IIA-1)=
     &           Ctot(NPCSF+IIAB+IIA-1)+ Work(ipAuxGaTi+IIA-1)
          IF( NTEST.GE.30 ) then
            write(6,*)'Ctot '
            call wrtmat(Ctot,MXPDIM,1,MXPDIM,1)
          end if
        end do
        call cwtime(C_oper2,W_oper2)
        C_oper = C_oper + C_oper2 - C_oper1
        W_oper = W_oper + W_oper2 - W_oper1

        IIAB = IIAB + NCSFA
      End Do
*     ^ End of the loop over iAlpha
      IF( NTEST.GE.30 ) THEN
        call cwtime(C_AlphaLoop2,W_AlphaLoop2)
        write(6,*) 'Total time needed to get_Cm in Alpha Loop'
        write(6,*) 'CPU timing : ', C_AlphaLoop2 - C_AlphaLoop1
        write(6,*) 'W. timing  : ', W_AlphaLoop2 - W_AlphaLoop1

        write(6,*) 'Total time to read H_AB :'
        write(6,*) 'CPU timing : ', C_ComputeH_AB
        write(6,*) 'W. timing  : ', W_ComputeH_AB

        write(6,*) 'Total time to read H_BB :'
        write(6,*) 'CPU timing : ', C_ComputeH_BB
        write(6,*) 'W. timing  : ', W_ComputeH_BB

        write(6,*) 'Total time to calculate (ddot+dscal+daxpy) :'
        write(6,*) 'CPU timing : ', C_Oper
        write(6,*) 'W. timing  : ', W_Oper
      END IF

      call cwtime(C_last1,W_last1)
      IIAB = 1
      Do Mindex = NPCNF+1,NCONF
C       Call FZero(Work(KLAUXD),MXCSFC*MXCSFC)
        IF( NTEST.GE.30 ) write(6,*)'Mindex last do loop = ', Mindex
        CALL GETCNF_LUCIA(Work(KACONF),IATYP,IPCNF(Mindex),
     &                    ICONF,IREFSM,NEL)
        NCSFA = NCSFTP(IATYP)
        IF( NTEST.GE.30 ) write(6,*)'NCSFA = ', NCSFA
        CALL CNHCN(Work(KACONF),IATYP,Work(KACONF),
     &             IATYP,Work(KLAUXD),Work(KLFREE),
     &             NAEL,NBEL,ECORE,ONEBOD,IPRODT,DTOC,NACTOB,TUVX,
     &             NTEST,ExFac,IREOTS)
        IF( NTEST.GE.30 ) then
          write(6,*)'Alpha_Alpha elements in BB-block'
          call wrtmat(Work(KLAUXD),MXCSFC,MXCSFC,MXCSFC,MXCSFC)
        END IF
        DO IIA = 1, NCSFA
          ILAI = IIA*IIA
*         Work(ipAuxD+IIA-1) = Work(KLAUXD+ILAI-1)
          IF( NTEST.GE.30 ) then
            write(6,*)'Work(KLAUXD+ILAI-1)', Work(KLAUXD+ILAI-1)
          END IF
*         Work(ipAuxD+IIA-1) = 1/(EnFin-Work(KLAUXD+ILAI-1))
          Ctot(NPCSF+IIAB+IIA-1) = Ctot(NPCSF+IIAB+IIA-1)/
     &                                  (EnFin-Work(KLAUXD+ILAI-1))
        END DO
        IIAB = IIAB +NCSFA
      End Do
      IF( NTEST.GE.30 ) then
        call cwtime(C_last2,W_last2)
        write(6,*) 'Time last iteration over Hbb :'
        write(6,*) 'CPU timing : ', C_last2 - C_last1
        write(6,*) 'W. timing  : ', W_last2 - W_last1
      END IF

      call dcopy_(NPCSF,Cn,1,Ctot,1)
      IF( NTEST.GE.30 ) then
        write(6,*)'final Ctot vector'
        call wrtmat(Ctot,MXPDIM,1,MXPDIM,1)
      end if
      CALL GETMEM('EXHSCR','FREE','REAL',LW2,MXXWS)
      call getmem('AuxBB'  ,'FREE','REAL',ipAuxBB,  MXCSFC*IBblockV)
      call getmem('AuxVer' ,'FREE','REAL',ipAuxV,   MXCSFC*NPCSF)
      call getmem('AuxGaTi','FREE','REAL',ipAuxGaTi,MXCSFC)
      call getmem('AuxGa'  ,'FREE','REAL',ipAuxGa,  MXCSFC)
      call getmem('AuxDia' ,'FREE','REAL', ipAuxD,  MXCSFC)

      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(lrootsPlit)
        CALL Unused_real_array(DIAG)
      END IF
      END
