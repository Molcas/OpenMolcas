!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
Module MCLR_Data

! Stuff from Pointers.fh
Integer ipMat(8,8),ipMatLT(8,8),ipCM(8),ipMC(8,8),ipmatba(8,8),ipMO(8,8,8),iADMO(8,8,8)
Integer nDens,nDensLT,nCMO,nscrtch,ipCI,nDens2,ipf0,ndensc,na(8),nb(8)
Integer nmba,nna,n1dens,n2dens,nconf1,nacpar,nacpr2

! Stuff from spin_mclr.fh
real*8     rms,rbetaa,rbetas,ralpha

! Stuff from genop.fh
!. Type of operator in action
!        I12        :         2:Both on and two electron integrals
!        IST        :
!        Square     :         Integrals square/triangular
Logical square
Integer I12,IST

!Stuff from machine.fh
integer  nrec

!Stuff from csfsd.fh
Integer i1,iAnders,lconf,lldet,iAllo

!Stuff from incdia.fh
integer ipdia

! Print flags from MCLR_Data.fh
Integer IPRSTR,IPRCIX,IPRORB,IPRDIA,IPRXT,IPRANA

!Stuff from sa.fh
Logical Fancy_Preconditioner,SA,esterr,isNAC,override
Logical isMECIMSPD
Integer istate,irlxroot,NACstates(2),NSSA(2)

!stuff from MCLR_Data.fh
#include "LenIn.fh"
Integer, Parameter:: LENIN6=LENIN+6
Integer lDisp(8)
Integer, Parameter,Private:: Mxdccc=500
Character ChDisp(Mxdccc*3)*(LENIN+6)
Character(LEN=8) SwLbl(Mxdccc)
integer   dspvec(mxdccc),nhess
Private LenIN,LENIN6

!
! Stuff from MCLR_Data.fh
!        iRefSm        :        Reference symmetry
!        iRefML        :        ML value
!        iRefPA        :        Parity
!        iRefL         :        L value
!        MS2           :        2*MS
!        MULTS         :        Spin multiplicity
!        nRoot         :        Number of roots
!        IDC           :        Ms combinations
!        PSSIGN        :        Ms combination PS factor
!        PLSIGN        :        Ms combination PL factor
!        IntSel        :        Internal space
!        iAlign        :        Not in use
!        Ethers        :        Threshold (energy cont)
!        Cthres        :        Threshold (coeff)
!        NGenSym       :        Number of reference symmetries
!        iGenSym       :        Reference symmetries
!        InvSym        :        Ger/UnGer inv sym
!        iKram         :        Kramer symmetry
!. CSTATE
Integer IREFSM,IREFML,IREFPA,IREFL,MS2,MULTS,NROOT,IDC,INTSEL,IALIGN,   &
        NGENSYM,IGENSYM(100),INVSYM,IKRAM
Real*8 PSSIGN,PLSIGN, ETHRES,CTHRES

! Stuff from Files_mclr.fh
Character(LEN=8) FnOne,FnJob,FnTwo,FnMol,FnRlx,FnMck,FnTemp,           &
                 FnHlf1,FnHlf2,FNHlf3,FNHlf4,FnPT2,                    &
                 FNTrI1,FNTrI2,FNTRi3,FNTRI4,FNTRI5,FNCSF2SD,FnMOTRA,  &
                 FnQDAT
Integer LuOne,LuJob,LuTwo,LuMol,LuRlx,LuMck,LuTemp,                    &
        LUHLF1,LUHlf2,LUHlf3,LUHLF4,LuPt2,                      &
        LUTrI1,LUTrI2,LUTRi3,LUTRI4,LUTRI5,LUCSF2SD,LuMOTRA,    &
        LuQDAT
Logical Nofile

! Stuff from cicisp_mclr.fh
#include "detdim.fh"
Private MXPIRR,MXPOBS,MXPR4T,MXINKA,MXPORB,MXPXOT,MXPXST,MXPSHL, &
        MXPL,MXPXT,MXPICI,MXPSTT,MXPCSM,MXPCTP,MXCNSM,MXPWRD, &
        MXNMS,MTYP,MXPNGAS,MXPNSMST,MXPPTSPC
!
! icisps **
! smost ***
!     niCISp                :        Number of CI-spaces (1) *
!     iAStFI                :       Alpha string type for internal CI space *
!     iBStFI                :       Beta string type for internal CI space (iBZTP)*
!     iActi                :        Active/Inactive space  (1)*
!     MnR1IC                :         Minimum number of electrons in RAS1*
!     MxR1IC                :        Maximum    -              "             - RAS1*
!     MnR3IC                :       Minimum    -             "       - RAS3*
!     MxR3IC                :        Maximum    -              "             - RAS3*
!     iZCI                :        Internal zero order space *
!     iRCI                :        Number of zero order space *
!     nElCI                :           Number of electrons per CI space *
!     nAElCI                :           Number of alpha electrons per CI space *
!     nBElCI                :           Number of beta electrons per CI space *
!     xispsm                :       Number of det. for each (symmetry,CI space) **
!     ismost                :       Symmetry operator ASym=ismost(BSym,iTOTSM) ***(istead of ieor)
!     MXSB                 :       Largest symmetry block **
!     MXSOOB                :       Largest block        **
!
Integer NICISP,IASTFI(MXPICI),IBSTFI(MXPICI),IACTI(MXPICI),MNR1IC(MXPICI),MXR1IC(MXPICI),    &
        MNR3IC(MXPICI),MXR3IC(MXPICI),IZCI,IRCI(3,7,7),NELCI(MXPICI),NAELCI(MXPICI),NBELCI(MXPICI), &
        ISMOST(MXPCSM,MXPCSM),MXSB,MXSOOB,ldet,lcsf
Real*8 XISPSM(MXPCSM,MXPICI)

! Stuff from crun_mclr.fh
!
!        MAXIT        :        Max iterations (NOT IN USE)
!        iRestr        :         Not in use
!       Intimp        :        Import mode for integrals
!         MxP1        :        Not in use
!         MxP2        :        Not in use
!        MxQ        :         Not in use
!              Incore        :        Integrals in core
!        MxCIV        :        Not in use
!        iCIStr        :        Storing mode for CI-vector
!       NoCSF        :        Transform to CSF base
!        iDiag        :        PICO
!        NoInt        :        Not in use
!        NoPair        :        Not in use
!        iH1SO        :        Not in use
!        iDelMS        :        Not in use
!       iSOSym        :        Not in use (symmetry of spin orbit)
!        inideg        :        Not in use
!        IFINSD        :        Not in use
!        iprerel        :        Not in use
!       nomspa  :        Not in use
!       ipredzn        :        Not in use
!       isupsp        :         not in use
!       iredstp :       Not in use
!       idmpin  :        Not in use
!       Maxrit  :         Not in use
!       Maxi    :         Max number of N str. treated simult.(not in use)y
!       MaxK        :        Largest number of inner resolution strings simult
!        NoPart         :        partitioning of strings in alpha-alpha + beta-beta loops
!        NOHSOO        :        No spin-other-orbit operator included
!        IDENMT        :        Calculate one-body density matrix
Integer MAXIT,IRESTR,INTIMP,MXP1,MXP2,MXQ,INCORE,MXCIV,       &
        ICISTR,NOCSF,IDIAG,NOINT,NOPAIR,IH1SO,IDELMS,         &
        ISOSYM(3),INIDEG,IFINSD,IPREREL,NOMSPA,IPREDZN,       &
        ISUPSP,IREDSTP,IDMPIN,MAXRIT,MAXI,MAXK,NOPART,        &
        NOHSOO,IDENMT

! Stuff from crun_mclr.fh
!
!        NINOB                :        Inactive orbitals         (0)
!        NACOB                :        Active orbiatls
!        NDEOB                :        Deleted                        (0)
!        NOCOB                :        Occupied
!        NTOOB                :        Total
!        NORB1                :        RAS1
!        NORB2                :        RAS2
!        NORB3                :        RAS3
!        NORB4                :        RAS4                        (0)
!        NOSPIR                :        ?                        (1)
!        IOSPIR(1,is)        :        ?                        iS
!        NINOBS                :        Inactive/symmetry        (0)
!        NR0OBS                :        Ras0/symmetry                (0)
!        NRSOBS                :        RASX/symmetry
!        NR4OBS                :        RAS4/symmetry                (0)
!        NACOBS                :        RAS/symmetry            (i)
!        NOCOBS                :        Occ/symmetry                (i)
!        NTOOBS                :        Orb/symmetry                (i)
!        NDEOBS                :        Deleted/symmetry        (0)
!        NRs4To                :         Not in use
!         IREOTS                 :         Reordering array type     => symmetry
!          IREOST                 :         Reordering array symmetry => type
!          ISMFTO          :         Symmetry array for type ordered orbitals
!        ITPFSO                :        Not in use
!          IBSO                   :         First orb of given symmetry (symmetry ordered)
!          NTSOB                  :         Number of active orb of give RAS type and SYM
!          IBTSOB                 :         Offset for active orb of given RAS type and SYM
!          ITSOB                  :        Orbitals of given RAS type and sym
!        NOBPTS                 :         Number of orbitals per subtype and symmetry
!          IOBPTS                 :         Off sets for orb of given subtype and symmetry
!        ITOBS                :        "Increaser"
!
!
!  ITFSO  : Type array for symmetry ordered orbitals( not activated )

!
!
Integer NINOB,NACOB,NDEOB,NOCOB,NTOOB,                            &
        NORB0,NORB1,NORB2,NORB3,NORB4,                      &
        NOSPIR(MXPIRR),IOSPIR(MXPOBS,MXPIRR),               &
        NINOBS(MXPOBS),NR0OBS(1,MXPOBS),NRSOBS(MXPOBS,3),   &
        NR4OBS(MXPOBS,MXPR4T),NACOBS(MXPOBS),NOCOBS(MXPOBS),&
        NTOOBS(MXPOBS),NDEOBS(MXPOBS),                      &
        NRS4TO(MXPR4T),                                     &
        IREOTS(MXPORB),IREOST(MXPORB),ISMFTO(MXPORB),       &
        ITPFSO(MXPORB),                                     &
        IBSO(MXPOBS),                                       &
        NTSOB(3,MXPOBS),IBTSOB(3,MXPOBS),ITSOB(MXPORB),     &
        NOBPTS(6+MXPR4T,MXPOBS),IOBPTS(6+MXPR4T,MXPOBS),    &
        ITOOBS(MXPOBS),ITPFTO(MXPORB),ISMFSO(MXPORB),       &
        NOBPT(6+MXPR4T)

! Stuff from spinfo_mclr.fh
!             MULTSP                        : Spin multiplicity
!        MS2P                        : 2*MS
!             MINOP                        : Minum open orbitals
!        MAXOP                        : Maximum open orbitals
!        NTYP                        : Maxop-MinOp+1
!        NDPCNT(MXPCTP)                : Number of det. / conf type
!        NCPCNT(MXPCTP)                : Number of csf / conf type
!        NCNATS(MXPCTP,MXPCSM)        :
!        NDTASM(MXPCSM)                :       Combinations
!        NCSASM(MXPCSM)                :        CSF
!        NCNASM(MXPCSM)                :i Configurations
!
Integer       MULTSP,MS2P,                                        &
              MINOP,MAXOP,NTYP,NDPCNT(MXPCTP),NCPCNT(MXPCTP),     &
              NCNATS(MXPCTP,MXPCSM),NDTASM(MXPCSM),NCSASM(MXPCSM),&
              NCNASM(MXPCSM)

End Module MCLR_Data
