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
Private MXPIRR,MXPOBS,MXPR4T,MXINKA,MXPORB,MXPXOT,MXPXST,MXPSHL, &
        MXPL,MXPXT,MXPICI,MXPSTT,MXPCSM,MXPCTP,MXCNSM,MXPWRD, &
        MXNMS,MTYP,MXPNGAS,MXPNSMST,MXPPTSPC

End Module MCLR_Data
