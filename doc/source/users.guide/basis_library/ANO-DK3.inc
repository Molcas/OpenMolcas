Relativistic minimal valence ANO basis sets -- ANO-DK3
.......................................................

A highly accurate minimal valence relativistic basis set developed for the 103 elements from
H to Lr :cite:`ANO-DK3`.
The orbital exponents were optimized by minimizing the atomic self-consistent field (SCF) energy
with relativistic third-order Douglas--Kroll approximation (the atomic ground states were:
H(:math:`^2S`), He(:math:`^1S`), Li(:math:`^2S`), Be(:math:`^1S`), B(:math:`^2P`), C(:math:`^3P`), N(:math:`^4S`), O(:math:`^3P`), F(:math:`^2P`),
Ne(:math:`^1S`), Na(:math:`^2S`), Mg(:math:`^1S`), Al(:math:`^2P`), Si(:math:`^3P`), P(:math:`^4S`), S(:math:`^3P`), Cl(:math:`^2P`), Ar(:math:`^1S`),
K(:math:`^2S`), Ca(:math:`^1S`), Sc(:math:`^2D`), Ti(:math:`^3F`), V(:math:`^4F`), Cr(:math:`^7S`), Mn(:math:`^6S`), Fe(:math:`^5D`), Co(:math:`^4F`),
Ni(:math:`^3F`), Cu(:math:`^2S`), Zn(:math:`^1S`), Ga(:math:`^2P`), Ge(:math:`^3P`), As(:math:`^4S`), Se(:math:`^3P`), Br(:math:`^2P`), Kr(:math:`^1S`),
Rb(:math:`^2S`), Sr(:math:`^1S`), Y(:math:`^2D`), Zr(:math:`^3F`), Nb(:math:`^6D`), Mo(:math:`^7S`), Tc(:math:`^6S`), Ru(:math:`^5F`), Rh(:math:`^4F`),
Pd(:math:`^1S`), Ag(:math:`^2S`), Cd(:math:`^1S`), In(:math:`^2P`), Sn(:math:`^3P`), Sb(:math:`^4S`), Te(:math:`^3P`), I(:math:`^2P`), Xe(:math:`^1S`),
Cs(:math:`^2S`), Ba(:math:`^1S`), La(:math:`^2D`), Ce(:math:`^1G`), Pr(:math:`^4I`), Nd(:math:`^5I`), Pm(:math:`^6H`), Sm(:math:`^7F`), Eu(:math:`^8S`),
Gd(:math:`^9D`), Tb(:math:`^6H`), Dy(:math:`^5I`), Ho(:math:`^4I`), Er(:math:`^3H`), Tm(:math:`^2F`), Yb(:math:`^1S`), Lu(:math:`^2D`), Hf(:math:`^3F`),
Ta(:math:`^4F`), W(:math:`^5D`), Re(:math:`^6S`), Os(:math:`^5D`), Ir(:math:`^4F`), Pt(:math:`^3D`), Au(:math:`^2S`), Hg(:math:`^1S`), Tl(:math:`^2P`),
Pb(:math:`^3P`), Bi(:math:`^4S`), Po(:math:`^3P`), At(:math:`^2P`), Rn(:math:`^1S`), Fr(:math:`^2S`), Ra(:math:`^1S`), Ac(:math:`^2D`), Th(:math:`^3F`),
Pa(:math:`^4J`), U(:math:`^5L`), Np(:math:`^6L`), Pu(:math:`^7F`), Am(:math:`^8S`), Cm(:math:`^9D`), Bk(:math:`^8H`), Cf(:math:`^5I`), Es(:math:`^4I`),
Fm(:math:`^3H`), Md(:math:`^2F`), No(:math:`^1S`), and Lr(:math:`^2D`)). The minimal valence basis set is of a quality better
than a typical 6-31G basis set. For more accurate calculations the basis need additional
polarization and diffuse functions. The latter can be accomplished by on-the-fly decontraction.

The primitive sets are:

* 6s for H--He.
* 12s for Li--Be.
* 12s8p for B--Ne.
* 16s11p for Na--Ar.
* 20s15p for K--Ca.
* 20s15p9d for Sc--Kr.
* 23s19p12d for Rb--Xe.
* 27s23p15d for Cs--La.
* 27s23p15d10f for Ce--Rn.
* 35s26p18d13f for Fr--Lr.

The maximum number of ANO-DK3's given in the library is:

* 1s for H--He.
* 2s for Li--Be.
* 2s1p for B--Ne.
* 3s1p for Na--Mg.
* 3s1p for Al--Ar.
* 4s2p for K--Ca.
* 4s2p1d for Sc--Zn.
* 4s3p1d for Ga--Kr.
* 5s3p1d for Rb--Rh.
* 4s3p2d for Pd.
* 5s3p2d for Ag--Cd.
* 5s4p2d for In--Xe.
* 6s4p2d for Cs--Ba.
* 6s4p3d for La.
* 6s4p3d1f for Ce.
* 6s4p2d1f for Pr--Eu.
* 6s4p3d1f for Gd.
* 6s4p2d1f for Tb--Yb.
* 6s4p3d1f for Lu--Hg.
* 6s5p3d1f for Tl--Rn.
* 7s5p3d1f for Fr--Ra.
* 7s5p4d1f for Ac--Th.
* 7s5p4d2f for Pa--Np.
* 7s5p3d2f for Pu--Am.
* 7s5p4d2f for Cm--Bk.
* 7s5p3d2f for Cf--No.
* 7s5p4d2f for Lr.

Note, however, the number of contracted functions can always be extended by automatic decontraction of the most
diffuse exponents in each shell by specifying more contracted functions as compared to what the basis set formally carries.

Below is a list of the basis sets available in this class. ::

  H.ANO-DK3.Tsuchiya.6s.1s.
  He.ANO-DK3.Tsuchiya.6s.1s.
  Li.ANO-DK3.Tsuchiya.12s.2s.
  Be.ANO-DK3.Tsuchiya.12s.2s.
  B.ANO-DK3.Tsuchiya.12s8p.2s1p.
  C.ANO-DK3.Tsuchiya.12s8p.2s1p.
  N.ANO-DK3.Tsuchiya.12s8p.2s1p.
  O.ANO-DK3.Tsuchiya.12s8p.2s1p.
  F.ANO-DK3.Tsuchiya.12s8p.2s1p.
  Ne.ANO-DK3.Tsuchiya.12s8p.2s1p.
  Na.ANO-DK3.Tsuchiya.16s11p.3s1p.
  Mg.ANO-DK3.Tsuchiya.16s11p.3s1p.
  Al.ANO-DK3.Tsuchiya.16s11p.3s2p.
  Si.ANO-DK3.Tsuchiya.16s11p.3s2p.
  P.ANO-DK3.Tsuchiya.16s11p.3s2p.
  S.ANO-DK3.Tsuchiya.16s11p.3s2p.
  Cl.ANO-DK3.Tsuchiya.16s11p.3s2p.
  Ar.ANO-DK3.Tsuchiya.16s11p.3s2p.
  K.ANO-DK3.Tsuchiya.20s15p.4s2p.
  Ca.ANO-DK3.Tsuchiya.20s15p.4s2p.
  Sc.ANO-DK3.Tsuchiya.20s15p9d.4s2p1d.
  Ti.ANO-DK3.Tsuchiya.20s15p9d.4s2p1d.
  V.ANO-DK3.Tsuchiya.20s15p9d.4s2p1d.
  Cr.ANO-DK3.Tsuchiya.20s15p9d.4s2p1d.
  Mn.ANO-DK3.Tsuchiya.20s15p9d.4s2p1d.
  Fe.ANO-DK3.Tsuchiya.20s15p9d.4s2p1d.
  Co.ANO-DK3.Tsuchiya.20s15p9d.4s2p1d.
  Ni.ANO-DK3.Tsuchiya.20s15p9d.4s2p1d.
  Cu.ANO-DK3.Tsuchiya.20s15p9d.4s2p1d.
  Zn.ANO-DK3.Tsuchiya.20s15p9d.4s2p1d.
  Ga.ANO-DK3.Tsuchiya.20s15p9d.4s3p1d.
  Ge.ANO-DK3.Tsuchiya.20s15p9d.4s3p1d.
  As.ANO-DK3.Tsuchiya.20s15p9d.4s3p1d.
  Se.ANO-DK3.Tsuchiya.20s15p9d.4s3p1d.
  Br.ANO-DK3.Tsuchiya.20s15p9d.4s3p1d.
  Kr.ANO-DK3.Tsuchiya.20s15p9d.4s3p1d.
  Rb.ANO-DK3.Tsuchiya.23s19p12d.5s3p1d.
  Sr.ANO-DK3.Tsuchiya.23s19p12d.5s3p1d.
  Y.ANO-DK3.Tsuchiya.23s19p12d.5s3p2d.
  Zr.ANO-DK3.Tsuchiya.23s19p12d.5s3p2d.
  Nb.ANO-DK3.Tsuchiya.23s19p12d.5s3p2d.
  Mo.ANO-DK3.Tsuchiya.23s19p12d.5s3p2d.
  Tc.ANO-DK3.Tsuchiya.23s19p12d.5s3p2d.
  Ru.ANO-DK3.Tsuchiya.23s19p12d.5s3p2d.
  Rh.ANO-DK3.Tsuchiya.23s19p12d.5s3p2d.
  Pd.ANO-DK3.Tsuchiya.23s19p12d.4s3p2d.
  Ag.ANO-DK3.Tsuchiya.23s19p12d.5s3p2d.
  Cd.ANO-DK3.Tsuchiya.23s19p12d.5s3p2d.
  In.ANO-DK3.Tsuchiya.23s19p12d.5s4p2d.
  Sn.ANO-DK3.Tsuchiya.23s19p12d.5s4p2d.
  Sb.ANO-DK3.Tsuchiya.23s19p12d.5s4p2d.
  Te.ANO-DK3.Tsuchiya.23s19p12d.5s4p2d.
  I.ANO-DK3.Tsuchiya.23s19p12d.5s4p2d.
  Xe.ANO-DK3.Tsuchiya.23s19p12d.5s4p2d.
  Cs.ANO-DK3.Tsuchiya.27s23p15d.6s4p2d.
  Ba.ANO-DK3.Tsuchiya.27s23p15d.6s4p2d.
  La.ANO-DK3.Tsuchiya.27s23p15d.6s4p3d.
  Ce.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p3d1f.
  Pr.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p2d1f.
  Nd.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p2d1f.
  Pm.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p2d1f.
  Sm.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p2d1f.
  Eu.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p2d1f.
  Gd.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p3d1f.
  Tb.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p2d1f.
  Dy.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p2d1f.
  Ho.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p2d1f.
  Er.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p2d1f.
  Tm.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p2d1f.
  Yb.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p2d1f.
  Lu.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p3d1f.
  Hf.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p3d1f.
  Ta.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p3d1f.
  W.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p3d1f.
  Re.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p3d1f.
  Os.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p3d1f.
  Ir.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p3d1f.
  Pt.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p3d1f.
  Au.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p3d1f.
  Hg.ANO-DK3.Tsuchiya.27s23p15d10f.6s4p3d1f.
  Tl.ANO-DK3.Tsuchiya.27s23p15d10f.6s5p3d1f.
  Pb.ANO-DK3.Tsuchiya.27s23p15d10f.6s5p3d1f.
  Bi.ANO-DK3.Tsuchiya.27s23p15d10f.6s5p3d1f.
  Po.ANO-DK3.Tsuchiya.27s23p15d10f.6s5p3d1f.
  At.ANO-DK3.Tsuchiya.27s23p15d10f.6s5p3d1f.
  Rn.ANO-DK3.Tsuchiya.27s23p15d10f.6s5p3d1f.
  Fr.ANO-DK3.Tsuchiya.35s26p18d13f.7s5p3d1f.
  Ra.ANO-DK3.Tsuchiya.35s26p18d13f.7s5p3d1f.
  Ac.ANO-DK3.Tsuchiya.35s26p18d13f.7s5p4d1f.
  Th.ANO-DK3.Tsuchiya.35s26p18d13f.7s5p4d1f.
  Pa.ANO-DK3.Tsuchiya.35s26p18d13f.7s5p4d2f.
  U.ANO-DK3.Tsuchiya.35s26p18d13f.7s5p4d2f.
  Np.ANO-DK3.Tsuchiya.35s26p18d13f.7s5p4d2f.
  Pu.ANO-DK3.Tsuchiya.35s26p18d13f.7s5p3d2f.
  Am.ANO-DK3.Tsuchiya.35s26p18d13f.7s5p3d2f.
  Cm.ANO-DK3.Tsuchiya.35s26p18d13f.7s5p4d2f.
  Bk.ANO-DK3.Tsuchiya.35s26p18d13f.7s5p4d2f.
  Cf.ANO-DK3.Tsuchiya.35s26p18d13f.7s5p3d2f.
  Es.ANO-DK3.Tsuchiya.35s26p18d13f.7s5p3d2f.
  Fm.ANO-DK3.Tsuchiya.35s26p18d13f.7s5p3d2f.
  Md.ANO-DK3.Tsuchiya.35s26p18d13f.7s5p3d2f.
  No.ANO-DK3.Tsuchiya.35s26p18d13f.7s5p3d2f.
  Lr.ANO-DK3.Tsuchiya.35s26p18d13f.7s5p4d2f.
