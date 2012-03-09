#include "Mosnbsim3SOI5T1T.h"

//Static members
const unsigned Mosnbsim3SOI5T1T::n_par = 206;

//Element information
ItemInfo Mosnbsim3SOI5T1T::einfo =
{
  "mosnbsim3soi5t1t",
  "BSIM3SOI Mosfet model version 3.2.4 with temperature",
  "Ramya Mohan",
  DEFAULT_ADDRESS"transistor>mosfet,electrothermal",
  "2003_00_00"
};

//Parameters
ParmInfo Mosnbsim3SOI5T1T::pinfo[] =
{
  {"l", "Length (m)", TR_DOUBLE, false},
  {"w", "Width (m)", TR_DOUBLE, false},
  {"dtoxcv", "---", TR_DOUBLE, false},
  {"llc", "Length reduction parameter for CV", TR_DOUBLE, false},
  {"lwc", "Length reduction parameter for CV", TR_DOUBLE, false},
  {"lwlc", "Length reduction parameter for CV", TR_DOUBLE, false},
  {"wlc", "Width reduction parameter for CV", TR_DOUBLE, false},
  {"wwc", "Width reduction parameter for CV", TR_DOUBLE, false},
  {"wwlc", "Width reduction parameter for CV", TR_DOUBLE, false},
  {"tsi", "Silicon film thickness (m)", TR_DOUBLE, false},
  {"tox", "Gate oxide thickness (m)", TR_DOUBLE, false},
  {"toxref", "--Target Oxide Thickness", TR_DOUBLE, false},
  {"tbox", "Buried Oxide thickness (m)", TR_DOUBLE, false},
  {"tnom", "Parameter measurement temperature (K)", TR_DOUBLE, false},
  {"rbody", "Intrinsic body contact sheet resistance (ohm/square)", TR_DOUBLE, false},
  {"rbsh", "Extrinsic body contact sheet resistance (ohm/square)", TR_DOUBLE, false},
  {"rsh", "S-D sheet resistancde (ohm/square)", TR_DOUBLE, false},
  {"rhalo", "Body halo sheet resistance (ohm/m)", TR_DOUBLE, false},
  {"wint", "Width offset fitting parameter from I-V without bias (m)", TR_DOUBLE, false},
  {"lint", "Length offset fitting parameter from I-V without bias (m)", TR_DOUBLE, false},
  {"wth0", "-----", TR_DOUBLE, false},
  {"ll", "Length reduction parameter", TR_DOUBLE, false},
  {"wl", "Width reduction parameter", TR_DOUBLE, false},
  {"lln", "Length reduction parameter", TR_DOUBLE, false},
  {"wln", "Width reduction parameter", TR_DOUBLE, false},
  {"lw", "Length reduction parameter", TR_DOUBLE, false},
  {"ww", "Width reduction parameter", TR_DOUBLE, false},
  {"lwn", "Length reduction parameter", TR_DOUBLE, false},
  {"wwn", "Width reduction parameter", TR_DOUBLE, false},
  {"lwl", "Length reduction parameter", TR_DOUBLE, false},
  {"wwl", "Width reduction parameter", TR_DOUBLE, false},
  {"ln", "Electron/hole diffusion length (m)", TR_DOUBLE, false},
  {"xpart", "-----Channel charge partititoning", TR_DOUBLE, false},
  {"xj", "S/DJunction depth (m)", TR_DOUBLE, false},
  {"k1b", "k1b", TR_DOUBLE, false},
  {"k2b", "k2b", TR_DOUBLE, false},
  {"dk2b", "dk2b", TR_DOUBLE, false},
  {"vbsa", "vbsa", TR_DOUBLE, false},
  {"aigc", "---", TR_DOUBLE, false},
  {"bigc", "----", TR_DOUBLE, false},
  {"cigc", "-----", TR_DOUBLE, false},
  {"aigsd", "-----", TR_DOUBLE, false},
  {"bigsd", "-----", TR_DOUBLE, false},
  {"cigsd", "-----", TR_DOUBLE, false},
  {"nigc", "-------", TR_DOUBLE, false},
  {"poxedge", "-----", TR_DOUBLE, false},
  {"pigcd", "-------", TR_DOUBLE, false},
  {"vth0", "Threshold voltage @Vbs=0 for long and wide device", TR_DOUBLE, false},
  {"k1", "First order body effect coefficient (V^0.5)", TR_DOUBLE, false},
  {"k1w1", "First body effect width dependent parameter (m)", TR_DOUBLE, false},
  {"k1w2", "Second body effect width dependent parameter (m)", TR_DOUBLE, false},
  {"k2", "Second order body effect coefficient", TR_DOUBLE, false},
  {"k3", "Narrow width effect coefficient", TR_DOUBLE, false},
  {"k3b", "Body effect coefficient of k3 (1/V)", TR_DOUBLE, false},
  {"kb1", "Backgate body charge coefficient", TR_DOUBLE, false},
  {"w0", "Narrow width effect parameter (m)", TR_DOUBLE, false},
  {"nlx", "Lateral non-uniform doping parameter (m)", TR_DOUBLE, false},
  {"nch", "Channel doping concentration (1/cm^3)", TR_DOUBLE, false},
  {"nsub", "Substrate doping concentration (1/cm^3)", TR_DOUBLE, false},
  {"ngate", "Poly-gate doping concentration (1/cm^3)", TR_DOUBLE, false},
  {"dvt0", "First coefficient of short-channel effect on Vth", TR_DOUBLE, false},
  {"dvt1", "Second coefficient of short-channel effect on Vth", TR_DOUBLE, false},
  {"dvt2", "Body-bias coefficient of short-channel effect on Vth (1/V)", TR_DOUBLE, false},
  {"dvt0w", "First coefficient of narrow width effect on Vth for small channel length", TR_DOUBLE, false},
  {"dvt1w", "Second coefficient of narrow width effect on Vth for small channel length", TR_DOUBLE, false},
  {"dvt2w", "Body-bias coefficient of narrow width effect on Vth for small channel length (1/V)", TR_DOUBLE, false},
  {"eta0", "DIBL coefficeint subthreshold region", TR_DOUBLE, false},
  {"etab", "Body bias coefficeint for the subthreshold DIBL effect (1/V)", TR_DOUBLE, false},
  {"dsub", "DIBL coefficient in the subthreshold region", TR_DOUBLE, false},
  {"voff", "Offset voltage in the threshold region for large W and L (V)", TR_DOUBLE, false},
  {"nfactor", "Subthreshold swing factor", TR_DOUBLE, false},
  {"cdsc", "Drain/Source to channel coupling capacitance (F/m^2)", TR_DOUBLE, false},
  {"cdscb", "Body-bias sensitivity of cdsc (F/m^2)", TR_DOUBLE, false},
  {"cdscd", "Drain-bias sensitivity of cdsc (F/m^2)", TR_DOUBLE, false},
  {"cit", "Interface trap capacitance (F/m^2)", TR_DOUBLE, false},
  {"u0", "Mobility at Temp=Tnom (cm^2/V-sec)", TR_DOUBLE, false},
  {"ua", "First-oreder mobility degradation coefficient (m/V)", TR_DOUBLE, false},
  {"ub", "Second-order mobility degradation coefficient (m/V)^2", TR_DOUBLE, false},
  {"uc", "Body-effect of mobility degradation coefficient (1/V)", TR_DOUBLE, false},
  {"prwg", "Gate-bias effect coefficient of rdsw", TR_DOUBLE, false},
  {"prwb", "Body effect coefficient of rdsw (1/V)", TR_DOUBLE, false},
  {"wr", "Width offset from Weff for Rds calculation", TR_DOUBLE, false},
  {"rdsw", "Parasitic resistance per unit width (ohm-um)", TR_DOUBLE, false},
  {"a0", "Bulk charge effect coefficient for channel length", TR_DOUBLE, false},
  {"ags", "Gate bias coefficient of Abulk (1/V)", TR_DOUBLE, false},
  {"a1", "First non-saturation effect parameter (1/V)", TR_DOUBLE, false},
  {"a2", "Second non-saturation effect parameter", TR_DOUBLE, false},
  {"b0", "Bulk charge effect coefficient for channel width (m)", TR_DOUBLE, false},
  {"b1", "Bulk charge effect width offset (m)", TR_DOUBLE, false},
  {"vsat", "Saturation velocity at Temp=Tnom (m/sec)", TR_DOUBLE, false},
  {"keta", "Body-bias coefficient of bulk charge effect (1/V)", TR_DOUBLE, false},
  {"ketas", "Surface potential adjustment for bulk charge effect (V)", TR_DOUBLE, false},
  {"dwg", "Coefficient of Weff's gate dependence (m/V)", TR_DOUBLE, false},
  {"dwb", "Coefficient of Weff's substrate body bias dependence (m/V)^0.5", TR_DOUBLE, false},
  {"dwbc", "Width offset for body contact isolation edge (m)", TR_DOUBLE, false},
  {"pclm", "Channel length modulation parameter", TR_DOUBLE, false},
  {"pdibl1", "First output resistance DIBL effect correction parameter", TR_DOUBLE, false},
  {"pdibl2", "Second output resistance DIBL effect correction parameter", TR_DOUBLE, false},
  {"pdiblb", "Body-effect on drain induced barrier lowering", TR_DOUBLE, false},
  {"drout", "----listed in the pgm", TR_DOUBLE, false},
  {"pvag", "Gate dependence of Early voltage", TR_DOUBLE, false},
  {"delta", "Effective Vds parameter", TR_DOUBLE, false},
  {"alpha0", "The first parameter of impact ionization current (m/V)", TR_DOUBLE, false},
  {"beta0", "First Vds dependent parameter of impact ionization current (1/V)", TR_DOUBLE, false},
  {"beta1", "Second Vds dependent parameter of impact ionization current", TR_DOUBLE, false},
  {"beta2", "Third Vds dependent parameter of impact ionization current (V)", TR_DOUBLE, false},
  {"fbjtii", "---", TR_DOUBLE, false},
  {"vdsatii0", "Nominal drain saturation voltage at threshold for impact ionization current (V)", TR_DOUBLE, false},
  {"tii", "Temperature dependent parameter for impact ionization current", TR_DOUBLE, false},
  {"lii", "Channel length dependent parameter at threshold for impact ionization current", TR_DOUBLE, false},
  {"esatii", "Saturation channel electric field for impact ionization current (V/m)", TR_DOUBLE, false},
  {"sii0", "First Vgs dependent parameter for impact ionization current (V^-1)", TR_DOUBLE, false},
  {"sii1", "Second Vgs dependent parameter for impact ionization current (V^-1)", TR_DOUBLE, false},
  {"sii2", "Third Vgs dependent parameter for impact ionization current", TR_DOUBLE, false},
  {"siid", "Vds dependent parameter of drain voltage for impact ionization current (V^-1)", TR_DOUBLE, false},
  {"agidl", "GIDL constant (ohm^-1)", TR_DOUBLE, false},
  {"bgidl", "GIDL Exponential coefficient (V/m)", TR_DOUBLE, false},
  {"ngidl", "GIDL Vds enhancement coefficient (V)", TR_DOUBLE, false},
  {"ebg", "-----", TR_DOUBLE, false},
  {"vgb1", "-----", TR_DOUBLE, false},
  {"vgb2", "------", TR_DOUBLE, false},
  {"voxh", "----", TR_DOUBLE, false},
  {"deltavox", "-------", TR_DOUBLE, false},
  {"ntox", "--------", TR_DOUBLE, false},
  {"ntun", "Reverse tunneling non-ideality factor", TR_DOUBLE, false},
  {"ndiode", " Diode non-ideality factor", TR_DOUBLE, false},
  {"nrecf0", "Recombination non-ideality factor at forward bias", TR_DOUBLE, false},
  {"nrecr0", "Recombination non-ideality factor at reversed bias", TR_DOUBLE, false},
  {"isbjt", "BJT Injection saturation current (A/m^2)", TR_DOUBLE, false},
  {"isdif", "Body to source/drain injection saturation current (A/m^2)", TR_DOUBLE, false},
  {"isrec", "Recombination in depletion saturation current (A/m^2)", TR_DOUBLE, false},
  {"istun", "Reverse tunneling saturation current (A/m^2)", TR_DOUBLE, false},
  {"vrec0", "Voltage dependent parameter for recombination current (V)", TR_DOUBLE, false},
  {"vtun0", "Voltage dependent parameter for tunneling current (V)", TR_DOUBLE, false},
  {"nbjt", "Power coeffcient of channel length dependency for bipolar current", TR_DOUBLE, false},
  {"lbjt0", "Reference channel length for bipolar current (m)", TR_DOUBLE, false},
  {"vabjt", "Early Voltage for bipolar current (V)", TR_DOUBLE, false},
  {"aely", "Channel length dependency of early voltage for bipolar current (V/m)", TR_DOUBLE, false},
  {"ahli", "High ;evel injection parameter for bipolar current", TR_DOUBLE, false},
  {"vevb", "----", TR_DOUBLE, false},
  {"vecb", "-----", TR_DOUBLE, false},
  {"cjswg", "Source/Drain (gate side) sidewall junction capacitance per unit width (normalized to 100nm tsi) (F/m^2)",
  TR_DOUBLE, false},
  {"mjswg", "Source/Drain (gate side) sidewall junction capacitance grading coefficient (V)", TR_DOUBLE, false},
  {"pbswg", "Source/Drain (gate side) sidewall junction capacitance built in potential (V)", TR_DOUBLE, false},
  {"tt", "Diffusion capacitance transit time coefficient (sec)", TR_DOUBLE, false},
  {"ldif0", "ldif0", TR_DOUBLE, false},
  {"cgeo", "Gate substrate overlap capacitance per unit channel length (F/m)", TR_DOUBLE, false},
  {"cgso", "------", TR_DOUBLE, false},
  {"cgdo", "-------", TR_DOUBLE, false},
  {"dlc", "Length offset fitting parameter for gate charge (m)", TR_DOUBLE, false},
  {"dwc", "Width offset fitting parameter from C-V (m)", TR_DOUBLE, false},
  {"dlcb", "Length offset fitting parameter for body charge (m)", TR_DOUBLE, false},
  {"dlbg", "Length offset fitting parameter for backgate charge (m)", TR_DOUBLE, false},
  {"fbody", "Scaling factor for body charge", TR_DOUBLE, false},
  {"clc", "Constant term for the short channel model (m)", TR_DOUBLE, false},
  {"cle", "Exponential term for the short channel model", TR_DOUBLE, false},
  {"cf", "---in the pgm", TR_DOUBLE, false},
  {"csdmin", "-----", TR_DOUBLE, false},
  {"asd", "--------", TR_DOUBLE, false},
  {"csdesw", "S/D sidewall fringing capacitance per unit length (F/m)", TR_DOUBLE, false},
  {"vsdfb", "-------", TR_DOUBLE, false},
  {"vsdth", "----------", TR_DOUBLE, false},
  {"delvt", "Threshold voltage adjust for C-V (V)", TR_DOUBLE, false},
  {"acde", "---in the pghm", TR_DOUBLE, false},
  {"moin", "Coefficient for the gate-bias dependent surface potential V^0.5", TR_DOUBLE, false},
  {"ckappa", "Coefficient for lightly doped region overlap capacitance fringing field capacitance (F/m)", TR_DOUBLE, false},
  {"cgdl", "Light doped drain-gate region overlap capacitance (F/m)", TR_DOUBLE, false},
  {"cgsl", "Light doped source-gate region overlap capacitance (F/m)", TR_DOUBLE, false},
  {"ndif", "ndif", TR_DOUBLE, false},
  {"rth0", "------", TR_DOUBLE, false},
  {"cth0", "----------", TR_DOUBLE, false},
  {"tpbswg", "----------", TR_DOUBLE, false},
  {"tcjswg", "----------", TR_DOUBLE, false},
  {"kt1", "Temperature coefficient of Vth (V)", TR_DOUBLE, false},
  {"kt1l", "Channel length dependence of the temperature coefficient of Vth (V*m)", TR_DOUBLE, false},
  {"kt2", "Body-bias coefficient of the Vth temperature effect", TR_DOUBLE, false},
  {"ute", "Temperature coefficient of mobility", TR_DOUBLE, false},
  {"ua1", "Temperature coefficient for ua (m/V)", TR_DOUBLE, false},
  {"ub1", "Temperature coefficient for ub ((m/V)^2)", TR_DOUBLE, false},
  {"uc1", "Temperature coefficient for uc (1/V)", TR_DOUBLE, false},
  {"prt", "Temperature coefficient of rdsw (ohm-um)", TR_DOUBLE, false},
  {"at", "Temperature coefficient of vsat (m/sec)", TR_DOUBLE, false},
  {"ntrecf", "Temperature coefficient for Nrecf", TR_DOUBLE, false},
  {"ntrecr", "Temperature coefficient for Nrecr", TR_DOUBLE, false},
  {"xbjt", "xbjt", TR_DOUBLE, false},
  {"xdif", "xdif", TR_DOUBLE, false},
  {"xrec", "xrec", TR_DOUBLE, false},
  {"xtun", "xtun", TR_DOUBLE, false},
  {"dlcig", "---", TR_DOUBLE, false},
  {"nbc", "----", TR_DOUBLE, false},
  {"nseg", "----", TR_DOUBLE, false},
  {"pdbcp", "-----", TR_DOUBLE, false},
  {"psbcp", "-----", TR_DOUBLE, false},
  {"toxqm", "Effective oxide thickness considering quantum effect", TR_DOUBLE, false},
  {"type", "-----", TR_DOUBLE, false},
  {"toxm", "=----", TR_DOUBLE, false},
  {"xt1", "Doping depth", TR_DOUBLE, false},
  {"dvbd0", "dvbd0", TR_DOUBLE, false},
  {"dvbd1", "dvbd1", TR_DOUBLE, false},
  {"temp", "Circuit temperature", TR_DOUBLE, false},
  {"npeak","--",TR_DOUBLE, false},
  {"capMod","Capacitance model ",TR_DOUBLE, false},
  {"vbm", "Maximum body voltage", TR_DOUBLE, false},
  {"nofffd", "nofffd", TR_DOUBLE, false},
  {"vofffd", "vofffd", TR_DOUBLE, false},
  {"moinFD", "moinfd", TR_DOUBLE, false},
  {"shmod", "self heat - shmod", TR_DOUBLE, false}
};

Mosnbsim3SOI5T1T::Mosnbsim3SOI5T1T(const string& iname) : ADInterface(&einfo, pinfo, n_par, iname)
{
  //Set default parameter values
  paramvalue[0] = &(l = 0.25e-6);
  paramvalue[1] = &(w = 10.0e-6);
  paramvalue[2] = &(dtoxcv = 0);
  paramvalue[3] = &(llc = 0);
  paramvalue[4] = &(lwc = 0);
  paramvalue[5] = &(lwlc = 0);
  paramvalue[6] = &(wlc = 0);
  paramvalue[7] = &(wwc = 0);
  paramvalue[8] = &(wwlc = 0);
  paramvalue[9] = &(tsi = 1e-07);
  paramvalue[10] = &(tox = 5e-09);
  paramvalue[11] = &(toxref = 5e-09);
  paramvalue[12] = &(tbox = 5e-07);
  paramvalue[13] = &(tnom = 293.15);
  paramvalue[14] = &(rbody = 1);
  paramvalue[15] = &(rbsh = 0);
  paramvalue[16] = &(rsh = 0);
  paramvalue[17] = &(rhalo = 1e+015);
  paramvalue[18] = &(wint = 0);
  paramvalue[19] = &(lint = 0);
  paramvalue[20] = &(wth0 = 0);
  paramvalue[21] = &(ll = 0);
  paramvalue[22] = &(wl = 0);
  paramvalue[23] = &(lln = 1);
  paramvalue[24] = &(wln = 1);
  paramvalue[25] = &(lw = 0);
  paramvalue[26] = &(ww = 0);
  paramvalue[27] = &(lwn = 1);
  paramvalue[28] = &(wwn = 1);
  paramvalue[29] = &(lwl = 0);
  paramvalue[30] = &(wwl = 0);
  paramvalue[31] = &(ln = 2e-06);
  paramvalue[32] = &(xpart = 1);
  paramvalue[33] = &(xj = 1e-07);
  paramvalue[34] = &(k1b = 0);
  paramvalue[35] = &(k2b = 0);
  paramvalue[36] = &(dk2b = 0);
  paramvalue[37] = &(vbsa = 0.0);
  paramvalue[38] = &(aigc = 1);
  paramvalue[39] = &(bigc = 1);
  paramvalue[40] = &(cigc = 1);
  paramvalue[41] = &(aigsd = 1);
  paramvalue[42] = &(bigsd = 1);
  paramvalue[43] = &(cigsd = 1);
  paramvalue[44] = &(nigc = 1);
  paramvalue[45] = &(pioxedge = 1);
  paramvalue[46] = &(pigcd = 1);
  paramvalue[47] = &(vth0 = 0.53);
  paramvalue[48] = &(k1 = 0.56);
  paramvalue[49] = &(k1w1 = 0);
  paramvalue[50] = &(k1w2 = 0);
  paramvalue[51] = &(k2 = 0);
  paramvalue[52] = &(k3 = -2);
  paramvalue[53] = &(k3b = 0);
  paramvalue[54] = &(kb1 = 1);
  paramvalue[55] = &(w0 = 0);
  paramvalue[56] = &(nlx = 0);
  paramvalue[57] = &(nch = 8e+17);
  paramvalue[58] = &(nsub = 5e+15);
  paramvalue[59] = &(ngate = 2e+20);
  paramvalue[60] = &(dvt0 = 1);
  paramvalue[61] = &(dvt1 = 0.15);
  paramvalue[62] = &(dvt2 = 0);
  paramvalue[63] = &(dvt0w = 0);
  paramvalue[64] = &(dvt1w = 2000000);
  paramvalue[65] = &(dvt2w = -0.032);
  paramvalue[66] = &(eta0 = 0.5);
  paramvalue[67] = &(etab = 0);
  paramvalue[68] = &(dsub = 0.35);
  paramvalue[69] = &(voff = -0.15);
  paramvalue[70] = &(nfactor = 0.4);
  paramvalue[71] = &(cdsc = 0.005);
  paramvalue[72] = &(cdscb = -0.01);
  paramvalue[73] = &(cdscd = 0);
  paramvalue[74] = &(cit = 0);
  paramvalue[75] = &(u0 = 0.05);
  paramvalue[76] = &(ua = 0);
  paramvalue[77] = &(ub = 1.2e-18);
  paramvalue[78] = &(uc = 0);
  paramvalue[79] = &(prwg = 0);
  paramvalue[80] = &(prwb = 0);
  paramvalue[81] = &(wr = 1);
  paramvalue[82] = &(rdsw = 100);
  paramvalue[83] = &(a0 = 0);
  paramvalue[84] = &(ags = 0);
  paramvalue[85] = &(a1 = 0);
  paramvalue[86] = &(a2 = 0.99);
  paramvalue[87] = &(b0 = 0);
  paramvalue[88] = &(b1 = 0);
  paramvalue[89] = &(vsat = 80000);
  paramvalue[90] = &(keta = 0);
  paramvalue[91] = &(ketas = 0);
  paramvalue[92] = &(dwg = 0);
  paramvalue[93] = &(dwb = 0);
  paramvalue[94] = &(dwbc = 0);
  paramvalue[95] = &(pclm = 1);
  paramvalue[96] = &(pdibl1 = 0.1);
  paramvalue[97] = &(pdibl2 = 0);
  paramvalue[98] = &(pdiblb = 0);
  paramvalue[99] = &(drout = 0.4);
  paramvalue[100] = &(pvag = 0);
  paramvalue[101] = &(delta = 0.001);
  paramvalue[102] = &(alpha0 = 8e-09);
  paramvalue[103] = &(beta0 = 0);
  paramvalue[104] = &(beta1 = 0);
  paramvalue[105] = &(beta2 = 0.05);
  paramvalue[106] = &(fbjtii = 0);
  paramvalue[107] = &(vdsatii0 = 0.8);
  paramvalue[108] = &(tii = -0.2);
  paramvalue[109] = &(lii = 5e-08);
  paramvalue[110] = &(esatii = 1e+08);
  paramvalue[111] = &(sii0 = 0.5);
  paramvalue[112] = &(sii1 = 0);
  paramvalue[113] = &(sii2 = 0);
  paramvalue[114] = &(siid = 0);
  paramvalue[115] = &(agidl = 2e-09);
  paramvalue[116] = &(bgidl = 2e+09);
  paramvalue[117] = &(ngidl = 0.5);
  paramvalue[118] = &(ebg = 1.2);
  paramvalue[119] = &(vgb1 = 300);
  paramvalue[120] = &(vgb2 = 17);
  paramvalue[121] = &(voxh = 5);
  paramvalue[122] = &(deltavox = 0.005);
  paramvalue[123] = &(ntox = 1);
  paramvalue[124] = &(ntun = 3.6);
  paramvalue[125] = &(ndiode = 1);
  paramvalue[126] = &(nrecf0 = 1.8);
  paramvalue[127] = &(nrecr0 = 1);
  paramvalue[128] = &(isbjt = 3e-07);
  paramvalue[129] = &(isdif = 3e-08);
  paramvalue[130] = &(isrec = 0.0005);
  paramvalue[131] = &(istun = 1e-08);
  paramvalue[132] = &(vrec0 = 0.05);
  paramvalue[133] = &(vtun0 = 5);
  paramvalue[134] = &(nbjt = 1);
  paramvalue[135] = &(lbjt0 = 2e-07);
  paramvalue[136] = &(vabjt = 10);
  paramvalue[137] = &(aely = 0);
  paramvalue[138] = &(ahli = 1e-15);
  paramvalue[139] = &(vevb = 0.075);
  paramvalue[140] = &(vecb = 0.026);
  paramvalue[141] = &(cjswg = 5e-10);
  paramvalue[142] = &(mjswg = 0.5);
  paramvalue[143] = &(pbswg = 0.8);
  paramvalue[144] = &(tt = 5e-10);
  paramvalue[145] = &(ldif0 = 0.001);
  paramvalue[146] = &(cgeo = 0);
  paramvalue[147] = &(cgso = 6.5e-10);
  paramvalue[148] = &(cgdo = 6e-10);
  paramvalue[149] = &(dlc = 0);
  paramvalue[150] = &(dwc = 0);
  paramvalue[151] = &(dlcb = 0);
  paramvalue[152] = &(dlbg = 0);
  paramvalue[153] = &(fbody = 1);
  paramvalue[154] = &(clc = 1e-07);
  paramvalue[155] = &(cle = 0.6);
  paramvalue[156] = &(cf = 0);
  paramvalue[157] = &(csdmin = 2.5e-05);
  paramvalue[158] = &(asd = 0.5);
  paramvalue[159] = &(csdesw = 0);
  paramvalue[160] = &(vsdfb = -0.8);
  paramvalue[161] = &(vsdth = -0.3);
  paramvalue[162] = &(delvt = 0);
  paramvalue[163] = &(acde = 0);
  paramvalue[164] = &(moin = 15);
  paramvalue[165] = &(ckappa = 0.6);
  paramvalue[166] = &(cgdl = 0);
  paramvalue[167] = &(cgsl = 0);
  paramvalue[168] = &(ndif = -1);
  paramvalue[169] = &(rth0 = 0.09);
  paramvalue[170] = &(cth0 = 1e-05);
  paramvalue[171] = &(tpbswg = 0);
  paramvalue[172] = &(tcjswg = 0.0005);
  paramvalue[173] = &(kt1 = -0.2);
  paramvalue[174] = &(kt1l = 8e-09);
  paramvalue[175] = &(kt2 = -0.06);
  paramvalue[176] = &(ute = -1.5);
  paramvalue[177] = &(ua1 = 3e-10);
  paramvalue[178] = &(ub1 = -3e-18);
  paramvalue[179] = &(uc1 = -6e-11);
  paramvalue[180] = &(prt = 10);
  paramvalue[181] = &(at = 65000);
  paramvalue[182] = &(ntrecf = 0.1);
  paramvalue[183] = &(ntrecr = -1);
  paramvalue[184] = &(xbjt = 1e-20);
  paramvalue[185] = &(xdif = 1.6);
  paramvalue[186] = &(xrec = 0.8);
  paramvalue[187] = &(xtun = 6);
  paramvalue[188] = &(dlcig = lint);
  paramvalue[189] = &(nbc = 0);
  paramvalue[190] = &(nseg = 1);
  paramvalue[191] = &(pdbcp = 0);
  paramvalue[192] = &(psbcp = 0);
  paramvalue[193] = &(toxqm = tox);
  paramvalue[194] = &(type = 1);
  paramvalue[195] = &(toxm = tox);
  paramvalue[196] = &(xt1 = 1.55e-7);
  paramvalue[197] = &(dvbd0 = 0.0);
  paramvalue[198] = &(dvbd1 = 0.0);
  paramvalue[199] = &(temp = 300.15);
  paramvalue[200] = &(npeak = 5.8e+17);
  paramvalue[201] = &(capMod = 2.0);
  paramvalue[202] = &(vbm = 0.0);
  paramvalue[203] = &(nofffd = 1.0);
  paramvalue[204] = &(vofffd = 0.0);
  paramvalue[205] = &(moinFD = 1e6);
  paramvalue[206] = &(shmod = 1.0);

  //Set flags
  setFlags(NONLINEAR | MULTI_REF | TR_TIME_DOMAIN);
}

void Mosnbsim3SOI5T1T::init() throw(string&)
{
  //Set number of terminals
  setNumTerms(6);

  //Set number of state variables
  setNumberOfStates(4);

  DenseIntVector var(4);
  var[0] = 0; //Vds
  var[1] = 1; //Vgs
  var[2] = 2; //Ves
  var[3] = 3; //temp
  initializeAD(var, var);
}

void Mosnbsim3SOI5T1T::getLocalRefIdx(UnsignedVector& local_ref_vec,
			      TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0));
  term_list.push_back(getTerminal(1));
  term_list.push_back(getTerminal(2));
  term_list.push_back(getTerminal(3));// Local reference terminal
  term_list.push_back(getTerminal(4));
  term_list.push_back(getTerminal(5));// Local reference terminal

  local_ref_vec.push_back(3); // Local reference index
  local_ref_vec.push_back(5); // Local reference index
}


void Mosnbsim3SOI5T1T::eval(AD * x, AD * effort, AD * flow)
{
  AD Vgs_eff, Vbsmos, Vbp, Vbsh, Vbseff, Phis, sqrtPhis, Xdep, lt1, ltw;
  AD DeltVthtemp, DIBL_Sft;
  AD sqrtPhisExt, Vth, n, Vgst, VgstNVt, ExpArg, Vgsteff, Vgst2Vtm, Weff, Rds;
  AD Abulk, Abulk0, Denomi, ueff, WVCox, WVCoxRds, AbovVgst2Vtm;
  AD Vdsat, Vdseff, diffVds, Vasat, VACLM, VADIBL, Va, CoxWovL, beta, fgche1, fgche2;
  AD gche, Vfb, V3, Vfbeff, Qac0, Qsub0, AbulkCV, VdsatCV, V4, VdseffCV, qbulk;
  AD qinv, qsrc, qgate, qbody, qsub, qdrn;
  AD T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, TMP, TMP3, TMP4, TMP11, TMP12, TMP13;
  AD t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, tmp, tmp1, tmp2, tmp3, tmp4, vbsc, V0, T00, TT0;
  AD Esat, EsatL;
  AD cjsbs, dcjsbs_dT, cjdbs,  dcjdbs_dT,  DioMax,arg, PhiBSWG,MJSWG,dT3_dVb,  qjs, qjd;
  AD SDphi, qse, qde,qgd, qgs, qge;
  AD vfbb, vearly;
  AD sdt1, st1, st2, st3, st4, dt2, dt3, dt4;
  AD Ahli = ahli;
  AD Ua = ua;
  AD Ub = ub;
  AD Uc = uc;
  AD Vsdfb = vsdfb;
  AD K2 = k2;
  AD K1 = k1;
  AD Vth0 = vth0;
  AD Nch = nch;


  int soiMod = 1 ;
  int rgateMod = 0;
  int selfheat = ((shmod ==1) && (rth0 !=0.0));

#define EPSOX 3.453133e-11
#define EPSSI 1.03594e-10
#define PI 3.141592654
#define Charge_q 1.60219e-19
#define Kb 1.3806226e-23
#define KboQ 8.617087e-5  /*  Kb / q   */
#define Eg300 1.115   /*  energy gap at 300K  */
#define DELTA_1 0.02
#define DELTA_2 0.02
#define DELTA_3 0.02
/* Original is 0.02, for matching IBM model, change to 0.08 */
#define DELTA_3_SOI 0.08
#define DELTA_4 0.02
#define DELT_Vbseff  0.005
#define DELTA_VFB  0.02
#define OFF_Vbsitf 0.02   /* v3.1*/
#define CONST_2OV3 0.6666666666
#define MAX_EXPL 2.688117142e+43
#define MIN_EXPL 3.720075976e-44
#define EXPL_THRESHOLD 100.0
#define DEXP(A,B,C) {                                       \
        if (A > EXPL_THRESHOLD) {                            \
            B = MAX_EXPL*(1.0+(A)-EXPL_THRESHOLD);            \
            C = MAX_EXPL;                                 	   \
        } else if (A < -EXPL_THRESHOLD)  {                		\
            B = MIN_EXPL;                                    \
            C = 0;                                          \
        } else   {                                          \
            B = exp(A);                                     \
            C = B;                                          \
        }                                                   \
    }


#define DEXP1(D,E) {                                                        \
        if (D > EXPL_THRESHOLD) {                                              \
            E = MAX_EXPL*(1.0+(D)-EXPL_THRESHOLD);                              \
        } else if (D < -EXPL_THRESHOLD)  {                                                \
            E = MIN_EXPL;                                                      \
        } else   {                                                            \
            E = exp(D);                                                       \
        }                                                                     \
    }


  double factor1 = sqrt(EPSSI / EPSOX * tox);
  double cox = 3.453133e-11 / tox;
  double ldrn = l;
  double wdrn = w;

  t0 = pow(ldrn, lln);
  t1 = pow(wdrn, lwn);
  tmp1 = ll / t0 + lw / t1 + lwl / (t0 * t1);
  AD dl = lint + tmp1;
  tmp1 = llc / t0 + lwc / t1 + lwlc / (t0 * t1);
  AD Dlc = dlc + tmp1;
  AD DLCIG = dlcig + tmp1;
  t2 = pow(ldrn, wln);
  t3 = pow(wdrn, wwn);
  tmp2 = wl / t2 + ww / t3 + wwl / (t2 * t3);
  AD dw = wint + tmp2;
  tmp2 = wlc / t2 + wwc / t3 + wwlc / (t2 * t3);
  AD Dwc = dwc + tmp2;
  double nbc = 0.0;
  AD leff = l - 2.0 * dl;
  AD weff = w - nbc * dwbc - (2.0 - nbc) * dw;
  double pdbcp = 0;
  double psbcp = 0;

  AD wdiod = weff / nseg + pdbcp;
  AD wdios = weff / nseg + psbcp;
  AD leffCV = l - 2.0 * Dlc;
  AD weffCV = w - nbc * dwbc - (2.0 - nbc) * Dwc;
  AD wdiodCV = weffCV / nseg + pdbcp;
  AD wdiosCV = weffCV / nseg + psbcp;
  AD leffCVb = l - 2.0 * Dlc - dlcb;
  AD leffCVbg = leffCVb + 2 * dlbg;
  AD abulkCVfactor = 1.0 + pow((clc / leff), cle);
  AD uatemp = ua;
  AD ubtemp = ub;
  AD uctemp = uc;
  AD rds0denom = pow(weff * 1e6, wr);
  AD rth = rth0 / (weff + wth0) * nseg;
  AD cth = cth0 * (weff + wth0) / nseg;
  AD frbody = 1.0;
  AD Rbody = frbody * rbody * rhalo / (2 * rbody + rhalo * leff) * weff / nseg;
  AD oxideRatio = pow(toxref/toxqm, ntox) / toxqm / toxqm;
  AD GatesidewallJctPotential = 0.9;
  AD bodyJctGateSideGradingCoeff = 0.5;
  AD unitLengthGateSidewallJctCap = 8.0e-10;
  AD tt = 1e-12;
  AD sourceArea = 0.0;
  AD drainArea =0.0;
  AD Cbox = 3.453133e-11 / tbox; //---Check
  AD csbox = Cbox * sourceArea;
  AD cdbox = Cbox * drainArea;
  AD csmin = Cbox * sourceArea;
  AD cdmin = Cbox * drainArea;
  double Vtm0 = KboQ * tnom;
  double Eg0 = 1.16 - 7.02e-4 * tnom * tnom / (tnom + 1108.0);
  double ni0 = 1.45e10 * (tnom / 300.15) * sqrt(tnom / 300.15) * exp(21.5565981 - Eg0 / (2.0 * Vtm0));/* ni is in cm^-3 */

  AD Temp=0.0;
  AD TempRatio=0.0;
  if(selfheat)
  {
    Temp = x[1] +temp;
    TempRatio = Temp / tnom; /* calulcatedtemp/tnom */
  }
  else
  {
    Temp = temp;
    TempRatio = Temp / tnom; /* ckt temp/tnom */
  }
  AD Vtm = KboQ * Temp;
  AD Eg = 1.16 - 7.02e-4 * Temp * Temp / (Temp + 1108.0);
  t1 = ((7.02e-4 * t5) - t0 * (14.04e-4 * Temp)) / t0 / t0;
  t2 = 1.9230584e-4;
  t5 = sqrt(Temp);
  t3 = 1.45e10 * Temp * t5 * t2;
  t4 = exp(21.5565981 - Eg / (2.0 * Vtm));
  AD  ni = t3 * t4;
  t0 = log(1.0e20 * Nch / (ni * ni));
  AD vbi = Vtm * t0;
  AD phi = 2.0 * Vtm * log(npeak / ni);
  if (nsub > 0)
  {
    t0 = log(Nch / nsub);
    vfbb = -type * Vtm * t0;
  }
  else
  {
    t0 = log(-Nch * nsub / ni / ni);
    vfbb = -type * Vtm * t0;
  }

  AD sqrtPhi = sqrt(phi);
  AD Xdep0 = sqrt(2.0 * EPSSI / (Charge_q * npeak * 1.0e6)) * sqrtPhi;
  t3 = TempRatio - 1.0;
  t8 = 1/ tnom;
  t4 = Eg300 / Vtm * t3;
  t7 = xbjt * t4 /ndiode;
  DEXP1(t7,t0)

  if (xbjt == xdif)
  {
    t1 = t0;
  }
  else
  {
    t7 = xdif * t4 / ndiode;
    DEXP1(t7,t1);
  }
  t7 = xrec * t4 / nrecf0;
  DEXP1(t7,t2);
  /* high level injection */
  Ahli = ahli * t0;
  AD jbjt = isbjt * t0;
  AD jdif = isdif * t1;
  AD jrec = isrec * t2;
  t7 = xtun * t3;
  DEXP1(t7,t0);
  AD jtun = istun * t0;
  AD u0temp = u0 * pow(TempRatio, ute);
  AD vsattemp = vsat -at * t3;
  AD rds0 = (rdsw + prt * t3) / rds0denom;

  if(selfheat)
  {
    Ua = uatemp + ua1 * t3;
    Ub = ubtemp + ub1 * t3;
    Uc = uctemp + uc1 * t3;
  }
  else
  {
    Ua = ua;
    Ub = ub;
    Uc = uc;
  }

  AD TempRatioMinus1 = Temp / tnom - 1.0;
  AD vtm = Vtm;
  /* v2.2.2 bug fix */
  if (Vsdfb==0)
  {
    if (nsub > 0)
      Vsdfb = -type * (vtm*log(1e20 * nsub / ni /ni) - 0.3);
    else if (nsub < 0)
      Vsdfb = -type * (vtm*log(-1e20 / nsub) + 0.3);
  }

  AD gamma1 = 5.753e-12 * sqrt(Nch) / cox;
  AD gamma2 = 5.753e-12 * sqrt(nsub) / cox;
  AD cf = 2.0 * EPSOX / PI * log(1.0 + 0.4e-6 / tox);
  AD Cgdo = (cgdo + cf) * wdiodCV;
  AD Cgso = (cgso + cf) * wdiosCV;
  AD Cgeo = cgeo * leffCV;

  if ((Nch==0) && (gamma1!=0))
  {
    T0 = gamma1 * cox;
    Nch = 3.021E22 * TempRatioMinus1 * TempRatioMinus1;
  }

  if (u0 > 1.0)
    u0 = u0 / 1.0e4;

  AD phis3 = sqrtPhi * phi;
  AD sqrtXdep0 = sqrt(Xdep0);
  AD litl = sqrt(3.0 * xj * tox);
  AD cdep0 = sqrt(Charge_q * EPSSI * npeak * 1.0e6 / 2.0 / phi);

  /* v3.0 */
  AD vfbsd;
  if (ngate > 0.0)
  {
    vfbsd = Vtm0 * log(ngate / 1.0e20);
  }
  else
    vfbsd = 0.0;

  AD poxedge = 1.0;
  AD ToxRatio = exp(ntox * log(toxref /toxqm)) /toxqm /toxqm;
  AD ToxRatioEdge = exp(ntox * log(toxref / (toxqm * poxedge))) / toxqm / toxqm / poxedge / poxedge;
  AD  Aechvb = (type == 1.0) ? 4.97232e-7 : 3.42537e-7;
  AD Bechvb = (type == 1.0) ? 7.45669e11 : 1.16645e12;
  AD AechvbEdge = Aechvb * weff/nseg * DLCIG * ToxRatioEdge; /* v3.1 bug fix  DLCIG*/
  AD BechvbEdge = -Bechvb * toxqm * poxedge;
  Aechvb *= weff/nseg * leff * ToxRatio; /* v3.1 bug fix */
  Bechvb *= -toxqm;

  /* v3.0 */
  /*Fix this k1 k2 problem */
  /* figure out a way to find out if K1 , K2 are given are not */
  AD vbx;
  if (vbx==0)
    vbx = phi - 7.7348e-4 * Nch * xt1 * xt1;
  if (vbx > 0.0)
    vbx = -vbx;
  if (vbm > 0.0)
    vbm = -vbm;

  if (gamma1==0)
    gamma1 = 5.753e-12 * sqrt(Nch) / cox;
  if (gamma2==0)
    gamma2 = 5.753e-12 * sqrt(nsub) / cox;

  t0 = gamma1 - gamma2;
  t1 = sqrt(phi - vbx) - sqrtPhi;
  t2 = sqrt(phi * (phi - vbm)) - phi;
  K2 = t0 * t1 / (2.0 * t2 + vbm);
  K1 = gamma2 - 2.0 * K2 * sqrt(phi - vbm);

  if (K2 < 0.0)
  {
    t0 = 0.5 * K1 / K2;
    vbsc = 0.9 * (phi - t0 * t0);
    if (vbsc > -3.0)
      vbsc = -3.0;
    else if (vbsc < -30.0)
      vbsc = -30.0;
  }
  else
    vbsc = -30.0;

  if (vbsc > vbm)
    vbsc = vbm;

  if ((t0 = weff + k1w2) < 1e-8)
    t0 = 1e-8;

  K1 = 0.56;
  AD k1eff = K1 * (1 + k1w1/t0);
  AD vfb;

  if (Vth0!=0)
  {
    vfb = type * Vth0 - phi - k1eff * sqrtPhi;
  }
  else
  {
    vfb = -1.0;
    Vth0 = type * (vfb + phi + k1eff * sqrtPhi);
  }

  k1eff *= tox / toxm;
  K2 *= tox / toxm;
  t1 = sqrt(EPSSI / EPSOX * tox * Xdep0);
  t0 = exp(-0.5 * dsub * leff / t1);
  AD theta0vb0 = (t0 + 2.0 * t0 * t0);
  t0 = exp(-0.5 * drout * leff / t1);
  t2 = (t0 + 2.0 * t0 * t0);
  AD thetaRout = pdibl1 * t2 + pdibl2;

  if (((nsub > 0) && (type > 0)) ||((nsub < 0) && (type < 0)))
  {
    t0 = vsdth - Vsdfb;
    sdt1 = Vsdfb + asd * t0;
    t1 = csbox - csmin;
    t2 = t1 / t0 / t0;
    st2 = t2 / asd;
    st3 = t2 /( 1 - asd);
    st4 = t0 * t1 * (1 + asd) / 3- csmin * Vsdfb;

    t1 = cdbox - cdmin;
    t2 = t1 / t0 / t0;
    dt2 = t2 / asd;
    dt3 = t2 /( 1 - asd);
    dt4 = t0 * t1 * (1 + asd) / 3 - cdmin * Vsdfb;
  }
  else
  {
    t0 = Vsdfb - vsdth;
    sdt1 = vsdth + asd * t0;
    t1 = csmin - csbox;
    t2 = t1 / t0 / t0;
    st2 = t2 / asd;
    st3 = t2 /( 1 - asd);
    st4 = t0 * t1 * (1 + asd) / 3 - csbox * vsdth;

    t1 = cdmin - cdbox;
    t2 = t1 / t0 / t0;
    dt2 = t2 / asd;
    dt3 = t2 /( 1 - asd);
    dt4 = t0 * t1 * (1 + asd) / 3 - cdbox * vsdth;
  }

  /* v2.0 release */
  if (ln < 1e-15) ln = 1e-15;
  t0 = -0.5 * leff * leff / ln / ln;
  DEXP1(t0,t1);
  AD arfabjt = t1;

  t0 = lbjt0 * (1.0 / leff + 1.0 / ln);
  AD lratio = pow(t0,nbjt);
  AD lratiodif = 1.0 + ldif0 * pow(t0,ndif);

  if ((vearly = vabjt + aely * leff) < 1)
    vearly = 1;

  AD qsi = Charge_q * npeak * (1.0 + nlx / leff) * 1e6 * tsi;
  AD ldeb = sqrt(EPSSI * Vtm0 / (Charge_q * Nch * 1.0e6)) / 3.0;
  AD csi = 1.03594e-10 / tsi;
  AD Vesfb = x[2] - vfbb;

  //**************Mosnbsim3SOI5T1T*****************//

  /* Poly Gate Si DEpeletion Effect */
  t0 = vfb + phi;

  if((ngate > 1.e18) && (ngate < 1.e25) && (x[3] > t0))
  {
    t1 = 1.0e6 * Charge_q * EPSSI * ngate / (cox * cox);
    T4 = sqrt(1.0 + 2.0 * (x[3] - t0) / t1);
    T2 = t1 * (T4 -1.0);
    T3 = 0.5 * T2 * T2 / t1;
    T7 = 1.12 - T3 - 0.05;
    T6 = sqrt(T7 * T7 + 0.224);
    T5 = 1.12 - 0.5 * (T7 + T6);
    Vgs_eff = x[3] - T5;
  }
  else
  {
    Vgs_eff = x[3];
  }

  AD Leff = leff;
  V0 = vbi - phi;

  /* begin of v3.0 block addition */
  /* B/S built-in potential lowering calculation */

  //or 0 for FD module
  /* Prepare Vbs0 & Vbs0mos for VthFD calculation */
  t0 = - dvbd1 * leff / litl;
  t1 = dvbd0 * (exp(0.5 * t0) + 2 * exp(t0));
  t2 = t1 * (vbi - phi);
  t3 = 0.5 * qsi / csi;
  AD Vbs0t = phi - t3 + vbsa + t2;
  t0 = 1 + csi / Cbox;
  t3 = -dk2b * leff / litl;
  t5 = k2b * (exp(0.5 * t3) + 2 * exp(t3));
  t1 = (k1b - t5) / t0;
  T2 = t1 * Vesfb;
  t4 = 1.0 / (1 + Cbox / csi);

  AD Vbs0 = t4 * Vbs0t + T2;

  /* Zero field body potential cal. */
  //Vbs0 is AD
  T1 = Vbs0t - Vbs0 - 0.005;
  T2 = sqrt(T1 * T1 + (2.5e-5));
  T3 = 0.5 * (T1 + T2);
  T4 = T3 * csi / qsi;
  AD Vbs0mos = Vbs0 - 0.5 * T3 * T4;
  T5 = 0.5 * T4 * (1 + T1 / T2);

  /* Set the upper bound of Vbs0mos to be phi for square root calc. */
  t1 = phi - 0.02;
  T2 = t1 - Vbs0mos - 0.005;
  T3 = sqrt(T2 * T2 + 4.0 * 0.005);
  Vbs0mos = t1 - 0.5 * (T2 + T3);
  T4 = 0.5 * (1 + T2 / T3);

  /* VthFD calculation */
  AD Theta0, thetavth, Delt_vth, DeltVthw;
  Phis = phi - Vbs0mos;
  sqrtPhis = sqrt(Phis);
  Xdep = Xdep0 * sqrtPhis / sqrtPhi;
  T3 = sqrt(Xdep);
  T0 = dvt2 * Vbs0mos;
  if (T0 >= -0.5)
  {
    T1 = 1.0 + T0;
    t2 = dvt2;
  }
  else
  {
    T4 = 1.0 / (3.0 + 8.0 * T0);
    T1 = (1.0 + 3.0 * T0) * T4;
    T2 = dvt2 * T4 * T4;
  }

  lt1 = factor1 * T3 * T1;

  T0 = dvt2w * Vbs0mos;
  if (T0 >= -0.5)
  {
    T1 = 1.0 + T0;
    t2 = dvt2w;
  }
  else
  {
    T4 = 1.0 / (3.0 + 8.0 * T0);
    T1 = (1.0 + 3.0 * T0) * T4;
    T2 = dvt2w * T4 * T4;
  }
  ltw = factor1 * T3 * T1;

  T0 = -0.5 * dvt1 * Leff / lt1;

  if (T0 > -EXPL_THRESHOLD)
  {
    T1 = exp(T0);
    Theta0 = T1 * (1.0 + 2.0 * T1);
  }
  else
  {
    T1 = MIN_EXPL;
    Theta0 = T1 * (1.0 + 2.0 * T1);
  }
  thetavth = dvt0 * Theta0;
  Delt_vth = thetavth * V0;

  T0 = -0.5 * dvt1w * weff * Leff / ltw;
  if (T0 > -EXPL_THRESHOLD)
  {
    T1 = exp(T0);
    T2 = T1 * (1.0 + 2.0 * T1);
  }
  else
  {
    T1 = MIN_EXPL;
    T2 = T1 * (1.0 + 2.0 * T1);
  }

  T0 = dvt0w * T2;
  DeltVthw = T0 * V0; //DeltVthw is T2 in Bsim3 code

  t0 = sqrt(1.0 + nlx / Leff);
  T1 = (kt1 + kt1l / Leff + kt2 * Vbs0mos);
  DeltVthtemp = k1eff * (t0 - 1.0) * sqrtPhi + T1 * TempRatioMinus1;

  tmp2 = tox * phi / (weff + w0);

  T3 = eta0 + etab * Vbs0mos;
  if (-T3 + 1.0e-4 > 0.0)
    T3 = (2.0e-4 - T3)/ (3.0 - 2.0e4 * T3);

  DIBL_Sft = T3 * theta0vb0 * x[0];

  K2=0.0;
  AD VthFD = type * Vth0 + k1eff * (sqrtPhis - sqrtPhi) - K2 * Vbs0mos
             - Delt_vth - DeltVthw + (k3 + k3b * Vbs0mos) * tmp2
             + DeltVthtemp - DIBL_Sft; //--Not changed '1.0' to "type" parameter

  /* VtgseffFD calculation for PhiFD */

  AD VtgsFD, ExpVtgsFD, VtgseffFD;
  VtgsFD = VthFD - Vgs_eff;
  t10 = nofffd * Vtm;

  DEXP((VtgsFD - vofffd) / t10, ExpVtgsFD, T0);

  VtgseffFD = t10 * log(1.0 + ExpVtgsFD);
  T0 /= (1.0 + ExpVtgsFD);

  /* Surface potential modeling at strong inversion: PhiON */

  AD VgstFD, ExpVgstFD, VgsteffFD;
  VgstFD = Vgs_eff - VthFD;

  AD A,B,C;
  A = (VgstFD - vofffd) / t10;
  B =  ExpVgstFD;
  C =  T0;

  if (A > EXPL_THRESHOLD)
  {
    B = MAX_EXPL*(1.0 + (A)-EXPL_THRESHOLD);
    C = MAX_EXPL;
  }
  else if (A < -EXPL_THRESHOLD)
  {
    B = MIN_EXPL;
    T0 = 0;
  }
  else
  {
    B = exp(A);
    C = B;
  }

  ExpVgstFD = B;
  T0 = C;

  VgsteffFD = t10 * log(1.0 + ExpVgstFD);

  T0 /= (1.0 + ExpVgstFD);
  T1 = moinFD * k1eff * Vtm * Vtm;
  T2 = VgsteffFD + 2 * k1eff * sqrt(phi);
  T0 = 1 + VgsteffFD * T2 / T1;
  AD PhiON = phi + Vtm * log(T0);

  /* Surface potential from subthreshold to inversion: PhiFD */
  t0 = cox / (cox + 1.0 / (1.0 / csi + 1.0 / Cbox));
  AD PhiFD = PhiON - t0 * VtgseffFD;

  /* built-in potential lowering: Vbs0 */
  t0 = -dvbd1 * leff / litl;
  t1 = dvbd0 * (exp(0.5 * t0) + 2 * exp(t0));
  t2 = t1 * (vbi - phi);
  t3 = 0.5 * qsi / csi;
  Vbs0t = PhiFD - t3 + vbsa + t2;

  t0 = 1 + csi / Cbox;
  t3 = -dk2b * leff / litl;
  t5 = k2b * (exp(0.5 * t3) + 2 * exp(t3));
  t1 = (k1b - t5) / t0;
  T2 = t1 * Vesfb;
  t0 = 1.0 / (1 + Cbox / csi);
  Vbs0 = t0 * Vbs0t + T2;

  /*set lowerbound for Vbs to VBS0: Vbsitf */
  /* soiMod = 1 */

  T1 = 0.0 - (Vbs0 + OFF_Vbsitf) - 0.01; //Vbs=0.0
  T2 = sqrt(T1 * T1 + 0.0001);
  T3 = 0.5 * (1 + T1/T2);
  AD Vbsitf = (Vbs0 + OFF_Vbsitf) + 0.5 * (T1 + T2);

  /* Based on Vbsitf, calculate zero-field body potential for MOS: Vbsmos */
  T1 = Vbs0t - Vbsitf - 0.005;
  T2 = sqrt(T1 * T1 + (2.5e-5));
  T3 = 0.5 * (T1 + T2);
  T4 = T3 * csi / qsi;
  Vbsmos = Vbsitf - 0.5 * T3 * T4;
  /* end of v3.0 block edition */

  /* v3.0 modification */
  /* T2 is Vbsmos limited above Vbsc=-5 */
  T0 = Vbsmos + 5 - 0.001;
  T1 = sqrt(T0 * T0 - 0.004 * (-5));
  T2 = (-5) + 0.5 * (T0 + T1);

  /* Vbsh is T2 limited below 1.5 */
  t0 = 1.5;
  T1 = t0 - T2 - 0.002;
  T3 = sqrt(T1 * T1 + 0.008 * t0);
  Vbsh = t0 - 0.5 * (T1 + T3);

  /* Vbseff is Vbsh limited to 0.95*phi */
  t0 = 0.95 * phi;
  T1 = t0 - Vbsh - 0.002;
  T2 = sqrt(T1 * T1 + 0.008 * t0);
  Vbseff = t0 - 0.5 * (T1 + T2);
  /* END OF MODIFICATION */

  Phis = phi - Vbseff;
  sqrtPhis = sqrt(Phis);
  Xdep = Xdep0 * sqrtPhis / sqrtPhi;
  //***********************************************
  //From here I am going to add lines of Bsim3 code
  //***********************************************

  /* Vth Calculation */
  T3 = sqrt(Xdep);
  T0 = dvt2 * Vbseff;
  if (T0 >= - 0.5)
  {
    T1 = 1.0 + T0;
    T2 = dvt2 ;
  }
  else /* Added to avoid any discontinuity problems caused by dvt2 */
  {
    T4 = 1.0 / (3.0 + 8.0 * T0);
    T1 = (1.0 + 3.0 * T0) * T4;
    T2 = dvt2 * T4 * T4 ;
  }

  lt1 = factor1 * T3 * T1;
  T0 = dvt2w * Vbseff;
  if (T0 >= - 0.5)
  {
    T1 = 1.0 + T0;
    T2 = dvt2w;
  }
  else /* Added to avoid any discontinuity problems caused by dvt2w */
  {
    T4 = 1.0 / (3.0 + 8.0 * T0);
    T1 = (1.0 + 3.0 * T0) * T4;
    T2 = dvt2w * T4 * T4 ;
  }

  ltw = factor1 * T3 * T1;

  T0 = -0.5 * dvt1 * Leff / lt1;
  if (T0 > -EXPL_THRESHOLD)
  {
    T1 = exp(T0);
    Theta0 = T1 * (1.0 + 2.0 * T1);
  }
  else
  {
    T1 = MIN_EXPL;
    Theta0 = T1 * (1.0 + 2.0 * T1);
  }

  thetavth = dvt0 * Theta0;
  Delt_vth = thetavth * V0;

  T0 = -0.5 * dvt1w * weff * Leff / ltw;
  if (T0 > -EXPL_THRESHOLD)
  {
    T1 = exp(T0);
    T2 = T1 * (1.0 + 2.0 * T1);
  }
  else
  {
    T1 = MIN_EXPL;
    T2 = T1 * (1.0 + 2.0 * T1);
  }

  T0 = dvt0w * T2;
  DeltVthw = T0 * V0; //DeltVthw is T2 in Bsim3 code

  T0 = sqrt(1.0 + nlx / Leff);

  T1 = (kt1 + kt1l / Leff + kt2 * Vbseff);
  DeltVthtemp = k1eff * (T0 - 1.0) * sqrtPhi + T1 * TempRatioMinus1;
  tmp2 = tox * phi / (weff + w0);

  T3 = eta0 + etab * Vbseff;
  DIBL_Sft = T3 * theta0vb0 * x[0];

  T9 = 2.2361 / sqrtPhi;
  sqrtPhisExt = sqrtPhis - T9 * (Vbsh - Vbseff);
  Vth = 1.0 * Vth0 + k1eff * (sqrtPhisExt - sqrtPhi) - K2 * Vbseff
        - Delt_vth - DeltVthw + (k3 + k3b *Vbseff) * tmp2 + DeltVthtemp - DIBL_Sft;
  AD von = Vth;

  /* Calculate n */
  T2 = nfactor * EPSSI / Xdep;
  T3 = cdsc + cdscb * Vbseff + cdscd * x[0];
  T4 = (T2 + T3 * Theta0 + cit) / cox;

  if (T4 >= -0.5)
  {
    n = 1.0 + T4;
  }/* avoid  discontinuity problems caused by T4 */
  else
  {
    T0 = 1.0 / (3.0 + 8.0 * T4);
    n = (1.0 + 3.0 * T4) * T0;
    T0 *= T0;
  }

  /* Poly Gate Si effect is already discussed in the very beginning */

  /* Effective Vgst (Vgsteff) Calculation */
  Vgst = Vgs_eff - Vth;
  T10 = 2.0 * n * Vtm;
  VgstNVt = Vgst / T10;
  ExpArg = (2.0 * voff - Vgst) / T10;

  AD ExpVgst, dT2_dVg;
  if (VgstNVt > EXPL_THRESHOLD)
  {
    Vgsteff = Vgst;
  }
  else if (ExpArg > EXPL_THRESHOLD)
  {
    T0 = (Vgst - voff) / (n * Vtm);
    ExpVgst = exp(T0);
    Vgsteff = Vtm * cdep0 / cox * ExpVgst;
    T3 = Vgsteff / (n * Vtm) ;
  }
  else
  {
    ExpVgst = exp(VgstNVt);
    T1 = T10 * log(1.0 + ExpVgst);
    dT2_dVg = -cox / (Vtm * cdep0) * exp(ExpArg);
    T2 = 1.0 - T10 * dT2_dVg;
    Vgsteff = T1 / T2;
  }

  Vgst2Vtm = Vgsteff + 2.0 * Vtm;

  /* Calculate Effective Channel Geometry */
  T9 = sqrtPhis - sqrtPhi;
  Weff = weff - (2.0 - nbc) * (dwg * Vgsteff + dwb * T9);
  if (Weff < 2.0e-8) /* to avoid the discontinuity problem due to Weff*/
  {
    T0 = 1.0 / (6.0e-8 - 2.0 * Weff);
    Weff = 2.0e-8 * (4.0e-8 - Weff) * T0;
  }
  T0 = prwg * Vgsteff + prwb * T9;
  if (T0 >= -0.9)
  {
    Rds = rds0 * (1.0 + T0);
  }
  else
   /* to avoid the discontinuity problem due to prwg and prwb*/
  {
    T1 = 1.0 / (17.0 + 20.0 * T0);
    Rds = rds0 * (0.8 + T0) * T1;
  }

  /* Calculate Abulk */
  if (a0 == 0.0)
  Abulk0 = Abulk = 1.0;
  else
  {
    T10 = keta * Vbsh;
    if (T10 >= -0.9)
    {
      T11 = 1.0 / (1.0 + T10);
    }
    else
    { /* added to avoid the problems caused by Keta */
      T12 = 1.0 / (0.8 + T10);
      T11 = (17.0 + 20.0 * T10) * T12;
    }

    T10 = phi + ketas;
    T13 = (Vbsh * T11) / T10;

    if(T13 < 0.96)
    {
      T14 = 1/ sqrt(1-T13);
      T10 = 0.5 * T14 / (1-T13);
    }
    else
    {
      T11 = 1.0 / (1.0 - 1.043406 * T13);
      T14 = (6.00167 - 6.26044 * T13) * T11;
      T10 = 0.001742 * T11 * T11;
    }

    t10 = 0.5 * k1eff / sqrt(phi + ketas);
    T1 = t10 * T14;

    T9 = sqrt(xj * Xdep);
    TMP11 = Leff + 2.0 * T9;
    T5 = Leff / TMP11;
    TMP12 = a0 * T5;
    tmp3 = weff + b1;
    tmp4 = b0 / tmp3;
    T2 = TMP12 + tmp4;
    T6 = T5 * T5;
    T7 = T5 * T6;

    Abulk0 = 1 + T1 * T2;

    T8 = ags * a0 * T7;
    Abulk = Abulk0 + (-T1 * T8) * Vgsteff;
  }

  if (Abulk0 < 0.01)
  {
    T9 = 1.0 / (3.0 - 200.0 * Abulk0);
    Abulk0 = (0.02 - Abulk0) * T9;
  }

  if (Abulk < 0.01)
  {
    T9 = 1.0 / (3.0 - 200.0 * Abulk);
    Abulk = (0.02 - Abulk) * T9;
  }

  /* Mobility calculation */
  //check the following
  // mobMod == 1 (default)
  T0 = Vgsteff + Vth + Vth;
  T2 = Ua + Uc * Vbseff;
  T3 = T0 / tox;
  T5 = T3 * (T2 + Ub * T3);

  if (T5 >= -0.8)
  {
    Denomi = 1.0 + T5;
  }
  else /* Added to avoid the discontinuity problem caused by ua and ub*/
  {
    T9 = 1.0 / (7.0 + 10.0 * T5);
    Denomi = (0.6 + T5) * T9;
    T9 *= T9;
  }
  ueff = u0temp / Denomi;

  /* Saturation Drain Voltage Vdsat*/
  WVCox = Weff * vsattemp * cox;
  WVCoxRds = WVCox * Rds;
  Esat = 2.0 * vsattemp / ueff;
  EsatL = Esat * Leff;

  //AD Lambda = a2; //CHECK
  /* Sqrt() */
  // a1 = a1;
  AD Lambda, dLambda_dVg;
  if (a1 == 0.0)
  {
    Lambda = a2;
    dLambda_dVg = 0.0;
  }

  AbovVgst2Vtm = Abulk / Vgst2Vtm;

  if ((Rds == 0.0) && (Lambda == 1.0))
  {
    T0 = 1.0 / (Abulk * EsatL + Vgst2Vtm);
    tmp1 = 0.0;
    T1 = T0 * T0;
    T2 = Vgst2Vtm * T0;
    T3 = EsatL * Vgst2Vtm;
    Vdsat = T3 * T0;
  }
  else
  {
    tmp1 = dLambda_dVg / (Lambda * Lambda);
    T9 = Abulk * WVCoxRds;
    T8 = Abulk * T9;
    T7 = Vgst2Vtm * T9;
    T6 = Vgst2Vtm * WVCoxRds;
    T0 = 2.0 * Abulk * (T9 - 1.0 + 1.0 / Lambda);
    T1 = Vgst2Vtm * (2.0 / Lambda - 1.0) + Abulk * EsatL + 3.0 * T7;
    T2 = Vgst2Vtm * (EsatL + 2.0 * T6);
    T3 = sqrt(T1 * T1 - 2.0 * T0 * T2);
    Vdsat = (T1 - T3) / T0;
  }

  /* Effective Vds (Vdseff) Calculation */
  T1 = Vdsat - x[0] - delta;
  T2 = sqrt(T1 * T1 + 4.0 * delta * Vdsat);
  T0 = T1 / T2;
  T3 = 2.0 * delta / T2;

  Vdseff = Vdsat - 0.5 * (T1 + T2);
  if (Vdseff > x[0])
    Vdseff = x[0];

  diffVds = x[0] - Vdseff;

  /* Calculate VAsat */
  TMP4 = 1.0 - 0.5 * Abulk * Vdsat / Vgst2Vtm;
  T9 = WVCoxRds * Vgsteff;//expanded
  T8 = T9 / Vgst2Vtm;
  T0 = EsatL + Vdsat + 2.0 * T9 * TMP4;
  T9 = WVCoxRds * Abulk;
  T1 = 2.0 / Lambda - 1.0 + T9;
  Vasat = T0 / T1;

  /* Calculate VACLM */
  if ((pclm > 0.0) && (diffVds > 1.0e-10))
  {
    T0 = 1.0 / (pclm * Abulk * litl);
    T2 = Vgsteff / EsatL;
    T1 = Leff * (Abulk + T2);
    T9 = T0 * T1;
    VACLM = T9 * diffVds;
  }
  else
  {
    VACLM = MAX_EXPL;
  }

  if (thetaRout > 0.0)
  {
    T8 = Abulk * Vdsat;
    T0 = Vgst2Vtm * T8;
    T1 = Vgst2Vtm + T8;

    T9 = T1 * T1;
    T2 = thetaRout;
    VADIBL = (Vgst2Vtm - T0 / T1) / T2;

    T7 = pdiblb * Vbseff;
    if (T7 >= -0.9)
    {
      T3 = 1.0 / (1.0 + T7);
      VADIBL *= T3;
    }
    else
    /* Added to avoid the discontinuity problem caused by pdiblcb */
    {
      T4 = 1.0 / (0.8 + T7);
      T3 = (17.0 + 20.0 * T7) * T4;
      VADIBL *= T3;
    }
  }
  else
  {
    VADIBL = MAX_EXPL;
  }

  /* Calculate VA */
  T8 = pvag / EsatL;
  T9 = T8 * Vgsteff;

  if (T9 > -0.9)
  {
    T0 = 1.0 + T9;
  }
  else /* Added to avoid the discontinuity problems caused by pvag */
  {
    T1 = 1.0 / (17.0 + 20.0 * T9);
    T0 = (0.8 + T9) * T1;
    T1 *= T1;
    T9 *= T1 / EsatL;
  }

  TMP11 = VACLM * VACLM;
  TMP12 = VADIBL * VADIBL;
  TMP13 = VACLM + VADIBL;

  T1 = VACLM * VADIBL / TMP13;
  Va = Vasat + T0 * T1;

  /* Calculate Ids */
  CoxWovL = cox * Weff / Leff;
  beta = ueff * CoxWovL;
  T0 = 1.0 - 0.5 * Abulk * Vdseff / Vgst2Vtm;
  fgche1 = Vgsteff * T0;
  T9 = Vdseff / EsatL;
  fgche2 = 1.0 + T9;
  gche = beta * fgche1 / fgche2;
  T0 = 1.0 + gche *Rds;
  T9 = Vdseff / T0;
  AD Idl = gche * T9;
  T9 = diffVds / Va;
  T0 = 1.0 + T9;
  AD Ids = Idl * T0 / nseg;

  // /****************************************/
  // /* all the extra currents */
  // /* soiMod != 2 v3.2*/
  //
  // /* Calculate GIDL current */
  //
  //
  AD Idgidl, Isgidl;
  T0 = 3 * tox;
  // /* For drain side */
  T1 = (x[0] - Vgs_eff - ngidl) / t0;
  //
  if ((agidl <= 0.0) || (bgidl <= 0.0) || (T1 <= 0.0))
    Idgidl = 0.0;
  else
  {
    T2 = bgidl / T1;
    if (T2 < EXPL_THRESHOLD)
    {
      Idgidl = wdiod * agidl * T1 * exp(-T2);
      T3 = Idgidl / T1 * (T2 + 1.0);
    }
    else
    {
      T3 = wdiod * agidl * MIN_EXPL;
      Idgidl = T3 * T1;
    }
  }

  // /* For source side */
  T1 = (- Vgs_eff - ngidl) / T0;

  if ((agidl <= 0.0) || (bgidl <= 0.0) || (T1 <= 0.0))
    Isgidl = 0.0;
  else
  {
    T2 = bgidl / T1;
    if (T2 < EXPL_THRESHOLD)
    {
      Isgidl = wdios * agidl * T1 * exp(-T2);
      T3 = Isgidl / T1 * (T2 + 1.0);
    }
    else
    {
      T3 = wdios * agidl * MIN_EXPL;
      Isgidl = T3 * T1;
    }
  }

  /* calculate diode and BJT current */
  AD ExpVbsNVtm, ExpVbdNVtm;
  AD WsTsi = wdios * tsi;
  AD WdTsi = wdiod * tsi;
  AD NVtm1 = Vtm * ndiode;
  AD E2ndFactor , ic, Ic;

  // T0 = Vbs / NVtm1;
  T0 = x[2] / NVtm1;
  DEXP(T0,ExpVbsNVtm,T1) ;
  T0 = (x[2] - x[0]) / NVtm1;
  DEXP(T0, ExpVbdNVtm, T1);

  /* Ibs1 / Ibd1 : diffusion current */
  AD Ibs1, Ibd1, Ibs2, Ibd2;
  jbjt = isbjt;
  jdif = isdif;
  jrec = isrec;
  jtun = istun;

  if (jdif == 0.0)
    Ibs1 = Ibd1 = 0.0;
  else
  {
    T0 = WsTsi * jdif;
    Ibs1 = T0 * (ExpVbsNVtm - 1.0);

    T0 = WdTsi * jdif;
    Ibd1 = T0 * (ExpVbdNVtm - 1.0);
  }

  /* Ibs2 : recombination/trap-assisted tunneling current */
  AD NVtmf = 0.026 * nrecf0 * (1 + ntrecf * (TempRatio - 1));
  AD NVtmr = 0.026 * nrecr0 * (1 + ntrecr * (TempRatio - 1));
  if (jrec == 0.0)
    Ibs2 = Ibd2 = 0.0;
  else
  {
    T0 = x[2] / NVtmf; //T0 = Vbs / NVtm1;
    DEXP(T0, T10, T2);
    T4 = 1 / NVtmf;
    if ((vrec0 - x[2]) < 1e-3)  //Vbs
    {
      T1 = 1e3;
      T0 = -x[2] / NVtmr * vrec0 * T1;   //Vbs
      T11 = -exp(T0);
    }
    else
    {
      T1 = 1 / (vrec0 - x[2]+x[0]);
      T0 = -x[2] / NVtmr * vrec0 * T1;
      DEXP(T0, T11, T2);
      T11 = -T11;
    }
    T3 = WsTsi * jrec;
    Ibs2 = T3 * (T10 + T11);

    // /* Ibd2 */
    T0 = (x[2] - x[0]) / NVtmf; //Vbd
    DEXP(T0, T10, T2);
    T4 = 1 / NVtmf;
    if ((vrec0 - (x[2] - x[0])) < 1e-3)
    {
      T1 = 1e3;
      T0 = -(x[2] - x[0]) / NVtmr * vrec0 * T1;
      T11 = -exp(T0);
    }
    else
    {
      T1 = 1 / (vrec0 - (x[2] - x[0]));
      T0 = (x[0] - x[2]) / NVtmr * vrec0 * T1;
      DEXP(T0, T11, T2);
      T11 = -T11;
    }
    T3 = WdTsi * jrec;
    Ibd2 = T3 * (T10 + T11);
  }

  // /* Ibs3/Ibd3 : recombination current in neutral body */
  //
  AD Ibs3, Ibd3, Ibs4, Ibd4, Ehlis, EhlisFactor, Ehlid, EhlidFactor, Ibsdif, Ibddif;
  AD Ien;
  //
  t0 = lbjt0 * (1.0 / leff + 1.0 / ln);
  lratio = pow(t0, nbjt);
  lratiodif = 1.0 + ldif0 * pow(t0, ndif);

  AD WTsi = weff / nseg * tsi;
  if (jbjt == 0.0)
    Ibs3 = Ibd3 = Ibsdif = Ibddif = 0.0;
  else
  {
    Ien = WTsi * jbjt * lratio;
    if ((Ehlis = Ahli * (ExpVbsNVtm - 1)) < 1e-5)
    {
      Ehlis = 0.0;
      EhlisFactor = 1.0;
    }
    else
    {
      EhlisFactor = 1.0 / sqrt(1 + Ehlis);
    }

    if ((Ehlid = Ahli * (ExpVbdNVtm - 1)) < 1e-5)
    {
      Ehlid = 0.0;
      EhlidFactor = 1.0;
    }
    else
    {
      EhlidFactor = 1.0 / sqrt(1 + Ehlid);
    }

    /* Ibjt(L) */
    if (ln < 1e-15)
    ln = 1e-15;
    t0 = -0.5 * leff * leff / ln / ln;
    // // //#define DEXP(A,B) {
    if (t0 > EXPL_THRESHOLD)
      t1 = MAX_EXPL*(1.0+(t0)-EXPL_THRESHOLD);
    else if (t0 < -EXPL_THRESHOLD)
      t1 = MIN_EXPL;
    else
      t1 = exp(t0);

    // DEXP(TempRatioMinus1, T11);
    AD arfabjt = t1;
    t0 = 1 - arfabjt;
    t1 = t0 * Ien;
    Ibs3 = t1 * (ExpVbsNVtm - 1) * EhlisFactor;
    Ibd3 = t1 * (ExpVbdNVtm - 1) * EhlidFactor;

    /* Effective diffusion current for capacitance calculation */
    AD Iendif = WTsi * jbjt * lratiodif;
    AD Ibsdif = Iendif * (ExpVbsNVtm - 1) * EhlisFactor;
    AD Ibddif = Iendif * (ExpVbdNVtm - 1) * EhlidFactor;

    /* Ic: Bjt colletor current */
    int bjtoff =0;

    if ((bjtoff == 1) || (x[0] == 0.0))
    {
      ic = Ic = 0.0;
    }
    else
    {
      //       /* second order effects */
      T0 = 1 + (x[2] + (x[2]-x[0])) / vearly;
      T1 = Ehlis + Ehlid;
      T3 = sqrt(T0 * T0 + 4 * T1);
      T2 = (T0 + T3) / 2.0;
      if (T2 < .1)
      {
        E2ndFactor = 10.0;
      }
      else
      {
        E2ndFactor = 1.0 / T2;
      }
      T0 = arfabjt * Ien;
      ic = Ic = T0 * (ExpVbsNVtm - ExpVbdNVtm) * E2ndFactor;
    }
  }

  // /* Ibs4/Ibd4 : tunneling currents */
  //
  AD NVtm2 = 0.026 * ntun;
  if (jtun == 0.0)
    Ibs4 = Ibd4 = 0.0;
  else
  {
    if ((vtun0 - x[2] < 1e-3))
    {
      T1 = 1e3;
      T0 = -x[2] / NVtm2 * vtun0 * T1;
      T1 = exp(T0);
      T3 = WsTsi * jtun;
      Ibs4 = T3 * (1 - T1);
    }
    else
    {
      T1 = 1 / (vtun0 - x[2]);
      T0 = -x[2] / NVtm2 * vtun0 * T1;
      DEXP(T0, T1, T2);
      T3 = WsTsi * jtun;
      Ibs4 = T3 * (1 - T1);
    }

    if ((vtun0 - (x[2] - x[0])) < 1e-3)
    {
      T1 = 1e3;
      T0 = - (x[2] - x[0]) / NVtm2 * vtun0 * T1;
      T1 = exp(T0);
      T3 = WdTsi * jtun;
      Ibd4 = T3 * (1 - T1);
    }
    else
    {
      T1 = 1 / (vtun0 - (x[2] - x[0]));
      T0 = - (x[2] - x[0]) / NVtm2 * vtun0 * T1;
      DEXP(T0, T1, T2);
      T3 = WdTsi * jtun;
      Ibd4 = T3 * (1 - T1);
    }
  }

  AD itun = - Ibd3 - Ibd4;
  AD Ibs = Ibs1 + Ibs2 + Ibs3 + Ibs4;
  AD Ibd = Ibd1 + Ibd2 + Ibd3 + Ibd4;

  /*Neglecting gate-tunneling from page 40 to page 47 */

  /* v3.0: gate-tunneling */
  double igbMod = 0;
  double igcMod = 0;
  AD Vgb, Voxacc, Voxdepinv, Vaux, vgs_eff, vgd_eff, Vox, Voxeff;
  AD Igc, Igcs, Igcd, Igs, Igb, Igd, Igb1, Igb2;
  AD VxNVt, ExpVxNVt;
  if ((igbMod != 0) || (igcMod != 0))
  {
    Vgb = Vgs_eff - x[2];

    /* Calculate Vox first */
    Vfb = type * Vth0 - phi - k1eff * sqrtPhi;

    T3 = Vfb - Vgs_eff + x[2] - DELTA_3;

    if (Vfb <= 0.0)
      T0 = sqrt(T3 * T3 - 4.0 * DELTA_3 * Vfb);
    else
      T0 = sqrt(T3 * T3 + 4.0 * DELTA_3 * Vfb);

    Vfbeff = Vfb - 0.5 * (T3 + T0);

    Voxacc = Vfb - Vfbeff;

    if (Voxacc < 0.0)
      Voxacc = 0.0;

    T0 = Vgs_eff - Vgsteff - Vfbeff - Vbseff;

    if (k1eff == 0.0)
      Voxdepinv =  0.0;
    else
    {
      if (T0 < 0.0)
        T1 = T0/k1eff;
      else
      {
        T1 = k1eff/2*(-1 + sqrt(1 + 4*T0/k1eff/k1eff));
        T2 = k1eff/2 *0.5/sqrt(1 + 4*T0/k1eff/k1eff) *4/k1eff/k1eff;
      }
      Voxdepinv = Vgs_eff - (T1*T1 + x[2]) - Vfb;
    }
  }

  /* gate-channel tunneling component */
  if (igcMod)
  {
    T0 = Vtm * nigc;
    VxNVt = (Vgs_eff - type * Vth0) / T0; /* Vth instead of Vth0 may be used */
    if (VxNVt > EXPL_THRESHOLD)
      Vaux = Vgs_eff - type * Vth0;
    else if (VxNVt < -EXPL_THRESHOLD)
      Vaux = T0 * log(1.0 + MIN_EXPL);
    else
    {
      ExpVxNVt = exp(VxNVt);
      Vaux = T0 * log(1.0 + ExpVxNVt);
    }

    T2 = Vgs_eff * Vaux;
    T11 = Aechvb;
    T12 = Bechvb;
    T3 = aigc * cigc - bigc;
    T4 = bigc * cigc;
    T5 = T12 * (aigc + T3 * Voxdepinv - T4 * Voxdepinv * Voxdepinv);

    if (T5 > EXPL_THRESHOLD)
      T6 = MAX_EXPL;
    else if (T5 < -EXPL_THRESHOLD)
      T6 = MIN_EXPL;
    else
      T6 = exp(T5);

    Igc = T11 * T2 * T6;

    T7 = -pigcd * x[0];
    T8 = T7 * T7 + 2.0e-4;

    if (T7 > EXPL_THRESHOLD)
      T9 = MAX_EXPL;
    else if (T7 < -EXPL_THRESHOLD)
      T9 = MIN_EXPL;
    else
      T9 = exp(T7);

    T0 = T8 * T8;
    T1 = T9 - 1.0 + 1.0e-4;
    T10 = (T1 - T7) / T8;

    Igcs = Igc * T10;

    T1 = T9 - 1.0 - 1.0e-4;
    T10 = (T7 * T9 - T1) / T8;

    Igcd = Igc * T10;
    Igcs = Igcs;
    Igcd = Igcd;

    T0 = x[3] - vfbsd;
    vgs_eff = sqrt(T0 * T0 + 1.0e-4);

    T2 = x[3] * vgs_eff;
    T11 = AechvbEdge;
    T12 = BechvbEdge;
    T3 = aigsd * cigsd - bigsd;
    T4 = bigsd * cigsd;
    T5 = T12 * (aigsd + T3 * vgs_eff- T4 * vgs_eff * vgs_eff);

    if (T5 > EXPL_THRESHOLD)
      T6 = MAX_EXPL;
    else if (T5 < -EXPL_THRESHOLD)
      T6 = MIN_EXPL;
    else
      T6 = exp(T5);

    Igs = T11 * T2 * T6;

    T0 = (x[3] - x[0]) - vfbsd;
    vgd_eff = sqrt(T0 * T0 + 1.0e-4);

    T2 = (x[3] - x[0]) * vgd_eff;

    T5 = T12 * (aigsd + T3 * vgd_eff - T4 * vgd_eff * vgd_eff);

    if (T5 > EXPL_THRESHOLD)
      T6 = MAX_EXPL;
    else if (T5 < -EXPL_THRESHOLD)
      T6 = MIN_EXPL;
    else
      T6 = exp(T5);

    Igd = T11 * T2 * T6;
  }
  else
  {
   Igs  = 0.0;
   Igd  = 0.0;
  }

  /* gate-body tunneling component */
  AD alphaGB1 = 0.35;
  AD betaGB1 = 0.03;
  AD alphaGB2 = 0.43;
  AD betaGB2 = 0.05;
  if ( soiMod != 2)  /* v3.2 */ /* v3.1: the Igb calculation is skipped for the ideal FD mode */
  {
    Vox = Voxdepinv;
    /* Voxeff is Vox limited below Voxh */
    T0 = voxh;
    T1 = T0 - Vox - deltavox;
    T3 = sqrt(T1 * T1 + 4* deltavox * T0);
    Voxeff = T0 - 0.5 * (T1 + T3);

    Vox = Voxeff;

    T0 = (Vox -  ebg)/ vevb;
    DEXP(T0, T1, T2);        /* T1=exp(T0), T2=dT1_dT0 */
    Vaux =  vevb * log(1 + T1);

    if ( vgb1 != 0)
    {
      T0 = 1 - Vox /  vgb1;
    }
    else
      T0 = 1;

    if (T0 < 0.01)
    {
      T0 = 0.01;
    }

    /* v2.2.3 bug fix */
    T1 = Leff * Weff * 3.7622e-7 * oxideRatio / nseg;
    T2 = -3.1051e10 *  toxqm;
    T3 = alphaGB1;
    T4 = betaGB1;
    T6 = T2*(T3 - T4 * Vox) / T0;

    DEXP(T6, T5, T7); /* T5=exp(T6), T7=dT5_dT6 */
    Igb1 = T1 * Vgb * Vaux * T5;

    Vox = Voxacc;
    /* Voxeff is Vox limited below Voxh */
    T0 =  voxh;
    T1 = T0 - Vox -  deltavox;
    T3 = sqrt(T1 * T1 + 4* deltavox * T0);
    Voxeff = T0 - 0.5 * (T1 + T3);
    Vox = Voxeff;
    T0 = (-Vgb+(Vfb))/ vecb;
    DEXP(T0, T1, T2); /* T1=exp(T0), T2=dT1_dT0 */
    Vaux =  vecb* log(1 + T1);
    if ( vgb2 != 0)
    {
      T0 = 1 - Vox /  vgb2;
    }
    else
      T0 = 1;

    if (T0 < 0.01)
    {
      T0 = 0.01;
    }

    /* v2.2.3 bug fix */
    T1 = Leff * Weff * 4.9758e-7  * oxideRatio / nseg;

    T2 = -2.357e10 *  toxqm;
    T3 = alphaGB2;
    T4 = betaGB2;

    T6 = T2*(T3 - T4 * Vox) / T0;
    DEXP(T6, T5, T7); /* T5=exp(T6), T7=dT5_dT6 */

    Igb2 = T1 * (x[3]-x[2]) * Vaux * T5;

    /* Igb1 dominates in inversion region, while Igb2 doninates in accumulation */
    /* v2.2.3 bug fix for residue at low Vgb */

    if (Vgb >= 0)
      Igb = Igb1;
    else
      Igb = Igb2;
  }
  else
  {
    Igb = 0.0;
  }

  /* end of gate-body tunneling */
  /* end of v3.0 gate-tunneling  */

  /* v3.1 on page 47 */
  /* Calculate substrate current Iii */
  AD Iii, VgsStep, Vdsatii, Ratio, Vdiff;
  AD Vdsatii0;

  if (alpha0 <= 0.0)
    Iii = 0.0;
  else
  {
    Vdsatii0 = vdsatii0 * (1 + tii * (TempRatio - 1.0)) - lii / Leff;
    /* Calculate VgsStep */
    T0 = esatii * Leff;
    T1 = sii0 * T0 / (1 + T0);
    T0 = 1 / (1 + sii1 * Vgsteff);
    T3 = T0 + sii2;
    T4 = Vgst * sii1 * T0 * T0;
    T2 = Vgst * T3;
    T3 = 1 / (1 + siid * x[0]);
    VgsStep = T1 * T2 * T3;
    Vdsatii = Vdsatii0 + VgsStep;
    Vdiff = x[0] - Vdsatii;

    T0 = beta2 + beta1 * Vdiff + beta0 * Vdiff * Vdiff;
    if (T0 < 1e-5)
      T0 = 1e-5;
    else
      T1 = beta1 + 2 * beta0 * Vdiff;

    if ((T0 < Vdiff / EXPL_THRESHOLD) && (Vdiff > 0.0))
      Ratio = alpha0 * MAX_EXPL;
    else if ((T0 < -Vdiff / EXPL_THRESHOLD) && (Vdiff < 0.0))
      Ratio = alpha0 * MIN_EXPL;
    else
      Ratio = alpha0 * exp(Vdiff / T0);


    /* Avoid too high ratio */
    if (Ratio > 10.0)
      Ratio = 10.0;
    T0 = Ids + fbjtii * Ic; //For the time being put Ic = 0
    Iii = Ratio * T0;
  }

  /***************************/
  /* C-V Model for capMod=2 */

  if (capMod == 2)
  {
    Vfb = Vth - phi - k1eff * sqrtPhis + delvt;
    V3 = Vfb - Vgs_eff + Vbseff - DELTA_3_SOI;
    if (Vfb <= 0.0)
    {
      T0 = sqrt(V3 * V3 - 4.0 * DELTA_3_SOI * Vfb);
      T2 = -DELTA_3_SOI / T0;
    }
    else
    {
      T0 = sqrt(V3 * V3 + 4.0 * DELTA_3_SOI * Vfb);
      T2 = DELTA_3_SOI / T0;
    }

    T1 = 0.5 * (1.0 + V3 / T0);
    Vfbeff = Vfb - 0.5 * (V3 + T0);

    AD CoxWL, CoxWLb;
    AD agbcp = 0.0;
    AD aebcp = 0.0;

    CoxWL  = cox * (weffCV / nseg * leffCV + agbcp);
    CoxWLb = fbody * cox * (weffCV / nseg * leffCVb + agbcp);

    Qac0 = CoxWLb * (Vfbeff - Vfb);
    t0 = 0.5 * K1;
    T3 = Vgs_eff - Vfbeff - Vbseff - Vgsteff;
    if (k1eff == 0.0)
    {
      T1 = 0.0;
      T2 = 0.0;
    }
    else if (T3 < 0.0)
    {
      T1 = t0 + T3 / k1eff;
      T2 = CoxWLb;
    }
    else
    {
      T1 = sqrt(t0 * t0 + T3);
      T2 = CoxWLb * T0 / T1;
    }

    Qsub0 = CoxWLb * K1 * (T1 - t0);
    AbulkCV = Abulk0 * abulkCVfactor;
    VdsatCV = Vgsteff / AbulkCV;
    V4 = VdsatCV - x[0] - DELTA_4;
    T0 = sqrt(V4 * V4 + 4.0 * DELTA_4 * VdsatCV);
    VdseffCV = VdsatCV - 0.5 * (V4 + T0);
    T1 = 0.5 * (1.0 + V4 / T0);
    T2 = DELTA_4 / T0;
    T3 = (1.0 - T1 - T2) / AbulkCV;
    T0 = AbulkCV * VdseffCV;
    T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1e-20);
    T2 = VdseffCV / T1;
    T3 = T0 * T2;
    T4 = (1.0 - 12.0 * T2 * T2 * AbulkCV);
    T5 = (6.0 * T0 * (4.0 * Vgsteff- T0) / (T1 * T1) - 0.5);
    T6 = 12.0 * T2 * T2 * Vgsteff;
    T7 = 1.0 - AbulkCV;
    qbulk = CoxWLb * T7 * (0.5 * VdseffCV - T3);
    T4 = -T7 * (T4 - 1.0);
    T5 = -T7 * T5;
    T6 = -(T7 * T6 + (0.5 * VdseffCV - T3));

    /* Total inversion charge */
    T0 = AbulkCV * VdseffCV;
    T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1e-20);
    T2 = T0 / T1;
    T3 = T0 * T2;

    T4 = (1.0 - 12.0 * T2 * T2);
    T7 = T2 * (2.0 + 6.0 * T2) - 0.5;

    T5 = T7 * AbulkCV;
    T6 = T7 * VdseffCV;

    qinv = CoxWL * (Vgsteff - 0.5 * T0 + T3);

    /* Inversion charge partitioning into S / D */
    if (xpart > 0.5)
    {   /* 0/100 Charge partition model */
      T1 = T1 + T1;
      qsrc = -CoxWL * (0.5 * Vgsteff + 0.25 * T0 - T0 * T0 / T1);
    }
    else if (xpart < 0.5)
    {   /* 40/60 Charge partition model */
      T1 = T1 / 12.0;
      T2 = 0.5 * CoxWL / (T1 * T1);
      T3 = Vgsteff * (2.0 * T0 * T0 / 3.0 + Vgsteff * (Vgsteff - 4.0 * T0 / 3.0))
           - 2.0 * T0 * T0 * T0 / 15.0;
      qsrc = -T2 * T3;
    }
    else
    {
      /* 50/50 Charge partition model */
      qsrc = - 0.5 * (qinv + qbulk);
    }

    /* Backgate charge */
    AD CboxWL = kb1 * fbody * Cbox * (weffCV / nseg * leffCVbg + aebcp);
    AD Qe1 = CboxWL * (Vesfb - 0.0);
    qgate = qinv + Qac0 + Qsub0;
    qbody = (qbulk - Qac0 - Qsub0 - Qe1);
    qsub = Qe1;
    qdrn = -(qgate + qsrc + qbody + qsub);
  } /* End of Capmod ==2 */

  /* C-V Model for capMod=3 */
  if(capMod==3)
  {
    AD dtoxcv = 0.0;
    AD Cox = 3.453133e-11 / (tox - dtoxcv);
    AD agbcp = 0.0;
    AD CoxWL  = cox * (weffCV / nseg *leffCV + agbcp);
    CoxWL = CoxWL * tox / (tox- dtoxcv);
    AD CoxWLb = CoxWLb * tox / (tox - dtoxcv);
    AD  Tox = 1.0e8 * (tox - dtoxcv);
    /* vfbzb calculation for capMod 3 */
    AD k1ox = K1;
    AD noff = 1.0;
    AD voffcv = 0.0;
    T0 = -0.5 * dvt1w * weff * leff / (factor1 * sqrt(Xdep0));
    if (T0 + EXPL_THRESHOLD > 0.0)
      T2 = exp(T0) * (1.0 + 2.0 * exp(t0));
    else
      T2 = MIN_EXPL * (1.0 + 2.0 * MIN_EXPL);
    T0 = dvt0w * T2;
    T2 = T0 * (vbi - phi);
    T0 = -0.5 * dvt1 * leff / (factor1 * sqrt(Xdep0));
    if (T0 + EXPL_THRESHOLD > 0.0)
      T3 = exp(T0) * (1.0 + 2.0 * exp(T0));
    else
      T3 = MIN_EXPL * (1.0 + 2.0 * MIN_EXPL);
    T3 = dvt0 * T3 * (vbi - phi);
    T4 = tox * phi / (weff + w0);
    T0 = sqrt(1.0 + nlx / leff);
    T5 = k1ox * (T0 - 1.0) * sqrtPhi + (kt1 + kt1l / leff) * (Temp/tnom - 1.0);
    T6 = Vth0 - T2 -T3 + k3 * T4 + T5;
    AD vfbzb = T6 - phi - K1 * sqrtPhi;
    AD tmp = sqrt(Xdep0);
    tmp1 = vbi - phi;
    tmp2 = factor1 * tmp;
    T0 = -0.5 * dvt1w * weff * leff / tmp2;
    if (T0 > -EXPL_THRESHOLD)
    {
      T1 = exp(T0);
      T2 = T1 * (1.0 + 2.0 * T1);
    }
    else
    {
      T1 = MIN_EXPL;
      T2 = T1 * (1.0 + 2.0 * T1);
    }
    T0 = dvt0w * T2;
    T2 = T0 * tmp1;

    T0 = -0.5 * dvt1 * leff / tmp2;
    if (T0 > -EXPL_THRESHOLD)
    {
      T1 = exp(T0);
      T3 = T1 * (1.0 + 2.0 * T1);
    }
    else
    {
      T1 = MIN_EXPL;
      T3 = T1 * (1.0 + 2.0 * T1);
    }
    T3 = dvt0 * T3 * tmp1;
    T4 = (tox - dtoxcv) * phi / (weff + w0);
    T0 = sqrt(1.0 + nlx / leff);
    T5 = k1eff * (T0 - 1.0) * sqrtPhi + (kt1 + kt1l / leff) * (TempRatio - 1.0);
    TMP3 = type * Vth0 - T2 - T3 + k3 * T4 + T5;
    vfbzb = TMP3 - phi - k1eff * sqrtPhi;
    /* End of vfbzb */

    //Calculation for VbseffCV
    AD VbseffCV;
    if (Vbseff > 0.0)
      VbseffCV = phi - Phis;
    else
      VbseffCV = Vbseff;
    //Calculation for VgsteffCV
    T0 = n * noff * Vtm0;
    T1 = (Vgs_eff - Vth) / T0;
    if (T1 > EXPL_THRESHOLD)
      Vgsteff = Vgs_eff - Vth - voffcv;
    else
      Vgsteff = T0 * log(1.0 + exp(T1));

    if (-T1 > EXPL_THRESHOLD)
      Vgsteff = T0 * log(1.0 + MIN_EXPL);
    else
      Vgsteff = T0 * log(1.0 + exp(T1));

    //Calculation for Vfbeff
    V3 = vfbzb - Vgs_eff + Vbseff - 0.02;
    if (vfbzb > 0.0)
      T0 = sqrt(V3 * V3 + 4.0 * 0.02 * vfbzb);
    else
      T0 = sqrt(V3 * V3 - 4.0 * 0.02 * vfbzb);
    Vfbeff = vfbzb - 0.5 * (V3 + T0);
    Tox = 1.0e8 * tox;
    T0 = (Vgs_eff - VbseffCV - vfbzb) / Tox;

    //Calculation for Tcen
    ldeb = sqrt(EPSSI * Vtm0/(Charge_q * Nch * 1.0e6)) / 3.0;
    AD acde = 1.0;
    AD Tcen;
    T1 = T0 * acde;
    if ((-EXPL_THRESHOLD < T1) && (T1 < EXPL_THRESHOLD))
      Tcen = ldeb * exp(T1);
    else if (T1 <= -EXPL_THRESHOLD)
      Tcen = ldeb * MIN_EXPL;
    else
      Tcen = ldeb * MAX_EXPL;
    V3 = ldeb - Tcen - (1.0e-3 * Tox);
    V4 = sqrt(V3 * V3 + 4.0 * (1.0e-3 * Tox) * ldeb);

    Tcen = ldeb - 0.5 * (V3 + V4);
    AD Ccen;
    Ccen = EPSSI / Tcen;
    AD Coxeff;
    Coxeff = Ccen * cox / (Ccen + cox);

    //Calculation for QoverlapCox
    AD CoxWLcen = cox * Weff * leff * Coxeff / cox;
    Qac0 = CoxWLcen * (Vfbeff - vfbzb);
    AD QovCox = Qac0 / Coxeff;

    T0 = 0.5 * k1ox;
    T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
    T1 = sqrt(T0 * T0 + T3);
    T2 = CoxWLcen * T0 / T1;

    if (-T3 > 0.0)
      T1 = (T0 + T3 / k1ox);
    else
      T1 = (sqrt(T0 * T0 + T3));

    if (-T3 > 0.0)
      T2 = CoxWLcen;
    else
      T2 = CoxWLcen * T0 / T1;

    Qsub0 = CoxWLcen * k1ox * (T1 - T0);
    QovCox = Qsub0 / Coxeff;

    //Calculation for Delta_phis
    if (k1ox > 0.0)
      t2 = moin * Vtm0 * k1ox * k1ox;
    else
      t2 = 0.25 * moin * Vtm0;

    if (k1ox > 0.0)
      t0 = k1ox * sqrt(phi);
    else
      t0 = 0.5 * sqrt(phi);

    AD DeltaPhi;
    T1 = 2.0 * t0 + Vgsteff;
    DeltaPhi = Vtm0 * log(1.0 + T1 * Vgsteff / t2);

    //The calculation for Tcen must be done once more
    T3 = 4.0 * (Vth - vfbzb - phi);
    Tox = Tox + Tox;
    if (T3 > 0.0)
      T0 = (Vgsteff + T3) / Tox;
    else
      T0 = (Vgsteff + 1.0e-20) / Tox;
    TMP = exp(0.7 * log(T0));
    T1 = 1.0 + TMP;
    T2 = 0.7 * TMP / (T0 * Tox);
    Tcen = 1.9e-9 / T1;
    Ccen = EPSSI / Tcen;
    Coxeff = Ccen * cox / (Ccen + cox);
    CoxWLcen = cox * Weff * leff * Coxeff / cox;
    AbulkCV = Abulk0 * (1.0 + pow((clc/leff),cle));
    VdsatCV = (Vgsteff - DeltaPhi) / AbulkCV;
    T0 = VdsatCV - x[0] - 0.02;
    T1 = sqrt(T0 * T0 + 4.0 * 0.02 * VdsatCV);

    if (T0 > 0.0)
      VdseffCV = VdsatCV - 0.5 * (T0 + T1);
    else
      VdseffCV = VdsatCV * (1.0 - 0.04/(T1-T0));

    if (x[0] == 0.0)
      VdseffCV = 0.0;

    T0 = AbulkCV * VdseffCV;
    T1 = Vgsteff - DeltaPhi;
    T2 = 12.0 * (T1 - 0.5 * T0 + 1.0e-20);
    T3 = T0 / T2;
    T4 = 1.0 - 12.0 * T3 * T3;
    T5 = AbulkCV * (6.0 * T0 * (4.0 * T1 - T0) / (T2 * T2) - 0.5);
    T6 = T5 * VdseffCV / AbulkCV;
    qgate = CoxWLcen * (T1 - T0 * (0.5 - T3));

    qbulk = CoxWLcen * (1.0 - AbulkCV) * (0.5*VdseffCV - T0*VdseffCV/T2);

    QovCox = qbulk / Coxeff;

    T2 = T2 / 12.0;
    T3 = 0.5 * CoxWLcen / (T2 * T2);
    T4 = T1 * (2.0 * T0 * T0 / 3.0 + T1 * (T1 - 4.0 * T0 / 3.0)) - 2.0 * T0 * T0 * T0 / 15.0;

    qsrc = -T3 * T4;

    /* Backgate charge */
    AD aebcp = 0.0;
    AD tbox = 3e-7;

    AD CboxWL = kb1 * fbody * Cbox * (weffCV / nseg * leffCVbg + aebcp);
    AD Qe1 = CboxWL * (Vesfb - 0.0);
    qgate = qgate + Qac0 + Qsub0 - qbulk;
    qbody = qbulk - Qac0 - Qsub0 - Qe1;
    qsub = Qe1;
    qdrn = -(qgate + qbody + qsub + qsrc);
  } /* End of Capmod==3 */

  /*
   *  COMPUTE EQUIVALENT DRAIN CURRENT SOURCE
   */
  if(soiMod < 2)
  {
    /*  Intrinsic S/D charge */
    PhiBSWG = GatesidewallJctPotential;
    MJSWG = bodyJctGateSideGradingCoeff;
    cjsbs = unitLengthGateSidewallJctCap * wdiosCV  * tsi / 1e-7;
    dcjsbs_dT = cjsbs * tcjswg;
    cjsbs += dcjsbs_dT * (Temp - tnom);
    cjdbs = unitLengthGateSidewallJctCap * wdiodCV * tsi / 1e-7;
    dcjdbs_dT = cjdbs * tcjswg;
    cjdbs += dcjdbs_dT * (Temp - tnom);
    AD DioMax = 0.9 * (PhiBSWG);

    AD argtemp=0.0;
    arg = 1.0 - x[2] / PhiBSWG;

    if (MJSWG == 0.5)
      dT3_dVb = 1.0 / sqrt(arg);
    else
      dT3_dVb = exp(-MJSWG * log(arg));
    T3 = (1.0 - arg * dT3_dVb) * PhiBSWG / (1.0 - MJSWG);

    qjs = cjsbs * T3 + tt * Ibsdif;
    arg = 1.0 - ((x[2]-x[0]) > DioMax ? DioMax : (x[2]-x[0])) / PhiBSWG;
    if (MJSWG == 0.5)
      dT3_dVb = 1.0 / sqrt(arg);
    else
      dT3_dVb = exp(-MJSWG * log(arg));
    T3 = (1.0 - arg * dT3_dVb) * PhiBSWG / (1.0 - MJSWG);
    if ((x[2]-x[0]) > DioMax)
      T3 += dT3_dVb * ((x[2]-x[0]) - DioMax);
    qjd = cjdbs * T3 + tt * Ibddif;
  }/* End of Intrinsic S/D charge */

  qdrn -= qjd;
  qbody += (qjs + qjd);
  qsrc = -(qgate + qbody + qdrn + qsub);

  /* Extrinsic Bottom S/D to substrate charge */
  T10 = 0.0;
  /* T10 is vse without type conversion */
  T11 = x[0] - T10;
  /* T11 is vde without type conversion */
  SDphi = 2.0*vtm*log(fabs(nsub) / ni);
  AD tmpp = sqrt(2.0 * EPSSI * SDphi / (Charge_q *  abs(nsub) * 1.0e6));
  AD tmp1p = EPSSI / tmpp;
  //csdmin = tmp1p * Cbox /(tmp1p + Cbox);
  if (csdmin != 0.0)
  {
    if ( ((nsub > 0) && (type > 0)) ||((nsub < 0) && (type < 0)))
    {
      if (Vsdfb > T10)
        qse = csbox * (T10 - Vsdfb);
      else if (T10 < sdt1)
      {
        T0 =  T10 - Vsdfb;
        T1 = T0 * T0;
        qse = T0 * (csbox - st2 / 3 * T1) ;
      }
      else if (T10 < vsdth)
      {
        T0 = T10 - vsdth;
        T1 = T0 * T0;
        qse = csmin * T10 + st4 + st3 / 3 * T0 * T1;
      }
      else
        qse = csmin * T10 + st4;
    }
    else
    {
      if (T10 < vsdth)
        qse = csmin * (T10 - vsdth);

      else if (T10 <sdt1)
      {
        T0 = T10 -vsdth;
        T1 = T0 * T0;
        qse = T0 * (csmin - st2 / 3 * T1) ;
      }
      else if (T10 < Vsdfb)
      {
        T0 = T10 - Vsdfb;
        T1 = T0 * T0;
        qse =csbox * T10 + st4 + st3 / 3 * T0 * T1;
      }
      else
        qse = csbox * T10 + st4;
    }
    if ( ((nsub > 0) && (type > 0)) ||((nsub < 0) && (type < 0)))
    {
      if (Vsdfb > T11)
        qde = cdbox * (T11 - Vsdfb);
      else if (T11 < sdt1)
      {
        T0 =  T10 - Vsdfb;
        T1 = T0 * T0;
        qde = T0 * (csbox - st2 / 3 * T1) ;
      }
      else if (T11 < vsdth)
      {
        T0 = T11 - vsdth;
        T1 = T0 * T0;
        qde = cdmin * T11 + dt4 + dt3 / 3 * T0 * T1;
      }
      else
        qde = cdmin * T11 + dt4;
    }
    else
    {
      if (T11 < vsdth)
        qde = cdmin * (T11 - vsdth);
      else if (T11 <sdt1)
      {
        T0 = T11 -vsdth;
        T1 = T0 * T0;
        qde = T0 * (cdmin - dt2 / 3 * T1) ;
      }
      else if (T11 < Vsdfb)
      {
        T0 = T11 - Vsdfb;
        T1 = T0 * T0;
        qde =cdbox * T11 + dt4 + dt3 / 3 * T0 * T1;
      }
      else
        qde = cdbox * T11 + st4;
    }
  }
  else
  {
    qse = csbox * T10;
    qde = cdbox * T11;
  } /*end of Extrinsic Bottom S/D charge */

  T0 = x[3]-x[2] + DELTA_1;
  T1 = sqrt(T0 * T0 + 4.0 * DELTA_1);
  T2 = 0.5 * (T0 - T1);

  /* v2.2.3 bug fix */
  T3 = wdiodCV * cgdl;

  /* v3.1 bug fix */
  T4 = sqrt(1.0 - 4.0 * T2 /ckappa);
  AD dcgdo = Cgdo + T3 - T3 * (1.0 - 1.0 / T4) * (0.5 - 0.5 * T0 / T1);
  AD qgdo = (dcgdo + T3) * (x[3] - x[2]) - T3 * (T2+ 0.5 * ckappa * (T4 - 1.0));
  T0 = x[3] + DELTA_1;
  T1 = sqrt(T0 * T0 + 4.0 * DELTA_1);
  T2 = 0.5 * (T0 - T1);
  T3 = wdiosCV * cgsl; /* v3.1 bug fix */

  T4 = sqrt(1.0 - 4.0 * T2 / ckappa);
  AD dcgso = Cgso + T3 - T3 * (1.0 - 1.0 / T4)* (0.5 - 0.5 * T0 / T1);
  AD qgso = (dcgso + T3) * x[3] - T3 * (T2 + 0.5 * ckappa * (T4 - 1.0));

  if (soiMod > 0)
  {
    /* v3.1 wanh added for RF */
    if (rgateMod == 3)
    {
      qgd = qgdo;
      qgs = qgso;
      qge = 0; /* v3.1 wanh changed */


      qgate += qge;
      qbody -= 0;
      qdrn += qde - qgd;
      qsub -= qse + qde;
      qsrc = -(qgate  + qbody + qdrn + qsub);
    }
    else
    {
      qgd = qgdo;
      qgs = qgso;
      qgate += qgd + qgs + qge;
      qdrn += qde - qgd;
      qsub -= qge + qse + qde;
      qsrc = -(qgate + qbody + qdrn + qsub);
    }
  }

  // Calculate dynamic current contributions due to charge
  // iqd is dynamic contribution to drain current
  // iqg is dynamic contribution to gate current
  // iqs is dynamic contribution to source current
  AD iqd, iqg, iqs;
  iqd = type * (qdrn.fastAccessDx(0)*x[3] + qdrn.fastAccessDx(1)*x[4] + qdrn.fastAccessDx(2)*x[5]);
  iqg = type * (qgate.fastAccessDx(0)*x[3] + qgate.fastAccessDx(1)*x[4] + qgate.fastAccessDx(2)*x[5]);
  iqs = type * (qsrc.fastAccessDx(0)*x[3] + qsrc.fastAccessDx(1)*x[4] + qsrc.fastAccessDx(2)*x[5]);

  AD Icjd = Ibd - Iii - Idgidl;
  AD Icjs = Ibs - Isgidl;

  flow[0] = Ids + Ic  -Icjd + Igd + iqd;
  flow[1] = Igd + Igs + Igb + iqd;
  flow[2] = -Ids - Ic + -Icjs + Igs - iqs;
  flow[3]= (-Ids - iqd) * (x[0]-x[2]);

  effort[0] = x[0] - x[2]; //Vdb
  effort[1] = x[1] - x[2]; //Vgb
  effort[2] = -x[2]; //Vsb
  effort[3] = x[1] + temp;
}
