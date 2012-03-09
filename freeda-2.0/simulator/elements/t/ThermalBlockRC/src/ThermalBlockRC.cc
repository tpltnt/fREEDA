#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "ThermalBlockRC.h"
#include <cstdio>

//Static members
// Number of netlist paramenters (remember to update!)

const unsigned ThermalBlockRC::n_par = 20;

// Element information
ItemInfo ThermalBlockRC::einfo = {
  "thermalblockrc",
  "Thermal Model for a block of material from layout parameters",
  DEFAULT_ADDRESS"category:thermal",
  "2005_04_07"
};


// Parameter information
// true means the paramenter is required and false means it 
// can be omitted in the netlist.
// See ../network/NetListItem.h for a list of types (TR_INT, TR_DOUBLE, etc.)

ParmInfo ThermalBlockRC::pinfo[] = {
  {"dn", "Density of metal in the North direction (metal/m2)", TR_DOUBLE, false},
  {"ds", "Density of metal in the South direction (metal/m2)", TR_DOUBLE, false},
  {"de", "Density of metal in the East direction (metal/m2)", TR_DOUBLE, false},
  {"dw", "Density of metal in the West direction (metal/m2)", TR_DOUBLE, false},
  {"dt", "Density of metal in the Top direction (metal/m2)", TR_DOUBLE, false},
  {"db", "Density of metal in the Bottom direction (metal/m2)", TR_DOUBLE, false},
  
  {"kmx", "Thermal conductivity of metal in the X direction (W/m.K)", TR_DOUBLE, true},
  {"kmy", "Thermal conductivity of metal in the Y direction (W/m.K)", TR_DOUBLE, true},
  {"kmz", "Thermal conductivity of metal in the Z direction (W/m.K)", TR_DOUBLE, false},
  {"kild","Thermal conductivity of active layer dielectric  (W/m.k)", TR_DOUBLE, false},
  {"kbulk", "Thermal conductivity of bulk material (W/m.K)", TR_DOUBLE, false},
  
  {"lx", "Length of block(metal) in X direction (m)", TR_DOUBLE, false},
  {"ly", "Length of block(metal) in Y direction (m)", TR_DOUBLE, false},
  
  {"habove", "Height of tier above the Silicon island in Z direction (m)", TR_DOUBLE, false},
  {"hbelow", "Height of tier below the Silicon island in Z direction (m)", TR_DOUBLE, false},
  
  {"dmx", "Distance between metal in X direction (m)", TR_DOUBLE, false},
  {"dmy", "Distance between metal in Y direction (m)", TR_DOUBLE, false},
  
  {"rho", "Density of bulk material (kg/m3)", TR_DOUBLE, false},
  {"cbulk", "Heat capacity of bulk material (J/kg.K)", TR_DOUBLE, false},
  {"dfactor", "Vertical Metal density factor ", TR_DOUBLE, false}
 
};

ThermalBlockRC::ThermalBlockRC(const string& iname): Element(&einfo, pinfo, n_par, iname)
{

  //setup parameters

  paramvalue[0] = &(dn=zero);
  paramvalue[1] = &(ds=zero);
  paramvalue[2] = &(de=zero);
  paramvalue[3] = &(dw=zero);
  paramvalue[4] = &(dt=zero);
  paramvalue[5] = &(db=zero);
  
  paramvalue[6] = &(kmx=zero);
  paramvalue[7] = &(kmy=zero);
  paramvalue[8] = &(kmz=zero);
  paramvalue[9] = &(kild=zero);
  paramvalue[10] = &(kbulk=zero);
  
  paramvalue[11] = &(lx=zero);
  paramvalue[12] = &(ly=zero);
  
  paramvalue[13] = &(habove=zero);
  paramvalue[14] = &(hbelow=zero);
  
  paramvalue[15]= &(dmx=zero);
  paramvalue[16] = &(dmy=zero);
  paramvalue[17] = &(rho=zero);
  paramvalue[18] = &(cbulk=zero);
  paramvalue[19] = &(dfactor=one);	
setNumTerms(7);

setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN);
};


void ThermalBlockRC::fillMNAM(FreqMNAM* mnam)
{
  double rn = zero,rs =zero , re =zero, rw =zero, rt = zero, rb =zero;
  double rn_die, rn_metal, rs_die, rs_metal, re_die, re_metal, rw_die, rw_metal, rt_die, rt_metal;
  double rb_die, rb_metal;
  double arean_die, arean_metal , areas_die, areas_metal,  areae_die, areae_metal, areaw_metal, areaw_die;
  double areat_metal, areat_die, areab_metal, areab_die;
  double cap = zero, int_g=zero;
  double kxeff = zero, kyeff = zero, kzeff = zero;
  string  myname;  
  char *cname;
  char rmessage[1024];
  int len=0, col=0, row=0, tier=0;
  FILE* fp1=NULL;
  cname = (char *) malloc(sizeof(char)*256);
  rn_die = rn_metal = rs_die = rs_metal = zero;
  re_die = re_metal = rw_die = rw_metal = rt_die = rt_metal = zero;
  rb_die = rb_metal = zero;
  
  dt = dt * dfactor;
  db = db * dfactor;
  
  arean_die = arean_metal = areas_die = areas_metal = zero;
  areae_die = areae_metal = areaw_metal = areaw_die = zero;
  areat_metal= areat_die = areab_metal = areab_die = zero;
  
  arean_metal = ( dn/100.0 ) * ly * (habove + hbelow);
  arean_die   = ((100.0-dn)/100) * ly * (habove + hbelow);
  
  areas_metal = ( ds/100.0 ) * ly * (habove + hbelow);
  areas_die   = ((100.0-ds)/100.0) * ly * (habove + hbelow);
  
  areae_metal = ( de/100.0 ) * lx * (habove + hbelow);
  areae_die   = ((100.0-de)/100.0) * lx * (habove + hbelow);
  
  areaw_metal = ( dw/100.0 ) * lx * (habove + hbelow);
  areaw_die   = ((100.0-dw)/100.0) * lx * (habove + hbelow);
  
  areat_metal = ( dt/100.0 ) * lx * ly;
  areat_die   = ((100-dt)/100.0) * lx * ly;
  
  areab_metal = ( db/100.0 ) * lx * ly;
  areab_die   = ((100.0-db)/100.0) * lx * ly;
  
  if(arean_metal)
  	rn_metal = lx / (2.0 * kmx * arean_metal);
  rn_die   = lx / (2.0 * kbulk * arean_die);
  
  if(areas_metal)
  	rs_metal = lx / (2.0 * kmx * areas_metal);
  rs_die   = lx / (2.0 * kbulk * areas_die);

  if(areae_metal)
	  re_metal = ly / (2.0 * kmy * areae_metal);
  re_die   = ly / (2.0 * kbulk * areae_die);
  
  if(areaw_metal)
	  rw_metal = ly / (2.0 * kmy * areaw_metal);
  rw_die   = ly / (2.0 * kbulk * areaw_die);
  
  if(areat_metal)
  	rt_metal = habove / ( kmz   * areat_metal);
  rt_die   = habove / ( kild * areat_die);
   
  if(areab_metal) 
  	rb_metal = hbelow / ( kmz * areab_metal);
  rb_die   = hbelow / (kild * areab_die);

  if(re_metal)
  	re = (re_die * re_metal) / (re_metal + re_die); // parallel resistance
  else
  	re = re_die;
  if(rw_metal)	
  	rw = (rw_die * rw_metal) / (rw_metal + rw_die); // parallel resistance
  else
  	rw = rw_die;
  if(rn_metal)
  	rn = (rn_die * rn_metal) / (rn_metal + rn_die); // parallel resistance
  else
  	rn = rn_die;
  if(rs_metal)
  	rs = (rs_die * rs_metal) / (rs_metal + rs_die); // parallel resistance
  else
  	rs= rs_die;
  if(rt_metal)
  	rt = (rt_die * rt_metal) / (rt_metal + rt_die); // parallel resistance
  else
  	rt= rt_die;
  if(rb_metal)
  	rb = (rb_die * rb_metal) / (rb_metal + rb_die); // parallel resistance
  else
  	rb = rb_die;
  	



  cap = rho*cbulk*hbelow*lx*ly;
  double_complex y(int_g, twopi*mnam->getFreq()*cap);
  myname = getInstanceName();
  
  
  cname = (char *)(myname.c_str());
  len = strlen(cname);
  col = cname[len-1] - 48;
  row = cname[len-2] - 48;
  tier = cname[len-3] - 48;
  
      
  kxeff = lx / (rn+rs);
  kxeff = kxeff / (ly*(habove+hbelow));
  
  kyeff = ly / (re+rw) ;
  kyeff = kyeff /(lx*(habove+hbelow));
  
  kzeff = (habove + hbelow) / (rt+rb) ;
  kzeff = kzeff / (lx*ly);
  
 fp1 = fopen("Thermal_values.txt","a");
 sprintf(rmessage,"%d\t%d\t%d\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\n",tier,row,col,kxeff,kyeff,kzeff,rn,rs,re,rw,rt,rb);
 fputs(rmessage,fp1);
 fclose(fp1);
 // sprintf(rmessage,"%0.4f\t%0.4f\t%0.4f\n",(rn+rs),(re+rw),(rt+rb));
 // fputs(rmessage,fp1);
  
//    cout<<" Re = "<<re<<" Rw = "<<rw ;
//    cout<<" Rn = "<<rn<<" Rs = " <<rs;
//    cout<<" Rt = "<<rt<<" Rb = "<<rb<<endl;
 /* cout<<"Capacitor ="<<cap<<endl;*/
  
  //cout<<tier<<"\t"<<row<<"\t"<<col<<"\t"<<(rn+rs)<<"\t"<<(re+rw)<<"\t"<<(rt+rb)<<endl;
  //cout<<tier<<"\t"<<row<<"\t"<<col<<"\t"<<kxeff<<"\t"<<kyeff<<"\t"<<kzeff<<endl;
  
  
  
  //Use superposition model

  mnam->setAdmittance(getTerminal(0)->getRC(), getTerminal(6)->getRC(),1.0/rn);
  mnam->setAdmittance(getTerminal(1)->getRC(), getTerminal(6)->getRC(),1.0/rs);
  mnam->setAdmittance(getTerminal(2)->getRC(), getTerminal(6)->getRC(),1.0/re);
  mnam->setAdmittance(getTerminal(3)->getRC(), getTerminal(6)->getRC(),1.0/rw);
  mnam->setAdmittance(getTerminal(4)->getRC(), getTerminal(6)->getRC(),1.0/rt);
  mnam->setAdmittance(getTerminal(5)->getRC(), getTerminal(6)->getRC(),1.0/rb);

 

  /*  mnam->setAdmittance(getTerminal(0)->getRC(), getTerminal(4)->getRC(),2.0/rx);
  mnam->setAdmittance(getTerminal(1)->getRC(), getTerminal(4)->getRC(),2.0/rx);
  mnam->setAdmittance(getTerminal(2)->getRC(), getTerminal(4)->getRC(),2.0/ry);
  mnam->setAdmittance(getTerminal(3)->getRC(), getTerminal(4)->getRC(),2.0/ry);
  mnam->setAdmittance(getTerminal(4)->getRC(), getTerminal(5)->getRC(),1/(rz_v+rz_x+rz_y));
  mnam->setAdmittance(getTerminal(5)->getRC(), getTerminal(6)->getRC(),2.0/rb);
  mnam->setAdmittance(getTerminal(5)->getRC(), getTerminal(7)->getRC(),2.0/rb); */
}

  

void ThermalBlockRC::fillMNAM(TimeMNAM* mnam)
{
  double rn = zero,rs =zero , re =zero, rw =zero, rt = zero, rb =zero;
  double rn_die, rn_metal, rs_die, rs_metal, re_die, re_metal, rw_die, rw_metal, rt_die, rt_metal;
  double rb_die, rb_metal;
  double arean_die, arean_metal , areas_die, areas_metal,  areae_die, areae_metal, areaw_metal, areaw_die;
  double areat_metal, areat_die, areab_metal, areab_die;
  double cap = zero, int_g=zero;
  double kxeff = zero, kyeff = zero, kzeff = zero;
  string  myname;  
  char *cname;
  char rmessage[1024];
  int len=0, col=0, row=0, tier=0;
  FILE* fp1=NULL;
  cname = (char *) malloc(sizeof(char)*256);
  rn_die = rn_metal = rs_die = rs_metal = zero;
  re_die = re_metal = rw_die = rw_metal = rt_die = rt_metal = zero;
  rb_die = rb_metal = zero;
  
  dt = dt * dfactor;
  db = db * dfactor;
  
  arean_die = arean_metal = areas_die = areas_metal = zero;
  areae_die = areae_metal = areaw_metal = areaw_die = zero;
  areat_metal= areat_die = areab_metal = areab_die = zero;
  
  arean_metal = ( dn/100.0 ) * ly * (habove + hbelow);
  arean_die   = ((100.0-dn)/100) * ly * (habove + hbelow);
  
  areas_metal = ( ds/100.0 ) * ly * (habove + hbelow);
  areas_die   = ((100.0-ds)/100.0) * ly * (habove + hbelow);
  
  areae_metal = ( de/100.0 ) * lx * (habove + hbelow);
  areae_die   = ((100.0-de)/100.0) * lx * (habove + hbelow);
  
  areaw_metal = ( dw/100.0 ) * lx * (habove + hbelow);
  areaw_die   = ((100.0-dw)/100.0) * lx * (habove + hbelow);
  
  areat_metal = ( dt/100.0 ) * lx * ly;
  areat_die   = ((100-dt)/100.0) * lx * ly;
  
  areab_metal = ( db/100.0 ) * lx * ly;
  areab_die   = ((100.0-db)/100.0) * lx * ly;
  
  if(arean_metal)
  	rn_metal = lx / (2.0 * kmx * arean_metal);
  rn_die   = lx / (2.0 * kbulk * arean_die);
  
  if(areas_metal)
  	rs_metal = lx / (2.0 * kmx * areas_metal);
  rs_die   = lx / (2.0 * kbulk * areas_die);

  if(areae_metal)
	  re_metal = ly / (2.0 * kmy * areae_metal);
  re_die   = ly / (2.0 * kbulk * areae_die);
  
  if(areaw_metal)
	  rw_metal = ly / (2.0 * kmy * areaw_metal);
  rw_die   = ly / (2.0 * kbulk * areaw_die);
  
  if(areat_metal)
  	rt_metal = habove / ( kmz   * areat_metal);
  rt_die   = habove / ( kild * areat_die);
   
  if(areab_metal) 
  	rb_metal = hbelow / ( kmz * areab_metal);
  rb_die   = hbelow / (kild * areab_die);

  if(re_metal)
  	re = (re_die * re_metal) / (re_metal + re_die); // parallel resistance
  else
  	re = re_die;
  if(rw_metal)	
  	rw = (rw_die * rw_metal) / (rw_metal + rw_die); // parallel resistance
  else
  	rw = rw_die;
  if(rn_metal)
  	rn = (rn_die * rn_metal) / (rn_metal + rn_die); // parallel resistance
  else
  	rn = rn_die;
  if(rs_metal)
  	rs = (rs_die * rs_metal) / (rs_metal + rs_die); // parallel resistance
  else
  	rs= rs_die;
  if(rt_metal)
  	rt = (rt_die * rt_metal) / (rt_metal + rt_die); // parallel resistance
  else
  	rt= rt_die;
  if(rb_metal)
  	rb = (rb_die * rb_metal) / (rb_metal + rb_die); // parallel resistance
  else
  	rb = rb_die;
  	



  cap = rho*cbulk*hbelow*lx*ly;
cap=10e-12;
 // double_complex y(int_g, twopi*mnam->getFreq()*cap);
  myname = getInstanceName();
  
  
  cname = (char *)(myname.c_str());
  len = strlen(cname);
  col = cname[len-1] - 48;
  row = cname[len-2] - 48;
  tier = cname[len-3] - 48;
  
      
  kxeff = lx / (rn+rs);
  kxeff = kxeff / (ly*(habove+hbelow));
  
  kyeff = ly / (re+rw) ;
  kyeff = kyeff /(lx*(habove+hbelow));
  
  kzeff = (habove + hbelow) / (rt+rb) ;
  kzeff = kzeff / (lx*ly);
  
 fp1 = fopen("Thermal_values.txt","a");
 sprintf(rmessage,"%d\t%d\t%d\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.8f\n",tier,row,col,kxeff,kyeff,kzeff,rn,rs,re,rw,rt,rb,cap);
 fputs(rmessage,fp1);
 fclose(fp1);
 // sprintf(rmessage,"%0.4f\t%0.4f\t%0.4f\n",(rn+rs),(re+rw),(rt+rb));
 // fputs(rmessage,fp1);
  
//    cout<<" Re = "<<re<<" Rw = "<<rw ;
//    cout<<" Rn = "<<rn<<" Rs = " <<rs;
//    cout<<" Rt = "<<rt<<" Rb = "<<rb<<endl;
 /* cout<<"Capacitor ="<<cap<<endl;*/
  
  //cout<<tier<<"\t"<<row<<"\t"<<col<<"\t"<<(rn+rs)<<"\t"<<(re+rw)<<"\t"<<(rt+rb)<<endl;
  //cout<<tier<<"\t"<<row<<"\t"<<col<<"\t"<<kxeff<<"\t"<<kyeff<<"\t"<<kzeff<<endl;
  
  
  
  //Use superposition model

  mnam->setMAdmittance(getTerminal(0)->getRC(), getTerminal(6)->getRC(),1.0/rn);
  mnam->setMAdmittance(getTerminal(1)->getRC(), getTerminal(6)->getRC(),1.0/rs);
  mnam->setMAdmittance(getTerminal(2)->getRC(), getTerminal(6)->getRC(),1.0/re);
  mnam->setMAdmittance(getTerminal(3)->getRC(), getTerminal(6)->getRC(),1.0/rw);
  mnam->setMAdmittance(getTerminal(4)->getRC(), getTerminal(6)->getRC(),1.0/rt);
  mnam->setMAdmittance(getTerminal(5)->getRC(), getTerminal(6)->getRC(),1.0/rb);
  mnam->setMpAdmittance(getTerminal(5)->getRC(), getTerminal(6)->getRC(),cap);

 

  /*  mnam->setAdmittance(getTerminal(0)->getRC(), getTerminal(4)->getRC(),2.0/rx);
  mnam->setAdmittance(getTerminal(1)->getRC(), getTerminal(4)->getRC(),2.0/rx);
  mnam->setAdmittance(getTerminal(2)->getRC(), getTerminal(4)->getRC(),2.0/ry);
  mnam->setAdmittance(getTerminal(3)->getRC(), getTerminal(4)->getRC(),2.0/ry);
  mnam->setAdmittance(getTerminal(4)->getRC(), getTerminal(5)->getRC(),1/(rz_v+rz_x+rz_y));
  mnam->setAdmittance(getTerminal(5)->getRC(), getTerminal(6)->getRC(),2.0/rb);
  mnam->setAdmittance(getTerminal(5)->getRC(), getTerminal(7)->getRC(),2.0/rb); */
}

  


/*void ThermalBlockRC::fillMNAM(TimeMNAM* mnam)
{
 
 double rn,rs, re, rw, rt, rb;
   double rn_die, rn_metal, rs_die, rs_metal, re_die, re_metal, rw_die, rw_metal, rt_die, rt_metal;
   double rb_die, rb_metal;
   double arean_die, arean_metal , areas_die, areas_metal,  areae_die, areae_metal, areaw_metal, areaw_die;
   double areat_metal, areat_die, areab_metal, areab_die;
   double cap, int_g=zero;
  
   dt = dt * dfactor;
   db = db * dfactor;
  
  areae_metal = ( de/100.0 ) * ly * (habove + hbelow);
   areae_die   = ((100.0-de)/100) * ly * (habove + hbelow);
   
   areaw_metal = ( dw/100.0 ) * ly * (habove + hbelow);
   areaw_die   = ((100.0-dw)/100.0) * ly * (habove + hbelow);
   
   arean_metal = ( dn/100.0 ) * lx * (habove + hbelow);
   arean_die   = ((100.0-dn)/100.0) * lx * (habove + hbelow);
   
   areas_metal = ( ds/100.0 ) * lx * (habove + hbelow);
   areas_die   = ((100.0-ds)/100.0) * lx * (habove + hbelow);
   
   areat_metal = ( dt/100.0 ) * lx * ly;
   areat_die   = ((100-dt)/100.0) * lx * ly;
   
   areab_metal = ( db/100.0 ) * lx * ly;
   areab_die   = ((100.0-db)/100.0) * lx * ly;
   
 
   re_metal = lx / (2.0 * kmx * areae_metal);
   re_die   = lx / (2.0 * kbulk * areae_die);
   
   rw_metal = lx / (2.0 * kmx * areaw_metal);
   rw_die   = lx / (2.0 * kbulk * areaw_die);
 
   rn_metal = ly / (2.0 * kmy * arean_metal);
   rn_die   = ly / (2.0 * kbulk * arean_die);
   
   rs_metal = ly / (2.0 * kmy * areas_metal);
   rs_die   = ly / (2.0 * kbulk * areas_die);
   
   rt_metal = habove / ( kmz   * areat_metal);
   rt_die   = habove / ( kbulk * areat_die);
     
   rb_metal = hbelow / (kmz * areab_metal);
   rb_die   = hbelow / (kild * areab_die);
 
   
   re = (re_die * re_metal) / (re_metal + re_die); // parallel resistance
   rw = (rw_die * rw_metal) / (rw_metal + rw_die); // parallel resistance
   rn = (rn_die * rn_metal) / (rn_metal + rn_die); // parallel resistance
   rs = (rs_die * rs_metal) / (rs_metal + rs_die); // parallel resistance
   rt = (rt_die * rt_metal) / (rt_metal + rt_die); // parallel resistance
   rb = (rb_die * rb_metal) / (rb_metal + rb_die); // parallel resistance
 
 
 
   cap = rho*cbulk*hbelow*lx*ly;
   //double_complex y(int_g, twopi*mnam->getFreq()*cap);
 
   cout<<"Re = "<<re<<" Rw = "<<rw ;
   cout<<" Rn = "<<rn<<" Rs = " <<rs;
   cout<<" Rt = "<<rt<<" Rb = "<<rb<<endl;
 
   cout<<"Capacitor ="<<cap<<endl;
   //Use superposition model
 
   mnam->setMAdmittance(getTerminal(0)->getRC(), getTerminal(6)->getRC(),1.0/rn);
   mnam->setMAdmittance(getTerminal(1)->getRC(), getTerminal(6)->getRC(),1.0/rs);
   mnam->setMAdmittance(getTerminal(2)->getRC(), getTerminal(6)->getRC(),1.0/re);
   mnam->setMAdmittance(getTerminal(3)->getRC(), getTerminal(6)->getRC(),1.0/rw);
   mnam->setMAdmittance(getTerminal(4)->getRC(), getTerminal(6)->getRC(),1.0/rt);
   mnam->setMAdmittance(getTerminal(5)->getRC(), getTerminal(6)->getRC(),1.0/rb);
   mnam->setMpAdmittance(getTerminal(5)->getRC(), getTerminal(6)->getRC(),cap);
 

    mnam->setAdmittance(getTerminal(0)->getRC(), getTerminal(4)->getRC(),2.0/rx);
  mnam->setAdmittance(getTerminal(1)->getRC(), getTerminal(4)->getRC(),2.0/rx);
  mnam->setAdmittance(getTerminal(2)->getRC(), getTerminal(4)->getRC(),2.0/ry);
  mnam->setAdmittance(getTerminal(3)->getRC(), getTerminal(4)->getRC(),2.0/ry);
  mnam->setAdmittance(getTerminal(4)->getRC(), getTerminal(5)->getRC(),1/(rz_v+rz_x+rz_y));
  mnam->setAdmittance(getTerminal(5)->getRC(), getTerminal(6)->getRC(),2.0/rb);
  mnam->setAdmittance(getTerminal(5)->getRC(), getTerminal(7)->getRC(),2.0/rb); 
  

}
*/
