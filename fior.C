#include <TMinuit.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TRandom3.h>


  float ne[11];
  float nee[11];
  float se[11];
  float ne2[5];
  float ne3[5];
  float nee2[5];
  float nee3[5];
  float se2[5];
  float ne2C[5];

  float s[10]={420.0238821667766,
	       459.34548824885945,
	       495.5012988341309,
	       527.325694462909,
	       556.7898631138137,
	       583.6103422220878,
	       607.7552734479062,
	       634.5241409759935,
	       656.2736709069407,
	       672.1717490765946};
  float z[10]={4.02719984301033, 
	      6.100021883891992,
	      8.077966749779366,
	      10.06461730975971,
	      12.121267074687356,
	      14.156059538013224,
	      16.09802494550964, 
	      18.100922234539482,
	      20.21015449815024, 
	      22.20830008426495};
  float l[10]={22447.861573005546,
	       24837.436279666035,
	       26846.660314655302,
	       28574.70888319828,
	       30016.766113536454,
	       31375.40027303347,
	       32510.711610145972,
	       33712.59708215327,
	       34483.00029022242,
	       35173.38196764031};


  float l2[5]={13515.733614963247,
	       16404.148081868607,
	       18215.143901460928,
	       19567.13027461969,
	       20106.094868989618};

  float l3[5]={10560.016658877503,
	       11446.723935797907,
	       11589.51639104171,
	       11544.376792991643,
	       11611.379350696265};

  float z2[5]={4.02719984301033, 
	       8.077966749779366,
	       12.121267074687356,
	       16.09802494550964, 
	       20.21015449815024};
  
  float s2[5]={420.0238821667766,
	       495.5012988341309,
	       556.7898631138137,
	       607.7552734479062,
	       656.2736709069407};


void FCN5(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {

  // Estrai i parametri dalla lista par
  Double_t p0 = par[0];    // Townsend   
  Double_t c = par[1];    // c                                                                                     
  Double_t p1 = par[2];    // p1                                                                                     
  
  // Calcola il chi-quadrato per il primo istogramma                                                                                        
  Double_t chi2_1 = 0;
  for (int i = 1; i <= 10; i++) {
    Double_t x = s[i];
    Double_t observed = ne[i];
    Double_t Vgem = 440;
    Double_t expected = TMath::Power(c*TMath::Exp(p0*Vgem),3) * TMath::Power(x, 3)/(TMath::Power(x,3)+(p1/Vgem)*TMath::Power(c*TMath::Exp(p0*Vgem),2)*(c*TMath::Exp(p0*Vgem)-1));
    if (expected > 0) {
      chi2_1 += TMath::Power((observed - expected), 2) / expected;
    }
  }
  
  // Calcola il chi-quadrato per il secondo istogramma                                                                                      
  Double_t chi2_2 = 0;
  for (int i = 1; i <= 5; i++) {
    Double_t x = s2[i];
    Double_t observed = ne2[i];
    Double_t Vgem = 430;
    Double_t expected = TMath::Power(c*TMath::Exp(p0*Vgem),3) * TMath::Power(x, 3)/(TMath::Power(x,3)+(p1/Vgem)*TMath::Power(c*TMath::Exp(p0*Vgem),2)*(c*TMath::Exp(p0*Vgem)-1));
    if (expected > 0) {
      chi2_2 += TMath::Power((observed - expected), 2) / expected;
    }
  }
  
  // Calcola il chi-quadrato per il terzo istogramma                         
  Double_t chi2_3 = 0;
  for (int i = 1; i <= 5; i++) {
    Double_t x = s2[i];
    Double_t observed = ne3[i];
    Double_t Vgem = 420;
    Double_t expected = TMath::Power(c*TMath::Exp(p0*Vgem),3) * TMath::Power(x, 3)/(TMath::Power(x,3)+(p1/Vgem)*TMath::Power(c*TMath::Exp(p0*Vgem),2)*(c*TMath::Exp(p0*Vgem)-1));
    if (expected > 0) {
      chi2_3 += TMath::Power((observed - expected), 2) / expected;
    }
  }
  // Somma i chi-quadrati dei tre dataset                                                                                                   
  f= chi2_1 + chi2_2 + chi2_3;
}




void fit(int flag = 5){
    
  for(int i=0; i<10; i++){
    // ne = nc/4 (per passare da conteggi a fotoni) *1/1e-0.03 (omega) *1/0.07(da fotoni a elettroni) 
    ne[i] = (3.50e3 * l[i]);
    // poi dividiamo per il numero di elettroni che arrivano sulla terza GEM 150e4
    ne[i] = ne[i]/168;
    //    ne[i] = ne[i]/1e4;
    //    ne[i] = pow(ne[i],1.0/3);
    nee[i] = 0.03*ne[i];
    se[i] = 10;
  };
  for(int i=0; i<5; i++){
    // ne2 = nc/4 (per passare da conteggi a fotoni) *1/1e-0.03 (omega) *1/0.07(da fotoni a elettroni) 
    ne2[i] = (3.50e3 * l2[i]);
    ne3[i] = (3.50e3 * l3[i]);
    // poi dividiamo per il numero di elettroni che arrivano sulla terza GEM 150e4
    ne2[i] = ne2[i]/168;
    ne3[i] = ne3[i]/168;
    //    ne2[i] = ne2[i]/0.64e4;
    //    ne[i] = pow(ne[i],1.0/3);
    nee2[i] = 0.03*ne2[i];
    nee3[i] = 0.03*ne3[i];
    se2[i] = 10;
  };
  
  TGraph *gz = new TGraph(10,z,l);
  TGraph *gs = new TGraph(10,s,l);
  TGraph *gn = new TGraph(10,s,ne);
  TGraphErrors *gne = new TGraphErrors(10,s,ne,se,nee);
  TGraphErrors *gne2 = new TGraphErrors(5,s2,ne2,se2,nee2);
  TGraphErrors *gne3 = new TGraphErrors(5,s2,ne3,se2,nee2);
  
  gz->SetName("gz");
  gz->SetMarkerStyle(24);
  gz->SetMarkerSize(0.8);
  gs->SetName("gs");
  gs->SetMarkerStyle(24);
  gs->SetMarkerSize(0.8);
  gn->SetName("gn");
  gn->SetMarkerStyle(24);
  gn->SetMarkerSize(0.8);
  gne->SetName("gne");
  gne->SetMarkerStyle(24);
  gne->SetMarkerSize(0.8);
  gne->SetMarkerColor(kRed);
  gne2->SetName("gne2");
  gne2->SetMarkerStyle(26);
  gne2->SetMarkerSize(0.8);
  gne2->SetMarkerColor(kBlue);
  TH2F *h = new TH2F("h","", 100, 300, 800, 100, 0, 1e6);
  TH2F *h1 = new TH2F("h","", 100, 350, 1e3, 100, 0, 1.1);

  if(flag==5){
    // Inizializza `TMinuit` con 2 parametri                                                                                             
    TMinuit minuit(3);
    minuit.SetFCN(FCN5);
    minuit.DefineParameter(0, "p0", 0.02, 0.0001, 0, 1); 
    minuit.DefineParameter(1, "c", 0.01, 0.0001, 0, 1); 
    minuit.DefineParameter(2, "p1", 33e3, 1, 0.1, 1e6);
    
    // Esegui il fit                                                                                                                     
    minuit.Migrad();
    
    Double_t p1, p0, c, p1e, p0e, ce;
    minuit.GetParameter(0, p0, p0e);
    minuit.GetParameter(1, c, ce);
    minuit.GetParameter(2, p1, p1e);
    
    h->Draw();
    gne->Draw("P");
    gne2->Draw("P");
    gne3->Draw("P");
    
    TF1 *fa = new TF1("fa", "TMath::Power([1]*TMath::Exp([0]*[3]),3) * TMath::Power(x, 3)/(TMath::Power(x,3)+([2]/[3])*TMath::Power([1]*TMath::Exp([0]*[3]),2)*([1]*TMath::Exp([0]*[3])-1))",0,1000);
    TF1 *fb = new TF1("fb","TMath::Power([1]*TMath::Exp([0]*[3]),3) * TMath::Power(x, 3)/(TMath::Power(x,3)+([2]/[3])*TMath::Power([1]*TMath::Exp([0]*[3]),2)*([1]*TMath::Exp([0]*[3])-1))",0,1000);
    TF1 *fc = new TF1("fc", "TMath::Power([1]*TMath::Exp([0]*[3]),3) * TMath::Power(x, 3)/(TMath::Power(x,3)+([2]/[3])*TMath::Power([1]*TMath::Exp([0]*[3]),2)*([1]*TMath::Exp([0]*[3])-1))",0,1000);
    TF1 *fd = new TF1("fd", "TMath::Power([1]*TMath::Exp([0]*[3]),3) * TMath::Power(x, 3)/(TMath::Power(x,3)+([2]/[3])*TMath::Power([1]*TMath::Exp([0]*[3]),2)*([1]*TMath::Exp([0]*[3])-1))",0,1000);
    
    fb->SetLineColor(kBlue);
    fb->SetParameter(0, p0);
    fb->SetParameter(1, c);
    fb->SetParameter(2, p1);
    fb->SetParameter(3, 440);
    fb->Draw("same");
    
    fc->SetParameter(0, p0);
    fc->SetParameter(1, c);
    fc->SetParameter(2, p1);
    fc->SetParameter(3, 430);
    fc->SetLineColor(kOrange+1);
    fc->Draw("same");
    
    fd->SetParameter(0, p0);
    fd->SetParameter(1, c);
    fd->SetParameter(2, p1);
    fd->SetParameter(3, 420);
    fd->SetLineColor(kRed);
    fd->Draw("same");
    
  }  
};
