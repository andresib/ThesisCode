#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>





Double_t fitfunc( Double_t* _x, Double_t* _par ){
  Double_t amp1 = _par[0];
  Double_t mean1= _par[1];
  Double_t sigma1 = _par[2];
  Double_t amp2 = _par[3];
  Double_t mean2= _par[4];
  Double_t sigma2 = _par[5];
  Double_t arg1 = (_x[0]-mean1)/sigma1;
  Double_t arg2 = (_x[0]-mean2)/sigma2;
  Double_t res1 = TMath::Exp(-0.5*arg1*arg1);
  Double_t res2 = TMath::Exp(-0.5*arg2*arg2);
  return amp1*res1+amp2*res2;//                                                                                                                                                                
}


void fit(){
 
  TFile f("fit_tf.root","update");

  Double_t amp[6];
  Double_t energy[6];
  const int parameter = 0;

  
  TH1D *tf1;
  f.GetObject("demo/h1_Difference_InEnergy_bJets_Gen_Reco_20;1",tf1);
  mean =tf1->GetMean();
  rms= tf1->GetRMS();
  cout<<mean<<" "<<rms<<endl;
  f1 = new TF1("f1",fitfunc,-20,20,6);
  f1->SetParameter(1, mean-mean/10);
  f1->SetParameter(2, rms);
  f1->SetParameter(4, mean+mean/10);
  f1->SetParameter(5, rms);
  tf1->Fit(f1);
  energy[0]=10;
  amp[0]= f1->GetParameter(0);


  TH1D *tf2;
  f.GetObject("demo/h1_Difference_InEnergy_bJets_Gen_Reco_40;1",tf2);
  mean =tf2->GetMean();
  rms= tf2->GetRMS();
  f2 = new TF1("f2",fitfunc,-20,20,6);
  f2->SetParameter(1, mean-mean/10);
  f2->SetParameter(2, rms);
  f2->SetParameter(4, mean+mean/10);
  f2->SetParameter(5, rms);
  tf2->Fit(f2);
  energy[1]=30;
  amp[1]= f2->GetParameter(0);



  TH1D *tf3;
  f.GetObject("demo/h1_Difference_InEnergy_bJets_Gen_Reco_60;1",tf3);
  mean =tf3->GetMean();
  rms= tf3->GetRMS();
  f3 = new TF1("f3",fitfunc,-20,20,6);
  f3->SetParameter(1, mean-mean/10);
  f3->SetParameter(2, rms);
  f3->SetParameter(4, mean+mean/10);
  f3->SetParameter(5, rms);
  tf3->Fit(f3);
  energy[2]=50;
  amp[2]= f3->GetParameter(0);



  TH1D *tf4;
  f.GetObject("demo/h1_Difference_InEnergy_bJets_Gen_Reco_80;1",tf4);
  mean =tf4->GetMean();
  rms= tf4->GetRMS();
  f4 = new TF1("f4",fitfunc,-20,20,6);
  f4->SetParameter(1, mean-mean/10);
  f4->SetParameter(2, rms);
  f4->SetParameter(4, mean+mean/10);
  f4->SetParameter(5, rms);
  tf4->Fit(f4);
  energy[3]=70;
  amp[3]= f4->GetParameter(0);

  TH1D *tf5;
  f.GetObject("demo/h1_Difference_InEnergy_bJets_Gen_Reco_100;1",tf5);
  mean =tf5->GetMean();
  rms= tf5->GetRMS();
  f5 = new TF1("f5",fitfunc,-20,20,6);
  f5->SetParameter(1, mean-mean/10);
  f5->SetParameter(2, rms);
  f5->SetParameter(4, mean+mean/10);
  f5->SetParameter(5, rms);
  tf5->Fit(f5);
  energy[4]=80;
  amp[4]= f5->GetParameter(0);


  TH1D *tf6;
  f.GetObject("demo/h1_Difference_InEnergy_bJets_Gen_Reco_120;1",tf6);
  mean =tf6->GetMean();
  rms= tf6->GetRMS();
  f6 = new TF1("f6",fitfunc,-20,20,6);
  f6->SetParameter(1, mean-mean/10);
  f6->SetParameter(2, rms);
  f6->SetParameter(4, mean+mean/10);
  f6->SetParameter(5, rms);
  tf6->Fit(f6);
  energy[5]=110;
  amp[5]= f6->GetParameter(0);



 


  TCanvas *c1 = new TCanvas("c1","Amplitud",200,10,600,400);
  
  amplitud = new TGraph(6,energy, amp);
  F1 = new TF1("F1","pol1",-1000,1000);
  amplitud->Fit(F1);
  amplitud->SetLineColor(kRed);
  amplitud->SetMarkerStyle(20);
  amplitud->SetMarkerSize(2.0);
  //hveff->SetMinimum(-0.01);
  //hveff->SetMaximum(110);
  //TAxis *axis = hveff->GetXaxis();
  //axis->SetLimits(8.5,9.9);
  amplitud->SetTitle("Amplitud");
  //hveff->GetXaxis()->SetTitle("HV_Eff(kV)");
  //hveff->GetYaxis()->SetTitle("Efficiency(%)");
  amplitud->Draw("AP");
  c1->SaveAs("amplitud0.png");
  c1->Clear();



  ///////////////////////

 
 
  TH1D *bjets;
  f.GetObject("demo/h1_discriminantbTotal;1",bjets);
  TH1D *cljets;
  f.GetObject("demo/h1_discriminantclTotal;1",cljets);
 

  TF1 *myfitcl = new TF1("myfitcl","pol8", 0 , 1);
  TF1 *myfitb = new TF1("myfitb","pol8", 0 , 1);
 

  cljets->Fit("myfitcl","R");
  bjets->Fit("myfitb","R");
 
  Double_t b,c;

  for (int i=0; i<9;i++)
    {
      b=myfitb->GetParameter(i);
      cout<<"pb"<<i<<": "<<b<<endl;
      cout<<endl;
      c=myfitc->GetParameter(i);
      cout<<"pc"<<i<<": "<<c<<endl;
      cout<<endl;
    }


  






  f.Write("",TObject::kOverwrite);
  f.Close();
  exit(0);
}
