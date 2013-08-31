#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>






void Discriminant(){
 
  TFile f("LikelihoodBtag.root");
 
  TFile fC("LikelihoodBtagChargino.root");

  TFile fS("top8M.root"); 

  
  TH1F *tf1;
  f.GetObject("MW Likelihood - RecoD: Inv. Mass t_ - b;6",tf1);
  //std::cout<<tf1->GetMean()<<std::endl;  

  TH1F *tf2;
  fC.GetObject("MW Likelihood - RecoD: Inv. Mass t_ - b;5",tf2);

  TH1D *tf25;
  fS.GetObject("demo/h_SAnal_Inv_Mass_t__bData0;1",tf25);
  tf25->Scale(1.22);
  
  TCanvas *c1 = new TCanvas("c1","Discriminant t_ - b (RecoD)",200,10,600,400);
  //std::cout<<tf1->GetMean()<<std::endl;

  tf1->Draw();
  tf2->SetLineColor(kRed);
  tf2->Draw("histsame"); 

  tf1->SetLineColor(kBlue);
  tf1->Draw("histsame");



  c1->SaveAs("Discriminant_t__b_RD.png");
  c1->Clear();

  TH1F *tf3;
  f.GetObject("MW Likelihood - RecoD: Inv. Mass t - b;6",tf3);
  //std::cout<<tf3->GetMean()<<std::endl;

  TH1F *tf4;
  fC.GetObject("MW Likelihood - RecoD: Inv. Mass t - b;5",tf4);

  TH1D *tf35;
  fS.GetObject("demo/h_SAnal_Inv_Mass_t_bData0;1",tf35);
  tf35->Scale(1.22);


  TCanvas *c1 = new TCanvas("c1","Discriminant t - b (RecoD)",200,10,600,400);
  //std::cout<<tf3->GetMean()<<std::endl;

  tf3->Draw();
  tf4->SetLineColor(kRed);
  tf4->Draw("histsame");

  tf35->SetLineColor(kBlue);
  tf35->Draw("histsame");


  c1->SaveAs("Discriminant_t_b_RD.png");
  c1->Clear();




  TH1F *tf11;
  f.GetObject("MW Likelihood - Reco: Inv. Mass t_ - b;6",tf11);
  //std::cout<<tf1->GetMean()<<std::endl;  

  TH1F *tf21;
  fC.GetObject("MW Likelihood - Reco: Inv. Mass t_ - b;5",tf21);

  
  

  
  TCanvas *c1 = new TCanvas("c1","Discriminant t_ - b (Reco)",200,10,600,400);
  //std::cout<<tf1->GetMean()<<std::endl;

  tf11->Draw();
  tf21->SetLineColor(kRed);
  tf21->Draw("histsame"); 

  c1->SaveAs("Discriminant_t__b_R.png");
  c1->Clear();

  TH1F *tf31;
  f.GetObject("MW Likelihood - Reco: Inv. Mass t - b;6",tf31);
  //std::cout<<tf3->GetMean()<<std::endl;

  TH1F *tf41;
  fC.GetObject("MW Likelihood - Reco: Inv. Mass t - b;5",tf41);

  TCanvas *c1 = new TCanvas("c1","Discriminant t - b (Reco)",200,10,600,400);
  //std::cout<<tf3->GetMean()<<std::endl;

  tf31->Draw();
  tf41->SetLineColor(kRed);
  tf41->Draw("histsame");

  c1->SaveAs("Discriminant_t_b_R.png");
  c1->Clear();





//


  TH1F *tf12;
  f.GetObject("MW Likelihood - Gen: Inv. Mass t_ - b;6",tf12);
  //std::cout<<tf1->GetMean()<<std::endl;  

  TH1F *tf22;
  fC.GetObject("MW Likelihood - Gen: Inv. Mass t_ - b;5",tf22);
  
  TCanvas *c1 = new TCanvas("c1","Discriminant t_ - b (Gen)",200,10,600,400);
  //std::cout<<tf1->GetMean()<<std::endl;

  tf12->Draw();
  tf22->SetLineColor(kRed);
  tf22->Draw("histsame"); 

  c1->SaveAs("Discriminant_t__b_G.png");
  c1->Clear();

  TH1F *tf32;
  f.GetObject("MW Likelihood - Gen: Inv. Mass t - b;6",tf32);
  //std::cout<<tf3->GetMean()<<std::endl;

  TH1F *tf42;
  fC.GetObject("MW Likelihood - Gen: Inv. Mass t - b;5",tf42);

  TCanvas *c1 = new TCanvas("c1","Discriminant t - b (Gen)",200,10,600,400);
  //std::cout<<tf3->GetMean()<<std::endl;

  tf32->Draw();
  tf42->SetLineColor(kRed);
  tf42->Draw("histsame");

  c1->SaveAs("Discriminant_t_b_G.png");
  c1->Clear();

//
















  TH1F *tf5;
  f.GetObject("MW likelihood (best) - RecoD;9",tf5);
 std::cout<<tf5->GetMean()<<std::endl;

 TH1F *tf6;
 fC.GetObject("MW likelihood (best) - RecoD;6",tf6);

 TCanvas *c1 = new TCanvas("c1","Discriminant Likelihood (RecoD)",200,10,600,400);
 std::cout<<tf5->GetMean()<<std::endl;

 tf5->Draw();
 tf6->SetLineColor(kRed);
 tf6->Draw("histsame");

 c1->SaveAs("DiscriminantLikelihood_RD.png");
 c1->Clear();



//

  TH1F *tf51;
  f.GetObject("MW likelihood (best) - Reco;9",tf51);
 //std::cout<<tf5->GetMean()<<std::endl;

 TH1F *tf61;
 fC.GetObject("MW likelihood (best) - Reco;6",tf61);

 TCanvas *c1 = new TCanvas("c1","Discriminant Likelihood (Reco)",200,10,600,400);
 //std::cout<<tf5->GetMean()<<std::endl;

 tf51->Draw();
 tf61->SetLineColor(kRed);
 tf61->Draw("histsame");

 c1->SaveAs("DiscriminantLikelihood_R.png");
 c1->Clear();

//

//

  TH1F *tf52;
  f.GetObject("MW likelihood (best) - Gen;9",tf52);
 //std::cout<<tf5->GetMean()<<std::endl;

 TH1F *tf62;
 fC.GetObject("MW likelihood (best) - Gen;6",tf62);

 TCanvas *c1 = new TCanvas("c1","Discriminant Likelihood (Gen)",200,10,600,400);
 //std::cout<<tf5->GetMean()<<std::endl;

 tf52->Draw();
 tf62->SetLineColor(kRed);
 tf62->Draw("histsame");

 c1->SaveAs("DiscriminantLikelihood_G.png");
 c1->Clear();

//

 TH1F *tf7;
 f.GetObject("MW Likelihood - RecoD: Inv. Mass t_;9",tf7);
 //std::cout<<tf7->GetMean()<<std::endl;

 TH1F *tf8;
 fC.GetObject("MW Likelihood - RecoD: Inv. Mass t_;6",tf8);

 TH1D *tf85;
 fS.GetObject("demo/h_SAnal_Inv_Mass_t__Data0;1",tf85);
 tf85->Scale(1.22);

 TCanvas *c1 = new TCanvas("c1","Discriminant t_ (RecoD)",200,10,600,400);
 //std::cout<<tf7->GetMean()<<std::endl;

 tf7->Draw();
 tf8->SetLineColor(kRed);
 tf8->Draw("histsame");

 tf85->SetLineColor(kBlue);
 tf85->Draw("histsame");

 c1->SaveAs("Discriminant_t__RD.png");
 c1->Clear();

 TH1F *tf9;
 f.GetObject("MW Likelihood - RecoD: Inv. Mass t;9",tf9);
 //std::cout<<tf7->GetMean()<<std::endl;

 TH1F *tf10;
 fC.GetObject("MW Likelihood - RecoD: Inv. Mass t;6",tf10);


 TH1D *tf105;
 fS.GetObject("demo/h_SAnal_Inv_Mass_t_Data0;1",tf105);
 tf105->Scale(1.22);

 TCanvas *c1 = new TCanvas("c1","Discriminant t (RecoD)",200,10,600,400);
 //std::cout<<tf9->GetMean()<<std::endl;

 tf9->Draw();
 tf10->SetLineColor(kRed);
 tf10->Draw("histsame");

 tf105->SetLineColor(kBlue);
 tf105->Draw("histsame");

 c1->SaveAs("Discriminant_t_RD.png");
 c1->Clear();

/////


 TH1F *tf71;
 f.GetObject("MW Likelihood - Reco: Inv. Mass t_;9",tf71);
 //std::cout<<tf7->GetMean()<<std::endl;

 TH1F *tf81;
 fC.GetObject("MW Likelihood - Reco: Inv. Mass t_;6",tf81);

 TCanvas *c1 = new TCanvas("c1","Discriminant t_ (Reco)",200,10,600,400);
 //std::cout<<tf7->GetMean()<<std::endl;

 tf71->Draw();
 tf81->SetLineColor(kRed);
 tf81->Draw("histsame");

 c1->SaveAs("Discriminant_t__R.png");
 c1->Clear();

 TH1F *tf91;
 f.GetObject("MW Likelihood - Reco: Inv. Mass t;9",tf91);
 //std::cout<<tf7->GetMean()<<std::endl;

 TH1F *tf101;
 fC.GetObject("MW Likelihood - Reco: Inv. Mass t;6",tf101);

 TCanvas *c1 = new TCanvas("c1","Discriminant t (Reco)",200,10,600,400);
 //std::cout<<tf9->GetMean()<<std::endl;

 tf91->Draw();
 tf101->SetLineColor(kRed);
 tf101->Draw("histsame");

 c1->SaveAs("Discriminant_t_R.png");
 c1->Clear();

///////J

 TH1F *tf72;
 f.GetObject("MW Likelihood - Gen: Inv. Mass t_;9",tf72);
 //std::cout<<tf7->GetMean()<<std::endl;

 TH1F *tf82;
 fC.GetObject("MW Likelihood - Gen: Inv. Mass t_;6",tf82);

 TCanvas *c1 = new TCanvas("c1","Discriminant t_ (Gen)",200,10,600,400);
 //std::cout<<tf7->GetMean()<<std::endl;

 tf72->Draw();
 tf82->SetLineColor(kRed);
 tf82->Draw("histsame");

 c1->SaveAs("Discriminant_t__G.png");
 c1->Clear();

 TH1F *tf92;
 f.GetObject("MW Likelihood - Gen: Inv. Mass t;9",tf92);
 //std::cout<<tf7->GetMean()<<std::endl;

 TH1F *tf102;
 fC.GetObject("MW Likelihood - Gen: Inv. Mass t;6",tf102);

 TCanvas *c1 = new TCanvas("c1","Discriminant t (Gen)",200,10,600,400);
 //std::cout<<tf9->GetMean()<<std::endl;

 tf92->Draw();
 tf102->SetLineColor(kRed);
 tf102->Draw("histsame");

 c1->SaveAs("Discriminant_t_G.png");
 c1->Clear();





//////
/*
 TH1F *tf11;
 f.GetObject("MW likelihood (best) - RecoD;9",tf11);
 std::cout<<tf11->GetMean()<<std::endl;

 TH1F *tf12;
 fC.GetObject("MW likelihood (best) - RecoD;6",tf12);
 //tf11->Scale(3500000);
 
 TCanvas *c1 = new TCanvas("c1","Discriminant Likelihood",200,10,600,400);
 std::cout<<tf11->GetMean()<<std::endl;

 tf11->Draw();
 tf12->SetLineColor(kRed);
 tf12->Draw("histsame");

 c1->SaveAs("DiscriminantLikelihood.png");
 c1->Clear();

*/

  //  f.Write("",TObject::kOverwrite);
  f.Close();
  exit(0);}

