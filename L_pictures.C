#include <sstream>
#include <string>
#include <fstream>

void sethist(TH1F* hist, int color=1, int rebin=1, int style=1, double scale=1)
{
  hist->SetLineColor(color);
  hist->Rebin(rebin);
  hist->SetLineStyle(style);
  hist->Scale(scale);
}

int L_pictures()
{
  const int fn=8;//number of files in "infile"
  std::ifstream infile("files_list_k.dat");//list of histograms
  std::string file;
  std::string directory="/lustre/nyx/hades/user/iciepal/Lambda1520_ic/";//directory comon for all files
  int n=0;

  
  float scale[]={
    300. ,
    130.*7.8e-5,
  };//ub

  TH1F *hL1520massDistZLpi0[fn];
  TH1F *hL1520massFTDistZLpi0[fn];
  TH1F *hL1520massFinalRLpi0[fn];
  TH1F *hL1520massFTFinalRLpi0[fn];
  TH1F *hL1520massFinalpi0[fn];
  TH1F *hL1520massFTFinalpi0[fn];
  
  TH1F *hL1520massFinalRLpi0_L[fn];
  TH1F *hL1520massDistZLRLpi0_L[fn];

  TH1F *hDLmassFinalRL_L[fn];
  TH1F *hDLmassDistZLRL_L[fn];
  TH1F *hDLmassDistZL[fn];
  TH1F *hDLmassFTDistZL[fn];
  
  //read file wth list of histograms
  if(!infile)
    {
      cout<<"Can't open the file"<<endl;
      return 1;
    }

  //loop over files names and clone histograms*********************************************************
  while (std::getline(infile,file))
    {
      cout<<"file number: "<<n<<" ";
      cout<<file<<endl;
      
      const char* file_name=(directory+file).c_str();
      TFile* hist_file=new TFile(file_name);
      if(hist_file->IsZombie())
	{
	  cout<<"can't open a file!"<<endl<<file_name<<endl;
	  n++;
	  continue;
	}
      hist_file->cd();
      
      hL1520massDistZLpi0[n]= (TH1F*)hist_file->Get("hL1520massDistZLpi0")->Clone();
      hL1520massFTDistZLpi0[n]= (TH1F*)hist_file->Get("hL1520massFTDistZLpi0")->Clone();
      hL1520massFinalRLpi0[n]= (TH1F*)hist_file->Get("hL1520massFinalRLpi0")->Clone();
      hL1520massFTFinalRLpi0[n]= (TH1F*)hist_file->Get("hL1520massFTFinalRLpi0")->Clone();
      hL1520massFinalpi0[n]= (TH1F*)hist_file->Get("hL1520massFinalpi0")->Clone();
      hL1520massFTFinalpi0[n]= (TH1F*)hist_file->Get("hL1520massFTFinalpi0")->Clone();
      
      hL1520massFinalRLpi0_L[n]=(TH1F*)hist_file->Get("hL1520massFinalRLpi0_L")->Clone();
      hL1520massDistZLRLpi0_L[n]=(TH1F*)hist_file->Get("hL1520massDistZLRLpi0_L")->Clone();

      hDLmassFinalRL_L[n]=(TH1F*)hist_file->Get("hDLmassFinalRL_L")->Clone();
      hDLmassDistZLRL_L[n]=(TH1F*)hist_file->Get("hDLmassDistZLRL_L")->Clone();
      hDLmassDistZL[n]= (TH1F*)hist_file->Get("hDLmassDistZL")->Clone();
      hDLmassFTDistZL[n]= (TH1F*)hist_file->Get("hDLmassFTDistZL")->Clone();
	
     n++;
    }
  //end of reading histograms*************************************************************************

  //Print resoults************************************************************************************
  for(int k=0; k<n;k++)
    {
      int bins=5;
      //set colors
          
      sethist(hL1520massFinalpi0[k],k+1,bins,1,scale[k]);
      sethist(hL1520massFTFinalpi0[k],k+1,bins,1,scale[k]);
      sethist(hL1520massFinalRLpi0_L[k],k+1,bins,2,scale[k]);
      sethist(hL1520massDistZLRLpi0_L[k],k+1,bins,2,scale[k]);
      sethist(hL1520massDistZLpi0[k],k+1,bins,1,scale[k]);
      sethist(hL1520massFTDistZLpi0[k],k+1,bins,1,scale[k]);
      
      sethist(hDLmassFinalRL_L[k],k+1,bins,2,scale[k]);
      sethist(hDLmassDistZLRL_L[k],k+1,bins,2,scale[k]);
      sethist(hDLmassDistZL[k],k+1,bins,1,scale[k]);
      sethist(hDLmassFTDistZL[k],k+1,bins,1,scale[k]);
      //sum FW with HADES
      
      hL1520massDistZLpi0[k]->Add(hL1520massFTDistZLpi0[k]);
      hL1520massFinalRLpi0[k]->Add(hL1520massFTFinalRLpi0[k]);
      hL1520massFinalpi0[k]->Add(hL1520massFTFinalpi0[k]);
      hDLmassDistZL[k]->Add(hDLmassFTDistZL[k]);
      
    }
  cout<<"all histograms set"<<endl;
  
  TLegend *legend = new TLegend(0.1,0.2,0.99,0.9);
  legend->AddEntry(hL1520massFinalpi0[1],"p K+ #Lambda(1520)[#Lambda(1115) e+ e-] 130#mub","l");
  legend->AddEntry(hL1520massFinalRLpi0_L[1],"p K+ #Lambda(1520)[#Lambda(1115) e+ e-] 130#mub - true signal","l");
  //legend->AddEntry(hL1520massDist[2],"p K+ #Lambda(1115) #pi^{0}  100#mub","l");
  //legend->AddEntry(hL1520massDist[3],"p K+ #Lambda(1115) 2#pi^{0} 20#mub","l");
  //legend->AddEntry(hL1520massDist[4],"p K+ #Lambda(1115) 3#pi^{0} 7#mub","l");

  //legend->AddEntry(hL1520massDist[0],"p p #pi^{+} #pi^{-} #pi^{0} 1840#mub","l");
  legend->AddEntry(hL1520massFinalpi0[0],"p p #pi^{+} #pi^{-} 2#pi^{0} 300#mub","l");
  // legend->AddEntry(hL1520massDist[7],"p n 2#pi^{+} #pi^{-} #pi^{0} 200#mub","l");
  //legend->AddEntry(hL1520massDist[6],"L1520 decays 130#mub","l");

  TCanvas *cPictures = new TCanvas("cPictures","cPictures");

  cPictures->Divide(2,2);
  double ymin=1e-4; //min value for y axis  
  
  cPictures->cd(1);
  gPad->SetLogy();
  hL1520massFinalRLpi0_L[0]->GetYaxis()->SetRangeUser(ymin,10e4);
  hL1520massFinalpi0[0]->GetYaxis()->SetRangeUser(ymin,10e4);
  for(int x=0;x<n;x++)
    {
      hL1520massFinalpi0[x]->Draw("same");
      hL1520massFinalRLpi0_L[x]->Draw("same");
    }
  
  cPictures->cd(2);
  gPad->SetLogy();
  hL1520massDistZLpi0[0]->GetYaxis()->SetRangeUser(ymin,10e4);
  hL1520massDistZLRLpi0_L[0]->GetYaxis()->SetRangeUser(ymin,10e4);
  for(int x=0;x<n;x++)
    {
      hL1520massDistZLpi0[x]->Draw("same");
      hL1520massDistZLRLpi0_L[x]->Draw("same");
    }

  cPictures->cd(3);
  gPad->SetLogy();
  hDLmassDistZL[0]->GetYaxis()->SetRangeUser(ymin,10e4);
  hDLmassDistZLRL_L[0]->GetYaxis()->SetRangeUser(ymin,10e4);
  for(int x=0;x<n;x++)
    {
      hDLmassDistZL[x]->Draw("same");
      hDLmassDistZLRL_L[x]->Draw("same");
    }

  cPictures->cd(4);
  legend->Draw();
  
  //end of printing results***************************************************************************
  return 1;
}
