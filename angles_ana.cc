#include "fwdet_res.h"

#include "hgeantfwdet.h"
#include "fwdetdef.h"
#include "hfwdetstrawcalsim.h"
#include "hfwdetcand.h"
#include "hfwdetcandsim.h"

#include "hparticlecandsim.h"
#include "hparticlecand.h"
#include "hparticletool.h"

#include "hloop.h"
#include "hcategory.h"

#include <TCanvas.h>
#include <TStyle.h>
#include <sstream>
#include <TGraph.h>
#include <TGraphErrors.h>

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

using namespace std;

Int_t   nentries;

//#include "CreateHistos.cxx"
//#include "WriteHistos.cxx"


    double trackDistance(HParticleCand* track1, HParticleCand*  track2)
    {
      double dist;
      HGeomVector base_1, base_2, dir_1, dir_2;
      HParticleTool p_tool;

      p_tool.calcSegVector(track1->getZ(),track1->getR(),TMath::DegToRad()*track1->getPhi(),TMath::DegToRad()*track1->getTheta(),base_1,dir_1);
      p_tool.calcSegVector(track2->getZ(),track2->getR(),TMath::DegToRad()*track2->getPhi(),TMath::DegToRad()*track2->getTheta(),base_2,dir_2);
      dist=p_tool.calculateMinimumDistance(base_1,dir_1,base_2,dir_2);
      return dist;
    }

    double trackDistance(HParticleCand* track1, HFwDetCand*  track2)
    {
      double dist;
      HGeomVector base_1, base_2, dir_1, dir_2;
      HParticleTool p_tool;

      p_tool.calcSegVector(track1->getZ(),track1->getR(),TMath::DegToRad()*track1->getPhi(),TMath::DegToRad()*track1->getTheta(),base_1,dir_1);

      base_2.setX(track2->getPointX());
      base_2.setY(track2->getPointY());
      base_2.setZ(track2->getPointZ());
      dir_2.setX(track2->getDirTx());
      dir_2.setY(track2->getDirTy());
      dir_2.setZ(1);//konwencja, tak jest ustawione w fwdetstrawvec
  
      dist=p_tool.calculateMinimumDistance(base_1,dir_1,base_2,dir_2);
      return dist;
    }

    HGeomVector trackVertex(HParticleCand* track1, HParticleCand*  track2)
    {
      HGeomVector ver;
      HGeomVector base_1, base_2, dir_1, dir_2;
      HParticleTool p_tool;

      p_tool.calcSegVector(track1->getZ(),track1->getR(),TMath::DegToRad()*track1->getPhi(),TMath::DegToRad()*track1->getTheta(),base_1,dir_1);
      p_tool.calcSegVector(track2->getZ(),track2->getR(),TMath::DegToRad()*track2->getPhi(),TMath::DegToRad()*track2->getTheta(),base_2,dir_2);
      ver=p_tool.calcVertexAnalytical(base_1,dir_1,base_2,dir_2);
      return ver;
    }

    HGeomVector trackVertex(HParticleCand* track1, HFwDetCand*  track2)
    {
      HGeomVector ver;
      HGeomVector base_1, base_2, dir_1, dir_2;
      HParticleTool p_tool;

      p_tool.calcSegVector(track1->getZ(),track1->getR(),TMath::DegToRad()*track1->getPhi(),TMath::DegToRad()*track1->getTheta(),base_1,dir_1);

      base_2.setX(track2->getPointX());
      base_2.setY(track2->getPointY());
      base_2.setZ(track2->getPointZ());
      dir_2.setX(track2->getDirTx());
      dir_2.setY(track2->getDirTy());
      dir_2.setZ(1);//konwencja, tak jest ustawione w fwdetstrawvec
  
      ver=p_tool.calcVertexAnalytical(base_1,dir_1,base_2,dir_2);
      return ver;
    }

    Int_t getMotherIndex(HGeantKine* particle)
    {
      Int_t trackID=particle->getTrack();
      HGeantKine* particleParent=particle->getParent(trackID);
      Int_t parentID=0;
      if(particleParent!=0)
	parentID=particleParent->getID();
      return parentID;
      
      //return particle->getGeneratorInfo1();
    }


    bool isLepton(HParticleCand* particle)
    {
      double delta=0.05;
      double mquality=particle->getRichMatchingQuality();

      double dphi=particle->getDeltaPhi();
      double dtheta=particle->getDeltaTheta();
      
          
      if(particle->isFlagBit(kIsUsed))
	{
	  if(mquality==-1 || mquality>5)
	    return false;
	  if(particle->getBeta()<(1-delta) || particle->getBeta()>(1+delta))
	    return false;
	  // if(dtheta<-0.4 || dtheta>0.4)
	  // return false;
	  // if(dphi*TMath::Sin(particle->getTheta())>0.4 || dphi*TMath::Sin(particle->getTheta())<-0.4)
	  // return false;
	}
      else
	return false;
     
      return true;
    }


//--------------------------------------------------------------------

Int_t fwdet_tests(HLoop * loop, const AnaParameters & anapars)
{
	    

	    if (!loop->setInput(""))
	      {                                                    // reading file structure
		std::cerr << "READBACK: ERROR : cannot read input !" << std::endl;
		std::exit(EXIT_FAILURE);
	      }
	    gStyle->SetOptStat(1);
	    gStyle->SetOptFit(1);

	    TStopwatch timer;
	    timer.Reset();
	    timer.Start();
	    


	    loop->printCategories(); 



    HCategory * fCatGeantKine = nullptr;
    fCatGeantKine = HCategoryManager::getCategory(catGeantKine, kTRUE, "catGeantKine");
    if (!fCatGeantKine)
    {
        cout << "No catGeantKine!" << endl;
        exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

    HCategory * fFwDetStrawCal = nullptr;
    fFwDetStrawCal = HCategoryManager::getCategory(catFwDetStrawCal, kTRUE, "catFwDetStrawCalSim");
    if (!fFwDetStrawCal)
    {
        cout << "No catFwDetStrawCal!" << endl;
        exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

    HCategory * fwDetCatSim = nullptr;
    fwDetCatSim = HCategoryManager::getCategory(catFwDetCand, kTRUE, "catFwDetCand");
    if (!fwDetCatSim)
    {
        cout << "No catFwDetCand!" << endl;
	//exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }

    HCategory * particleCatSim= nullptr;
    particleCatSim = HCategoryManager::getCategory(catParticleCand, kTRUE, "catParticleCand");
    if(!particleCatSim)
      {
	cout<< "No catParticleCandSim!"<<endl;
      }

    Int_t entries = loop->getEntries();
    //     //setting numbers of events regarding the input number of events by the user
    if (anapars.events < entries and anapars.events >= 0 ) entries = anapars.events;



    
    TFile * output_file = TFile::Open(anapars.outfile, "RECREATE");
    output_file->cd();
    
    cout << "NEW ROOT TREE , vertex_ana" << endl;
    Int_t evnb;
    
 //  if (fChain == 0) return;
    
    vector<HParticleCandSim*> lH, pH, pimH;
    vector<HParticleCandSim*> pFT;
 

  for (Long_t event=0; event<nentries; event++) 
    {

      loop->nextEvent(event); 
      if(!(event%100000))   cout<<"event no. "<<event<<endl;


     
      evnb=event;
      //fChain->GetEntry(event);
      //cout<< eventNum<<endl;
      lH.clear();
      pH.clear();  
      pFT.clear();  
      pimH.clear();
    	
    


        	// Geant pairs ana
	int hnum=particleCatSim->getEntries();
	int fnum=fwDetCatSim->getEntries();
	int knum=fCatGeantKine->getEntries();
	HParticleCandSim* ep=nullptr;
	HParticleCandSim* em=nullptr;

	HParticleCandSim* pionH=nullptr;
   
	HParticleCandSim* protonH=nullptr;
	HFwDetCandSim*  protonFT=nullptr;

	HParticleCandSim* partH=nullptr;
	HFwDetCandSim*  partFT=nullptr;


	HGeantKine* kine=nullptr;
	Int_t isDilepton=0;
	Int_t isLambda=0;
	Int_t isLambda1520=0;
	Int_t leptonparticle1=-1;
	Int_t leptonparticle2=-1;
	Int_t lambdaparticle1=-1;
	Int_t lambdaparticle2=-1;
	Int_t leptonparticle=-1;
	HParticleTool tool;
	
	HGeomVector vertexL;
	HGeomVector vertexDL;
	HGeomVector vertexL1520;
	HGeomVector dirL;
	HGeomVector dirDL;

	
	double min_dist_dl=10;
	double min_dist_l=50;
	double min_angle=4;



      
     if (hnum){

       //HADES
	for (int i=0;i<hnum;i++){
	 
	  //float th=lv->Theta()*180./TMath::Pi();;
	  //float ph=lv->Phi()*180./TMath::Pi();;	  
	  //float mom=lv->Rho()*1000;

	  // if(pid==9)thmom9_4pi->Fill(mom,th);
	  //if(pid==14)thmom14_4pi->Fill(mom,th);
	  //if(pid==2)thmom2_4pi->Fill(mom,th);
	  //if(pid==3)thmom3_4pi->Fill(mom,th);
	  //if(pid==9)thmom9_acc->Fill(mom,th,ae);
	  //if(pid==14)thmom14_acc->Fill(mom,th,ae);


	  
	  partH=HCategoryManager::getObject(partH, particleCatSim,i);
	
	  //lep=HCategoryManager::getObject(lep, particleCatSim,i);
	  //pion=HCategoryManager::getObject(pion, particleCatSim,i);
	  //protonFT=HCategoryManager::getObject(protonFT,fwDetCatSim,i);

	  if(partH->getRichMatchingQuality()!=-1){
	    //e+
	    if(partH->getGeantPID()==2) {
	      
	      partH->calc4vectorProperties(HPhysicsConstants::mass(partH->getGeantPID()));
	      
	      ep.push_back(partH);  
 
	      //h2IIleptonsInAcceptance->Fill(partH->getMomentum(),partH->getTheta());
	    }
	    //e-		  
	    if(partH->getGeantPID()==3){
	      em.push_back(partH);  
	    }
	  }


	  if(partH->getGeantPID()==9)// proton->getChi2()<10)
	    {
  
	      pimH.push_back(partH);  
	    }

	  if(partH->getGeantPID()==14)// proton->getChi2()<10)
	    {
	      pH.push_back(partH);  
	    }

	}		 
	 
	//FT
	if (fnum){

	  for (int i=0;i<fnum;i++){
	    partFT=HCategoryManager::getObject(partFT,fwDetCatSim,i);

	    if(partFT->getGeantPID()==14 || partFT->getGeantPID()==9){
	      pFT.push_back(partFT);  
	    }
	 

	  }

	}
	//TLorentzVector lv_lambda1;
	//vector<TLorentzVector> lv_pm;
     
	//**************
	//HADES-HADES
	//*************

	if (pH.size()>=1){
	  for (int k=0;k<pH.size();k++){
	  for (int j=0;j<pimH.size();j++){
	  	    

	      double lambdaM=(*pH+*pimH).M();
	      double lambdaD=trackDistance(ph,pimH);

	      hinvM_pmHpHAll->Fill(lambdaM);
	      
	 
	      if (lambdaM>1105 && lambdaM<1125)
		//hinvM_pmHpHSigw->Fill(lv_lambdaHH.M()*1000,wH);

	     
	     
	   
		/*
	    if (lv22.size() && lv3.size()){
	    //cout<<lv22.size()<<" "<< lv3.size()<<endl;
	  if (lv_lambdaHH.M()*1000>1105 && lv_lambdaHH.M()*1000<1125){  
        	for (int s=0;s<lv22.size();s++){
		  for (int ss=0;ss<lv3.size();ss++){	     


	  //***********************************

	  double invMdilLam=(lv22[s]+lv3[ss]+lv_lambdaHH).M();
	  double invMepem=(lv22[s]+lv3[ss]).M();
		  //cout<<invMdilLam<<endl;
		  //hinvMepem->Fill(invMepem);
	  double oa = lv22[s].Angle(lv3[ss].Vect())*180./TMath::Pi();
		 		 
		hinvM_epemHH->Fill(invMepem,gloW);

		  if(oa>5.){
		    hinvM_LamepemoaHH->Fill(invMdilLam,gloW);
		    //cout<<w1<<endl;

		    hinvM_epemoaHH->Fill(invMepem,gloW);
		    hinvM_epemoawHH->Fill(invMepem,w1*gloW);
	            hinvM_epemoawwHH->Fill(invMepem,w2*gloW);
		
		    hinvM_LamepemoawwHH->Fill(invMdilLam,w2*gloW);
		      
		  if(invMdilLam>1.4 && invMdilLam<1.7)hinvM_epemoawwHHL->Fill(invMepem,w2*gloW);	    
		      

		  }
		  }
		}*/
	  }
	    }
	    
	    
	  }
	}
	//cout<<nbHH<<endl;
    }


  /*
  if (pidH.size()>=1 && pidFT.size()>=1){
      for (int k=0;k<pidH.size();k++){

	for (int j=0;j<pidFT.size();j++){

	  //cout<<event<<" "<<pidH.size()<<" "<<pidFT.size()<<" "<<pidH[k]<<" "<<pidFT[j]<<endl;
	  if (pidFT[j]>-1 && pidH[k]>-1){

	    //cout<<pidFT.size()<<" "<<pidH.size()<<"xxx "<<nbdil2<<" "<<nbdil3<<endl;
	    //cout<<momFT[j].Mag()<<" "<<Tof<<endl;

	    float wH=0;
	    float ac=0, acer=0;
	    float eff=0, effer=0;
	    int pid11=0;
	    
	     
	    float thH=lvH[k].Theta()*180./TMath::Pi();;
	    float phH=lvH[k].Phi()*180./TMath::Pi();;	  
	    float momH1=lvH[k].Rho()*1000;
	    float phHa;	  
	
	    if (phH<0)phHa=360.-fabs(phH);

            //filter->GetAcc(pidH[k],phHa,thH,momH1,ac,acer);
            //if (pidH[k]==9)pid11=8;
            //else pid11=pidH[k];
            filter->GetAcc(pidH[k],phHa,thH,momH1,ac,acer);
   
	    //GetAcc1(pidH[k],phH,thH,momH1,ac,acer);
	    GetEff1(pidH[k],phH,thH,momH1,eff,effer);
	  
	    //cout<<pidH[k]<<" "<<ac<<" "<<eff<<endl;

	    lv_lambda1=lvH[k]+lvFT[j];
	    hinvM_pmFT->Fill(lv_lambda1.M()*1000);  

	    wH=ac*eff;
	   //cout<<pidH[k]<<" "<<ac<<" "<<eff<<" "<<wH<<endl;
	    //if (momH1>2000 && pidH[k]==9)cout<<pidH[k]<<" "<<momH1<<" "<<ac<<" "<<eff<<endl;
	    hinvM_pmFTw->Fill(lv_lambda1.M()*1000,wH); 

	    if (pinFT[j]==5 && pinH[k]==5)hinvM_pmFTsigw->Fill(lv_lambda1.M()*1000,wH);

	    //distance between track in FT and pi- in Hades

	    hVzFT->Fill(vertFT[j].Z());
	    TVector3 newVert=vertexNEW1(momFT[j],momH[k],crosspSTS1[j],vertH[k]);
	    hVz_newVert->Fill(newVert.Z());
	    double dist;
	    dist=TwoTracksMinDist2(momH[k], momFT[j], vertH[k],crosspSTS1[j]);	
   
	    // cout<<"p-pi- min. dist.: "<<dist_ppm<<endl;
	    //cout<<vertFT[j].X()<<" "<<vertFT[j].Y()<<" "<<vertFT[j].Z()<<" "<<endl;	   
	    
	    hDistpmFT->Fill(dist);
	    if (pinH[k]==5 && pinFT[j]==5){
	     //cout<<"--------->> "<<dist<<endl;
	      hDist_Lam->Fill(dist);
	     
 	    }

	    if (pinFT[j]==5 && pinH[k]==5){
	      hVzFT_Lam->Fill(vertFT[j].Z());
	      hVz_newVert_Lam->Fill(newVert.Z());
	    }	    

	    if (pinH[k]==0 && pinFT[j]==0) hVzFT_Primary->Fill(vertFT[j].Z());
	    
	    const float minVectDist=4.;
	    double mInvLam=lv_lambda1.M()*1000;

	    if (lv22.size() && lv3.size()) hinvM_epem0->Fill((lv22[0]+lv3[0]).M(),gloW);
	   
	    if(dist<minVectDist){
	    // if(dist<minVectDist && mInvLam<1125 && mInvLam>1105){
	    
	      hVzFT_Dist->Fill(vertFT[j].Z());
	    
	      
	      hinvM_pmFTdist->Fill(mInvLam); 


	      hinvM_pmFTdistw->Fill(mInvLam,wH); 
	      if (pinFT[j]==5 && pinH[k]==5)hinvM_pmFTdistSigw->Fill(mInvLam,wH);
	      //invMLam=mInvLam;
	      //lvL=lv_lambda1;
	      //flagLam=1;
	      if (vertH[k].Z()>30 && vertFT[j].Z()>30 ){
	     
	       
		hinvM_pmFTdistZ->Fill(mInvLam);
		hinvM_pmFTdistZw->Fill(mInvLam,wH);

		//if (pinFT[j]==5 && pinH[k]==5)hinvM_pmFTdistZsigw->Fill(mInvLam,wH);
		if (mInvLam>1105 && mInvLam<1125 )hinvM_pmFTdistZsigw->Fill(mInvLam,wH);
		
		
		flagLam=1;
		invMLam=mInvLam;
	
		if (lv22.size() && lv3.size()){
		  if(pinH[k]==5 && pinFT[j]==5)nbL++;

		  for (int s=0;s<lv22.size();s++){
		    for (int ss=0;ss<lv3.size();ss++){
		  
       

	  //***********************************
	  double thep=lv22[s].Theta()*180./TMath::Pi();;
	  double phep=lv22[s].Phi()*180./TMath::Pi();;	  
	  double momep=lv22[s].Rho()*1000;
	  double phepa;
	  float acep, aceper;
	  float efep, efeper;

	   if (phep<0)phepa=360.-fabs(phep);
	   filter->GetAcc(2,phepa,thep,momep,acep,aceper);
	   //filter->GetEff(2,phep,thep,momep,efep,efeper);
	   

	   //GetAcc1(2,phep,thep,momep,acep,aceper);
	   GetEff1(2,phep,thep,momep,efep,efeper);

	  //***********************

	  double them=lv3[ss].Theta()*180./TMath::Pi();;
	  double phem=lv3[ss].Phi()*180./TMath::Pi();;	  
	  double momem=lv3[ss].Rho()*1000;
	  double phema;
	  
	  float acem, acemer;
	  float efem, efemer;

	  if (phem<0)phema=360.-fabs(phem);    
          filter->GetAcc(3,phema,them,momem,acem,acemer);
	  //filter->GetEff(3,phem,them,momem,efem,efemer);
	   	 
	    //GetAcc1(3,phem,them,momem,acem,acemer);
	    GetEff1(3,phem,them,momem,efem,efemer);

	  //***********************************
	 
	  double wela=acem;
	  double welef=efem;
	 
	  double wpoa=acep;
	  double wpoef=efep;
	  
	  double w1=wela*welef*wpoa*wpoef;
	  double w2=w1*wH;
	  //cout<<"2::: "<<w1<<" "<<wH<<" "<<w2<<endl;
	  //cout<<wela<<" "<<welef<<" "<<wpoa<<" "<<wpoef<<" "<<w1<<" "<<wH<<endl;

	  
		  
		  double invMepem=(lv22[s]+lv3[ss]).M();
		  double invMdilLam=(lv22[s]+lv3[ss]+lv_lambda1).M();
		  //cout<<invMdilLam<<endl;
		  //hinvMepem->Fill(invMepem);
		  double oa = lv22[s].Angle(lv3[ss].Vect())*180./TMath::Pi();
		  hOAepem->Fill(oa);
		  //hinvMLamDil->Fill(invMdilLam);
		
		  //lvep=lv22;
		  //lvem=lv3;
		  hinvM_Lamepem->Fill(invMdilLam,gloW);
		  hinvM_Lamepemw->Fill(invMdilLam,w2*gloW);
		  hinvM_epem->Fill(invMepem,gloW);
		  hinvM_epemw->Fill(invMepem,w1*gloW);
                  hinvM_epemww->Fill(invMepem,w2*gloW);

		  if(oa>5.){
		    hinvM_Lamepemoa->Fill(invMdilLam,gloW);
		    hinvM_Lamepemoaw->Fill(invMdilLam,w1*gloW);
		    
		    //cout<<w1<<endl;

		    hinvM_epemoa->Fill(invMepem,gloW);
		    hinvM_epemoaw->Fill(invMepem,w1*gloW);
		    hinvM_epemoaww->Fill(invMepem,w2*gloW);

		    if (invMdilLam>1.4&& invMdilLam<1.7){
		    hinvM_epemoawwL->Fill(invMepem,w2*gloW);
                    //if(pidH[k]==72) {
                    //hLam1520->Fill(invMdilLam*1000.); cout<<"!!"<<endl;} 
                      }
		    hinvM_Lamepemoaww->Fill(invMdilLam,w2*gloW);
		      //nbL1520sig+=w2;

		    //hinvMLamDiloa1->Fill(invMdilLam);
		    
		   

		  }
	if(invMLam<1120 && invMLam>1110){
       
	  //hinvMLamDil->Fill(invMdilLam);
	  //if(oa>5.)hinvMLamDiloa->Fill(invMdilLam);
	 
	  //cout<<invMdilLam<<" "<<endl;
	  //cout<<invMLam<<" "<<invMepem<<endl;
	}
      
		  }//s
 }//ss
 
		}



	      }
	      hVXvsVZdist->Fill(vertFT[j].Z(),vertFT[j].X());
	    }
	    
	    
	    
	  }
	  
	}
      }
      */
      
      //double oa;
      /*
      if(nbdil2==1 && nbdil3==1 ){

	double invMepem=(lv22+lv3).M()*1000.;
	double invMdilLam=(lv22+lv3+lvL).M()*1000.;
	hinvMepem->Fill(invMepem);
	double oa = lv22.Angle(lv3.Vect())*180./TMath::Pi();
	hOAepem->Fill(oa);

      if(oa>5. && flagLam)hinvMLamDiloa1->Fill(invMdilLam);

	if(invMLam<1120 && invMLam>1110 && flagLam){
       
	  hinvMLamDil->Fill(invMdilLam);
	  if(oa>5.)hinvMLamDiloa->Fill(invMdilLam);
	  //cout<<invMdilLam<<" "<<endl;
	  //cout<<invMLam<<" "<<invMepem<<endl;
	}
      
     
      }
      */

  //  WriteHistos();

}
      



