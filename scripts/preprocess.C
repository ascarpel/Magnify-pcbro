/*
Script to merge histograms from six APAs
Wenqiang||||||||||||||||G|u (wgu@bnl.gov)
*/

#include <iostream>
#include <cassert>
#include <string>
#include <regex>
#include "TFile.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TCollection.h"
#include "TKey.h"
#include "TClass.h"
#include "TString.h"
#include "TROOT.h"
#include "TTree.h"

/*
 * Get File Name from a Path with or without extension
 * e.g., std::string filePath = "/home/user/sample/temp/data.csv";
 * std::string name = getFileName(filePath);
 * assert(name == "data.csv");
 */
std::string getFileName(const std::string& fullPath, const std::string extension=".root")
{
  const size_t lastSlashIndex = fullPath.find_last_of("/\\");
  return fullPath.substr(lastSlashIndex + 1, fullPath.size() - lastSlashIndex - extension.size() - 1);
}


double getMedian( std::vector<double> a ) {

  double median =0.0;

  //not the quickest, but the simplest
  std::sort( a.begin(), a.end() );
  size_t n = a.size();

  //even case is simple
  if( n % 2 !=0 ) {
    median =  a[n/2];
  }

  //odd is the mean between the two central values
  median = (a[ (n-1) / 2 ] + a[n/2]) / 2.0;

  return  median;
}


void MergeByTag(TFile* f1, TH2* hall, const char* tag="hu_orig", bool set_baseline=false, double scale=1){

     //TH1I* hbase = new TH1I("baseline","baseline", 4096, 0, 4096);
     std::vector<double> waveform;

     // Get the histogram by key ( there should be only one )
     TH2* h = (TH2*)f1->Get( tag );

    int nX = h->GetNbinsX();
    int nY = h->GetNbinsY();

    for(int i=1; i<=nX; i++) {

      double X = h->GetXaxis()->GetBinCenter(i);

      if(set_baseline) waveform.clear();

      for(int j=1; j<=nY; j++){
        double Y = h->GetYaxis()->GetBinCenter(j);
        double content = h->GetBinContent(i,j);
        if(set_baseline) { waveform.push_back( content ); }
        else{
          int bin1 = hall->FindBin(X,Y);
          hall->SetBinContent(bin1, content);
        }
      }

      if(set_baseline){
        float baseline = getMedian( waveform );
        for(int j=1; j<=nY; j++) {
          double Y = h->GetYaxis()->GetBinCenter(j);
          double content = h->GetBinContent(i,j);
          int bin1 = hall->FindBin(X,Y);
          hall->SetBinContent(bin1, content - baseline);
        }

      }
    }

    std::cout << "Merging " << h->GetName() << " to " << hall->GetName() << std::endl;

    //hbase->Delete();
}

void Merge1DByTag(TFile* f1, TH1* hall, const char* tag="hu_threshold"){

  TH1 *h = (TH1*)f1->Get(tag);

  int nX = h->GetNbinsX();

  for(int i=1; i<=nX; i++) {
    double X = h->GetXaxis()->GetBinCenter(i);
    double content = h->GetBinContent(i);
    int bin1 = hall->FindBin(X);
    hall->SetBinContent(bin1, h->GetBinContent(i));
  }

  std::cout << "Merging " << h->GetName() << " to " << hall->GetName() << std::endl;

}


void getChannelsByTag( TString tag, double &xmin, double &xmax, int &nbinsx ) {


  int nChannels = 64;
  int firstChannel = 64;

  if (tag.Contains("hv")) {
      firstChannel = 128;
  }
  else if (tag.Contains("hw")) {
      firstChannel = 0;
  }

  xmin = firstChannel-0.5;
  xmax = (firstChannel+nChannels)-0.5;
  nbinsx = xmax-xmin;

}

// merge histograms given in- and out-tags
void preprocess(
      std::string inPath = "data/magnify_5141_23468.root",
      std::string outDir = "data/",
      const char* intag= "orig",
      const char* outtag = "orig",
      std::string suffix="v2",
      bool set_baseline=false,
      const char* file_open_mode = "update"
) {

  double xmin=0, xmax=0, ymin=0, ymax=648;
  int nbinsx=0;  int nbinsy = ymax-ymin;

  std::string fileName = getFileName(inPath);
  std::string outPath = outDir + "/" + fileName + "-" + suffix + ".root";
  std::cout << "input file: " << inPath << std::endl;
  std::cout << "output file: " << outPath << std::endl;
  std::cout << "in tag: " << intag << std::endl;
  std::cout << "out tag: " << outtag << std::endl;


  TFile *f1 = TFile::Open(inPath.c_str());
  // Merge trees
  if(std::string(intag).size()>0 &&
     std::string(intag).find("tree:")==0){
    std::string tree_name = std::string(intag).substr(5, std::string(intag).size());
    TTree* treelist[6];
    int ntree=0;
    TFile* fout = new TFile(outPath.c_str(), file_open_mode);
    TList *list = new TList;
    for(int i=0; i<6; i++){
      std::string tree_name1 = tree_name + std::to_string(i);
      treelist[i] = (TTree*)f1->Get(tree_name1.c_str());
      if(treelist[i]){
        ntree ++;
        list->Add(treelist[i]);
      }
    }

    if (ntree>0){
      fout->cd();
      TTree *newtree = TTree::MergeTrees(list);
      std::cout << "\n No. of bad channel regions: "<< newtree->GetEntries() << "\n";
      if(std::string(outtag)=="") outtag = tree_name.c_str();
      newtree->SetName(outtag);
      newtree->SetTitle(outtag);
      newtree->Write();
    }
    else {
      std::cout << "\n No. of bad channels: 0 \n";
    }

    std::cout << "\n Now try: ./magnify.sh " + outPath + "\n";
    fout->Close();
    return;
  }

  // Merge histograms
  auto h = f1->Get(Form("hu_%s", intag));
  //if(h){
  //  std::cout << "No histogram with tag: " << Form("hu_%s", intag) << " found" << std::endl;
  //  return;
  //}

  if(std::string(outtag)=="threshold"){

    getChannelsByTag( Form("hu_%s", outtag), xmin, xmax, nbinsx );
    TH1I* hu = new TH1I(Form("hu_%s", outtag),Form("hu_%s", outtag), nbinsx, xmin, xmax);

    getChannelsByTag( Form("hv_%s", outtag), xmin, xmax, nbinsx );
    TH1I* hv = new TH1I(Form("hv_%s", outtag),Form("hv_%s", outtag), nbinsx, xmin, xmax);

    getChannelsByTag( Form("hw_%s", outtag), xmin, xmax, nbinsx );
    TH1I* hw = new TH1I(Form("hw_%s", outtag),Form("hw_%s", outtag), nbinsx, xmin, xmax);

    TFile* fout = new TFile(outPath.c_str(), file_open_mode);
    Merge1DByTag(f1, hu, Form("hu_%s", intag));
    Merge1DByTag(f1, hv, Form("hv_%s", intag));
    Merge1DByTag(f1, hw, Form("hw_%s", intag));
    fout->cd();
    hu->Write();
    hv->Write();
    hw->Write();
    fout->Close();
  }
  else if (std::string(outtag)=="orig"){

    getChannelsByTag( Form("hu_%s", outtag), xmin, xmax, nbinsx );
    TH2I* hu = new TH2I(Form("hu_%s", outtag),Form("hu_%s", outtag), nbinsx, xmin, xmax, nbinsy, ymin,ymax);

    getChannelsByTag( Form("hv_%s", outtag), xmin, xmax, nbinsx );
    TH2I* hv = new TH2I(Form("hv_%s", outtag),Form("hv_%s", outtag), nbinsx, xmin, xmax, nbinsy, ymin,ymax);

    getChannelsByTag( Form("hw_%s", outtag), xmin, xmax, nbinsx );
    TH2I* hw = new TH2I(Form("hw_%s", outtag),Form("hw_%s", outtag), nbinsx, xmin, xmax, nbinsy, ymin,ymax);

    TFile* fout = new TFile(outPath.c_str(), file_open_mode);
    MergeByTag(f1, hu, Form("hu_%s", intag));
    MergeByTag(f1, hv, Form("hv_%s", intag));
    MergeByTag(f1, hw, Form("hw_%s", intag));
    fout->cd();
    hu->Write();
    hv->Write();
    hw->Write();
    fout->Close();
  }
  else{

    getChannelsByTag( Form("hu_%s", outtag), xmin, xmax, nbinsx );
    TH2F* hu = new TH2F(Form("hu_%s", outtag),Form("hu_%s", outtag), nbinsx, xmin, xmax, nbinsy, ymin,ymax);

    getChannelsByTag( Form("hv_%s", outtag), xmin, xmax, nbinsx );
    TH2F* hv = new TH2F(Form("hv_%s", outtag),Form("hv_%s", outtag), nbinsx, xmin, xmax, nbinsy, ymin,ymax);

    getChannelsByTag( Form("hw_%s", outtag), xmin, xmax, nbinsx );
    TH2F* hw = new TH2F(Form("hw_%s", outtag),Form("hw_%s", outtag), nbinsx, xmin, xmax, nbinsy, ymin,ymax);

    TFile* fout = new TFile(outPath.c_str(), file_open_mode);
    MergeByTag(f1, hu, Form("hu_%s", intag), set_baseline);
    MergeByTag(f1, hv, Form("hv_%s", intag), set_baseline);
    MergeByTag(f1, hw, Form("hw_%s", intag), set_baseline);
    fout->cd();
    hu->Write();
    hv->Write();
    hw->Write();
    fout->Close();
   }

   f1->Close();

   std::cout << "\n Now try: ./magnify.sh " + outPath + "\n";
}
