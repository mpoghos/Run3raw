// test1

#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "RStringView.h"
#include <Rtypes.h>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>

#include "DetectorsRaw/RDHUtils.h"
#include "DetectorsRaw/RawFileReader.h"
#include "EMCALBase/Mapper.h"
#include "Framework/Logger.h"

#include "EMCALReconstruction/AltroDecoder.h"
#include "EMCALReconstruction/Bunch.h"
#include "EMCALReconstruction/CaloFitResults.h"
#include "EMCALReconstruction/CaloRawFitterStandard.h"

#include "EMCALReconstruction/AltroDecoder.h"
#include "EMCALReconstruction/Bunch.h"
#include "EMCALReconstruction/CaloFitResults.h"
#include "EMCALReconstruction/CaloRawFitterGamma2.h"
#include "EMCALReconstruction/CaloRawFitterStandard.h"
//#include "EMCALReconstruction/CaloRawFitterGamma2v2.h"

#endif

#include "Clusterizer3x3.h"

using namespace o2::emcal;

const int debug = 2;
const Bool_t kDoFit = kTRUE; // kFALSE;

bool IsGoodChannel(int col, int row);
Clusterizer3x3 clmaker;

void AnalData_TEDwithClusterizer() {
  const Int_t NoiseThreshold = 3;
  const float cellAmpMin = 5.;
  const float cellAmpMax = 800.;

  //  std::unique_ptr<Clusterizer3x3> clmaker =
  //      std::unique_ptr<Clusterizer3x3>(new Clusterizer3x3);
  clmaker = Clusterizer3x3();
  clmaker.Init();
  clmaker.setDebugLevel(debug);

  const char *aliceO2env = std::getenv("O2_ROOT");
  std::string inputDir = " ";
  o2::raw::RawFileReader reader;
  reader.setDefaultDataOrigin(o2::header::gDataOriginEMC);
  reader.setDefaultDataDescription(o2::header::gDataDescriptionRawData);
  reader.setDefaultReadoutCardType(o2::raw::RawFileReader::RORC);
  // TED_20211001
  // reader.addFile("data/readoutSMA3_2021_10_01_13_52_11.raw");
  // 1- injection
  //  reader.addFile("data/readoutSMA3_2021_10_01_17_08_00.raw");
  //  std::string outFileName = "injData";

  // reader.addFile("data/readoutSMA3_2021_10_01_22_29_17.raw"); // no injection

  // 2- no_injection_L3on
  // reader.addFile("data/readoutSMA3_2021_10_01_23_16_35.raw"); //
  // std::string outFileName = "noInjection";

  // 3- no_injection_L3on_rnd10hz
  //  reader.addFile("data/readoutSMA3_2021_10_02_09_53_51.raw");
  //  std::string outFileName = "rnd10hz";

  // 4 - no injection, L3 on , rnd10hz trig, L3 will be switched on in 15mins
  //  reader.addFile("data/readoutSMA3_2021_10_04_06_51_53.raw");
  //  std::string outFileName = "noInjection_L3on_10hz_2";

  // 5 - no injection, L3 off, rnd10hz trig, short run because of FLP upgrade
  //  reader.addFile("data/readoutSMA3_2021_10_04_09_06_49.raw");
  //  std::string outFileName = "noInjection_L3off_10hz_2";

  // TED_20211008
  // 1 - injection, L3 on
  // reader.addFile("data/readoutSMA3_2021_10_09__13_01_33__.raw"); //run503701
  reader.addFile("data/readoutSMA3_2021_10_09__12_02_49__.raw"); // run503695

  //  reader.addFile("data/TED_20211008/readoutSMA3_2021_10_08_12_04_51.raw");
  //  //

  //  reader.addFile("data/TED_20211008/readoutSMA3_2021_10_08_12_04_51.raw");
  //  reader.addFile("data/TED_20211008/readoutSMA3_2021_10_08_17_18_12.raw");
  // reader.addFile("data/TED_20211008/readoutSMA3_2021_10_08_09_13_54.raw");
  std::string outFileName = "TED_20211008_InjectionHigh";

  auto hAmpVsTime =
      new TH2F("hAmpVsTime", "hAmpVsTime", 150, 0, 1500, 200, 0, 1024);
  hAmpVsTime->GetXaxis()->SetTitle("time [ns]");
  hAmpVsTime->GetYaxis()->SetTitle("amp [ADC]");

  auto hColVsRow_Amp = new TH2F("hColVsRow_Amp", "hColVsRow_Amp", 48, -0.5,
                                47.5, 24, -0.5, 23.5);
  auto hColVsRow_Acc = new TH2F("hColVsRow_Acc", "hColVsRow_Acc", 48, -0.5,
                                47.5, 24, -0.5, 23.5);
  auto hADC_min = new TH1F("hADC_min", "hADC_min", 100, -0.5, 99.5);
  auto hBC = new TH1I("hBC", "hBC", 3564, -0.5, 3563.5);
  auto hBC_good = new TH1I("hBC_good", "hBC_good", 3564, -0.5, 3563.5);
  hBC_good->SetLineColor(2);

  auto hBCgoodVsNclus = new TH2I("hBCgoodVsNclus", "hBCgoodVsNclus", 3564, -0.5,
                                 3563.5, 10, -0.5, 9.5);
  auto hBCgoodVsLeadClusMult =
      new TH2I("hBCgoodVsLeadClusMult", "hBCgoodVsLeadClusMult", 3564, -0.5,
               3563.5, 10, -0.5, 9.5);
  auto hBCgoodVsLeadClusE = new TH2I("hBCgoodVsLeadClusE", "hBCgoodVsLeadClusE",
                                     3564, -0.5, 3563.5, 102, -0.5, 1023.5);

  auto hAmpVsTime_clusLead =
      new TH2F("hAmpVsTime_clusLead", "hAmpVsTime_clusLead", 150, 0, 1500, 1024,
               0, 1024);
  hAmpVsTime_clusLead->GetXaxis()->SetTitle("clus time [ns]");
  hAmpVsTime_clusLead->GetYaxis()->SetTitle("clus amp [ADC]");

  auto hAmpVsMult_clusLead =
      new TH2F("hAmpVsMult_clusLead", "hAmpVsMult_clusLead", 10, -0.5, 9.5, 200,
               0, 1024);
  hAmpVsMult_clusLead->GetXaxis()->SetTitle("clus mult");
  hAmpVsMult_clusLead->GetYaxis()->SetTitle("clus amp [ADC]");

  auto hColVsRow_clusLead_Amp =
      new TH2F("hColVsRow_clusLead_Amp", "hColVsRow_clusLead_Amp", 48, -0.5,
               47.5, 24, -0.5, 23.5);
  auto hColVsRow_clusLead_Acc =
      new TH2F("hColVsRow_clusLead_Acc", "hColVsRow_clusLead_Acc", 48, -0.5,
               47.5, 24, -0.5, 23.5);

  TH1F *hTimeClusBC[4];
  int histColor[4] = {1, 2, 4, 3};
  for (int i = 0; i < 4; i++) {
    hTimeClusBC[i] = new TH1F(Form("hTimeClusBC_%d", i),
                              Form("hTimeClusBC_%d", i), 150, 0, 1500);
    hTimeClusBC[i]->SetLineColor(histColor[i]);
  }

  TH2F *hNcellMultVsBC_ddl0 =
      new TH2F("hNcellMultVsBC_ddl0", "hNcellMultVsBC_ddl0", 3564, -0.5, 3563.5,
               50, -0.5, 399.5);
  TH2F *hNcellMultVsBC_ddl1 =
      new TH2F("hNcellMultVsBC_ddl1", "hNcellMultVsBC_ddl1", 3564, -0.5, 3563.5,
               50, -0.5, 399.5);

  reader.init();

  // define the standard raw fitter
  // o2::emcal::CaloRawFitterStandard RawFitter;
  o2::emcal::CaloRawFitterGamma2 RawFitter;
  RawFitter.setAmpCut(NoiseThreshold);
  RawFitter.setL1Phase(0.);

  int EventCounter = 0;

  std::unique_ptr<o2::emcal::MappingHandler> mMapper = nullptr;
  if (!mMapper) {
    mMapper = std::unique_ptr<o2::emcal::MappingHandler>(
        new o2::emcal::MappingHandler);
  }
  if (!mMapper) {
    LOG(ERROR) << "Failed to initialize mapper";
  }

  while (1) {
    int tfID = reader.getNextTFToRead();
    if (tfID >= reader.getNTimeFrames()) {
      std::cerr << "nothing left to read after " << tfID << " TFs read";
      break;
    }
    std::vector<char> dataBuffer; // where to put extracted data
    for (int il = 0; il < reader.getNLinks(); il++) {
      auto &link = reader.getLink(il);

      if (debug > 2)
        std::cout << "Decoding link " << il << std::endl;

      auto sz = link.getNextTFSize(); // size in bytes needed for the next TF of
                                      // this link
      dataBuffer.resize(sz);
      link.readNextTF(dataBuffer.data());

      // Parse
      o2::emcal::RawReaderMemory parser(dataBuffer);

      bool isGoodEvent = false;
      bool isLargePulseFound = false;
      int BC_good = -1;
      int FEE_ID_good = -1;
      int Nceels_good = 0;

      EventCounter++;

      while (parser.hasNext()) {
        parser.next();
        //        std::cout << "new page" << std::endl;

        auto &header = parser.getRawHeader();
        auto triggerBC = o2::raw::RDHUtils::getTriggerBC(header);
        auto triggerOrbit = o2::raw::RDHUtils::getTriggerOrbit(header);
        auto feeID = o2::raw::RDHUtils::getFEEID(header);
        auto triggerbits = o2::raw::RDHUtils::getTriggerType(header);

        bool kIsCalib = triggerbits & 0x1 << 6;
        bool kIsPHYSICS = triggerbits & 0x1 << 4;

        BC_good = triggerBC;
        FEE_ID_good = feeID;

        clmaker.setColRange(0, 48);
        if (feeID % 2 == 0) {
          clmaker.setRowRange(0, 7);
        } else {
          clmaker.setRowRange(14, 23);
        }

        if (triggerBC == 3208 && triggerOrbit == 45210564)
          continue;

        if (triggerBC == 764 && triggerOrbit == 157133537)
          continue;

        //        if (triggerBC < 340 || triggerBC > 356)
        //          continue;

        hBC->Fill(triggerBC);

        if (debug > 1) {
          std::cout << "FEEID: " << feeID << std::endl;
          std::cout << "next page: (TriggerBits=" << triggerbits << ",  ";
          if (kIsCalib)
            cout << "CALIB type, ";
          if (kIsPHYSICS)
            cout << "PHYSICS type,  ";
          std::cout << "Orbit=" << triggerOrbit << ", BC=" << triggerBC << ")"
                    << std::endl;
        }

        if (o2::raw::RDHUtils::getFEEID(parser.getRawHeader()) >= 40)
          continue;

        o2::emcal::AltroDecoder decoder(parser);
        decoder.decode();

        int nchannels = 0;
        o2::emcal::Mapper map = mMapper->getMappingForDDL(feeID);
        int iSM = feeID / 2;

        // Loop over all the channels
        for (auto &chan : decoder.getChannels()) {

          int row, col;
          ChannelType_t chantype;
          try {
            row = map.getRow(chan.getHardwareAddress());
            col = map.getColumn(chan.getHardwareAddress());
            chantype = map.getChannelType(chan.getHardwareAddress());
          } catch (Mapper::AddressNotFoundException &ex) {
            std::cerr << "FEE_ID=" << feeID << " : " << ex.what() << std::endl;
            continue;
          };

          if (chantype != 0) // select only HG
            continue;

          if (!IsGoodChannel(col, row))
            continue;

          nchannels++;

          if (debug > 1)
            std::cout << "SM" << iSM << "/DDL" << feeID % 2
                      << "/HW=" << chan.getHardwareAddress() << "/FEC"
                      << chan.getFECIndex() << "/CF" << chantype << " ("
                      << setw(2) << col << "," << setw(2) << row << "): ";

          int ADCmin = 1000;
          int ADCmax = 0;

          for (auto &bunch : chan.getBunches()) {
            for (auto const e : bunch.getADC()) {
              if (debug > 1)
                std::cout << setw(4) << e << " ";

              if (ADCmin > e)
                ADCmin = e;

              if (ADCmax < e)
                ADCmax = e;
            }

            hADC_min->Fill(ADCmin);

            if (ADCmax < cellAmpMin) {
              if (debug > 1)
                std::cout << std::endl;
              continue;
            }

            if (ADCmax > 15)
              Nceels_good++;

            o2::emcal::CaloFitResults fitResults;
            if (kDoFit) {
              try {
                fitResults = RawFitter.evaluate(chan.getBunches(), 0, 0);
              } catch (o2::emcal::CaloRawFitter::RawFitterError_t &fiterror) {
                //      std::cerr << "Error processing raw fit: " <<
                //      o2::emcal::CaloRawFitter::createErrorMessage(fiterror)
                //      << std::endl;
                if (debug > 1)
                  std::cout << std::endl;
                continue;
              }
            }

            if (fitResults.getAmp() < cellAmpMin) {
              if (debug > 1)
                std::cout << std::endl;
              continue;
            }

            isGoodEvent = true;

            if (fitResults.getAmp() > cellAmpMax) {
              isLargePulseFound = true;
            }

            float Ecell_fit = fitResults.getAmp();
            float Tcell_fit = fitResults.getTime();
            Tcell_fit -= ((triggerBC % 4) * 25);

            if (TMath::Abs(Tcell_fit - 570) > 100)
              continue;

            clmaker.setCell(row, col, Ecell_fit, Tcell_fit);

            if (debug > 1) {
              std::cout << "| t=" << Tcell_fit << " ns,  A= " << Ecell_fit
                        << " ADC " << std::endl;
            }

            if (((ADCmax - ADCmin) > 50) && kIsPHYSICS && debug >= 0 && 0) {
              std::cout << "ERROR -- ";
              std::cout << "FEEID: " << feeID << std::endl;
              std::cout << "next page: (TriggerBits=" << triggerbits << ",  ";
              if (kIsCalib)
                cout << "CALIB type, ";
              if (kIsPHYSICS)
                cout << "PHYSICS type,  ";
              std::cout << "Orbit=" << triggerOrbit << ", BC=" << triggerBC
                        << ")" << std::endl;

              std::cout << "SM" << iSM << "/FEC" << chan.getFECIndex() << "/CF"
                        << chantype << " (" << setw(2) << col << "," << setw(2)
                        << row << "): ";

              for (auto const e : bunch.getADC()) {
                std::cout << setw(4) << e << " ";
              }
              std::cout << std::endl;
            }

            if (kDoFit) {
              hAmpVsTime->Fill(Tcell_fit, Ecell_fit);
              hColVsRow_Amp->Fill(col, row, Ecell_fit);
              hColVsRow_Acc->Fill(col, row, 1);
            }
          }
        }
        if (debug > 1)
          std::cout << "channels found : " << nchannels << std::endl;

        if (!isGoodEvent || isLargePulseFound)
          continue;

        hBC_good->Fill(BC_good);
        clmaker.Process();

        if (clmaker.getNclusters() > 0) {
          hAmpVsTime_clusLead->Fill(clmaker.getClusterT(0),
                                    clmaker.getClusterE(0));
          hColVsRow_clusLead_Amp->Fill(clmaker.getClusterCol(0),
                                       clmaker.getClusterRow(0),
                                       clmaker.getClusterE(0));
          hColVsRow_clusLead_Acc->Fill(clmaker.getClusterCol(0),
                                       clmaker.getClusterRow(0));
          hAmpVsMult_clusLead->Fill(clmaker.getClusterMult(0),
                                    clmaker.getClusterE(0));

          hTimeClusBC[BC_good % 4]->Fill(clmaker.getClusterT(0));

          hBCgoodVsNclus->Fill(BC_good, clmaker.getNclusters());
          hBCgoodVsLeadClusMult->Fill(BC_good, clmaker.getClusterMult(0));
          hBCgoodVsLeadClusE->Fill(BC_good, clmaker.getClusterE(0));

          if (FEE_ID_good % 2 == 0)
            hNcellMultVsBC_ddl0->Fill(BC_good, Nceels_good);
          else
            hNcellMultVsBC_ddl1->Fill(BC_good, Nceels_good);
        }
        clmaker.Reset();
      }

    } // DDL loop
    reader.setNextTFToRead(++tfID);
  }

  std::cout << "EventCounter = " << EventCounter << std::endl;

  printf("------ plotting ------ \n");

  auto *cAcc = new TCanvas("cAcc", "cAcc", 10, 10, 900, 300);
  cAcc->Divide(3, 1);
  cAcc->cd(1)->SetLogy();
  hADC_min->Draw("colz");
  cAcc->cd(2);
  hColVsRow_Acc->Draw("colz");
  cAcc->cd(3);
  hColVsRow_Amp->Draw("colz");

  auto *cAmpVsTime = new TCanvas("cAmpVsTime", "cAmpVsTime", 10, 320, 900, 300);
  cAmpVsTime->Divide(3, 1);
  cAmpVsTime->cd(1);
  hAmpVsTime->Draw("colz");
  cAmpVsTime->cd(2);
  hAmpVsTime->ProjectionX()->Draw();
  cAmpVsTime->cd(3)->SetLogy();
  hAmpVsTime->ProjectionY()->Draw();

  TCanvas *cLHC = new TCanvas("cLHC", "cLHC", 10, 410, 1200, 400);
  hBC->Draw("LF2");
  hBC_good->Draw("same LF2");

  auto cClusterLead =
      new TCanvas("cClusterLead", "cClusterLead", 10, 10, 900, 600);
  cClusterLead->Divide(3, 2);
  cClusterLead->cd(1);
  hColVsRow_clusLead_Amp->Draw("colz");
  cClusterLead->cd(2);
  hColVsRow_clusLead_Acc->Draw("colz");
  cClusterLead->cd(3);
  hAmpVsMult_clusLead->Draw("colz");
  cClusterLead->cd(4);
  hAmpVsTime_clusLead->Draw("colz");
  cClusterLead->cd(5);
  hAmpVsTime_clusLead->ProjectionX()->Draw();
  cClusterLead->cd(6);
  hAmpVsTime_clusLead->ProjectionY()->Draw();

  THStack *hs_ClusterLeadingTime_BC =
      new THStack("hs_ClusterLeadingTime_BC", "hs_ClusterLeadingTime_BC");
  for (int i = 0; i < 4; i++)
    hs_ClusterLeadingTime_BC->Add(hTimeClusBC[i]);

  auto cClusterLeadingTime_BC = new TCanvas(
      "cClusterLeadingTime_BC", "cClusterLeadingTime_BC", 10, 10, 300, 300);
  hs_ClusterLeadingTime_BC->Draw("nostack");

  TFile *fout = TFile::Open(Form("%s.root", outFileName.data()), "recreate");
  hAmpVsTime_clusLead->Write(
      Form("%s_%s", hAmpVsTime_clusLead->GetName(), outFileName.data()));

  auto cBCgoodVsCluster =
      new TCanvas("cBCgoodVsCluster", "cBCgoodVsCluster", 10, 10, 300, 900);
  cBCgoodVsCluster->Divide(1, 3);
  cBCgoodVsCluster->cd(1);
  hBCgoodVsNclus->Draw("colz");
  cBCgoodVsCluster->cd(2);
  hBCgoodVsLeadClusMult->Draw("colz");
  cBCgoodVsCluster->cd(3);
  hBCgoodVsLeadClusE->Draw("colz");

  auto cNcellMultVsBC =
      new TCanvas("cNcellMultVsBC", "cNcellMultVsBC", 10, 10, 500, 500);
  cNcellMultVsBC->Divide(1, 2);
  cNcellMultVsBC->cd(1);
  hBC_good->GetXaxis()->SetTitle("BC");
  hBC_good->GetYaxis()->SetTitle("Nevents");

  hBC_good->Draw("LF2");
  //  hNcellMultVsBC_ddl0->DrawCopy("colz");
  //  cNcellMultVsBC->cd(2)->SetLogz();
  //  hNcellMultVsBC_ddl1->DrawCopy("colz");

  TLatex *text = new TLatex();
  text->SetNDC(kTRUE);
  text->SetTextSize(0.04);
  text->DrawLatex(0.7, 0.91, "TED shots 09/10/21");
  text->DrawLatex(0.7, 0.86, "run #503695 ");
  text->DrawLatex(0.7, 0.81, "SMA3");

  cNcellMultVsBC->cd(2); //->SetLogz();
  hNcellMultVsBC_ddl0->Add(hNcellMultVsBC_ddl1);
  //  hNcellMultVsBC_ddl0->GetXaxis()->SetRangeUser(339.5 - 10, 356.5 + 10);
  hNcellMultVsBC_ddl0->GetXaxis()->SetTitle("BC");
  //  hNcellMultVsBC_ddl0->GetYaxis()->SetRangeUser(0, 400);
  hNcellMultVsBC_ddl0->GetYaxis()->SetTitle("Ncells with ADC > 15");
  hNcellMultVsBC_ddl0->DrawCopy("colz");

  return;
}

bool IsGoodChannel(int col, int row) {

  if (col == 4 && row == 23)
    return false;
  if (col == 11 && row == 23)
    return false;
  if (col == 29 && row == 7)
    return false;
  if (col == 31 && row == 7)
    return false;
  if (col == 39 && row == 16)
    return false;
  if (col == 45 && row == 0)
    return false;

  return true;
}
