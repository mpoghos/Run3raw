#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "RStringView.h"
#include <Rtypes.h>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include <TCanvas.h>
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

/*
#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "RStringView.h"
#include "TFile.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTree.h"
#include <Rtypes.h>
#include <array>
#include <fstream>
#include <iostream>
#include <vector>

//#include "EMCALReconstruction/RawHeaderStream.h"
#endif
*/

using namespace o2::emcal;

const int debug = -3;

const Bool_t kDoFit = kFALSE;

void AnalData() {
  const Int_t NoiseThreshold = 3;

  // Use the RawReaderFile to read the raw data file
  const char *aliceO2env = std::getenv("O2_ROOT");
  std::string inputDir = " ";

  o2::raw::RawFileReader reader;
  reader.setDefaultDataOrigin(o2::header::gDataOriginEMC);
  reader.setDefaultDataDescription(o2::header::gDataDescriptionRawData);
  reader.setDefaultReadoutCardType(o2::raw::RawFileReader::RORC);
  //  reader.addFile("../public/readoutSMA10_1FEC.raw");
  //  reader.addFile("../public/readout_ForSTFB_BothLinks.raw");
  //  reader.addFile("../public/readoutSMA10_aaa.raw");
  //  reader.addFile("/eos/home-p/poghos/P2_data/SRUdata/readoutSMA9_0.raw");
  //  reader.addFile("data/readoutSMC5_ALIECS_1.raw");
  //  reader.addFile("data/readoutSMC5_1.raw");
  //  reader.addFile("/home/flp/data/ForAnal/readoutSMC345_PHYSandCALIB_LowStat.raw");
  //  reader.addFile("/home/flp/data/readoutSMC3.raw");

  //  reader.addFile("/home/emc/readoutSMC9.raw");
  //  reader.addFile("/home/emc/data/readoutSMA910.raw");
  //  reader.addFile("/home/flp/data/readoutDDL_23.raw");
  //  reader.addFile("/home/flp/data/readoutEMCAL_SMA4_calib.raw");
  //  reader.addFile("/home/flp/data/readoutSMA3_run0.raw");
  reader.addFile("/home/flp/data/readoutEMCAL.raw");
  //  reader.addFile("/home/flp/data/TED_20211008/readoutSMA3_2021_10_08_08_46_42.raw");

  //  reader.addFile("/home/flp/data/SMA1_FEC36_channel/readoutDDL_1103.raw");
  //  reader.addFile("/home/flp/data/SMA1_FEC36/readoutDDL_1103.raw");

  //  reader.addFile("/home/emc/readoutSMAC9.raw");
  //  reader.addFile("/tmp/data.raw");

  //  reader.addFile("/tmp/readoutVS.raw");
  //  reader.addFile("/home/flp/data/readoutSMC4.raw");
  //  reader.addFile("readoutSMC5_1_LEDwithPed.raw");
  //  reader.addFile("readoutSMC5_1_PHYSwoPed.raw");

  // auto hAmpSM = new TH1F("hAmpSM", "hAmpSM", 200, 0, 1024);
  auto hAmpVsTime_SM =
      new TH2F("hAmpVsTime_SM", "hAmpVsTime_SM", 150, 0, 1500, 200, 0, 1024);
  hAmpVsTime_SM->GetXaxis()->SetTitle("time [ns]");
  hAmpVsTime_SM->GetYaxis()->SetTitle("amp [ADC]");

  auto hColVsRow_Amp = new TH2F("hColVsRow_Amp", "hColVsRow_Amp", 2 * 48, -0.5,
                                48 + 47.5, 10 * 24, -0.5, 239.5);

  TH2F *hColVsRow[2];
  TH1F *hADC_min[20][2];
  TH1F *hADC_max[20][2];
  TH1F *hAmp_PHYS[20];
  TH1F *hAmp_CAL[20];

  for (int iG = 0; iG < 2; iG++) {
    hColVsRow[iG] = new TH2F(Form("hColVsRow_%d", iG), Form("hColVsRow_%d", iG),
                             2 * 48, -0.5, 48 + 47.5, 10 * 24, -0.5, 239.5);
  }

  for (int iSM = 0; iSM < 20; iSM++) {
    hAmp_PHYS[iSM] = new TH1F(Form("hAmp_PHYS_SM%02d", iSM),
                              Form("hAmp_PHYS_SM%02d", iSM), 500, 0, 500);
    hAmp_CAL[iSM] = new TH1F(Form("hAmp_CAL_SM%02d", iSM),
                             Form("hAmp_CAL_SM%02d", iSM), 500, 0, 500);
    hAmp_PHYS[iSM]->SetLineColor(1);
    hAmp_CAL[iSM]->SetLineColor(2);

    for (int iG = 0; iG < 2; iG++) {
      hADC_min[iSM][iG] =
          new TH1F(Form("hADC_min_SM%02d_%s", iSM, iG == 0 ? "HG" : "LG"),
                   Form("hADC_min_SM%002d_%s", iSM, iG == 0 ? "HG" : "LG"), 100,
                   -0.5, 99.5);
      hADC_max[iSM][iG] =
          new TH1F(Form("hADC_max_SM%02d_%s", iSM, iG == 0 ? "HG" : "LG"),
                   Form("hADC_max_SM%002d_%s", iSM, iG == 0 ? "HG" : "LG"), 500,
                   -0.5, 499.5);
    }
  }

  auto *hBC = new TH1I("hBC", "hBC", 3564, -0.5, 3563.5);
  auto *hBCvsNChannel =
      new TH1I("hBCvsNChannel", "hBCvsNChannel", 3564, -0.5, 3563.5);

  reader.init();

  // return;

  // define the standard raw fitter
  // o2::emcal::CaloRawFitterStandard RawFitter;
  o2::emcal::CaloRawFitterGamma2 RawFitter;
  RawFitter.setAmpCut(NoiseThreshold);
  RawFitter.setL1Phase(0.);

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

      // parser.getRawHeader();

      // std::cout << "  " <<parser.getRawHeader().getLink() << std::endl;

      // o2::header::RDHAny RawReaderMemory::decodeRawHeader(const void*
      // payloadwords)

      while (parser.hasNext()) {
        parser.next();

        auto &header = parser.getRawHeader();
        auto triggerBC = o2::raw::RDHUtils::getTriggerBC(header);
        auto triggerOrbit = o2::raw::RDHUtils::getTriggerOrbit(header);
        auto feeID = o2::raw::RDHUtils::getFEEID(header);
        auto triggerbits = o2::raw::RDHUtils::getTriggerType(header);

        bool kIsCalib = triggerbits & 0x1 << 6;
        bool kIsPHYSICS = triggerbits & 0x1 << 4;

        // if(feeID!=0)
        // continue;

        // if(triggerBC!=0)
        // continue;

        kIsCalib = false;

        //      if(!kIsCalib)
        // continue;

        hBC->Fill(triggerBC);

        //		if(!kIsCalib )
        //		  continue;

        //		if(!kIsPHYSICS )
        //		  continue;

        if (debug > 1 || kIsCalib)
        //      if(debug > 1)
        {
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

        //        std::cout<<rawreader.getRawHeader()<<std::endl;

        // use the altro decoder to decode the raw data, and extract the RCU
        // trailer
        o2::emcal::AltroDecoder decoder(parser);
        decoder.decode();

        //        std::cout << decoder.getRCUTrailer() << std::endl;

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
            std::cerr << "FEE_ID=" << feeID << " : ";
            std::cerr << ex.what() << std::endl;
            //	  std::cerr << "FEE_ID"<<feeID  << "/HW=" <<
            //chan.getHardwareAddress() <<"\n";
            continue;
          };

          if (!(chantype == 0 || chantype == 1))
            // if(chantype!=2)
            // if(chantype!=0)
            continue;

          o2::emcal::CaloFitResults fitResults;
          if (kDoFit) {
            // fitResults = RawFitter.evaluate(chan.getBunches());

            try {
              fitResults = RawFitter.evaluate(chan.getBunches());
              //       std::cout << "Fit done" << std::endl;
            } catch (o2::emcal::CaloRawFitter::RawFitterError_t &fiterror) {
              //      std::cerr << "Error processing raw fit: " <<
              //      o2::emcal::CaloRawFitter::createErrorMessage(fiterror) <<
              //      std::endl;
              continue;
            }
          }
          // fitResults = RawFitter.evaluate(chan.getBunches(), 0, 0);
          // print the fit output

          nchannels++;

          //	std::cout << "HW=" << chan.getHardwareAddress() <<  " (FEC"<<
          //chan.getFECIndex() << "): ";

          if (debug > 1 || kIsCalib)
            std::cout << setw(4) << "SM" << iSM << "/DDL" << setw(2)
                      << feeID % 2 << "/HW=" << setw(4)
                      << chan.getHardwareAddress() << "/FEC"
                      << chan.getFECIndex() << "/CF" << chantype << " ("
                      << setw(2) << col << "," << setw(2) << row << "): ";

          int ADCmin = 1000;
          int ADCmax = 0;

          for (auto &bunch : chan.getBunches()) {
            //		std::cout << "BunchLength=" << (int) bunch.getBunchLength() <<
            //"  StartTiime=" <<  (int) bunch.getStartTime() << "  :";
            for (auto const e : bunch.getADC()) {
              if (debug > 1 || kIsCalib)
                std::cout << setw(4) << e << " ";

              if (ADCmin > e)
                ADCmin = e;

              if (ADCmax < e)
                ADCmax = e;
            }
            //		std::cout << std::endl;

            if (((ADCmax - ADCmin) > 50) && 1) {

              std::cout << setw(4) << "SM" << iSM << "/FEC"
                        << chan.getFECIndex() << "/CF" << chantype << " ("
                        << setw(2) << col << "," << setw(2) << row << "): ";

              for (auto const e : bunch.getADC()) {
                std::cout << setw(4) << e << " ";
              }

              /*
                      std::cout << "ERROR -- " ;
                      std::cout << "FEEID: " << feeID  << std::endl;
                      std::cout << "next page: (TriggerBits="<< triggerbits <<
                 ",  "; if(kIsCalib ) cout << "CALIB type, "; if(kIsPHYSICS)
                              cout << "PHYSICS type,  ";
                      std::cout << "Orbit=" << triggerOrbit  << ", BC=" <<
                 triggerBC <<")"<< std::endl;
              */

              std::cout << "| BC=" << triggerBC;

              std::cout << std::endl;

              hBCvsNChannel->Fill(triggerBC);
            }
          }

          if ((debug > 1 || kIsCalib)) {
            if ((fitResults.getAmp() > NoiseThreshold) && kDoFit)
              std::cout << "| t=" << fitResults.getTime()
                        << " ns,  A= " << fitResults.getAmp() << " ADC "
                        << std::endl;
            else
              std::cout << std::endl;
          }

          int colGlob = col;
          int rowGlob = row;

          if (iSM % 2 == 1)
            col += 48;
          row += (iSM / 2) * 24;

          if (chantype < 2) {
            hADC_min[iSM][chantype]->Fill(ADCmin);
            hADC_max[iSM][chantype]->Fill(ADCmax);
            if (kDoFit) {
              hAmpVsTime_SM->Fill(fitResults.getTime(), fitResults.getAmp());
              hColVsRow_Amp->Fill(col, row, fitResults.getAmp());
            }
            hColVsRow[chantype]->Fill(col, row);

            if (kIsCalib)
              hAmp_CAL[iSM]->Fill(ADCmax - ADCmin);
            //		    hAmp_CAL[iSM]->Fill(fitResults.getAmp());

            if (kIsPHYSICS)
              hAmp_PHYS[iSM]->Fill(ADCmax - ADCmin);
            //		    hAmp_PHYS[iSM]->Fill(fitResults.getAmp());
          }
        }

        if (debug > 1)
          std::cout << "channels found : " << nchannels << std::endl;
      }
    }
    reader.setNextTFToRead(++tfID);
  }

  printf("------ plotting ------ \n");

  auto *cColVsRow = new TCanvas("cColVsRow", "cColVsRow", 10, 10, 600, 300);
  cColVsRow->Divide(2, 1);
  cColVsRow->cd(1);
  hColVsRow[0]->Draw("colz");
  cColVsRow->cd(2);
  hColVsRow[1]->Draw("colz");

  auto *cADCminHG = new TCanvas("cADCminHG", "cADCminHG", 10, 10, 400, 1200);
  auto *cADCminLG = new TCanvas("cADCminLG", "cADCminLG", 10, 10, 400, 1200);
  cADCminHG->Divide(2, 10);
  cADCminLG->Divide(2, 10);

  auto *cADCmaxHG = new TCanvas("cADCmaxHG", "cADCmaxHG", 10, 10, 400, 1200);
  auto *cADCmaxLG = new TCanvas("cADCmaxLG", "cADCmaxLG", 10, 10, 400, 1200);
  cADCmaxHG->Divide(2, 10);
  cADCmaxLG->Divide(2, 10);

  for (int ic = 0; ic < 20; ic++) {
    cADCminHG->cd(ic + 1);
    hADC_min[ic][0]->Draw();
    cADCminLG->cd(ic + 1);
    hADC_min[ic][1]->Draw();

    cADCmaxHG->cd(ic + 1);
    hADC_max[ic][0]->Draw();
    cADCmaxLG->cd(ic + 1);
    hADC_max[ic][1]->Draw();
  }

  TCanvas *cLHC = new TCanvas("cLHC", "cLHC", 10, 410, 400, 400);
  hBC->Draw("LF2");

  TCanvas *cBCvsNChannel =
      new TCanvas("cBCvsNChannel", "cBCvsNChannel", 10, 410, 400, 400);
  hBCvsNChannel->Divide(hBC);
  hBCvsNChannel->Draw("LF2");

  return;
  auto *cAmpVsTime_SM =
      new TCanvas("cAmpVsTime_SM", "cAmpVsTime_SM", 10, 10, 900, 300);
  cAmpVsTime_SM->Divide(3, 1);
  cAmpVsTime_SM->cd(1);
  hAmpVsTime_SM->Draw("colz");
  cAmpVsTime_SM->cd(2);
  hAmpVsTime_SM->ProjectionX()->Draw();
  cAmpVsTime_SM->cd(3);
  hAmpVsTime_SM->ProjectionY()->Draw();

  cADCminHG->Divide(4, 1);
  cADCminHG->cd(1);
  hADC_min[13][0]->Draw();
  cADCminHG->cd(2);
  hADC_min[15][0]->Draw();
  cADCminHG->cd(3);
  hADC_min[17][0]->Draw();
  cADCminHG->cd(4);
  hADC_min[19][0]->Draw();
  // cADCminHG->cd(3);
  //  hADC_min[11][0]->Draw();

  THStack *hs[20];
  for (Int_t iSM = 0; iSM < 20; iSM++) {
    hs[iSM] = new THStack(Form("hs_iSM%d", iSM), Form("hs_iSM%d", iSM));
    hs[iSM]->Add(hAmp_PHYS[iSM]);
    hs[iSM]->Add(hAmp_CAL[iSM]);
  }

  auto *cAmp = new TCanvas("cAmp", "cAmp", 10, 10, 900, 300);
  cAmp->Divide(4, 1);
  cAmp->cd(1)->SetLogy();
  //  hAmp_PHYS[7]->Draw();
  //  hAmp_CAL[7]->Draw("same");
  hs[13]->Draw("nostack");
  cAmp->cd(2)->SetLogy();
  //  hAmp_PHYS[9]->Draw();
  //  hAmp_CAL[9]->Draw("same");
  hs[15]->Draw("nostack");
  // cAmp->cd(3)->SetLogy();
  //  hAmp_PHYS[11]->Draw();
  //  hAmp_CAL[11]->Draw("same");
  //  hs[11]->Draw("nostack");
  cAmp->cd(3)->SetLogy();
  hs[17]->Draw("nostack");
  cAmp->cd(4)->SetLogy();
  hs[19]->Draw("nostack");
}
