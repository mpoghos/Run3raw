

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

const int debug = -2;

const Bool_t kDoFit = kTRUE; // kFALSE;

bool IsGoodChannel(int col, int row);

void AnalData_TED() {
  const Int_t NoiseThreshold = 3;
  const float cellAmpMin = 5.;

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
  //  reader.addFile("/home/flp/data/readoutSMA3_2021_10_01_11_26_55.raw");
  //  reader.addFile("/home/flp/data/TED_20211008/readoutSMA3_2021_10_08__08_46_42__.raw");
  reader.addFile(
      "/home/flp/data/TED_20211017/readoutSMA3_2021_10_16__17_56_45__.raw");
  //  reader.addFile("/home/flp/data/readoutEMCAL.raw");

  //  reader.addFile("/home/flp/data/SMA1_FEC36_channel/readoutDDL_1103.raw");
  //  reader.addFile("/home/flp/data/SMA1_FEC36/readoutDDL_1103.raw");

  //  reader.addFile("/home/emc/readoutSMAC9.raw");
  //  reader.addFile("/tmp/data.raw");

  //  reader.addFile("/tmp/readoutVS.raw");
  //  reader.addFile("/home/flp/data/readoutSMC4.raw");
  //  reader.addFile("readoutSMC5_1_LEDwithPed.raw");
  //  reader.addFile("readoutSMC5_1_PHYSwoPed.raw");

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

      bool isGoodEvent = false;
      int BC_good = -1;

      while (parser.hasNext()) {
        parser.next();

        auto &header = parser.getRawHeader();
        auto triggerBC = o2::raw::RDHUtils::getTriggerBC(header);
        auto triggerOrbit = o2::raw::RDHUtils::getTriggerOrbit(header);
        auto feeID = o2::raw::RDHUtils::getFEEID(header);
        auto triggerbits = o2::raw::RDHUtils::getTriggerType(header);

        bool kIsCalib = triggerbits & 0x1 << 6;
        bool kIsPHYSICS = triggerbits & 0x1 << 4;

        BC_good = triggerBC;

        if (triggerBC == 3208 && triggerOrbit == 45210564)
          continue;

        hBC->Fill(triggerBC);

        if (debug > 1 || feeID == 0) {
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
            // chan.getHardwareAddress() <<"\n";
            continue;
          };

          if (!(chantype == 0 || chantype == 1))
            // if(chantype!=2)
            //	if(chantype!=0)
            continue;

          if (!IsGoodChannel(col, row))
            continue;

          nchannels++;

          //	std::cout << "HW=" << chan.getHardwareAddress() <<  " (FEC"<<
          // chan.getFECIndex() << "): ";

          if (debug > 1)
            std::cout << "SM" << iSM << "/DDL" << feeID % 2
                      << "/HW=" << chan.getHardwareAddress() << "/FEC"
                      << chan.getFECIndex() << "/CF" << chantype << " ("
                      << setw(2) << col << "," << setw(2) << row << "): ";

          int ADCmin = 1000;
          int ADCmax = 0;

          for (auto &bunch : chan.getBunches()) {
            //		std::cout << "BunchLength=" << (int)
            //bunch.getBunchLength()
            //<< "  StartTiime=" <<  (int) bunch.getStartTime() << "  :";
            for (auto const e : bunch.getADC()) {
              if (debug > 1)
                std::cout << setw(4) << e << " ";

              if (ADCmin > e)
                ADCmin = e;

              if (ADCmax < e)
                ADCmax = e;
            }
            //		std::cout << std::endl;

            hADC_min->Fill(ADCmin);

            //  }

            if (ADCmax < cellAmpMin) {
              if (debug > 1)
                std::cout << std::endl;
              continue;
            }
            o2::emcal::CaloFitResults fitResults;
            if (kDoFit) {
              try {
                fitResults = RawFitter.evaluate(chan.getBunches());
                //       std::cout << "Fit done" << std::endl;
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

            if (debug > 1) {
              std::cout << "| t=" << fitResults.getTime()
                        << " ns,  A= " << fitResults.getAmp() << " ADC "
                        << std::endl;
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

            int colGlob = col;
            int rowGlob = row;

            // if(iSM%2==1) col+=48;
            // row+=(iSM/2)*24;

            if (chantype == 0) {
              //		hADC_min[iSM][chantype]->Fill(ADCmin);
              if (kDoFit) {
                hAmpVsTime->Fill(fitResults.getTime(), fitResults.getAmp());
                hColVsRow_Amp->Fill(col, row, fitResults.getAmp());
                hColVsRow_Acc->Fill(col, row, 1);
              }
            }
          }
        }
        if (debug > 1)
          std::cout << "channels found : " << nchannels << std::endl;
      }

      if (isGoodEvent)
        hBC_good->Fill(BC_good);

    } // DDL loop
    reader.setNextTFToRead(++tfID);
  }

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
  if (col == 45 && row == 0)
    return false;

  return true;
}
