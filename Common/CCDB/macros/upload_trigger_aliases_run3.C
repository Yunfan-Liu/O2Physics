// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "DataFormatsCTP/Configuration.h"
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "TObjArray.h"
#include "TriggerAliases.h"
#include "TTree.h"
#include "TString.h"
#include <fstream>
#include <map>
#include <string>
#include <vector>
using std::map;
using std::string;

void createDefaultAliases(map<int, TString>& mAliases)
{
  mAliases[kEMC7] = "CTVXEMC-B-NOPF-EMC";
  mAliases[kDMC7] = "CTVXDMC-B-NOPF-EMC";
  mAliases[kTVXinTRD] = "CMTVX-B-NOPF-TRD,minbias_TVX";
  mAliases[kTVXinEMC] = "C0TVX-B-NOPF-EMC,minbias_TVX_L0,CMTVXTSC-B-NOPF-EMC,CMTVXTCE-B-NOPF-EMC";
  mAliases[kTVXinPHOS] = "C0TVX-B-NOPF-PHSCPV,minbias_TVX_L0,CMTVXTSC-B-NOPF-PHSCPV,CMTVXTSC-B-NOPF-PHSCPV";
  mAliases[kTVXinHMP] = "C0TVX-B-NOPF-HMP,minbias_TVX_L0,CMTVXTSC-B-NOPF-HMP";
  mAliases[kPHOS] = "CTVXPH0-B-NOPF-PHSCPV,mb_PH0_TVX,CPH0SC-B-NOPF-PHSCPV,CPH0CE-B-NOPF-PHSCPV";
}

void createPbPbAliases(map<int, TString>& mAliases)
{
  mAliases[kTVXinTRD] = "CMTVXTSC-B-NOPF-TRD,CMTVXTCE-B-NOPF-TRD";
  mAliases[kTVXinEMC] = "CMTVXTSC-B-NOPF-EMC,CMTVXTCE-B-NOPF-EMC,C0TVXTSC-B-NOPF-EMC,C0TVXTCE-B-NOPF-EMC";
  mAliases[kTVXinPHOS] = "CMTVXTSC-B-NOPF-PHSCPV,CMTVXTCE-B-NOPF-PHSCPV,C0TVXTSC-B-NOPF-PHSCPV,C0TVXTCE-B-NOPF-PHSCPV";
  mAliases[kTVXinHMP] = "CMTVXTSC-B-NOPF-HMP,CMTVXTCE-B-NOPF-HMP";
  mAliases[kPHOS] = "CPH0SC-B-NOPF-PHSCPV,CPH0CE-B-NOPF-PHSCPV";
}

void upload_trigger_aliases_run3()
{
  map<int, TString> mAliases;
  createDefaultAliases(mAliases);
  // createPbPbAliases(mAliases);

  TObjArray* classNames[kNaliases];
  for (auto& al : mAliases) {
    classNames[al.first] = (al.second).Tokenize(",");
  }

  o2::ccdb::CcdbApi ccdb;
  ccdb.init("http://alice-ccdb.cern.ch");
  // ccdb.init("http://ccdb-test.cern.ch:8080");
  // ccdb.truncate("EventSelection/TriggerAliases");

  map<string, string> metadata;

  // read list of runs from text file
  std::ifstream f("runs.txt");
  std::vector<int> runs;
  int r = 0;
  while (f >> r) {
    runs.push_back(r);
  }

  if (0) {
    ULong64_t sor = 1672531200000;
    ULong64_t eor = 1893456000000;
    TriggerAliases* aliases = new TriggerAliases();
    metadata["runNumber"] = "default";
    ccdb.storeAsTFileAny(aliases, "EventSelection/TriggerAliases", metadata, sor, eor);
  }

  for (auto& run : runs) {
    LOGP(info, "run = {}", run);

    //    if (run != 535983)
    //      continue; // already filled
    //    if (run <= 539580)
    //      continue; // already filled
    if (run < 519903)
      continue; // no CTP info
    if (run == 527349)
      continue; // no CTP info
    if (run == 527963)
      continue; // no CTP info
    if (run == 528537)
      continue; // no CTP info
    if (run == 528543)
      continue; // no CTP info

    // read SOR and EOR timestamps from RCT CCDB via utility function
    auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdb, run);
    auto eor = soreor.second;
    auto sor = soreor.first;
    auto ts = sor;

    // read CTP config
    metadata["runNumber"] = Form("%d", run);
    auto ctpcfg = ccdb.retrieveFromTFileAny<o2::ctp::CTPConfiguration>("CTP/Config/Config", metadata, ts);
    if (!ctpcfg)
      continue;

    if (run == 529414) { // adding tolerance to sor for this run
      sor = 1668809980000;
    }

    std::vector<o2::ctp::CTPClass> classes = ctpcfg->getCTPClasses();
    ctpcfg->printConfigString();
    // create trigger aliases
    TriggerAliases* aliases = new TriggerAliases();
    for (auto& al : mAliases) {
      int aliasId = al.first;
      LOGP(debug, "alias = {}", al.second.Data());
      for (const auto& className : *(classNames[aliasId])) {
        TString sname = className->GetName();
        LOGP(debug, " className = {}", sname.Data());
        sname.ToUpper();
        for (const auto& cl : classes) {
          TString clname = cl.name;
          clname.ToUpper();
          if (clname != sname) {
            continue;
          }
          int classId = cl.getIndex();
          TString cluster = TString(cl.cluster->name);
          cluster.ToUpper();
          if (aliasId == kTVXinTRD && cluster != "TRD") { // workaround for configs with ambiguous class names
            continue;
          }
          if (aliasId == kTVXinEMC && cluster != "EMC") { // workaround for configs with ambiguous class names
            continue;
          }
          if (aliasId == kTVXinPHOS && cluster != "PHSCPV") { // workaround for configs with ambiguous class names
            continue;
          }
          if (aliasId == kTVXinHMP && cluster != "HMP") { // workaround for configs with ambiguous class names
            continue;
          }
          LOGP(debug, " class index = {}, name = {}, cluster = {}", classId, cl.name, cl.cluster->name);
          aliases->AddClassIdToAlias(aliasId, classId);
          break;
        }
      }
    }
    aliases->Print();
    ccdb.storeAsTFileAny(aliases, "EventSelection/TriggerAliases", metadata, sor - 1000, eor + 10000); // adding tolerance of 10s to eor
  }
}
