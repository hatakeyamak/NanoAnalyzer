import ROOT

# Enable multi-threading
#ROOT.ROOT.EnableImplicitMT()
ROOT.ROOT.DisableImplicitMT()

# Load NanoAOD file
rdf = ROOT.RDataFrame("Events", "/cms/data/hatake/ana/DarkHiggs/CMSSW_15_0_10/src/NanoAnalyzer/DarkHiggs/data/PrivSample_DarkHiggsToWW_Zp2000_s300_Chi200_custom_RunIII2024Summer24NanoAODv15_job1.root")

print("start")
#-------------------------------------------------------------------------------
# 1. Define a mask for GenParticles that are electrons (pdgId 11)
# GenPart_pdgId is a vector branch, Define applies this per-event
rdf_with_mask = rdf.Define("GenElecMask", "abs(GenPart_pdgId) == 11")

# 2. Filter events: Keep only events where ANY electron passes the mask
rdf_filtered = rdf_with_mask.Filter("ROOT::VecOps::Any(GenElecMask)")

# 3. (Optional) Filter by pT: Keep events with electrons > 20 GeV
rdf_final = rdf_filtered.Filter("ROOT::VecOps::Any(GenPart_pt[GenElecMask] > 20)")

#-------------------------------------------------------------------------------
# Define quarks
rdf = rdf.Define(
    "GenPart_WParentQuark",
    """
    ROOT::VecOps::RVec<int> motherIdx = GenPart_genPartIdxMother;
    ROOT::VecOps::RVec<int> pdgId = GenPart_pdgId;
    ROOT::VecOps::RVec<bool> isWQuark(pdgId.size(), false);
    
    for (size_t i = 0; i < pdgId.size(); ++i) {
        int mothIdx = motherIdx[i];
        if (mothIdx >= 0 && mothIdx < pdgId.size()) {
            if (std::abs(pdgId[i]) <= 5 && std::abs(pdgId[mothIdx]) == 24) {
                isWQuark[i] = true;
            }
        }
    }
    return isWQuark;
    """
)

# Count quarks
rdf = rdf.Define("nWQuarks", "Sum(GenPart_WParentQuark)")

# Define a subset based on a condition
rdf = rdf.Define("WQuarkMask", "GenPart_WParentQuark == 1") \
         .Define("WQuark_pt", "GenPart_pt[WQuarkMask]") \
         .Define("WQuark_eta", "GenPart_eta[WQuarkMask]") \
         .Define("WQuark_phi", "GenPart_phi[WQuarkMask]") \
         .Define("WQuark_mass", "GenPart_mass[WQuarkMask]")

# This function handles a variable number of entries (particles) per event
# --- ComputeTotalInvMass ---
ROOT.gInterpreter.Declare("""
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include "Math/LorentzVector.h"

double ComputeTotalInvMass(ROOT::RVec<double> pt, ROOT::RVec<double> eta, 
                           ROOT::RVec<double> phi, ROOT::RVec<double> mass) {
    if (pt.size() < 2) return 0.0; // Need at least two particles
    
    ROOT::Math::PtEtaPhiMVector totalVector(0,0,0,0);
    for (size_t i = 0; i < pt.size(); ++i) {
        ROOT::Math::PtEtaPhiMVector p(pt[i], eta[i], phi[i], mass[i]);
        totalVector += p;
    }
    return totalVector.M();
}
""")

# Declare the C++ helper function
# Find closest GenPart (status 1 or 23, prompt) to each Jet
# We assume GenPart_pt, Jet_pt etc., are available
rdf = rdf.Define("JetGenMatchedIdx", 
    """
    #include "ROOT/RVec.hxx"
    #include "Math/Vector4D.h"
    #include "Math/LorentzVector.h"

    ROOT::RVec<int> matchedIdx(Jet_pt.size(), -1);
    for (size_t i = 0; i < Jet_pt.size(); ++i) {
        float minDR = 0.4; // Max DeltaR
        int bestGen = -1;
        for (size_t j = 0; j < GenPart_pt.size(); ++j) {
            // Apply GenPart selection (e.g., status 1)
            if (GenPart_WParentQuark[j]) { 
                float dr = ROOT::VecOps::DeltaR(Jet_eta[i], GenPart_eta[j], 
                                              Jet_phi[i], GenPart_phi[j]);
                if (dr < minDR) {
                    minDR = dr;
                    bestGen = j;
                }
            }
        }
        matchedIdx[i] = bestGen;
    }
    return matchedIdx;
    """)

# Define a subset based on a condition
rdf = rdf.Define("MatchJetMask", "JetGenMatchedIdx > 0") \
         .Define("MatchJet_pt", "Jet_pt[MatchJetMask]") \
         .Define("MatchJet_eta", "Jet_eta[MatchJetMask]") \
         .Define("MatchJet_phi", "Jet_phi[MatchJetMask]") \
         .Define("MatchJet_mass", "Jet_mass[MatchJetMask]")

# WQuark invariant mass
rdf = rdf.Define("WQuark_inv_mass", 
                 "ComputeTotalInvMass(WQuark_pt, WQuark_eta, WQuark_phi, WQuark_mass)")

# Define leptons
rdf = rdf.Define(
    "GenPart_WParentLepton",
    """
    ROOT::VecOps::RVec<int> motherIdx = GenPart_genPartIdxMother;
    ROOT::VecOps::RVec<int> pdgId = GenPart_pdgId;
    ROOT::VecOps::RVec<bool> isWLepton(pdgId.size(), false);
    
    for (size_t i = 0; i < pdgId.size(); ++i) {
        int mothIdx = motherIdx[i];
        if (mothIdx >= 0 && mothIdx < pdgId.size()) {
            if (std::abs(pdgId[i]) >= 11 && std::abs(pdgId[i]) <= 16 &&
                std::abs(pdgId[mothIdx]) == 24) {
                isWLepton[i] = true;
            }
        }
    }
    return isWLepton;
    """
)

# Count leptons
rdf = rdf.Define("nWLeptons", "Sum(GenPart_WParentLepton)")

# Total
rdf = rdf.Define("nWDecayProducts", "nWQuarks + nWLeptons")

# AK4 invariant mass
rdf = rdf.Define("AK4_inv_mass", 
                 "ComputeTotalInvMass(Jet_pt, Jet_eta, Jet_phi, Jet_mass)")
rdf = rdf.Define("MatchJet_inv_mass", 
                 "ComputeTotalInvMass(MatchJet_pt, MatchJet_eta, MatchJet_phi, MatchJet_mass)")

# Make histogram
h_nQuarks = rdf.Histo1D(("nQuarks", "Number of W-parent Quarks;N;Events", 5, 0, 5), "nWQuarks")
h_nWDecayProducts = rdf.Histo1D(
    ("h_nWDecayProducts", "W decay products;N;Events", 10, 0, 10),
    "nWDecayProducts"
)
h_nFatJet = rdf.Histo1D(("nFatJet","nFatJet;nFatJet;Events", 5, 0, 5), "nFatJet")

# Draw to show
c1 = ROOT.TCanvas()
h_nQuarks.Draw()
h_nWDecayProducts.Draw()
c1.SaveAs("nWQuarks.png")
c1.SaveAs("nWDecayProducts.png")

#-------------------------------------------------------------------------------
# 1. Both W decays hadronically
rdf_FullHad = rdf.Filter("nWQuarks == 4", "Events with four WQuarks")

h_nFatJet_FullHad = rdf_FullHad.Histo1D(("nFatJet_FullHad","nFatJet_FullHad;nFatJet;Events", 5, 0, 5), "nFatJet")
h_nQuarks_FullHad = rdf_FullHad.Histo1D(("nQuarks_FullHad", "Number of W-parent Quarks;NWQuarks;Events", 5, 0, 5), "nWQuarks")
h_WQuark_inv_mass_FullHad = rdf_FullHad.Histo1D(("WQuark_inv_mass_FullHad","WQuark_inv_mass;WQuark_inv_mass;Events", 100, 0, 500), "WQuark_inv_mass")
h_AK4_inv_mass_FullHad = rdf_FullHad.Histo1D(("AK4_inv_mass","AK4_inv_mass;AK4_inv_mass;Events", 100, 0, 500), "AK4_inv_mass")
h_MatchJet_inv_mass_FullHad = rdf_FullHad.Histo1D(("MatchJet_inv_mass","MatchJet_inv_mass;MatchJet_inv_mass;Events", 100, 0, 500), "MatchJet_inv_mass")

#rdf_FullHad.Display(["WQuark_inv_mass", "AK4_inv_mass", "MatchJet_inv_mass", "Jet_pt", "JetGenMatchedIdx"],10).Print()

#-------------------------------------------------------------------------------
# 2. Basic event selection
rdf_Baseline = rdf_FullHad.Filter("nFatJet >=1 && FatJet_pt[0]>200 && abs(FatJet_eta[0])<2.4 && PuppiMET_pt>250.", "Events with four WQuarks")

# Declare the C++ helper function
# Find closest GenPart (status 1 or 23, prompt) to each Jet
# We assume GenPart_pt, Jet_pt etc., are available
rdf_Baseline = rdf_Baseline.Define("FatJetMatchedIdx", 
    """
    #include "ROOT/RVec.hxx"
    #include "Math/Vector4D.h"
    #include "Math/LorentzVector.h"

    ROOT::RVec<float> deltaRToFatJet(Jet_pt.size(), -1);
    for (size_t i = 0; i < Jet_pt.size(); ++i) {
        float minDR = 100.; // Max DeltaR
        float dr = ROOT::VecOps::DeltaR(Jet_eta[i], FatJet_eta[0], 
                                        Jet_phi[i], FatJet_phi[0]);
        deltaRToFatJet[i] = dr;
    }
    return deltaRToFatJet;
    """)

# Define a subset based on a condition
#rdf = rdf.Define("MatchJetMask", "JetGenMatchedIdx > 0") \
#         .Define("MatchJet_pt", "Jet_pt[MatchJetMask]") \
#         .Define("MatchJet_eta", "Jet_eta[MatchJetMask]") \
#         .Define("MatchJet_phi", "Jet_phi[MatchJetMask]") \
#         .Define("MatchJet_mass", "Jet_mass[MatchJetMask]")

display2 = rdf_Baseline.Display(["WQuark_inv_mass", "AK4_inv_mass", "MatchJet_inv_mass", 
    "FatJet_pt", "FatJet_eta", "FatJet_phi", "FatJet_mass", 
    "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass",
    "JetGenMatchedIdx", "FatJetMatchedIdx"],10)
print(display2.AsString())

#-------------------------------------------------------------------------------
# Specifically inspect nFatJet==2 events
rdf_FullHad_2FatJet = rdf_FullHad.Filter("nFatJet == 2", "Events with 2 FatJet")

#rdf.Display({"nFatJet"}, 10).Print()
rdf_FullHad_2FatJet = rdf_FullHad_2FatJet.Define("WMask", "abs(GenPart_pdgId) == 24") \
             .Define("GenPart_pt_W", "GenPart_pt[WMask]") \
             .Define("GenPart_phi_W", "GenPart_phi[WMask]") \
             .Define("GenPart_eta_W", "GenPart_eta[WMask]")

columns = ["FatJet_pt","FatJet_eta","FatJet_phi","GenPart_pt","GenPart_eta","GenPart_phi","GenPart_pdgId","GenPart_genPartIdxMother"]

#rdf_limited = rdf_FullHad_2FatJet.Range(10)
#rdf_FullHad_2FatJet.Display(columns,10,20).Print()
display = rdf_FullHad_2FatJet.Display(columns,10)
print(display.AsString())

#-------------------------------------------------------------------------------

# Define C++ function to find specific GenParticles
# Here we find GenPart_pt of prompt leptons (status 1)
ROOT.gInterpreter.Declare("""
#include "ROOT/RVec.hxx"
#include <vector>

ROOT::RVecF getPromptLeptonPts(const ROOT::RVecF& pt, const ROOT::RVecI& status, const ROOT::RVecI& pdgId) {
    // Select prompt leptons (status 1, abs(pdgId) is 11 or 13)
    auto mask = (status == 1) && (abs(pdgId) == 11 || abs(pdgId) == 13);
    return pt[mask]; // Returns RVec of pts
}
""")

# Apply the function to the DataFrame
rdf_with_pts = rdf.Define("PromptLepton_pt", "getPromptLeptonPts(GenPart_pt, GenPart_status, GenPart_pdgId)")

# Create a histogram of the result
hist = rdf_with_pts.Histo1D(("PromptLepPt", "Prompt Lepton Pt;p_{T} [GeV];Events", 50, 0, 200), "PromptLepton_pt")

hist.Draw()

#-------------------------------------------------------------------------------

# Example 1: Select only GenParticles with status 1 and pt > 10 GeV
# GenPart_pt is typically a ROOT::RVec<float>
# rdf_gen = rdf.Filter("GenPart_status == 1 && GenPart_pt > 10.0")

# Example 2: Define a new column (e.g., Rapidity) using a lambda function
# The function acts on the whole RVec for each event
#rdf_with_rapidity = rdf_gen.Define("GenPart_rapidity",
#                                 "0.5 * log((sqrt(GenPart_mass*GenPart_mass + GenPart_pt*GenPart_pt*cosh(GenPart_eta)*cosh(GenPart_eta)) + GenPart_pt*sinh(GenPart_eta)) / (sqrt(GenPart_mass*GenPart_mass + GenPart_pt*GenPart_pt*cosh(GenPart_eta)*cosh(GenPart_eta)) - GenPart_pt*sinh(GenPart_eta)))")

# Example 3: Histogramming PT of specific GenParticles
#hist_gen_pt = rdf_with_rapidity.Histo1D(("gen_pt", "GenPart Pt;Pt [GeV];Counts", 100, 0, 100), "GenPart_pt")

#-------------------------------------------------------------------------------

# Save the histogram
outfile = ROOT.TFile("DarkHiggs_RD_histos.root", "RECREATE")
hist.Write()
#---
h_nQuarks.Write()
h_nWDecayProducts.Write()
h_nFatJet.Write()
#---
h_nFatJet_FullHad.Write()
h_nQuarks_FullHad.Write()
h_WQuark_inv_mass_FullHad.Write()
h_AK4_inv_mass_FullHad.Write()
h_MatchJet_inv_mass_FullHad.Write()
#---
outfile.Close()
