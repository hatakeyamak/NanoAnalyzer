import ROOT

# Enable multi-threading
ROOT.ROOT.EnableImplicitMT()

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
# 1. Define quark selection: |pdgId| <= 5 and parent is W (24)
df_with_Ws = rdf.Define(
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

# 2. Count quarks per event
df_with_Ws = df_with_Ws.Define("nWQuarks", "Sum(GenPart_WParentQuark)")

# 3. Create histogram of count
h_nQuarks = df_with_Ws.Histo1D(("nQuarks", "Number of W-parent Quarks;N;Events", 5, 0, 5), "nWQuarks")

# Draw to show
c1 = ROOT.TCanvas()
h_nQuarks.Draw()
c1.SaveAs("nWQuarks.png")

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

# Make histogram
h_nWDecayProducts = rdf.Histo1D(
    ("h_nWDecayProducts", "W decay products;N;Events", 10, 0, 10),
    "nWDecayProducts"
)
h_nWDecayProducts.Draw()
c1.SaveAs("nWDecayProducts.png")
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
outfile = ROOT.TFile("output.root", "RECREATE")
hist.Write()
h_nQuarks.Write()
h_nWDecayProducts.Write()
outfile.Close()
