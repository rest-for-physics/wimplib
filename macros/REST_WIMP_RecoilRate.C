

#include <TRestWimpSensitivity.h>

void REST_WIMP_RecoilRate(const std::string& rmlFile, const double wimpMass = 1,
                          const double crossSection = 1E-45) {
    TRestWimpSensitivity WS(rmlFile.c_str());
    WS.PrintMetadata();

    auto recoilRate = WS.GetRecoilSpectra(wimpMass, crossSection);

    std::stringstream ss;
    ss << "WimpSensitivity_WimpMass_"<<wimpMass<<"_CrossSect_"<<crossSection;

    TCanvas* can = new TCanvas(ss.str().c_str(), ss.str().c_str());
    TLegend* leg = new TLegend(0.1, 0.7, 0.3, 0.9);

    int color = 1;

    for (const auto& [element, hist] : recoilRate) {
        hist->GetXaxis()->SetTitle("Energy (keV)");
        hist->GetYaxis()->SetTitle("Recoil rate (c/kev/day)");
        hist->SetLineColor(color);
        hist->SetLineWidth(3);
        leg->AddEntry(hist, element.c_str());
        can->cd();
        hist->Draw("SAME");
        color++;
    }

    can->cd();
    leg->Draw();
}
