

#include <TRestWimpSensitivity.h>

void REST_WIMP_QuenchingFactor(const std::string& rmlFile) {
    TRestWimpSensitivity WS(rmlFile.c_str());
    WS.PrintMetadata();
    WS.CalculateQuenchingFactor();

    auto quenchingFactor = WS.GetQuenchingFactor();

    TCanvas* can = new TCanvas("Quenching Factor", "Quenching Factor");
    TLegend* leg = new TLegend(0.1, 0.7, 0.3, 0.9);

    int color = 1;

    for (const auto& [element, hist] : quenchingFactor) {
        hist->GetXaxis()->SetTitle("Energy (keV)");
        hist->GetYaxis()->SetTitle("Quenching Factor");
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
