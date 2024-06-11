#include <TCanvas.h>
#include <TGraph.h>
#include <TFile.h>
#include <TRandom.h>

void createScatterPlot() {
    // 创建一个画布
    TCanvas *c1 = new TCanvas("c1", "Scatter Plot", 800, 600);

    // 生成随机数据点
    const int n = 100;
    double x[n], y[n];
    TRandom rand;
    for (int i = 0; i < n; ++i) {
        x[i] = rand.Uniform(0, 10);
        y[i] = rand.Uniform(0, 10);
    }

    // 创建散点图
    TGraph *graph = new TGraph(n, x, y);
    graph->SetTitle("Scatter Plot;X;Y");
    graph->SetMarkerStyle(20);

    // 绘制散点图
    graph->Draw("AP");

    // 保存画布为 ROOT 文件
    TFile *file = new TFile("scatterPlot.root", "RECREATE");
    c1->Write();
    file->Close();

    // 释放内存
    delete c1;
    delete graph;
    delete file;
}