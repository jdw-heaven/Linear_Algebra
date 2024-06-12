#include <TCanvas.h>
#include <TGraph.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <TAxis.h>

void DrawGraph() {
    // 打开文件
    std::ifstream file("/home/heaven/Desktop/doc/Linear_Algebra/lanczos/draw/data/data12_13_2.txt");
    double x, y;
    TGraph *graph = new TGraph();

    if (!file.is_open()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    // 读取数据
    while (file >> x >> y) {
        graph->SetPoint(graph->GetN(), x, y);
    }
    graph->GetXaxis()->SetTitle("Delta");
    graph->GetYaxis()->SetTitle("Total Ground Energy");
    // 关闭文件
    file.close();

    // 设置标记样式和大小
    graph->SetMarkerStyle(20); // 设置标记样式，20 是一个常用的圆形标记
    graph->SetMarkerSize(0.2); // 设置标记大小，值越大点越大

    // 创建画布
    TCanvas *c1 = new TCanvas("c1", "Graph Example", 800, 600);
    // 设置网格样式
    c1->SetGrid(); // 开启网格
    c1->SetGridx(); // 网格水平线
    c1->SetGridy(); // 网格垂直线
    graph->SetTitle("Lanczos alogrithm, num = 12, error = 1e-13");
    // 绘制图形
    graph->Draw("ALP");


    // 保存画布为 ROOT 文件
    TFile *fileRoot = new TFile("/home/heaven/Desktop/doc/Linear_Algebra/lanczos/draw/graph.root", "RECREATE");
    c1->Write();
    fileRoot->Close();

    // 显示画布
    c1->Update();

    // 释放内存
    delete fileRoot;
    delete graph;
    delete c1;
}

int draw() {
    DrawGraph();
    return 0;
}