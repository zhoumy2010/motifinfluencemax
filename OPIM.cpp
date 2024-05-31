/**
* @file OPIM.cpp
* @brief This project implements the OPIM and OPIM-C for the following paper:
* Jing Tang, Xueyan Tang, Xiaokui Xiao, Junsong Yuan, "Online Processing Algorithms for Influence Maximization," in Proc. ACM SIGMOD, 2018.
*
* @author Jing Tang (Nanyang Technological University)
*
* Copyright (C) 2018 Jing Tang and Nanyang Technological University. All rights reserved.
*
*/

#include <iomanip>
#include "stdafx.h"
#include "SFMT/dSFMT/dSFMT.c"
#include "alg.cpp"

int main(int argc, char *argv[]) {
    for (int k = 0; k < 1; k++) {

        // Randomize the seed for generating random numbers
        dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
        const TArgument Arg(argc, argv);
        const std::string infilename = Arg._dir + "/" + Arg._graphname;
        std::vector<size_t> vecInDegRev;
//    Graph forward_graph;
        if (Arg._func == 0 || Arg._func == 2 || true) {
            // Format the graph
            vecInDegRev = GraphBase::format_graph(infilename, Arg._mode);
            if (Arg._func == 0) return 1;
        }


        std::cout << "---The Begin of " << Arg._outFileName << "---\n";
        Timer mainTimer("main");

        // Load the reverse graph
        Graph graph = GraphBase::load_graph(infilename, true, Arg._probDist, Arg._probEdge);

        int sum = 0;
        for (int i = 0; i < graph.size(); i++) {
            sum += graph[i].size();
        }

        int model = 0;
        if (Arg._model == LT) {
            model = 1;
            std::cout << "----------------LT------------------\n";
        } else {

            std::cout << "----------------IC------------------\n";
        }

        // 获取当前时间点
        auto now = std::chrono::system_clock::now();
        // 转换为 time_t 类型
        std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);
        // 使用 std::put_time 格式化日期和时间，替换空格为下划线
        std::tm *now_tm = std::localtime(&now_time_t);
        std::ostringstream oss;
        oss << std::put_time(now_tm, "%Y_%m_%d_%H_%M_%S");
        // 输出当前日期和时间的字符串
        std::string currentDateTime = oss.str();

        std::string fileName =
                Arg._graphname + "_" + Arg._groupSide + "_" + std::to_string(Arg._groupActiveNum) + "_" +
                std::to_string(model) + currentDateTime +
                ".txt";
        std::ofstream outputFile(fileName);
        if (!outputFile.is_open()) {
            std::cerr << "无法打开文件!" << std::endl;
            return 1;  // 返回错误码
        }

        outputFile << "Arg._func = " << Arg._func << std::endl;
        outputFile << "Arg._seedsize = " << Arg._seedsize << std::endl;
        outputFile << "Arg._samplesize = " << Arg._samplesize << std::endl;
        outputFile << "Arg._groupNum = " << Arg._groupNum << std::endl;
        outputFile << "Arg._groupSide = " << Arg._groupSide << std::endl;
        outputFile << "Arg._groupActiveNum = " << Arg._groupActiveNum << std::endl;
        outputFile << "Arg._probEdge = " << Arg._probEdge << std::endl;
        outputFile << "Arg._eps = " << Arg._eps << std::endl;
        outputFile << "Arg._delta = " << Arg._delta << std::endl;
        outputFile << "Arg._model = " << Arg._model << std::endl;
        outputFile << "Arg._graphname = " << Arg._graphname << std::endl;
        outputFile << "Arg._mode = " << Arg._mode << std::endl;
        outputFile << "Arg._dir = " << Arg._dir << std::endl;
        outputFile << "Arg._resultFolder = " << Arg._resultFolder << std::endl;
        outputFile << "Arg._algName = " << Arg._algName << std::endl;
        outputFile << "Arg._probDist = " << Arg._probDist << std::endl;
        outputFile << "Arg._outFileName = " << Arg._outFileName << std::endl;


        std::cout << "---The end of count_group_weight ---\n";

        if (Arg._model == LT) {
            // Normalize the propagation probabilities in accumulation format for LT cascade model for quickly generating RR sets
            to_normal_accum_prob(graph);
        }

        int mode = 2; // Default is to use the minimum upper bound among all the rounds
        if (Arg._mode == "0" || Arg._mode == "vanilla") {
            mode = 0;
        } else if (Arg._mode == "1" || Arg._mode == "last") {
            mode = 1;
        }
        auto delta = Arg._delta;

        int groupActiveNum = Arg._groupActiveNum;


//    for (Group &group: groups) {
//        GraphBase::count_group_weight(group, graph);
//    }
        std::cout << "  ==>Group loaded end!\n";
        for (auto i = 10; i <= 200; i += 10) {
            std::cout << "group init start: " << mainTimer.get_total_time() << '\n';
            if (delta < 0) delta = 1.0 / graph.size();
            Groups groups = GroupBase::init_group(graph, Arg._groupSide, Arg._groupNum);
            std::cout << "group init end: " << mainTimer.get_total_time() << '\n';
            outputFile << "\nseeds num : " << i << std::endl;
            // Initialize a result object to record the results
            TResult tRes;
            TResult tRes2;
            TResult tRes3;
            TAlg tAlg(graph, tRes);
            TAlg tAlg2(graph, tRes2);
            TAlg tAlg3(graph, tRes3);
            tAlg.set_cascade_model(Arg._model); // Set propagation model
            tAlg2.set_cascade_model(Arg._model); // Set propagation model
            tAlg3.set_cascade_model(Arg._model); // Set propagation model
            if (Arg._algName == "opim-c" || Arg._algName == "OPIM-C") {

                tAlg.opimc(i, Arg._eps, delta, graph, groups, outputFile, mode, model, groupActiveNum);
                tAlg2.LGMOIM(i, Arg._eps, vecInDegRev, graph, groups, outputFile, model, groupActiveNum);
                tAlg3.GIA(i, 0.3, 1.0 / (graph.size() * 1.0), graph, groups, outputFile, model, groupActiveNum);
            }
            TIO::write_result(Arg._outFileName, tRes, Arg._resultFolder);
            TIO::write_order_seeds(Arg._outFileName, tRes, Arg._resultFolder);
            std::cout << "---The End of " << Arg._outFileName << "---\n";
        }
        outputFile.close();
        return 0;
    }
}